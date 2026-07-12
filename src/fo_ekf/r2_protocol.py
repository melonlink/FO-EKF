"""Machine validation for the R2 finite-record data-to-disk protocol.

The validator checks protocol completeness, one-sided numerical propagation,
and evidence semantics.  A PASS certifies internal consistency with the
declared evidence; it does not prove that an external calibration artifact is
truthful or that a hash was published before data access.  It never interprets
an R2 failure as rejection of the fractional surrogate; model rejection belongs
to the downstream R1 joint-feasibility calculation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import date, datetime, time
from enum import Enum
from fractions import Fraction
from pathlib import Path
from typing import Any

import tomllib


class R2ProtocolStatus(str, Enum):
    """Allowed outcomes of protocol validation."""

    PASS_DETERMINISTIC = "PASS_DETERMINISTIC"
    PASS_PROBABILISTIC = "PASS_PROBABILISTIC"
    NOT_CERTIFIABLE = "NOT_CERTIFIABLE"
    EXCLUDE_PROTOCOL = "EXCLUDE_PROTOCOL"


@dataclass(frozen=True)
class R2ProtocolValidation:
    """Validation result with stable machine-readable reason codes."""

    status: R2ProtocolStatus
    reasons: tuple[str, ...]
    canonical_sha256: str

    @property
    def passed(self) -> bool:
        return self.status in {
            R2ProtocolStatus.PASS_DETERMINISTIC,
            R2ProtocolStatus.PASS_PROBABILISTIC,
        }

    def as_dict(self) -> dict[str, Any]:
        return {
            "status": self.status.value,
            "reasons": list(self.reasons),
            "canonical_sha256": self.canonical_sha256,
        }


BOUND_NAMES = (
    "history",
    "dwell",
    "hrv",
    "phase_anchor",
    "window_leakage",
    "sampling",
    "delay",
    "filter_initialization",
    "measurement_noise",
    "model_residual",
)

DETERMINISTIC_SOURCE_KINDS = (
    "sensor_spec",
    "independent_bounded_calibration",
    "proved_engineering_envelope",
    "observable_window_geometry",
)

PROBABILISTIC_SOURCE_KINDS = (
    "gaussian_model",
    "bootstrap",
    "empirical_quantile",
    "confidence_interval",
    "prediction_interval",
)

COVERAGE_SCOPE = "simultaneous_all_selected_rate_harmonic_disks"

_RADIUS_FIELDS = (
    ("history", "radius_history"),
    ("dwell", "radius_dwell"),
    ("hrv", "radius_hrv"),
    ("phase_anchor", "radius_phase_anchor"),
    ("window_leakage", "radius_window_leakage"),
    ("sampling", "radius_sampling"),
    ("sampling", "radius_quadrature"),
    ("delay", "radius_delay"),
    ("filter_initialization", "radius_filter_initialization"),
    ("measurement_noise", "radius_measurement_noise"),
    ("model_residual", "radius_model_residual"),
)

_RADIUS_SOURCE_FIELDS = (
    ("history", "radius_history", "computed_w_norm_bound", True),
    ("dwell", "radius_dwell", "computed_w_norm_bound", True),
    ("hrv", "radius_hrv", "computed_w_norm_bound", True),
    ("phase_anchor", "radius_phase_anchor", "computed_w_norm_bound", True),
    (
        "window_leakage",
        "radius_window_leakage",
        "computed_coefficient_radius",
        False,
    ),
    ("sampling", "radius_sampling", "computed_w_norm_bound", True),
    ("sampling", "radius_quadrature", "quadrature_radius", False),
    ("delay", "radius_delay", "computed_coefficient_radius", False),
    (
        "filter_initialization",
        "radius_filter_initialization",
        "computed_coefficient_radius",
        False,
    ),
    (
        "measurement_noise",
        "radius_measurement_noise",
        "computed_w_norm_bound",
        True,
    ),
    ("model_residual", "radius_model_residual", "computed_w_norm_bound", True),
)

_MISSING = object()
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_GIT_COMMIT_RE = re.compile(r"^(?:[0-9a-f]{40}|[0-9a-f]{64})$")


def load_r2_protocol(path: str | Path) -> dict[str, Any]:
    """Load an R2 TOML document using the standard-library parser."""

    with Path(path).open("rb") as stream:
        return tomllib.load(stream)


def _canonical_object(
    value: Any,
    *,
    excluded_keys: frozenset[str] = frozenset(),
) -> Any:
    if isinstance(value, Mapping):
        return {
            str(key): _canonical_object(item, excluded_keys=excluded_keys)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
            if str(key) not in excluded_keys
        }
    if isinstance(value, (list, tuple)):
        return [_canonical_object(item, excluded_keys=excluded_keys) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        if math.isnan(value):
            marker = "nan"
        elif value > 0.0:
            marker = "+inf"
        else:
            marker = "-inf"
        return {"__nonfinite_float__": marker}
    if isinstance(value, (datetime, date, time)):
        return {"__toml_datetime__": value.isoformat()}
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"unsupported canonical value: {type(value).__name__}")


def canonical_protocol_sha256(config: Mapping[str, Any]) -> str:
    """Hash only the preregistered protocol, excluding run-time results."""

    protocol_only = {
        str(key): value for key, value in config.items() if str(key) != "window_result"
    }

    payload = json.dumps(
        _canonical_object(protocol_only, excluded_keys=frozenset({"config_sha256"})),
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def canonical_window_result_sha256(result: Mapping[str, Any]) -> str:
    """Hash one result row, excluding only its self-referential result hash."""

    payload = json.dumps(
        _canonical_object(result, excluded_keys=frozenset({"result_sha256"})),
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _get(config: Mapping[str, Any], *path: str, default: Any = _MISSING) -> Any:
    current: Any = config
    for part in path:
        if not isinstance(current, Mapping) or part not in current:
            return default
        current = current[part]
    return current


def _table(
    config: Mapping[str, Any],
    name: str,
    reasons: list[str],
) -> Mapping[str, Any]:
    value = config.get(name, _MISSING)
    if not isinstance(value, Mapping):
        reasons.append(f"{name}:missing_table")
        return {}
    return value


def _is_nonempty(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip()) and value.strip().upper() != "UNSET"


def _is_finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _is_finite_nonnegative(value: Any) -> bool:
    return _is_finite_number(value) and value >= 0


def _upward(value: float) -> float:
    """Move a positive binary float outward by one representable value."""

    if value == 0.0:
        return 0.0
    return math.nextafter(value, math.inf)


def _upward_sum(values: Sequence[float]) -> float:
    return _upward(math.fsum(values))


def _upward_product(left: float, right: float) -> float:
    if left == 0.0 or right == 0.0:
        return 0.0
    return _upward(left * right)


def _upward_ratio(numerator: float, denominator: float) -> float:
    return _upward(numerator / denominator)


def _downward(value: float) -> float:
    """Move a positive binary float inward toward a safe lower bound."""

    if value == 0.0:
        return 0.0
    return math.nextafter(value, -math.inf)


def _downward_product(left: float, right: float) -> float:
    if left == 0.0 or right == 0.0:
        return 0.0
    return _downward(left * right)


def _downward_ratio(numerator: float, denominator: float) -> float:
    return _downward(numerator / denominator)


def _fraction_lower_float(value: Fraction) -> float:
    """Convert an exact binary-rational expression to a containing lower float."""

    candidate = float(value)
    if Fraction.from_float(candidate) > value:
        return math.nextafter(candidate, -math.inf)
    return candidate


def _require_true(table: Mapping[str, Any], key: str, prefix: str, reasons: list[str]) -> None:
    if table.get(key, _MISSING) is not True:
        reasons.append(f"{prefix}.{key}:must_be_true")


def _require_false(
    table: Mapping[str, Any],
    key: str,
    prefix: str,
    reasons: list[str],
) -> None:
    if table.get(key, _MISSING) is not False:
        reasons.append(f"{prefix}.{key}:must_be_false")


def _require_nonempty(
    table: Mapping[str, Any],
    key: str,
    prefix: str,
    reasons: list[str],
) -> None:
    if not _is_nonempty(table.get(key, _MISSING)):
        reasons.append(f"{prefix}.{key}:missing")


def _require_sha256(
    table: Mapping[str, Any],
    key: str,
    prefix: str,
    reasons: list[str],
) -> None:
    value = table.get(key, _MISSING)
    if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
        reasons.append(f"{prefix}.{key}:invalid_sha256")


def _validate_interval(
    value: Any,
    *,
    lower_open: float,
    upper_closed: float | None,
    prefix: str,
    reasons: list[str],
) -> None:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or not all(_is_finite_number(item) for item in value)
    ):
        reasons.append(f"{prefix}:invalid_interval")
        return
    lower, upper = value
    if lower <= lower_open or upper < lower:
        reasons.append(f"{prefix}:invalid_interval")
    if upper_closed is not None and upper > upper_closed:
        reasons.append(f"{prefix}:invalid_interval")


def _iter_numeric(
    value: Any,
    path: tuple[str, ...] = (),
) -> Sequence[tuple[tuple[str, ...], int | float]]:
    found: list[tuple[tuple[str, ...], int | float]] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            child_path = (*path, str(key))
            if str(key) == "explicit_omitted_harmonics":
                continue
            found.extend(_iter_numeric(item, child_path))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            found.extend(_iter_numeric(item, (*path, str(index))))
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        found.append((path, value))
    return found


def _validate_exclusions(config: Mapping[str, Any]) -> list[str]:
    reasons: list[str] = []
    preprocessing = _get(config, "preprocessing", default={})
    window = _get(config, "window_design", default={})
    estimator = _get(config, "estimator", default={})

    if isinstance(preprocessing, Mapping):
        if preprocessing.get("all_steps_causal", _MISSING) is False:
            reasons.append("preprocessing.all_steps_causal:explicitly_false")
        if preprocessing.get("uses_filtfilt_or_centered_smoothing", _MISSING) is True:
            reasons.append("preprocessing:acausal_filter_enabled")
    if isinstance(window, Mapping):
        if window.get("per_window_amplitude_normalization", _MISSING) is True:
            reasons.append("window_design:per_window_amplitude_normalization")
        if window.get("per_window_free_phase_alignment", _MISSING) is True:
            reasons.append("window_design:per_window_free_phase_alignment")
    if isinstance(estimator, Mapping):
        if estimator.get("uses_samples_after_right_anchor", _MISSING) is True:
            reasons.append("estimator:uses_future_samples")
    return reasons


def _validate_freeze_and_hash(
    config: Mapping[str, Any],
    protocol: Mapping[str, Any],
    canonical_sha256: str,
    reasons: list[str],
) -> None:
    _require_nonempty(protocol, "protocol_id", "protocol", reasons)
    _require_nonempty(protocol, "version", "protocol", reasons)
    _require_nonempty(protocol, "title", "protocol", reasons)
    _require_nonempty(protocol, "owner", "protocol", reasons)
    _require_true(protocol, "frozen", "protocol", reasons)
    _require_true(protocol, "frozen_before_target_ecg_access", "protocol", reasons)
    frozen_at = protocol.get("frozen_at_utc", _MISSING)
    if not _is_nonempty(frozen_at):
        reasons.append("protocol.frozen_at_utc:missing")
    else:
        try:
            parsed = datetime.fromisoformat(frozen_at.replace("Z", "+00:00"))
        except ValueError:
            reasons.append("protocol.frozen_at_utc:invalid")
        else:
            if parsed.tzinfo is None or parsed.utcoffset().total_seconds() != 0:
                reasons.append("protocol.frozen_at_utc:not_utc")
    git_commit = protocol.get("git_commit", _MISSING)
    if not isinstance(git_commit, str) or _GIT_COMMIT_RE.fullmatch(git_commit) is None:
        reasons.append("protocol.git_commit:invalid")

    stored_hash = protocol.get("config_sha256", _MISSING)
    if not isinstance(stored_hash, str) or _SHA256_RE.fullmatch(stored_hash) is None:
        reasons.append("protocol.config_sha256:invalid_sha256")
    elif stored_hash != canonical_sha256:
        reasons.append("protocol.config_sha256:mismatch")

    statuses = protocol.get("allowed_result_statuses", _MISSING)
    allowed = {status.value for status in R2ProtocolStatus}
    if (
        not isinstance(statuses, list)
        or not all(isinstance(status, str) for status in statuses)
        or len(statuses) != len(allowed)
        or set(statuses) != allowed
    ):
        reasons.append("protocol.allowed_result_statuses:invalid")
    if isinstance(statuses, list) and "REJECT_MODEL" in statuses:
        reasons.append("protocol.allowed_result_statuses:reject_forbidden")


def _validate_data_boundary(table: Mapping[str, Any], reasons: list[str]) -> None:
    for key, expected in (
        ("raw_data_root_env", "ECG_DATA_DIR"),
        ("processed_data_root_env", "ECG_PROCESSED_DIR"),
        ("output_root_env", "ECG_OUTPUT_DIR"),
    ):
        if table.get(key, _MISSING) != expected:
            reasons.append(f"data_boundary.{key}:invalid")
    _require_false(table, "copy_raw_ecg_into_repository", "data_boundary", reasons)
    _require_true(table, "patient_level_split_required", "data_boundary", reasons)
    _require_true(
        table,
        "calibration_disjoint_from_identification",
        "data_boundary",
        reasons,
    )
    split_ids: list[str] = []
    for key in ("calibration_split_id", "identification_split_id", "validation_split_id"):
        _require_nonempty(table, key, "data_boundary", reasons)
        value = table.get(key)
        if _is_nonempty(value):
            split_ids.append(value)
    if len(split_ids) == 3 and len(set(split_ids)) != 3:
        reasons.append("data_boundary:split_ids_not_distinct")
    _require_sha256(table, "split_manifest_sha256", "data_boundary", reasons)


def _validate_units(table: Mapping[str, Any], reasons: list[str]) -> None:
    _require_nonempty(table, "signal_unit", "units", reasons)
    for key, expected in (
        ("physical_time_unit", "s"),
        ("phase_unit", "rad"),
        ("frequency_unit", "rad/s"),
        ("dimensionless_time_symbol", "xi"),
    ):
        if table.get(key, _MISSING) != expected:
            reasons.append(f"units.{key}:invalid")
    _require_true(table, "tau_star_provided", "units", reasons)
    tau_star = table.get("tau_star_s", _MISSING)
    if not _is_finite_number(tau_star) or tau_star <= 0:
        reasons.append("units.tau_star_s:must_be_positive_finite")
    _require_true(table, "lambda_dimensionless", "units", reasons)
    _require_true(table, "q_z_epsilon_same_signal_unit", "units", reasons)


def _validate_model(
    table: Mapping[str, Any],
    harmonics: list[int] | None,
    signal_unit: Any,
    reasons: list[str],
) -> None:
    _require_true(table, "alpha_interval_provided", "model", reasons)
    _validate_interval(
        table.get("alpha_interval", _MISSING),
        lower_open=0.0,
        upper_closed=1.0,
        prefix="model.alpha_interval",
        reasons=reasons,
    )
    _require_true(table, "lambda_interval_provided", "model", reasons)
    _validate_interval(
        table.get("lambda_interval", _MISSING),
        lower_open=0.0,
        upper_closed=None,
        prefix="model.lambda_interval",
        reasons=reasons,
    )
    for key in (
        "principal_fractional_branch",
        "common_alpha_across_rates",
        "common_positive_lambda_or_preregistered_drift",
        "common_q_or_preregistered_restitution",
        "q_bounds_independent_of_target_fit",
    ):
        _require_true(table, key, "model", reasons)
    _require_nonempty(table, "q_magnitude_bounds_source", "model", reasons)
    q_bounds_unit = table.get("q_bounds_unit", _MISSING)
    if q_bounds_unit != "same_as_signal" and q_bounds_unit != signal_unit:
        reasons.append("model.q_bounds_unit:unresolved")

    lowers = table.get("q_abs_lower_by_harmonic", _MISSING)
    uppers = table.get("q_abs_upper_by_harmonic", _MISSING)
    if (
        harmonics is None
        or not isinstance(lowers, list)
        or not isinstance(uppers, list)
        or len(lowers) != len(harmonics)
        or len(uppers) != len(harmonics)
    ):
        reasons.append("model.q_bounds:wrong_length")
        return
    for harmonic, lower, upper in zip(harmonics, lowers, uppers, strict=True):
        valid = (
            _is_finite_nonnegative(lower)
            and _is_finite_number(upper)
            and upper > 0
            and lower <= upper
        )
        if harmonic != 0:
            valid = valid and lower > 0
        if not valid:
            reasons.append(f"model.q_bounds[{harmonic}]:invalid")


def _validate_estimator(
    table: Mapping[str, Any],
    reasons: list[str],
) -> list[int] | None:
    if table.get("kind", _MISSING) != "delayed_causal_weighted_least_squares":
        reasons.append("estimator.kind:invalid")
    if table.get("phase_model", _MISSING) != "constant_rate_between_window_anchors":
        reasons.append("estimator.phase_model:invalid")
    _require_true(table, "output_only_after_right_r_peak", "estimator", reasons)
    _require_false(table, "uses_samples_after_right_anchor", "estimator", reasons)
    _require_true(table, "retained_harmonics_provided", "estimator", reasons)
    _require_true(table, "weights_sum_to_one", "estimator", reasons)
    _require_true(table, "gram_positive_definite_required", "estimator", reasons)
    _require_true(table, "gram_thresholds_provided", "estimator", reasons)
    _require_true(table, "coefficient_level_filter_correction", "estimator", reasons)
    _require_true(table, "require_nonzero_nominal_transfer", "estimator", reasons)
    _require_nonempty(table, "weight_rule", "estimator", reasons)

    min_eigenvalue = table.get("gram_min_eigenvalue_threshold", _MISSING)
    max_condition = table.get("gram_condition_number_max", _MISSING)
    if not _is_finite_number(min_eigenvalue) or min_eigenvalue <= 0:
        reasons.append("estimator.gram_min_eigenvalue_threshold:invalid")
    if not _is_finite_number(max_condition) or max_condition < 1:
        reasons.append("estimator.gram_condition_number_max:invalid")

    harmonics = table.get("retained_harmonics", _MISSING)
    if (
        not isinstance(harmonics, list)
        or not harmonics
        or any(not isinstance(item, int) or isinstance(item, bool) for item in harmonics)
        or len(set(harmonics)) != len(harmonics)
        or not any(item != 0 for item in harmonics)
    ):
        reasons.append("estimator.retained_harmonics:invalid")
        return None

    representation = table.get("signal_representation", _MISSING)
    if representation not in {"real_two_sided_fourier", "delayed_causal_analytic"}:
        reasons.append("estimator.signal_representation:invalid")
    elif representation == "real_two_sided_fourier":
        _require_true(table, "include_dc_if_real_signal", "estimator", reasons)
        _require_true(table, "include_negative_harmonics_if_real_signal", "estimator", reasons)
        if 0 not in harmonics or any(-item not in harmonics for item in harmonics):
            reasons.append("estimator.retained_harmonics:not_two_sided")
    return harmonics


def _validate_window(table: Mapping[str, Any], reasons: list[str]) -> None:
    if table.get("left_anchor", _MISSING) != "R_peak":
        reasons.append("window_design.left_anchor:must_be_R_peak")
    if table.get("right_anchor", _MISSING) != "R_peak":
        reasons.append("window_design.right_anchor:must_be_R_peak")
    _require_false(
        table,
        "fixed_seconds_window_allowed",
        "window_design",
        reasons,
    )
    for key in ("same_subject_required", "same_lead_required", "same_physical_gain_required"):
        _require_true(table, key, "window_design", reasons)
    _require_false(
        table,
        "window_selection_uses_target_fit_quality",
        "window_design",
        reasons,
    )
    _require_true(
        table,
        "complete_rr_intervals_provided",
        "window_design",
        reasons,
    )
    intervals = table.get("complete_rr_intervals", _MISSING)
    if not isinstance(intervals, int) or isinstance(intervals, bool) or intervals <= 0:
        reasons.append("window_design.complete_rr_intervals:invalid")
    for flag, key in (
        ("burn_in_provided", "burn_in_s_min"),
        ("dwell_provided", "dwell_s_min"),
        ("phase_deviation_threshold_provided", "maximum_observed_rr_phase_deviation_rad"),
    ):
        _require_true(table, flag, "window_design", reasons)
        if not _is_finite_nonnegative(table.get(key, _MISSING)):
            reasons.append(f"window_design.{key}:invalid")


def _validate_preprocessing(table: Mapping[str, Any], reasons: list[str]) -> None:
    _require_true(table, "all_steps_causal", "preprocessing", reasons)
    _require_false(
        table,
        "uses_filtfilt_or_centered_smoothing",
        "preprocessing",
        reasons,
    )
    _require_true(
        table,
        "same_chain_and_state_policy_across_windows",
        "preprocessing",
        reasons,
    )
    _require_true(
        table,
        "relative_ecg_rpeak_timing_calibrated",
        "preprocessing",
        reasons,
    )
    _require_nonempty(table, "pipeline_manifest", "preprocessing", reasons)
    _require_sha256(table, "pipeline_manifest_sha256", "preprocessing", reasons)
    if table.get("filter_state_policy", _MISSING) not in {
        "continuous_from_record_start",
        "restored_certified_state",
        "zero_state_with_tail_bound",
    }:
        reasons.append("preprocessing.filter_state_policy:invalid")
    if table.get("identity_front_end", _MISSING) is not True:
        _require_true(
            table,
            "nominal_transfer_calibrated",
            "preprocessing",
            reasons,
        )
        _require_nonempty(table, "nominal_transfer_source", "preprocessing", reasons)


def _validate_evidence_policy(table: Mapping[str, Any], reasons: list[str]) -> None:
    for key in (
        "target_fit_residuals_may_set_bounds",
        "post_hoc_bound_tuning_allowed",
        "post_hoc_window_exclusion_allowed",
    ):
        _require_false(table, key, "evidence_policy", reasons)
    _require_true(
        table,
        "probabilistic_terms_require_failure_probability",
        "evidence_policy",
        reasons,
    )
    if table.get("bound_numeric_scope", _MISSING) != "uniform_all_windows_satisfying_frozen_design":
        reasons.append("evidence_policy.bound_numeric_scope:invalid")
    deterministic = table.get("deterministic_source_kinds", _MISSING)
    probabilistic = table.get("probabilistic_source_kinds", _MISSING)
    if (
        not isinstance(deterministic, list)
        or not all(isinstance(item, str) for item in deterministic)
        or len(deterministic) != len(DETERMINISTIC_SOURCE_KINDS)
        or set(deterministic) != set(DETERMINISTIC_SOURCE_KINDS)
    ):
        reasons.append("evidence_policy.deterministic_source_kinds:not_fixed_allowlist")
    if (
        not isinstance(probabilistic, list)
        or not all(isinstance(item, str) for item in probabilistic)
        or len(probabilistic) != len(PROBABILISTIC_SOURCE_KINDS)
        or set(probabilistic) != set(PROBABILISTIC_SOURCE_KINDS)
    ):
        reasons.append("evidence_policy.probabilistic_source_kinds:not_fixed_allowlist")
    if table.get("joint_probability_method", _MISSING) not in {
        "deterministic_zero",
        "union_bound",
        "directly_calibrated_joint_bound",
    }:
        reasons.append("evidence_policy.joint_probability_method:invalid")


def _validate_audit_outputs(table: Mapping[str, Any], reasons: list[str]) -> None:
    _require_false(table, "write_raw_ecg", "audit_outputs", reasons)
    for key in (
        "write_window_manifest",
        "write_bound_breakdown",
        "write_gram_diagnostics",
        "write_filter_transfer_diagnostics",
        "write_probability_accounting",
        "write_input_and_result_hashes",
        "preserve_excluded_and_not_certifiable_rows",
    ):
        _require_true(table, key, "audit_outputs", reasons)


def _validate_bounds(
    bounds: Mapping[str, Any],
    signal_unit: Any,
    requested_semantics: Any,
    calibration_split_id: Any,
    identification_split_id: Any,
    validation_split_id: Any,
    reasons: list[str],
) -> list[float]:
    deterministic_sources = frozenset(DETERMINISTIC_SOURCE_KINDS)
    probabilistic_sources = frozenset(PROBABILISTIC_SOURCE_KINDS)
    probabilities: list[float] = []
    semantics_seen: list[str] = []

    for name in BOUND_NAMES:
        prefix = f"bounds.{name}"
        table = bounds.get(name, _MISSING)
        if not isinstance(table, Mapping):
            reasons.append(f"{prefix}:missing_table")
            probabilities.append(math.nan)
            continue
        for key in ("required", "applicability_decided", "applicable", "provided"):
            _require_true(table, key, prefix, reasons)
        _require_nonempty(table, "source_reference", prefix, reasons)
        _require_nonempty(table, "calibration_method", prefix, reasons)
        _require_true(table, "independent_of_identification_fit", prefix, reasons)
        _require_false(table, "uses_target_fit_residuals", prefix, reasons)

        if "calibration_split_id" in table:
            bound_split = table.get("calibration_split_id", _MISSING)
            if (
                _is_nonempty(identification_split_id) and bound_split == identification_split_id
            ) or (_is_nonempty(validation_split_id) and bound_split == validation_split_id):
                reasons.append(f"{prefix}.calibration_split_id:data_leakage")
            if bound_split != "not_applicable" and not (
                _is_nonempty(calibration_split_id) and bound_split == calibration_split_id
            ):
                reasons.append(f"{prefix}.calibration_split_id:invalid")

        semantics = table.get("semantics", _MISSING)
        source_kind = table.get("source_kind", _MISSING)
        if semantics not in {"deterministic", "probabilistic"}:
            reasons.append(f"{prefix}.semantics:invalid")
        else:
            semantics_seen.append(semantics)
            valid_sources = (
                deterministic_sources if semantics == "deterministic" else probabilistic_sources
            )
            if not isinstance(source_kind, str) or source_kind not in valid_sources:
                reasons.append(f"{prefix}.source_kind:incompatible_with_semantics")

        unit = table.get("radius_unit", _MISSING)
        if unit != "same_as_signal" and unit != signal_unit:
            reasons.append(f"{prefix}.radius_unit:unresolved")

        for numeric_path, value in _iter_numeric(table):
            if not _is_finite_number(value):
                reasons.append(f"{prefix}.{'.'.join(numeric_path)}:nonfinite")
            elif value < 0:
                reasons.append(f"{prefix}.{'.'.join(numeric_path)}:negative")

        probability = table.get("failure_probability", _MISSING)
        if not _is_finite_nonnegative(probability) or probability >= 1:
            reasons.append(f"{prefix}.failure_probability:invalid")
            probabilities.append(math.nan)
        else:
            probabilities.append(float(probability))
            if semantics == "deterministic" and probability != 0:
                reasons.append(f"{prefix}.failure_probability:deterministic_must_be_zero")
            if semantics == "probabilistic" and probability <= 0:
                reasons.append(f"{prefix}.failure_probability:probabilistic_must_be_positive")

    if requested_semantics == "deterministic":
        if any(item != "deterministic" for item in semantics_seen) or len(semantics_seen) != len(
            BOUND_NAMES
        ):
            reasons.append("protocol.requested_output_semantics:contains_probabilistic_bound")
    elif requested_semantics == "probabilistic":
        if "probabilistic" not in semantics_seen:
            reasons.append("protocol.requested_output_semantics:no_probabilistic_bound")
    else:
        reasons.append("protocol.requested_output_semantics:invalid")
    return probabilities


def _validate_bound_specific_gates(
    bounds: Mapping[str, Any],
    reasons: list[str],
) -> None:
    history = bounds.get("history", {})
    if isinstance(history, Mapping) and history.get("input_and_state_unit") != "same_as_signal":
        reasons.append("bounds.history.input_and_state_unit:unresolved")

    dwell = bounds.get("dwell", {})
    if isinstance(dwell, Mapping) and dwell.get("input_unit") != "same_as_signal":
        reasons.append("bounds.dwell.input_unit:unresolved")

    hrv = bounds.get("hrv", {})
    if isinstance(hrv, Mapping):
        _require_true(
            hrv,
            "phase_defined_operationally_from_r_peaks",
            "bounds.hrv",
            reasons,
        )
        _require_true(hrv, "r_peak_timing_error_provided", "bounds.hrv", reasons)
        _require_true(hrv, "q_harmonic_envelope_provided", "bounds.hrv", reasons)
        if (
            hrv.get("observed_rr_phase_deviation_rule", _MISSING)
            != "sup_unwrapped_difference_from_nominal_rate"
        ):
            reasons.append("bounds.hrv.observed_rr_phase_deviation_rule:invalid")
        if hrv.get("latent_physiological_phase_claimed", _MISSING) is True:
            _require_true(hrv, "latent_phase_bound_provided", "bounds.hrv", reasons)

    phase = bounds.get("phase_anchor", {})
    if isinstance(phase, Mapping):
        _require_true(
            phase,
            "common_radian_gauge_absorbed_into_q",
            "bounds.phase_anchor",
            reasons,
        )
        _require_true(
            phase,
            "phase_error_disjoint_from_hrv_term",
            "bounds.phase_anchor",
            reasons,
        )

    leakage = bounds.get("window_leakage", {})
    if isinstance(leakage, Mapping):
        _require_true(
            leakage,
            "noninteger_phase_length_retained",
            "bounds.window_leakage",
            reasons,
        )
        for key in (
            "alias_factors_computed",
            "harmonic_tail_bound_provided",
            "front_end_transfer_envelope_provided",
        ):
            _require_true(leakage, key, "bounds.window_leakage", reasons)

    sampling = bounds.get("sampling", {})
    if isinstance(sampling, Mapping):
        _require_true(sampling, "sample_rate_provided", "bounds.sampling", reasons)
        sample_rate = sampling.get("sample_rate_hz", _MISSING)
        if not _is_finite_number(sample_rate) or sample_rate <= 0:
            reasons.append("bounds.sampling.sample_rate_hz:must_be_positive_finite")
        _require_true(sampling, "anti_alias_certified", "bounds.sampling", reasons)
        _require_true(
            sampling,
            "first_derivative_bound_provided",
            "bounds.sampling",
            reasons,
        )
        implementation = sampling.get("implementation", _MISSING)
        if implementation not in {"discrete_wls", "trapezoidal_continuous_lockin"}:
            reasons.append("bounds.sampling.implementation:invalid")
        if implementation == "trapezoidal_continuous_lockin":
            _require_true(
                sampling,
                "second_derivative_bound_provided",
                "bounds.sampling",
                reasons,
            )
        if sampling.get("missing_sample_policy", _MISSING) != "reject_window":
            reasons.append("bounds.sampling.missing_sample_policy:invalid")

    delay = bounds.get("delay", {})
    if isinstance(delay, Mapping):
        _require_true(
            delay,
            "ecg_to_rpeak_channel_alignment_included",
            "bounds.delay",
            reasons,
        )
        for key in (
            "nominal_transfer_magnitude_lower_bound",
            "target_harmonic_abs_upper_bound",
        ):
            value = delay.get(key, _MISSING)
            if not _is_finite_number(value) or value <= 0:
                reasons.append(f"bounds.delay.{key}:must_be_positive_finite")

    measurement = bounds.get("measurement_noise", {})
    if isinstance(measurement, Mapping):
        bound_type = measurement.get("bound_type", _MISSING)
        field_by_type = {
            "weighted_l2": "weighted_l2_bound",
            "continuous_l2": "continuous_l2_bound_signal_sqrt_s",
            "linfinity": "linfinity_bound",
        }
        if bound_type not in field_by_type:
            reasons.append("bounds.measurement_noise.bound_type:invalid")
        elif not _is_finite_nonnegative(measurement.get(field_by_type[bound_type], _MISSING)):
            reasons.append(f"bounds.measurement_noise.{field_by_type[bound_type]}:invalid")

    model = bounds.get("model_residual", {})
    if isinstance(model, Mapping):
        location = model.get("residual_location", _MISSING)
        field_by_location = {
            "output": "output_weighted_l2_bound",
            "dynamics_input": "dynamics_input_sup_bound",
        }
        if location not in field_by_location:
            reasons.append("bounds.model_residual.residual_location:invalid")
        elif not _is_finite_nonnegative(model.get(field_by_location[location], _MISSING)):
            reasons.append(f"bounds.model_residual.{field_by_location[location]}:invalid")


def _validate_thresholds(
    table: Mapping[str, Any],
    requested_semantics: Any,
    window_design: Mapping[str, Any],
    reasons: list[str],
) -> None:
    for flag in (
        "thresholds_provided",
        "target_lower_bound_provided",
        "minimum_energy_snr_provided",
        "probability_threshold_provided",
        "dwell_thresholds_provided",
    ):
        _require_true(table, flag, "decision_thresholds", reasons)
    for key in (
        "deterministic_total_radius_max",
        "deterministic_relative_radius_max",
        "target_harmonic_abs_lower_bound",
        "minimum_energy_snr",
    ):
        value = table.get(key, _MISSING)
        if not _is_finite_number(value) or value <= 0:
            reasons.append(f"decision_thresholds.{key}:must_be_positive_finite")
    for key in ("minimum_burn_in_s", "minimum_dwell_s"):
        if not _is_finite_nonnegative(table.get(key, _MISSING)):
            reasons.append(f"decision_thresholds.{key}:invalid")
    for threshold_key, design_key in (
        ("minimum_burn_in_s", "burn_in_s_min"),
        ("minimum_dwell_s", "dwell_s_min"),
    ):
        threshold_value = table.get(threshold_key, _MISSING)
        design_value = window_design.get(design_key, _MISSING)
        if (
            _is_finite_nonnegative(threshold_value)
            and _is_finite_nonnegative(design_value)
            and threshold_value < design_value
        ):
            reasons.append(f"decision_thresholds.{threshold_key}:below_window_design")
    max_probability = table.get("maximum_joint_failure_probability", _MISSING)
    if not _is_finite_nonnegative(max_probability) or max_probability >= 1:
        reasons.append("decision_thresholds.maximum_joint_failure_probability:invalid")
    elif requested_semantics == "probabilistic" and max_probability <= 0:
        reasons.append("decision_thresholds.maximum_joint_failure_probability:must_be_positive")


def _validate_probability_accounting(
    table: Mapping[str, Any],
    expected_probabilities: list[float],
    requested_semantics: Any,
    maximum_probability: Any,
    evidence_method: Any,
    calibration_split_id: Any,
    reasons: list[str],
) -> None:
    _require_true(
        table,
        "joint_failure_probability_provided",
        "probability_accounting",
        reasons,
    )
    _require_true(
        table,
        "deterministic_terms_have_zero_failure_probability",
        "probability_accounting",
        reasons,
    )
    _require_true(
        table,
        "probabilistic_result_label_required",
        "probability_accounting",
        reasons,
    )
    if table.get("coverage_scope", _MISSING) != COVERAGE_SCOPE:
        reasons.append("probability_accounting.coverage_scope:invalid")
    method = table.get("simultaneous_coverage_method", _MISSING)
    if method not in {"deterministic_zero", "union_bound", "directly_calibrated_joint_bound"}:
        reasons.append("probability_accounting.simultaneous_coverage_method:invalid")
    if method != evidence_method:
        reasons.append("probability_accounting.simultaneous_coverage_method:not_preregistered")
    if method == "directly_calibrated_joint_bound":
        _require_nonempty(
            table,
            "direct_joint_reference",
            "probability_accounting",
            reasons,
        )
        _require_true(
            table,
            "direct_joint_independent_of_identification_fit",
            "probability_accounting",
            reasons,
        )
        _require_sha256(
            table,
            "direct_joint_evidence_sha256",
            "probability_accounting",
            reasons,
        )
        if table.get("direct_joint_calibration_split_id", _MISSING) != calibration_split_id:
            reasons.append("probability_accounting.direct_joint_calibration_split_id:invalid")

    reported = table.get("term_failure_probabilities", _MISSING)
    if (
        not isinstance(reported, list)
        or len(reported) != len(BOUND_NAMES)
        or any(not _is_finite_nonnegative(item) or item >= 1 for item in reported)
    ):
        reasons.append("probability_accounting.term_failure_probabilities:invalid")
    elif any(
        not math.isfinite(expected) or actual != expected
        for actual, expected in zip(reported, expected_probabilities, strict=True)
    ):
        reasons.append("probability_accounting.term_failure_probabilities:mismatch")

    joint = table.get("joint_failure_probability", _MISSING)
    if not _is_finite_nonnegative(joint) or joint >= 1:
        reasons.append("probability_accounting.joint_failure_probability:invalid")
        return
    if _is_finite_nonnegative(maximum_probability) and joint > maximum_probability:
        reasons.append("probability_accounting.joint_failure_probability:above_threshold")

    if requested_semantics == "deterministic":
        if joint != 0 or any(probability != 0 for probability in expected_probabilities):
            reasons.append("probability_accounting:deterministic_probability_nonzero")
        if method != "deterministic_zero":
            reasons.append("probability_accounting:deterministic_method_required")
    elif method == "union_bound" and all(math.isfinite(item) for item in expected_probabilities):
        union_bound = _upward_sum(expected_probabilities)
        if union_bound >= 1 or joint < union_bound:
            reasons.append("probability_accounting:invalid_union_bound")


def _validate_decision_semantics(table: Mapping[str, Any], reasons: list[str]) -> None:
    expected_rules = {
        "missing_required_bound": "NOT_CERTIFIABLE",
        "bound_uses_target_fit_residuals": "NOT_CERTIFIABLE",
        "nonfinite_bound": "NOT_CERTIFIABLE",
        "mixed_or_unresolved_units": "NOT_CERTIFIABLE",
        "probabilistic_term_in_deterministic_request": "NOT_CERTIFIABLE",
        "unknown_harmonic_tail": "NOT_CERTIFIABLE",
        "unknown_antialias_or_derivative_bound": "NOT_CERTIFIABLE",
        "unknown_relative_delay": "NOT_CERTIFIABLE",
        "gram_singular_or_ill_conditioned": "NOT_CERTIFIABLE",
        "total_radius_above_threshold": "NOT_CERTIFIABLE",
        "acausal_preprocessing": "EXCLUDE_PROTOCOL",
        "per_window_amplitude_normalization": "EXCLUDE_PROTOCOL",
        "per_window_free_phase_alignment": "EXCLUDE_PROTOCOL",
    }
    for key, expected in expected_rules.items():
        if table.get(key, _MISSING) != expected:
            reasons.append(f"decision_rules.{key}:invalid")
    _require_false(
        table,
        "r2_alone_may_return_reject_model",
        "decision_rules",
        reasons,
    )
    _require_true(
        table,
        "r1_empty_joint_set_required_for_model_rejection",
        "decision_rules",
        reasons,
    )


def _window_numeric(
    table: Mapping[str, Any],
    key: str,
    prefix: str,
    reasons: list[str],
    *,
    positive: bool = False,
) -> float | None:
    value = table.get(key, _MISSING)
    valid = _is_finite_number(value) and (value > 0 if positive else value >= 0)
    if not valid:
        qualifier = "positive_finite" if positive else "finite_nonnegative"
        reasons.append(f"{prefix}.{key}:must_be_{qualifier}")
        return None
    return float(value)


def _window_complex_magnitude(
    table: Mapping[str, Any],
    real_key: str,
    imag_key: str,
    prefix: str,
    reasons: list[str],
) -> float | None:
    real = table.get(real_key, _MISSING)
    imag = table.get(imag_key, _MISSING)
    if not _is_finite_number(real) or not _is_finite_number(imag):
        reasons.append(f"{prefix}.{real_key}_{imag_key}:nonfinite")
        return None
    return math.hypot(real, imag)


def _measurement_snr(
    measurement: Mapping[str, Any],
    transfer_magnitude_lower: float,
    target_lower: float,
    duration_lower: float,
) -> float | None:
    bound_type = measurement.get("bound_type")
    if bound_type == "continuous_l2":
        noise = measurement.get("continuous_l2_bound_signal_sqrt_s")
        scale = _downward(math.sqrt(duration_lower))
    elif bound_type == "weighted_l2":
        noise = measurement.get("weighted_l2_bound")
        scale = 1.0
    elif bound_type == "linfinity":
        noise = measurement.get("linfinity_bound")
        scale = 1.0
    else:
        return None
    if not _is_finite_nonnegative(noise):
        return None
    if noise == 0:
        return math.inf
    signal_floor = _downward_product(
        _downward_product(transfer_magnitude_lower, target_lower),
        scale,
    )
    return _downward_ratio(signal_floor, float(noise))


def _validate_computed_bound_floors(
    bounds: Mapping[str, Any],
    *,
    duration: float | None,
    nominal_rate: float | None,
    target_harmonic: int | None,
    reasons: list[str],
) -> None:
    """Check primitive bounds cannot be hidden behind smaller computed fields."""

    measurement = bounds.get("measurement_noise", {})
    if isinstance(measurement, Mapping):
        computed = _window_numeric(
            measurement,
            "computed_w_norm_bound",
            "bounds.measurement_noise",
            reasons,
        )
        bound_type = measurement.get("bound_type", _MISSING)
        raw_field = {
            "weighted_l2": "weighted_l2_bound",
            "continuous_l2": "continuous_l2_bound_signal_sqrt_s",
            "linfinity": "linfinity_bound",
        }.get(bound_type)
        raw = (
            _window_numeric(
                measurement,
                raw_field,
                "bounds.measurement_noise",
                reasons,
            )
            if raw_field is not None
            else None
        )
        required = None
        if raw is not None:
            if bound_type == "continuous_l2" and duration is not None:
                sqrt_lower = math.nextafter(math.sqrt(duration), -math.inf)
                required = _upward_ratio(_upward(raw), sqrt_lower)
            elif bound_type in {"weighted_l2", "linfinity"}:
                required = _upward(raw)
        if computed is not None and required is not None and computed < required:
            reasons.append("bounds.measurement_noise.computed_w_norm_bound:below_raw_floor")

    model = bounds.get("model_residual", {})
    if isinstance(model, Mapping) and model.get("residual_location") == "output":
        computed = _window_numeric(
            model,
            "computed_w_norm_bound",
            "bounds.model_residual",
            reasons,
        )
        raw = _window_numeric(
            model,
            "output_weighted_l2_bound",
            "bounds.model_residual",
            reasons,
        )
        if computed is not None and raw is not None and computed < _upward(raw):
            reasons.append("bounds.model_residual.computed_w_norm_bound:below_raw_floor")

    sampling = bounds.get("sampling", {})
    if isinstance(sampling, Mapping):
        fields = {
            key: _window_numeric(sampling, key, "bounds.sampling", reasons)
            for key in (
                "anti_alias_residual_bound",
                "adc_quantization_step",
                "timestamp_jitter_s",
                "r_peak_endpoint_interpolation_bound",
                "first_derivative_bound_signal_per_s",
                "computed_w_norm_bound",
            )
        }
        if all(value is not None for value in fields.values()):
            anti_alias = fields["anti_alias_residual_bound"]
            adc_step = fields["adc_quantization_step"]
            jitter = fields["timestamp_jitter_s"]
            endpoint = fields["r_peak_endpoint_interpolation_bound"]
            derivative = fields["first_derivative_bound_signal_per_s"]
            computed = fields["computed_w_norm_bound"]
            assert anti_alias is not None
            assert adc_step is not None
            assert jitter is not None
            assert endpoint is not None
            assert derivative is not None
            assert computed is not None
            required = _upward_sum(
                [
                    anti_alias,
                    _upward_ratio(adc_step, 2.0),
                    _upward_product(derivative, jitter),
                    endpoint,
                ]
            )
            if computed < required:
                reasons.append("bounds.sampling.computed_w_norm_bound:below_raw_floor")

    delay = bounds.get("delay", {})
    if (
        isinstance(delay, Mapping)
        and nominal_rate is not None
        and target_harmonic is not None
        and target_harmonic != 0
    ):
        fields = {
            key: _window_numeric(delay, key, "bounds.delay", reasons)
            for key in (
                "relative_amplitude_error_bound",
                "relative_group_delay_error_s",
                "residual_phase_calibration_error_rad",
                "target_harmonic_abs_upper_bound",
                "computed_coefficient_radius",
            )
        }
        if all(value is not None for value in fields.values()):
            amplitude = fields["relative_amplitude_error_bound"]
            delay_s = fields["relative_group_delay_error_s"]
            phase_zero = fields["residual_phase_calibration_error_rad"]
            target_upper = fields["target_harmonic_abs_upper_bound"]
            computed = fields["computed_coefficient_radius"]
            assert amplitude is not None
            assert delay_s is not None
            assert phase_zero is not None
            assert target_upper is not None
            assert computed is not None
            harmonic_rate = _upward_product(float(abs(target_harmonic)), nominal_rate)
            phase = _upward_sum([_upward_product(harmonic_rate, delay_s), phase_zero])
            if phase >= math.pi:
                chord = 2.0
            else:
                chord = _upward_product(2.0, _upward(math.sin(phase / 2.0)))
            one_plus_amplitude = _upward_sum([1.0, amplitude])
            relative_error = _upward_sum([amplitude, _upward_product(one_plus_amplitude, chord)])
            required = _upward_product(target_upper, relative_error)
            if computed < required:
                reasons.append("bounds.delay.computed_coefficient_radius:below_raw_floor")


def _validate_window_results(
    config: Mapping[str, Any],
    bounds: Mapping[str, Any],
    estimator: Mapping[str, Any],
    window_design: Mapping[str, Any],
    preprocessing: Mapping[str, Any],
    thresholds: Mapping[str, Any],
    probability: Mapping[str, Any],
    requested_semantics: Any,
    harmonics: list[int] | None,
    reasons: list[str],
) -> None:
    results = config.get("window_result", _MISSING)
    if not isinstance(results, list) or not results:
        reasons.append("window_result:missing")
        return

    expected_status = (
        R2ProtocolStatus.PASS_DETERMINISTIC.value
        if requested_semantics == "deterministic"
        else R2ProtocolStatus.PASS_PROBABILISTIC.value
    )
    min_eigenvalue = estimator.get("gram_min_eigenvalue_threshold", _MISSING)
    max_condition = estimator.get("gram_condition_number_max", _MISSING)
    planned_intervals = window_design.get("complete_rr_intervals", _MISSING)
    target_lower = thresholds.get("target_harmonic_abs_lower_bound", _MISSING)
    total_threshold = thresholds.get("deterministic_total_radius_max", _MISSING)
    relative_threshold = thresholds.get("deterministic_relative_radius_max", _MISSING)
    minimum_burn_in = thresholds.get("minimum_burn_in_s", _MISSING)
    minimum_dwell = thresholds.get("minimum_dwell_s", _MISSING)
    minimum_snr = thresholds.get("minimum_energy_snr", _MISSING)
    expected_probability = probability.get("joint_failure_probability", _MISSING)
    measurement = bounds.get("measurement_noise", {})
    dwell = bounds.get("dwell", {})
    delay = bounds.get("delay", {})
    protocol_sha256 = canonical_protocol_sha256(config)
    joint_identity: tuple[str, str, str] | None = None
    window_harmonics: set[tuple[str, int]] = set()

    semantics_by_bound = {
        name: (bounds[name].get("semantics") if isinstance(bounds.get(name), Mapping) else None)
        for name in BOUND_NAMES
    }

    for index, result in enumerate(results):
        prefix = f"window_result[{index}]"
        if not isinstance(result, Mapping):
            reasons.append(f"{prefix}:not_table")
            continue
        for key in (
            "window_id",
            "rate_id",
            "subject_pseudonym",
            "lead_id",
            "physical_gain_id",
            "record_manifest_id",
        ):
            _require_nonempty(result, key, prefix, reasons)

        window_id = result.get("window_id", _MISSING)

        identity_values = tuple(
            result.get(key, _MISSING)
            for key in ("subject_pseudonym", "lead_id", "physical_gain_id")
        )
        if all(_is_nonempty(value) for value in identity_values):
            current_identity = identity_values
            if joint_identity is None:
                joint_identity = current_identity
            elif current_identity != joint_identity:
                reasons.append(f"{prefix}:joint_identity_mismatch")

        if result.get("protocol_sha256", _MISSING) != protocol_sha256:
            reasons.append(f"{prefix}.protocol_sha256:mismatch")
        stored_result_hash = result.get("result_sha256", _MISSING)
        if (
            not isinstance(stored_result_hash, str)
            or _SHA256_RE.fullmatch(stored_result_hash) is None
        ):
            reasons.append(f"{prefix}.result_sha256:invalid_sha256")
        elif stored_result_hash != canonical_window_result_sha256(result):
            reasons.append(f"{prefix}.result_sha256:mismatch")

        left = _window_numeric(result, "left_r_peak_time_s", prefix, reasons)
        right = _window_numeric(result, "right_r_peak_time_s", prefix, reasons)
        duration = _window_numeric(
            result,
            "window_duration_s",
            prefix,
            reasons,
            positive=True,
        )
        rate = _window_numeric(
            result,
            "nominal_rate_rad_s",
            prefix,
            reasons,
            positive=True,
        )
        intervals = result.get("complete_rr_intervals", _MISSING)
        if not isinstance(intervals, int) or isinstance(intervals, bool) or intervals <= 0:
            reasons.append(f"{prefix}.complete_rr_intervals:invalid")
            intervals = None
        elif isinstance(planned_intervals, int) and intervals != planned_intervals:
            reasons.append(f"{prefix}.complete_rr_intervals:plan_mismatch")

        sample_count = result.get("sample_count", _MISSING)
        minimum_samples = len(harmonics) if harmonics is not None else 1
        if (
            not isinstance(sample_count, int)
            or isinstance(sample_count, bool)
            or sample_count < minimum_samples
        ):
            reasons.append(f"{prefix}.sample_count:insufficient")
        sampling = bounds.get("sampling", {})
        sample_rate = (
            sampling.get("sample_rate_hz", _MISSING) if isinstance(sampling, Mapping) else _MISSING
        )
        if (
            isinstance(sample_count, int)
            and not isinstance(sample_count, bool)
            and duration is not None
            and _is_finite_number(sample_rate)
            and abs(sample_count - sample_rate * duration) > 1.0
        ):
            reasons.append(f"{prefix}.sample_count:rate_duration_mismatch")

        if left is not None and right is not None and right <= left:
            reasons.append(f"{prefix}.r_peak_times:not_increasing")
        if left is not None and right is not None and duration is not None:
            if not math.isclose(right - left, duration, rel_tol=1.0e-9, abs_tol=1.0e-12):
                reasons.append(f"{prefix}.window_duration_s:anchor_mismatch")
        certified_duration_fraction = (
            Fraction.from_float(duration) if duration is not None else None
        )
        if left is not None and right is not None and right > left:
            anchor_duration_fraction = Fraction.from_float(right) - Fraction.from_float(left)
            certified_duration_fraction = (
                anchor_duration_fraction
                if certified_duration_fraction is None
                else min(certified_duration_fraction, anchor_duration_fraction)
            )
        certified_duration_lower = (
            _fraction_lower_float(certified_duration_fraction)
            if certified_duration_fraction is not None
            else None
        )
        certified_rate_upper = None
        if duration is not None and intervals is not None and rate is not None:
            expected_rate = 2.0 * math.pi * intervals / duration
            if not math.isclose(rate, expected_rate, rel_tol=1.0e-9, abs_tol=1.0e-12):
                reasons.append(f"{prefix}.nominal_rate_rad_s:geometry_mismatch")
        if certified_duration_lower is not None and intervals is not None:
            pi_upper = math.nextafter(math.pi, math.inf)
            angular_cycles_upper = _upward_product(
                _upward_product(2.0, pi_upper),
                float(intervals),
            )
            certified_rate_upper = _upward_ratio(
                angular_cycles_upper,
                certified_duration_lower,
            )

        gram_min = _window_numeric(
            result,
            "gram_min_eigenvalue",
            prefix,
            reasons,
            positive=True,
        )
        gram_condition = _window_numeric(
            result,
            "gram_condition_number",
            prefix,
            reasons,
            positive=True,
        )
        wls_row_norm_upper = _window_numeric(
            result,
            "wls_row_norm_upper_bound",
            prefix,
            reasons,
            positive=True,
        )
        _require_true(
            result,
            "gram_bounds_outward_certified",
            prefix,
            reasons,
        )
        _require_sha256(result, "gram_certificate_sha256", prefix, reasons)
        if gram_min is not None and _is_finite_number(min_eigenvalue) and gram_min < min_eigenvalue:
            reasons.append(f"{prefix}.gram_min_eigenvalue:below_threshold")
        if gram_min is not None and wls_row_norm_upper is not None:
            sqrt_lower = math.nextafter(math.sqrt(gram_min), -math.inf)
            row_norm_floor = _upward_ratio(1.0, sqrt_lower)
            if wls_row_norm_upper < row_norm_floor:
                reasons.append(f"{prefix}.wls_row_norm_upper_bound:below_gram_floor")
        if gram_condition is not None:
            if gram_condition < 1:
                reasons.append(f"{prefix}.gram_condition_number:below_one")
            if _is_finite_number(max_condition) and gram_condition > max_condition:
                reasons.append(f"{prefix}.gram_condition_number:above_threshold")

        target_harmonic = result.get("target_harmonic", _MISSING)
        if (
            not isinstance(target_harmonic, int)
            or isinstance(target_harmonic, bool)
            or target_harmonic == 0
            or harmonics is None
            or target_harmonic not in harmonics
        ):
            reasons.append(f"{prefix}.target_harmonic:invalid")
            valid_target_harmonic = None
        else:
            valid_target_harmonic = target_harmonic
            if _is_nonempty(window_id):
                result_key = (window_id, target_harmonic)
                if result_key in window_harmonics:
                    reasons.append(f"{prefix}:duplicate_window_harmonic")
                window_harmonics.add(result_key)

        _validate_computed_bound_floors(
            bounds,
            duration=certified_duration_lower,
            nominal_rate=certified_rate_upper,
            target_harmonic=valid_target_harmonic,
            reasons=reasons,
        )

        transfer = _window_complex_magnitude(
            result,
            "nominal_transfer_real",
            "nominal_transfer_imag",
            prefix,
            reasons,
        )
        transfer_real = result.get("nominal_transfer_real", _MISSING)
        transfer_imag = result.get("nominal_transfer_imag", _MISSING)
        if preprocessing.get("identity_front_end", _MISSING) is True and (
            transfer_real != 1.0 or transfer_imag != 0.0
        ):
            reasons.append(f"{prefix}.nominal_transfer:identity_must_equal_one")
        _window_complex_magnitude(
            result,
            "z_tilde_real",
            "z_tilde_imag",
            prefix,
            reasons,
        )
        if transfer is not None:
            if transfer <= 0:
                reasons.append(f"{prefix}.nominal_transfer:zero")
            transfer_lower = (
                delay.get("nominal_transfer_magnitude_lower_bound", _MISSING)
                if isinstance(delay, Mapping)
                else _MISSING
            )
            if preprocessing.get("identity_front_end", _MISSING) is True and transfer_lower != 1.0:
                reasons.append(
                    "bounds.delay.nominal_transfer_magnitude_lower_bound:identity_must_equal_one"
                )
            if _is_finite_number(transfer_lower) and transfer < transfer_lower:
                reasons.append(f"{prefix}.nominal_transfer:below_calibrated_lower_bound")

        radius_values: dict[str, float] = {}
        for _, field in _RADIUS_FIELDS:
            value = _window_numeric(result, field, prefix, reasons)
            if value is not None:
                radius_values[field] = value
        all_radii_present = len(radius_values) == len(_RADIUS_FIELDS)

        transfer_lower = (
            delay.get("nominal_transfer_magnitude_lower_bound", _MISSING)
            if isinstance(delay, Mapping)
            else _MISSING
        )
        amplification = None
        if wls_row_norm_upper is not None and _is_finite_number(transfer_lower):
            if transfer_lower > 0:
                amplification = _upward_ratio(wls_row_norm_upper, float(transfer_lower))

        if all_radii_present:
            for (
                bound_name,
                result_field,
                source_field,
                needs_amplification,
            ) in _RADIUS_SOURCE_FIELDS:
                bound_table = bounds.get(bound_name, {})
                source_value = (
                    _window_numeric(
                        bound_table,
                        source_field,
                        f"bounds.{bound_name}",
                        reasons,
                    )
                    if isinstance(bound_table, Mapping)
                    else None
                )
                if source_value is None:
                    continue
                if needs_amplification:
                    if amplification is None:
                        continue
                    radius_floor = _upward_product(amplification, source_value)
                else:
                    radius_floor = _upward(source_value)
                if radius_values[result_field] < radius_floor:
                    reasons.append(f"{prefix}.{result_field}:below_recomputed_floor")

        radius_sum = _upward_sum(list(radius_values.values())) if all_radii_present else None
        deterministic_sum = (
            _upward_sum(
                [
                    radius_values[field]
                    for bound_name, field in _RADIUS_FIELDS
                    if semantics_by_bound.get(bound_name) == "deterministic"
                ]
            )
            if all_radii_present
            else None
        )

        total_deterministic = _window_numeric(
            result,
            "radius_total_deterministic",
            prefix,
            reasons,
        )
        total_probability = _window_numeric(
            result,
            "radius_total_probabilistic",
            prefix,
            reasons,
        )
        if deterministic_sum is not None and total_deterministic is not None:
            if total_deterministic < deterministic_sum:
                reasons.append(f"{prefix}.radius_total_deterministic:below_sum")
        if requested_semantics == "deterministic":
            active_total = total_deterministic
            if total_probability is not None and total_probability != 0:
                reasons.append(f"{prefix}.radius_total_probabilistic:must_be_zero")
        else:
            active_total = total_probability
            if radius_sum is not None and total_probability is not None:
                if total_probability < radius_sum:
                    reasons.append(f"{prefix}.radius_total_probabilistic:below_sum")

        if active_total is not None:
            if _is_finite_number(total_threshold) and active_total > total_threshold:
                reasons.append(f"{prefix}:total_radius_above_threshold")
            if _is_finite_number(target_lower) and target_lower > 0:
                relative_radius = _upward_ratio(active_total, float(target_lower))
                if _is_finite_number(relative_threshold) and relative_radius > relative_threshold:
                    reasons.append(f"{prefix}:relative_radius_above_threshold")

        burn_in = dwell.get("burn_in_s", _MISSING) if isinstance(dwell, Mapping) else _MISSING
        certified_duration = (
            dwell.get("window_duration_s", _MISSING) if isinstance(dwell, Mapping) else _MISSING
        )
        if (
            duration is not None
            and _is_finite_number(certified_duration)
            and not math.isclose(
                duration,
                certified_duration,
                rel_tol=1.0e-9,
                abs_tol=1.0e-12,
            )
        ):
            reasons.append(f"{prefix}.window_duration_s:dwell_bound_mismatch")
        burn_ok = (
            _is_finite_nonnegative(burn_in)
            and _is_finite_nonnegative(minimum_burn_in)
            and burn_in >= minimum_burn_in
        )
        if (
            certified_duration_fraction is not None
            and _is_finite_nonnegative(minimum_dwell)
            and _is_finite_nonnegative(burn_in)
        ):
            dwell_exact = Fraction.from_float(float(burn_in)) + certified_duration_fraction
            burn_ok = burn_ok and dwell_exact >= Fraction.from_float(float(minimum_dwell))
        if not burn_ok:
            reasons.append(f"{prefix}:minimum_burn_in_or_dwell_not_met")
        if result.get("minimum_burn_in_pass", _MISSING) is not burn_ok:
            reasons.append(f"{prefix}.minimum_burn_in_pass:inconsistent")

        snr_ok = False
        snr_transfer_lower = (
            delay.get("nominal_transfer_magnitude_lower_bound", _MISSING)
            if isinstance(delay, Mapping)
            else _MISSING
        )
        if (
            isinstance(measurement, Mapping)
            and _is_finite_number(snr_transfer_lower)
            and snr_transfer_lower > 0
            and certified_duration_lower is not None
            and _is_finite_number(target_lower)
            and target_lower > 0
            and _is_finite_number(minimum_snr)
        ):
            computed_snr = _measurement_snr(
                measurement,
                float(snr_transfer_lower),
                float(target_lower),
                certified_duration_lower,
            )
            snr_ok = computed_snr is not None and computed_snr >= minimum_snr
        if not snr_ok:
            reasons.append(f"{prefix}:minimum_snr_not_met")
        if result.get("minimum_snr_pass", _MISSING) is not snr_ok:
            reasons.append(f"{prefix}.minimum_snr_pass:inconsistent")

        result_probability = _window_numeric(
            result,
            "joint_failure_probability",
            prefix,
            reasons,
        )
        if (
            result_probability is not None
            and _is_finite_number(expected_probability)
            and result_probability != expected_probability
        ):
            reasons.append(f"{prefix}.joint_failure_probability:mismatch")

        if result.get("result_status", _MISSING) != expected_status:
            reasons.append(f"{prefix}.result_status:inconsistent")
        reason_codes = result.get("reason_codes", _MISSING)
        if not isinstance(reason_codes, list) or reason_codes:
            reasons.append(f"{prefix}.reason_codes:pass_requires_empty_list")


def validate_r2_protocol(config: Mapping[str, Any]) -> R2ProtocolValidation:
    """Validate one parsed protocol with exclusion taking strict precedence."""

    if not isinstance(config, Mapping):
        empty_hash = canonical_protocol_sha256({})
        return R2ProtocolValidation(
            R2ProtocolStatus.NOT_CERTIFIABLE,
            ("root:not_mapping",),
            empty_hash,
        )

    canonical_sha256 = canonical_protocol_sha256(config)
    exclusions = _validate_exclusions(config)
    if exclusions:
        return R2ProtocolValidation(
            R2ProtocolStatus.EXCLUDE_PROTOCOL,
            tuple(sorted(set(exclusions))),
            canonical_sha256,
        )

    reasons: list[str] = []
    protocol = _table(config, "protocol", reasons)
    data_boundary = _table(config, "data_boundary", reasons)
    units = _table(config, "units", reasons)
    model = _table(config, "model", reasons)
    estimator = _table(config, "estimator", reasons)
    window = _table(config, "window_design", reasons)
    preprocessing = _table(config, "preprocessing", reasons)
    evidence = _table(config, "evidence_policy", reasons)
    bounds = _table(config, "bounds", reasons)
    thresholds = _table(config, "decision_thresholds", reasons)
    decision_rules = _table(config, "decision_rules", reasons)
    probability = _table(config, "probability_accounting", reasons)
    audit_outputs = _table(config, "audit_outputs", reasons)

    _validate_freeze_and_hash(config, protocol, canonical_sha256, reasons)
    _validate_data_boundary(data_boundary, reasons)
    _validate_units(units, reasons)
    harmonics = _validate_estimator(estimator, reasons)
    _validate_model(model, harmonics, units.get("signal_unit", _MISSING), reasons)
    _validate_window(window, reasons)
    _validate_preprocessing(preprocessing, reasons)
    _validate_evidence_policy(evidence, reasons)
    _validate_audit_outputs(audit_outputs, reasons)
    requested_semantics = protocol.get("requested_output_semantics", _MISSING)
    expected_probabilities = _validate_bounds(
        bounds,
        units.get("signal_unit", _MISSING),
        requested_semantics,
        data_boundary.get("calibration_split_id", _MISSING),
        data_boundary.get("identification_split_id", _MISSING),
        data_boundary.get("validation_split_id", _MISSING),
        reasons,
    )
    _validate_bound_specific_gates(bounds, reasons)
    _validate_thresholds(thresholds, requested_semantics, window, reasons)
    _validate_probability_accounting(
        probability,
        expected_probabilities,
        requested_semantics,
        thresholds.get("maximum_joint_failure_probability", _MISSING),
        evidence.get("joint_probability_method", _MISSING),
        data_boundary.get("calibration_split_id", _MISSING),
        reasons,
    )
    _validate_decision_semantics(decision_rules, reasons)
    _validate_window_results(
        config,
        bounds,
        estimator,
        window,
        preprocessing,
        thresholds,
        probability,
        requested_semantics,
        harmonics,
        reasons,
    )

    if reasons:
        return R2ProtocolValidation(
            R2ProtocolStatus.NOT_CERTIFIABLE,
            tuple(sorted(set(reasons))),
            canonical_sha256,
        )
    status = (
        R2ProtocolStatus.PASS_DETERMINISTIC
        if requested_semantics == "deterministic"
        else R2ProtocolStatus.PASS_PROBABILISTIC
    )
    return R2ProtocolValidation(status, (), canonical_sha256)


def validate_r2_protocol_file(path: str | Path) -> R2ProtocolValidation:
    """Load and validate one TOML file without turning parse errors into crashes."""

    try:
        config = load_r2_protocol(path)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        return R2ProtocolValidation(
            R2ProtocolStatus.NOT_CERTIFIABLE,
            (f"toml:{type(exc).__name__}",),
            canonical_protocol_sha256({}),
        )
    return validate_r2_protocol(config)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("protocol", type=Path)
    args = parser.parse_args(argv)
    result = validate_r2_protocol_file(args.protocol)
    print(json.dumps(result.as_dict(), ensure_ascii=False, sort_keys=True))
    return 0 if result.passed else 2


if __name__ == "__main__":
    sys.exit(main())
