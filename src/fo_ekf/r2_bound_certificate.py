"""Replayable deterministic R2 component-bound certificates.

This module closes only the numerical propagation for six declared R2
components: the WLS Gram factor, finite-record history, previous-rate dwell,
sampling/interpolation, optional trapezoidal quadrature, and relative front-end
delay.  Primitive engineering envelopes remain external evidence.  A
certificate from this module is therefore *not* an R2 data-to-disk PASS.

Every issued bound is recomputed with Arb ball arithmetic from an exact
binary-rational input manifest.  Missing provenance, an uncertified Gram
matrix, unsupported inputs, resource exhaustion, and unresolved ball signs
all produce ``UNKNOWN``.  No outcome in this module can reject a model.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from typing import Any

import flint
from flint import acb, arb, ctx

from .certified_set import (
    _CERTIFICATION_LOCK,
    DEFAULT_PRECISION_SCHEDULE,
    _validate_precision_schedule,
)

SCHEMA = "fo-ekf.r2-bound-engine.v2"
ALGORITHM = "arb-mode-separated-r2-bounds-v2"
BOUND_NUMERIC_SCOPE = "uniform_all_windows_satisfying_frozen_design"
MAX_SAMPLES = 4096
MAX_HARMONICS = 64
MAX_CERTIFICATE_BYTES = 4_000_000
MAX_REPLAY_PRECISION = 4096

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_GROUP_SOURCE_KINDS = {
    "history": frozenset(
        {
            "independent_bounded_calibration",
            "proved_engineering_envelope",
        }
    ),
    "dwell": frozenset(
        {
            "independent_bounded_calibration",
            "proved_engineering_envelope",
        }
    ),
    "sampling": frozenset(
        {
            "sensor_spec",
            "independent_bounded_calibration",
            "proved_engineering_envelope",
        }
    ),
    "delay": frozenset(
        {
            "sensor_spec",
            "independent_bounded_calibration",
            "proved_engineering_envelope",
        }
    ),
}
_ALL_SOURCE_KINDS = frozenset(
    {
        "sensor_spec",
        "independent_bounded_calibration",
        "proved_engineering_envelope",
        "observable_window_geometry",
    }
)


class R2BoundRelation(str, Enum):
    """The only claims expressible by the component-bound engine."""

    CERTIFIED_BOUND = "CERTIFIED_BOUND"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class EvidenceReference:
    """Frozen provenance for one group of primitive deterministic bounds."""

    source_kind: str
    source_reference: str
    evidence_sha256: str
    calibration_split_id: str
    applicability_domain: str
    signal_unit: str
    time_unit: str
    independent_of_identification_fit: bool
    uses_target_fit_residuals: bool = False
    numeric_scope: str = BOUND_NUMERIC_SCOPE


@dataclass(frozen=True)
class HistoryBoundInput:
    """Primitive bounds for finite-lower-limit and pre-record history."""

    record_origin_s: float
    rate_segment_start_s: float
    initial_state_abs_bound: float
    new_input_sup_bound: float
    pre_record_history_coefficient: float
    evidence: EvidenceReference
    preconvolved_output_envelope: tuple[float, ...] | None = None


@dataclass(frozen=True)
class DwellBoundInput:
    """Primitive input envelopes for previous-rate memory."""

    previous_input_sup_bound: float
    new_input_sup_bound: float
    evidence: EvidenceReference
    preconvolved_output_envelope: tuple[float, ...] | None = None


@dataclass(frozen=True)
class SamplingBoundInput:
    """Per-sample acquisition bounds and an optional quadrature contract."""

    timestamp_error_s: tuple[float, ...]
    interpolation_error_bound: tuple[float, ...]
    anti_alias_error_bound: tuple[float, ...]
    adc_quantization_step: float
    first_derivative_bound_signal_per_s: float
    implementation: str
    panel_second_derivative_bound: tuple[float, ...]
    evidence: EvidenceReference


@dataclass(frozen=True)
class DelayBoundInput:
    """Primitive relative front-end amplitude, timing, and phase bounds."""

    relative_amplitude_error_bound: float
    relative_group_delay_error_s: float
    residual_phase_calibration_error_rad: float
    target_harmonic_abs_upper_bound: float
    nominal_transfer_magnitude_lower_bound: float
    ecg_to_rpeak_channel_alignment_included: bool
    evidence: EvidenceReference


@dataclass(frozen=True)
class R2BoundRequest:
    """All frozen inputs required for one window/harmonic component proof."""

    window_id: str
    protocol_sha256: str
    record_manifest_sha256: str
    estimator_manifest_sha256: str
    calibration_split_id: str
    identification_split_id: str
    frozen_before_target_ecg_access: bool
    signal_unit: str
    dynamic_front_end_mode: str
    tau_star_s: float
    alpha_interval: tuple[float, float]
    lambda_interval: tuple[float, float]
    left_r_peak_time_s: float
    right_r_peak_time_s: float
    complete_rr_intervals: int
    sample_times_s: tuple[float, ...]
    raw_weights: tuple[float, ...]
    retained_harmonics: tuple[int, ...]
    target_harmonic: int
    gram_min_eigenvalue_threshold: float
    gram_condition_number_max: float
    history: HistoryBoundInput | None
    dwell: DwellBoundInput | None
    sampling: SamplingBoundInput | None
    delay: DelayBoundInput | None


@dataclass(frozen=True)
class R2BoundCertificate:
    """A sealed component certificate or a sealed fail-closed UNKNOWN."""

    relation: R2BoundRelation
    reason: str
    precision_bits: int
    attempted_precisions: tuple[int, ...]
    input_sha256: str
    certificate_sha256: str
    certificate_json: str

    def to_json(self) -> str:
        return self.certificate_json


@dataclass(frozen=True)
class R2BoundReplayResult:
    """Independent replay verdict; malformed material never certifies."""

    valid: bool
    relation: R2BoundRelation
    reason: str


def _fraction(value: float) -> Fraction:
    return Fraction.from_float(float(value))


def _rat(value: Fraction) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _arb_fraction(value: Fraction) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _environment_manifest() -> dict[str, str]:
    return {
        "python_flint": flint.__version__,
        "flint": flint.__FLINT_VERSION__,
    }


def _reject_duplicate_object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate JSON object key")
        value[key] = item
    return value


def _evidence_manifest(value: EvidenceReference) -> dict[str, Any]:
    return {
        "source_kind": value.source_kind,
        "source_reference": value.source_reference,
        "evidence_sha256": value.evidence_sha256,
        "calibration_split_id": value.calibration_split_id,
        "applicability_domain": value.applicability_domain,
        "signal_unit": value.signal_unit,
        "time_unit": value.time_unit,
        "independent_of_identification_fit": value.independent_of_identification_fit,
        "uses_target_fit_residuals": value.uses_target_fit_residuals,
        "numeric_scope": value.numeric_scope,
    }


def _float_manifest(value: float) -> str:
    return _rat(_fraction(value))


def _request_manifest(request: R2BoundRequest) -> dict[str, Any]:
    history = request.history
    dwell = request.dwell
    sampling = request.sampling
    delay = request.delay
    return {
        "model": "r2-finite-record-component-bounds-v2",
        "window_id": request.window_id,
        "protocol_sha256": request.protocol_sha256,
        "record_manifest_sha256": request.record_manifest_sha256,
        "estimator_manifest_sha256": request.estimator_manifest_sha256,
        "calibration_split_id": request.calibration_split_id,
        "identification_split_id": request.identification_split_id,
        "frozen_before_target_ecg_access": request.frozen_before_target_ecg_access,
        "signal_unit": request.signal_unit,
        "dynamic_front_end_mode": request.dynamic_front_end_mode,
        "tau_star_s": _float_manifest(request.tau_star_s),
        "alpha_interval": [_float_manifest(value) for value in request.alpha_interval],
        "lambda_interval": [_float_manifest(value) for value in request.lambda_interval],
        "geometry": {
            "left_r_peak_time_s": _float_manifest(request.left_r_peak_time_s),
            "right_r_peak_time_s": _float_manifest(request.right_r_peak_time_s),
            "complete_rr_intervals": request.complete_rr_intervals,
            "sample_times_s": [_float_manifest(value) for value in request.sample_times_s],
            "raw_weights": [_float_manifest(value) for value in request.raw_weights],
            "retained_harmonics": list(request.retained_harmonics),
            "target_harmonic": request.target_harmonic,
        },
        "gram_thresholds": {
            "minimum_eigenvalue": _float_manifest(request.gram_min_eigenvalue_threshold),
            "maximum_condition_number": _float_manifest(request.gram_condition_number_max),
        },
        "history": None
        if history is None
        else {
            "record_origin_s": _float_manifest(history.record_origin_s),
            "rate_segment_start_s": _float_manifest(history.rate_segment_start_s),
            "initial_state_abs_bound": _float_manifest(history.initial_state_abs_bound),
            "new_input_sup_bound": _float_manifest(history.new_input_sup_bound),
            "pre_record_history_coefficient": _float_manifest(
                history.pre_record_history_coefficient
            ),
            "preconvolved_output_envelope": None
            if history.preconvolved_output_envelope is None
            else [_float_manifest(value) for value in history.preconvolved_output_envelope],
            "evidence": _evidence_manifest(history.evidence),
        },
        "dwell": None
        if dwell is None
        else {
            "previous_input_sup_bound": _float_manifest(dwell.previous_input_sup_bound),
            "new_input_sup_bound": _float_manifest(dwell.new_input_sup_bound),
            "preconvolved_output_envelope": None
            if dwell.preconvolved_output_envelope is None
            else [_float_manifest(value) for value in dwell.preconvolved_output_envelope],
            "evidence": _evidence_manifest(dwell.evidence),
        },
        "sampling": None
        if sampling is None
        else {
            "timestamp_error_s": [_float_manifest(value) for value in sampling.timestamp_error_s],
            "interpolation_error_bound": [
                _float_manifest(value) for value in sampling.interpolation_error_bound
            ],
            "anti_alias_error_bound": [
                _float_manifest(value) for value in sampling.anti_alias_error_bound
            ],
            "adc_quantization_step": _float_manifest(sampling.adc_quantization_step),
            "first_derivative_bound_signal_per_s": _float_manifest(
                sampling.first_derivative_bound_signal_per_s
            ),
            "implementation": sampling.implementation,
            "panel_second_derivative_bound": [
                _float_manifest(value) for value in sampling.panel_second_derivative_bound
            ],
            "evidence": _evidence_manifest(sampling.evidence),
        },
        "delay": None
        if delay is None
        else {
            "relative_amplitude_error_bound": _float_manifest(delay.relative_amplitude_error_bound),
            "relative_group_delay_error_s": _float_manifest(delay.relative_group_delay_error_s),
            "residual_phase_calibration_error_rad": _float_manifest(
                delay.residual_phase_calibration_error_rad
            ),
            "target_harmonic_abs_upper_bound": _float_manifest(
                delay.target_harmonic_abs_upper_bound
            ),
            "nominal_transfer_magnitude_lower_bound": _float_manifest(
                delay.nominal_transfer_magnitude_lower_bound
            ),
            "ecg_to_rpeak_channel_alignment_included": (
                delay.ecg_to_rpeak_channel_alignment_included
            ),
            "evidence": _evidence_manifest(delay.evidence),
        },
    }


def _safe_manifest(request: Any) -> dict[str, Any]:
    try:
        return _request_manifest(request)
    except (AttributeError, OverflowError, TypeError, ValueError):
        # UNKNOWN documents need a stable input binding even for malformed
        # requests.  The repr is inert JSON text and carries no proof authority.
        return {
            "model": "r2-finite-record-component-bounds-v2",
            "malformed_request_repr": repr(request),
        }


def _finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _nonnegative(value: Any) -> bool:
    return _finite_number(value) and value >= 0


def _positive(value: Any) -> bool:
    return _finite_number(value) and value > 0


def _valid_sha(value: Any) -> bool:
    return isinstance(value, str) and _SHA256_RE.fullmatch(value) is not None


def _normalized_weight_fractions(values: tuple[float, ...]) -> tuple[Fraction, ...]:
    raw = tuple(_fraction(value) for value in values)
    total = sum(raw, Fraction(0))
    return tuple(value / total for value in raw)


def _composite_trapezoid_weights(
    sample_times_s: tuple[float, ...],
) -> tuple[Fraction, ...]:
    times = tuple(_fraction(value) for value in sample_times_s)
    panels = tuple(right - left for left, right in zip(times[:-1], times[1:], strict=True))
    duration = times[-1] - times[0]
    numerators = (
        (panels[0],)
        + tuple(panels[index - 1] + panels[index] for index in range(1, len(times) - 1))
        + (panels[-1],)
    )
    return tuple(value / (2 * duration) for value in numerators)


def _validate_evidence(
    evidence: Any,
    request: R2BoundRequest,
    name: str,
) -> str | None:
    if not isinstance(evidence, EvidenceReference):
        return f"{name}_evidence_missing"
    if evidence.source_kind not in _ALL_SOURCE_KINDS:
        return f"{name}_source_kind_invalid"
    if evidence.source_kind not in _GROUP_SOURCE_KINDS[name]:
        return f"{name}_source_kind_not_allowed"
    for value, suffix in (
        (evidence.source_reference, "source_reference_missing"),
        (evidence.applicability_domain, "applicability_domain_missing"),
    ):
        if not isinstance(value, str) or not value.strip() or value.strip().upper() == "UNSET":
            return f"{name}_{suffix}"
    if evidence.signal_unit != request.signal_unit:
        return f"{name}_signal_unit_mismatch"
    if evidence.time_unit != "s":
        return f"{name}_time_unit_must_be_s"
    if not _valid_sha(evidence.evidence_sha256):
        return f"{name}_evidence_sha256_invalid"
    if evidence.calibration_split_id not in {
        request.calibration_split_id,
        "not_applicable",
    }:
        return f"{name}_calibration_split_mismatch"
    if evidence.calibration_split_id == request.identification_split_id:
        return f"{name}_calibration_identification_leakage"
    if evidence.independent_of_identification_fit is not True:
        return f"{name}_not_independent_of_identification_fit"
    if evidence.uses_target_fit_residuals is not False:
        return f"{name}_uses_target_fit_residuals"
    if evidence.numeric_scope != BOUND_NUMERIC_SCOPE:
        return f"{name}_numeric_scope_invalid"
    return None


def _validate_request(request: Any) -> str | None:
    if not isinstance(request, R2BoundRequest):
        return "request_type_invalid"
    if not isinstance(request.window_id, str) or not request.window_id.strip():
        return "window_id_missing"
    for name in ("protocol_sha256", "record_manifest_sha256", "estimator_manifest_sha256"):
        if not _valid_sha(getattr(request, name, None)):
            return f"{name}_invalid"
    if (
        not isinstance(request.calibration_split_id, str)
        or not request.calibration_split_id.strip()
        or not isinstance(request.identification_split_id, str)
        or not request.identification_split_id.strip()
        or request.calibration_split_id == request.identification_split_id
    ):
        return "data_split_contract_invalid"
    if request.frozen_before_target_ecg_access is not True:
        return "protocol_not_frozen_before_target_access"
    if not isinstance(request.signal_unit, str) or not request.signal_unit.strip():
        return "signal_unit_missing"
    if request.dynamic_front_end_mode not in {
        "identity",
        "preconvolved_output_envelopes",
    }:
        return "dynamic_front_end_mode_invalid"
    if not _positive(request.tau_star_s):
        return "tau_star_invalid"
    if (
        not isinstance(request.alpha_interval, tuple)
        or len(request.alpha_interval) != 2
        or not all(_finite_number(value) for value in request.alpha_interval)
        or not 0 < request.alpha_interval[0] <= request.alpha_interval[1] <= 1
    ):
        return "alpha_interval_invalid"
    if (
        not isinstance(request.lambda_interval, tuple)
        or len(request.lambda_interval) != 2
        or not all(_finite_number(value) for value in request.lambda_interval)
        or not 0 < request.lambda_interval[0] <= request.lambda_interval[1]
    ):
        return "lambda_interval_invalid"
    if not _finite_number(request.left_r_peak_time_s) or not _finite_number(
        request.right_r_peak_time_s
    ):
        return "window_anchor_nonfinite"
    if request.right_r_peak_time_s <= request.left_r_peak_time_s:
        return "window_anchor_order_invalid"
    if (
        not isinstance(request.complete_rr_intervals, int)
        or isinstance(request.complete_rr_intervals, bool)
        or request.complete_rr_intervals <= 0
    ):
        return "complete_rr_intervals_invalid"
    sample_times = request.sample_times_s
    weights = request.raw_weights
    if (
        not isinstance(sample_times, tuple)
        or not 2 <= len(sample_times) <= MAX_SAMPLES
        or not all(_finite_number(value) for value in sample_times)
        or any(right <= left for left, right in zip(sample_times, sample_times[1:], strict=False))
        or sample_times[0] < request.left_r_peak_time_s
        or sample_times[-1] > request.right_r_peak_time_s
    ):
        return "sample_geometry_invalid"
    if (
        not isinstance(weights, tuple)
        or len(weights) != len(sample_times)
        or not all(_positive(value) for value in weights)
    ):
        return "raw_weights_invalid"
    harmonics = request.retained_harmonics
    if (
        not isinstance(harmonics, tuple)
        or not 1 <= len(harmonics) <= MAX_HARMONICS
        or len(harmonics) > len(sample_times)
        or any(not isinstance(value, int) or isinstance(value, bool) for value in harmonics)
        or len(set(harmonics)) != len(harmonics)
    ):
        return "retained_harmonics_invalid"
    if (
        not isinstance(request.target_harmonic, int)
        or isinstance(request.target_harmonic, bool)
        or request.target_harmonic == 0
        or request.target_harmonic not in harmonics
    ):
        return "target_harmonic_invalid"
    if not _positive(request.gram_min_eigenvalue_threshold):
        return "gram_min_threshold_invalid"
    if (
        not _finite_number(request.gram_condition_number_max)
        or request.gram_condition_number_max < 1
    ):
        return "gram_condition_threshold_invalid"

    history = request.history
    if history is None:
        return "history_bound_missing"
    if not isinstance(history, HistoryBoundInput):
        return "history_bound_invalid"
    if not all(
        _nonnegative(value)
        for value in (
            history.initial_state_abs_bound,
            history.new_input_sup_bound,
            history.pre_record_history_coefficient,
        )
    ):
        return "history_magnitude_invalid"
    if (
        not _finite_number(history.record_origin_s)
        or not _finite_number(history.rate_segment_start_s)
        or history.rate_segment_start_s < history.record_origin_s
        or history.rate_segment_start_s > request.left_r_peak_time_s
    ):
        return "history_time_geometry_invalid"
    reason = _validate_evidence(history.evidence, request, "history")
    if reason is not None:
        return reason

    dwell = request.dwell
    if dwell is None:
        return "dwell_bound_missing"
    if not isinstance(dwell, DwellBoundInput):
        return "dwell_bound_invalid"
    if not all(
        _nonnegative(value) for value in (dwell.previous_input_sup_bound, dwell.new_input_sup_bound)
    ):
        return "dwell_magnitude_invalid"
    reason = _validate_evidence(dwell.evidence, request, "dwell")
    if reason is not None:
        return reason

    sampling = request.sampling
    if sampling is None:
        return "sampling_bound_missing"
    if not isinstance(sampling, SamplingBoundInput):
        return "sampling_bound_invalid"
    for values, name in (
        (sampling.timestamp_error_s, "timestamp_error"),
        (sampling.interpolation_error_bound, "interpolation_error"),
        (sampling.anti_alias_error_bound, "anti_alias_error"),
    ):
        if (
            not isinstance(values, tuple)
            or len(values) != len(sample_times)
            or not all(_nonnegative(value) for value in values)
        ):
            return f"sampling_{name}_invalid"
    if not _nonnegative(sampling.adc_quantization_step) or not _nonnegative(
        sampling.first_derivative_bound_signal_per_s
    ):
        return "sampling_primitive_invalid"
    if sampling.implementation not in {
        "discrete_wls",
        "trapezoidal_continuous_lockin",
    }:
        return "sampling_implementation_unsupported"
    if sampling.implementation == "discrete_wls":
        if sampling.panel_second_derivative_bound:
            return "discrete_wls_must_not_declare_quadrature_panels"
    elif (
        sample_times[0] != request.left_r_peak_time_s
        or sample_times[-1] != request.right_r_peak_time_s
        or len(sampling.panel_second_derivative_bound) != len(sample_times) - 1
        or not all(_nonnegative(value) for value in sampling.panel_second_derivative_bound)
    ):
        return "quadrature_panel_contract_invalid"
    if sampling.implementation == "trapezoidal_continuous_lockin" and _normalized_weight_fractions(
        weights
    ) != _composite_trapezoid_weights(sample_times):
        return "trapezoidal_weights_must_match_composite_rule"
    reason = _validate_evidence(sampling.evidence, request, "sampling")
    if reason is not None:
        return reason

    delay = request.delay
    if delay is None:
        return "delay_bound_missing"
    if not isinstance(delay, DelayBoundInput):
        return "delay_bound_invalid"
    if not all(
        _nonnegative(value)
        for value in (
            delay.relative_amplitude_error_bound,
            delay.relative_group_delay_error_s,
            delay.residual_phase_calibration_error_rad,
            delay.target_harmonic_abs_upper_bound,
        )
    ) or not _positive(delay.nominal_transfer_magnitude_lower_bound):
        return "delay_primitive_invalid"
    if delay.ecg_to_rpeak_channel_alignment_included is not True:
        return "delay_channel_alignment_missing"
    reason = _validate_evidence(delay.evidence, request, "delay")
    if reason is not None:
        return reason
    output_envelopes = (
        history.preconvolved_output_envelope,
        dwell.preconvolved_output_envelope,
    )
    if request.dynamic_front_end_mode == "identity":
        if any(value is not None for value in output_envelopes):
            return "identity_mode_forbids_preconvolved_output_envelopes"
        if (
            delay.relative_amplitude_error_bound != 0
            or delay.relative_group_delay_error_s != 0
            or delay.residual_phase_calibration_error_rad != 0
        ):
            return "identity_mode_requires_zero_front_end_error"
        if delay.nominal_transfer_magnitude_lower_bound != 1:
            return "identity_mode_requires_unit_transfer"
    else:
        for name, envelope in zip(
            ("history", "dwell"),
            output_envelopes,
            strict=True,
        ):
            if (
                not isinstance(envelope, tuple)
                or len(envelope) != len(sample_times)
                or not all(_nonnegative(value) for value in envelope)
            ):
                return f"{name}_preconvolved_output_envelope_invalid"
    return None


def _outer_dyadic(value: arb) -> tuple[Fraction, Fraction] | None:
    """Return finite binary-rational endpoints containing one Arb ball."""

    if not value.is_finite():
        return None
    try:
        lower_float = math.nextafter(float(value.lower()), -math.inf)
        upper_float = math.nextafter(float(value.upper()), math.inf)
    except (OverflowError, ValueError):
        return None
    if not (math.isfinite(lower_float) and math.isfinite(upper_float)):
        return None
    for _ in range(64):
        lower = Fraction.from_float(lower_float)
        if value >= _arb_fraction(lower):
            break
        lower_float = math.nextafter(lower_float, -math.inf)
    else:
        return None
    for _ in range(64):
        upper = Fraction.from_float(upper_float)
        if value <= _arb_fraction(upper):
            break
        upper_float = math.nextafter(upper_float, math.inf)
    else:
        return None
    return lower, upper


def _upper_claim(value: arb) -> Fraction | None:
    enclosure = _outer_dyadic(value)
    return None if enclosure is None else enclosure[1]


def _weighted_norm_upper(
    weights: tuple[Fraction, ...],
    bounds: tuple[Fraction, ...],
) -> Fraction | None:
    squared = sum(
        (weight * bound * bound for weight, bound in zip(weights, bounds, strict=True)),
        Fraction(0),
    )
    return _upper_claim(_arb_fraction(squared).sqrt())


def _relaxation_upper(time: Fraction, alpha: Fraction, damping: Fraction) -> Fraction | None:
    """Uniform upper bound for E_alpha(-damping*time**alpha).

    Simon's optimal hyperbolic inequality and Gamma(1+alpha) <= 1 on
    alpha in [0,1] give E_alpha(-x) <= 1/(1+x).  The caller chooses the
    alpha endpoint that minimizes ``time**alpha`` over the frozen interval.
    """

    if time == 0:
        return Fraction(1)
    powered = _arb_fraction(time) ** _arb_fraction(alpha)
    enclosure = _outer_dyadic(powered)
    if enclosure is None or enclosure[0] <= 0:
        return None
    denominator_lower = Fraction(1) + damping * enclosure[0]
    return Fraction(1, 1) / denominator_lower


def _discrete_wls_gram_proof(
    request: R2BoundRequest,
    weights: tuple[Fraction, ...],
) -> tuple[dict[str, Any] | None, Fraction | None, str | None]:
    times = tuple(_fraction(value) for value in request.sample_times_s)
    left = _fraction(request.left_r_peak_time_s)
    duration = _fraction(request.right_r_peak_time_s) - left
    turns = Fraction(2 * request.complete_rr_intervals, 1) / duration
    harmonics = request.retained_harmonics
    phases = tuple(arb.pi() * _arb_fraction(turns * (time - left)) for time in times)

    gram: list[list[acb]] = []
    for row_index, row_harmonic in enumerate(harmonics):
        row: list[acb] = []
        for column_index, column_harmonic in enumerate(harmonics):
            if row_index == column_index:
                row.append(acb(1))
                continue
            delta = column_harmonic - row_harmonic
            value = acb(0)
            for weight, phase in zip(weights, phases, strict=True):
                value += _arb_fraction(weight) * acb(0, delta * phase).exp()
            row.append(value)
        gram.append(row)

    rows: list[dict[str, Any]] = []
    lower_claims: list[Fraction] = []
    upper_claims: list[Fraction] = []
    for index, row in enumerate(gram):
        radius = arb(0)
        for column, value in enumerate(row):
            if column != index:
                radius += value.abs_upper()
        radius_upper = _upper_claim(radius)
        if radius_upper is None:
            return None, None, "gram_radius_nonfinite"
        lower = Fraction(1) - radius_upper
        upper = Fraction(1) + radius_upper
        lower_claims.append(lower)
        upper_claims.append(upper)
        rows.append(
            {
                "harmonic": harmonics[index],
                "offdiagonal_radius_upper": _rat(radius_upper),
                "eigenvalue_interval": [_rat(lower), _rat(upper)],
            }
        )

    lambda_lower = min(lower_claims)
    lambda_upper = max(upper_claims)
    if lambda_lower <= 0:
        return None, None, "gram_not_positive_by_gershgorin"
    threshold = _fraction(request.gram_min_eigenvalue_threshold)
    if lambda_lower < threshold:
        return None, None, "gram_min_below_threshold"
    condition = lambda_upper / lambda_lower
    if condition > _fraction(request.gram_condition_number_max):
        return None, None, "gram_condition_above_threshold"
    row_norm = _upper_claim(arb(1) / _arb_fraction(lambda_lower).sqrt())
    if row_norm is None:
        return None, None, "gram_row_norm_nonfinite"
    transfer_lower = _fraction(request.delay.nominal_transfer_magnitude_lower_bound)  # type: ignore[union-attr]
    amplification = row_norm / transfer_lower
    return (
        {
            "method": "hermitian_gershgorin_v1",
            "rows": rows,
            "lambda_min_lower": _rat(lambda_lower),
            "lambda_max_upper": _rat(lambda_upper),
            "condition_number_upper": _rat(condition),
            "wls_row_norm_upper": _rat(row_norm),
            "nominal_transfer_magnitude_lower": _rat(transfer_lower),
            "coefficient_amplification_upper": _rat(amplification),
        },
        amplification,
        None,
    )


def _continuous_lockin_gram_proof(
    request: R2BoundRequest,
) -> tuple[dict[str, Any] | None, Fraction | None, str | None]:
    """Use continuous full-cycle Fourier orthogonality, not a sampled Gram."""

    unit = Fraction(1)
    if _fraction(request.gram_min_eigenvalue_threshold) > unit:
        return None, None, "gram_min_below_threshold"
    if _fraction(request.gram_condition_number_max) < unit:
        return None, None, "gram_condition_above_threshold"
    delay = request.delay
    assert delay is not None
    transfer_lower = _fraction(delay.nominal_transfer_magnitude_lower_bound)
    amplification = unit / transfer_lower
    return (
        {
            "method": "continuous_full_cycle_orthogonality_v1",
            "normalization": "one_over_window_duration",
            "complete_rr_intervals": request.complete_rr_intervals,
            "harmonics": list(request.retained_harmonics),
            "lambda_min_lower": "1",
            "lambda_max_upper": "1",
            "condition_number_upper": "1",
            "wls_row_norm_upper": "1",
            "nominal_transfer_magnitude_lower": _rat(transfer_lower),
            "coefficient_amplification_upper": _rat(amplification),
        },
        amplification,
        None,
    )


def _gram_proof(
    request: R2BoundRequest,
    weights: tuple[Fraction, ...],
) -> tuple[dict[str, Any] | None, Fraction | None, str | None]:
    sampling = request.sampling
    assert sampling is not None
    if sampling.implementation == "trapezoidal_continuous_lockin":
        return _continuous_lockin_gram_proof(request)
    return _discrete_wls_gram_proof(request, weights)


def _dynamic_proof(
    request: R2BoundRequest,
    weights: tuple[Fraction, ...],
    amplification: Fraction,
) -> tuple[dict[str, Any] | None, str | None]:
    history = request.history
    dwell = request.dwell
    assert history is not None and dwell is not None
    tau = _fraction(request.tau_star_s)
    origin = _fraction(history.record_origin_s)
    switch = _fraction(history.rate_segment_start_s)
    segment_origin = (switch - origin) / tau
    sample_since_switch = tuple(
        (_fraction(value) - switch) / tau for value in request.sample_times_s
    )
    alpha_lower = _fraction(request.alpha_interval[0])
    alpha_upper = _fraction(request.alpha_interval[1])
    damping_lower = _fraction(request.lambda_interval[0])

    history_coefficient = (
        _fraction(history.initial_state_abs_bound)
        + _fraction(history.new_input_sup_bound) / damping_lower
        + _fraction(history.pre_record_history_coefficient)
    )
    dwell_coefficient = (
        _fraction(dwell.previous_input_sup_bound) + _fraction(dwell.new_input_sup_bound)
    ) / damping_lower
    history_input_values: list[Fraction] = []
    dwell_input_values: list[Fraction] = []
    history_relaxation: list[Fraction] = []
    dwell_relaxation: list[Fraction] = []
    for since_switch in sample_since_switch:
        history_time = segment_origin + since_switch
        history_alpha = alpha_upper if history_time <= 1 else alpha_lower
        dwell_alpha = alpha_upper if since_switch <= 1 else alpha_lower
        history_r = _relaxation_upper(history_time, history_alpha, damping_lower)
        dwell_r = _relaxation_upper(since_switch, dwell_alpha, damping_lower)
        if history_r is None or dwell_r is None:
            return None, "relaxation_envelope_undecided"
        history_relaxation.append(history_r)
        dwell_relaxation.append(dwell_r)
        history_input_values.append(history_coefficient * history_r)
        dwell_input_values.append(dwell_coefficient * dwell_r)

    if request.dynamic_front_end_mode == "identity":
        history_output_values = tuple(history_input_values)
        dwell_output_values = tuple(dwell_input_values)
        output_source = "identity_front_end_from_input_envelope"
    else:
        assert history.preconvolved_output_envelope is not None
        assert dwell.preconvolved_output_envelope is not None
        history_output_values = tuple(
            _fraction(value) for value in history.preconvolved_output_envelope
        )
        dwell_output_values = tuple(
            _fraction(value) for value in dwell.preconvolved_output_envelope
        )
        output_source = "independent_preconvolved_output_envelope"

    history_norm = _weighted_norm_upper(weights, history_output_values)
    dwell_norm = _weighted_norm_upper(weights, dwell_output_values)
    if history_norm is None or dwell_norm is None:
        return None, "dynamic_weighted_norm_nonfinite"
    history_radius = amplification * history_norm
    dwell_radius = amplification * dwell_norm
    return (
        {
            "dynamic_front_end_mode": request.dynamic_front_end_mode,
            "input_envelope_diagnostic": "simon-hyperbolic-gamma-unit-v1",
            "output_bound_source": output_source,
            "history": {
                "coefficient": _rat(history_coefficient),
                "relaxation_upper_by_sample": [_rat(value) for value in history_relaxation],
                "input_envelope_diagnostic_by_sample": [
                    _rat(value) for value in history_input_values
                ],
                "front_end_output_envelope_by_sample": [
                    _rat(value) for value in history_output_values
                ],
                "weighted_residual_upper": _rat(history_norm),
                "coefficient_radius_upper": _rat(history_radius),
            },
            "dwell": {
                "coefficient": _rat(dwell_coefficient),
                "relaxation_upper_by_sample": [_rat(value) for value in dwell_relaxation],
                "input_envelope_diagnostic_by_sample": [
                    _rat(value) for value in dwell_input_values
                ],
                "front_end_output_envelope_by_sample": [
                    _rat(value) for value in dwell_output_values
                ],
                "weighted_residual_upper": _rat(dwell_norm),
                "coefficient_radius_upper": _rat(dwell_radius),
            },
        },
        None,
    )


def _sampling_proof(
    request: R2BoundRequest,
    weights: tuple[Fraction, ...],
    amplification: Fraction,
) -> tuple[dict[str, Any] | None, str | None]:
    sampling = request.sampling
    delay = request.delay
    assert sampling is not None and delay is not None
    adc_half = _fraction(sampling.adc_quantization_step) / 2
    derivative = _fraction(sampling.first_derivative_bound_signal_per_s)
    per_sample = tuple(
        _fraction(anti_alias)
        + derivative * _fraction(time_error)
        + adc_half
        + _fraction(interpolation)
        for time_error, interpolation, anti_alias in zip(
            sampling.timestamp_error_s,
            sampling.interpolation_error_bound,
            sampling.anti_alias_error_bound,
            strict=True,
        )
    )
    sampling_norm = _weighted_norm_upper(weights, per_sample)
    if sampling_norm is None:
        return None, "sampling_weighted_norm_nonfinite"
    sampling_radius = amplification * sampling_norm

    if sampling.implementation == "discrete_wls":
        quadrature_radius = Fraction(0)
    else:
        times = tuple(_fraction(value) for value in request.sample_times_s)
        numerator = sum(
            (
                (right - left) ** 3 * _fraction(second_derivative)
                for left, right, second_derivative in zip(
                    times[:-1],
                    times[1:],
                    sampling.panel_second_derivative_bound,
                    strict=True,
                )
            ),
            Fraction(0),
        )
        duration = _fraction(request.right_r_peak_time_s) - _fraction(request.left_r_peak_time_s)
        transfer_lower = _fraction(delay.nominal_transfer_magnitude_lower_bound)
        quadrature_radius = numerator / (12 * duration * transfer_lower)
    return (
        {
            "implementation": sampling.implementation,
            "per_sample_output_error_upper": [_rat(value) for value in per_sample],
            "weighted_residual_upper": _rat(sampling_norm),
            "coefficient_radius_upper": _rat(sampling_radius),
            "quadrature_coefficient_radius_upper": _rat(quadrature_radius),
        },
        None,
    )


def _delay_proof(request: R2BoundRequest) -> tuple[dict[str, Any] | None, str | None]:
    delay = request.delay
    assert delay is not None
    duration = _fraction(request.right_r_peak_time_s) - _fraction(request.left_r_peak_time_s)
    omega = arb(2) * arb.pi() * request.complete_rr_intervals / _arb_fraction(duration)
    omega_upper = _upper_claim(omega)
    if omega_upper is None:
        return None, "nominal_rate_nonfinite"
    theta = abs(request.target_harmonic) * omega_upper * _fraction(
        delay.relative_group_delay_error_s
    ) + _fraction(delay.residual_phase_calibration_error_rad)
    theta_ball = _arb_fraction(theta)
    if theta == 0:
        chord = Fraction(0)
        chord_method = "exact_zero"
    elif theta_ball < arb.pi():
        chord = _upper_claim(arb(2) * (theta_ball / 2).sin())
        if chord is None:
            return None, "delay_chord_nonfinite"
        chord_method = "rigorous_sine"
    else:
        # The chord distance is globally at most two.  This branch is also a
        # safe fallback if theta is too close to pi for a strict comparison.
        chord = Fraction(2)
        chord_method = "global_chord_cap"
    amplitude = _fraction(delay.relative_amplitude_error_bound)
    relative_error = amplitude + (Fraction(1) + amplitude) * chord
    radius = _fraction(delay.target_harmonic_abs_upper_bound) * relative_error
    return (
        {
            "nominal_rate_rad_s_upper": _rat(omega_upper),
            "phase_error_rad_upper": _rat(theta),
            "chord_method": chord_method,
            "chord_upper": _rat(chord),
            "relative_transfer_error_upper": _rat(relative_error),
            "coefficient_radius_upper": _rat(radius),
        },
        None,
    )


def _compute_proof(
    request: R2BoundRequest,
    precision: int,
) -> tuple[dict[str, Any] | None, str | None]:
    weights = _normalized_weight_fractions(request.raw_weights)
    gram, amplification, reason = _gram_proof(request, weights)
    if reason is not None or gram is None or amplification is None:
        return None, reason or "gram_unknown"
    dynamic, reason = _dynamic_proof(request, weights, amplification)
    if reason is not None or dynamic is None:
        return None, reason or "dynamic_unknown"
    sampling, reason = _sampling_proof(request, weights, amplification)
    if reason is not None or sampling is None:
        return None, reason or "sampling_unknown"
    delay, reason = _delay_proof(request)
    if reason is not None or delay is None:
        return None, reason or "delay_unknown"

    partial_total = sum(
        (
            Fraction(dynamic["history"]["coefficient_radius_upper"]),
            Fraction(dynamic["dwell"]["coefficient_radius_upper"]),
            Fraction(sampling["coefficient_radius_upper"]),
            Fraction(sampling["quadrature_coefficient_radius_upper"]),
            Fraction(delay["coefficient_radius_upper"]),
        ),
        Fraction(0),
    )
    return (
        {
            "precision_bits": precision,
            "normalized_weights": [_rat(value) for value in weights],
            "gram": gram,
            "dynamic": dynamic,
            "sampling": sampling,
            "delay": delay,
            "partial_radius_upper": _rat(partial_total),
            "omitted_r2_components": [
                "hrv",
                "phase_anchor",
                "window_leakage_and_harmonic_tail",
                "filter_initialization",
                "measurement_noise",
                "model_residual",
            ],
        },
        None,
    )


def _seal(payload: dict[str, Any]) -> tuple[str, str]:
    digest = _sha256_json(payload)
    document = dict(payload)
    document["certificate_sha256"] = digest
    return _canonical_json(document), digest


def _unknown_certificate(
    manifest: dict[str, Any],
    reason: str,
    attempted: tuple[int, ...],
    precision: int,
) -> R2BoundCertificate:
    input_hash = _sha256_json(manifest)
    payload = {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "input_manifest": manifest,
        "input_sha256": input_hash,
        "relation": R2BoundRelation.UNKNOWN.value,
        "reason": reason,
        "precision_bits": precision,
        "attempted_precisions": list(attempted),
        "environment": _environment_manifest(),
    }
    certificate_json, digest = _seal(payload)
    return R2BoundCertificate(
        R2BoundRelation.UNKNOWN,
        reason,
        precision,
        attempted,
        input_hash,
        digest,
        certificate_json,
    )


def certify_r2_bounds(
    request: Any,
    *,
    precision_schedule: tuple[int, ...] = DEFAULT_PRECISION_SCHEDULE,
) -> R2BoundCertificate:
    """Issue a replayable component certificate or a fail-closed UNKNOWN."""

    schedule = _validate_precision_schedule(precision_schedule)
    manifest = _safe_manifest(request)
    if any(precision > MAX_REPLAY_PRECISION for precision in schedule):
        return _unknown_certificate(manifest, "precision_resource_limit", (), 0)
    invalid_reason = _validate_request(request)
    if invalid_reason is not None:
        return _unknown_certificate(manifest, invalid_reason, (), 0)

    attempted: list[int] = []
    final_reason = "precision_resource_limit"
    with _CERTIFICATION_LOCK:
        for precision in schedule:
            attempted.append(precision)
            with ctx.workprec(precision):
                try:
                    proof, reason = _compute_proof(request, precision)
                except (ArithmeticError, OverflowError, ValueError, ZeroDivisionError):
                    proof, reason = None, "ball_arithmetic_failure"
            if proof is None:
                final_reason = reason or "bound_undecided"
                continue
            input_hash = _sha256_json(manifest)
            payload = {
                "schema": SCHEMA,
                "algorithm": ALGORITHM,
                "input_manifest": manifest,
                "input_sha256": input_hash,
                "relation": R2BoundRelation.CERTIFIED_BOUND.value,
                "reason": "certified_all_requested_components",
                "precision_bits": precision,
                "attempted_precisions": list(attempted),
                "environment": _environment_manifest(),
                "proof": proof,
            }
            certificate_json, digest = _seal(payload)
            if len(certificate_json.encode("utf-8")) > MAX_CERTIFICATE_BYTES:
                return _unknown_certificate(
                    manifest,
                    "certificate_resource_bytes",
                    tuple(attempted),
                    precision,
                )
            return R2BoundCertificate(
                R2BoundRelation.CERTIFIED_BOUND,
                "certified_all_requested_components",
                precision,
                tuple(attempted),
                input_hash,
                digest,
                certificate_json,
            )
    return _unknown_certificate(
        manifest,
        final_reason,
        tuple(attempted),
        attempted[-1],
    )


def replay_r2_bound_certificate(
    request: R2BoundRequest,
    certificate: R2BoundCertificate | str,
) -> R2BoundReplayResult:
    """Rebuild the input manifest and independently recompute an issued proof."""

    certificate_json = (
        certificate.to_json() if isinstance(certificate, R2BoundCertificate) else certificate
    )
    if (
        not isinstance(certificate_json, str)
        or len(certificate_json.encode("utf-8")) > MAX_CERTIFICATE_BYTES
    ):
        return R2BoundReplayResult(False, R2BoundRelation.UNKNOWN, "invalid_certificate_envelope")
    try:
        try:
            document = json.loads(
                certificate_json,
                object_pairs_hook=_reject_duplicate_object_pairs,
            )
        except RecursionError as exc:
            raise ValueError("certificate nesting exceeds the replay limit") from exc
        if not isinstance(document, dict):
            raise ValueError("certificate must be an object")
        digest = document.get("certificate_sha256")
        if not _valid_sha(digest):
            raise ValueError("invalid digest")
        unsigned = dict(document)
        del unsigned["certificate_sha256"]
        if _sha256_json(unsigned) != digest:
            raise ValueError("digest mismatch")
        if document.get("schema") != SCHEMA or document.get("algorithm") != ALGORITHM:
            raise ValueError("schema mismatch")
        if document.get("environment") != _environment_manifest():
            raise ValueError("environment mismatch")
        relation = R2BoundRelation(document.get("relation"))
        manifest = _safe_manifest(request)
        input_hash = _sha256_json(manifest)
        if document.get("input_manifest") != manifest or document.get("input_sha256") != input_hash:
            raise ValueError("input mismatch")
        attempted_raw = document.get("attempted_precisions")
        if not isinstance(attempted_raw, list) or any(
            not isinstance(value, int) or isinstance(value, bool) for value in attempted_raw
        ):
            raise ValueError("attempt schedule invalid")
        attempted = tuple(attempted_raw)
        if attempted:
            _validate_precision_schedule(attempted)
        precision = document.get("precision_bits")
        if (
            not isinstance(precision, int)
            or isinstance(precision, bool)
            or precision < 0
            or precision > MAX_REPLAY_PRECISION
        ):
            raise ValueError("precision invalid")
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return R2BoundReplayResult(False, R2BoundRelation.UNKNOWN, "invalid_certificate_envelope")

    if relation is R2BoundRelation.UNKNOWN:
        return R2BoundReplayResult(True, R2BoundRelation.UNKNOWN, "replay_unknown")
    if (
        _validate_request(request) is not None
        or precision < 64
        or not attempted
        or attempted[-1] != precision
    ):
        return R2BoundReplayResult(False, R2BoundRelation.UNKNOWN, "invalid_certified_request")
    with _CERTIFICATION_LOCK:
        with ctx.workprec(precision):
            try:
                proof, reason = _compute_proof(request, precision)
            except (ArithmeticError, OverflowError, ValueError, ZeroDivisionError):
                proof, reason = None, "ball_arithmetic_failure"
    if proof is None or reason is not None:
        return R2BoundReplayResult(False, R2BoundRelation.UNKNOWN, "proof_no_longer_certifies")
    if document.get("proof") != proof:
        return R2BoundReplayResult(False, R2BoundRelation.UNKNOWN, "proof_mismatch")
    if document.get("reason") != "certified_all_requested_components":
        return R2BoundReplayResult(False, R2BoundRelation.UNKNOWN, "certified_reason_invalid")
    return R2BoundReplayResult(True, R2BoundRelation.CERTIFIED_BOUND, "replay_verified")


__all__ = [
    "ALGORITHM",
    "BOUND_NUMERIC_SCOPE",
    "SCHEMA",
    "DelayBoundInput",
    "DwellBoundInput",
    "EvidenceReference",
    "HistoryBoundInput",
    "R2BoundCertificate",
    "R2BoundRelation",
    "R2BoundReplayResult",
    "R2BoundRequest",
    "SamplingBoundInput",
    "certify_r2_bounds",
    "replay_r2_bound_certificate",
]
