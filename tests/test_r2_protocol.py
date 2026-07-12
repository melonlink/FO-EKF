from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import pytest

from fo_ekf.r2_protocol import (
    BOUND_NAMES,
    R2ProtocolStatus,
    _upward_product,
    _upward_ratio,
    _upward_sum,
    canonical_protocol_sha256,
    canonical_window_result_sha256,
    load_r2_protocol,
    validate_r2_protocol,
    validate_r2_protocol_file,
)

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "config" / "r2_protocol_template.toml"


def _up(value: float) -> float:
    return 0.0 if value == 0.0 else math.nextafter(value, math.inf)


def _up_sum(values: list[float]) -> float:
    return _up(math.fsum(values))


def _set_path(config: dict, path: tuple[str | int, ...], value: object) -> None:
    current: object = config
    for part in path[:-1]:
        current = current[part]  # type: ignore[index]
    current[path[-1]] = value  # type: ignore[index]


def _valid_deterministic_config() -> dict:
    config = deepcopy(load_r2_protocol(TEMPLATE))

    protocol = config["protocol"]
    protocol.update(
        {
            "protocol_id": "r2-unit-test",
            "frozen": True,
            "frozen_at_utc": "2026-07-11T00:00:00Z",
            "frozen_before_target_ecg_access": True,
            "git_commit": "a" * 40,
            "owner": "test",
            "requested_output_semantics": "deterministic",
        }
    )

    data = config["data_boundary"]
    data.update(
        {
            "calibration_split_id": "calibration-subjects",
            "identification_split_id": "identification-subjects",
            "validation_split_id": "validation-subjects",
            "calibration_disjoint_from_identification": True,
            "calibration_disjoint_from_validation": True,
            "identification_disjoint_from_validation": True,
            "subject_specific_identification_precedes_validation": True,
            "record_interval_hashes_required": True,
            "minimum_temporal_guard_s": 1.0,
            "temporal_guard_provided": True,
            "split_manifest_sha256": "b" * 64,
            "record_interval_manifest_sha256": "d" * 64,
        }
    )

    units = config["units"]
    units.update(
        {
            "signal_unit": "mV",
            "tau_star_s": 1.0,
            "tau_star_provided": True,
        }
    )

    estimator = config["estimator"]
    estimator.update(
        {
            "signal_representation": "real_two_sided_fourier",
            "retained_harmonics": [-2, -1, 0, 1, 2],
            "retained_harmonics_provided": True,
            "include_dc_if_real_signal": True,
            "include_negative_harmonics_if_real_signal": True,
            "weights_sum_to_one": True,
            "weight_rule": "normalized_trapezoidal",
            "gram_min_eigenvalue_threshold": 0.1,
            "gram_condition_number_max": 20.0,
            "gram_thresholds_provided": True,
        }
    )

    model = config["model"]
    model.update(
        {
            "alpha_interval": [0.4, 1.0],
            "alpha_interval_provided": True,
            "lambda_interval": [0.1, 5.0],
            "lambda_interval_provided": True,
            "q_magnitude_bounds_source": "independent calibration",
            "q_bounds_independent_of_target_fit": True,
            "q_abs_lower_by_harmonic": [0.05, 0.05, 0.0, 0.05, 0.05],
            "q_abs_upper_by_harmonic": [2.0, 2.0, 2.0, 2.0, 2.0],
        }
    )

    window = config["window_design"]
    window.update(
        {
            "complete_rr_intervals": 8,
            "complete_rr_intervals_provided": True,
            "burn_in_s_min": 10.0,
            "burn_in_provided": True,
            "dwell_s_min": 20.0,
            "dwell_provided": True,
            "maximum_observed_rr_phase_deviation_rad": 0.2,
            "phase_deviation_threshold_provided": True,
        }
    )

    config["preprocessing"].update(
        {
            "pipeline_manifest": "identity-front-end-v1",
            "pipeline_manifest_sha256": "c" * 64,
            "relative_ecg_rpeak_timing_calibrated": True,
        }
    )
    config["evidence_policy"]["joint_probability_method"] = "deterministic_zero"

    for name in BOUND_NAMES:
        bound = config["bounds"][name]
        bound.update(
            {
                "applicability_decided": True,
                "applicable": True,
                "provided": True,
                "semantics": "deterministic",
                "source_kind": "sensor_spec",
                "source_reference": f"independent-{name}-evidence",
                "calibration_method": f"frozen-{name}-upper-bound",
                "independent_of_identification_fit": True,
                "uses_target_fit_residuals": False,
                "failure_probability": 0.0,
            }
        )
        if "calibration_split_id" in bound:
            bound["calibration_split_id"] = data["calibration_split_id"]

    config["bounds"]["dwell"].update(
        {
            "burn_in_s": 12.0,
            "window_duration_s": 8.0,
        }
    )
    config["bounds"]["hrv"].update(
        {
            "r_peak_timing_error_s": 0.001,
            "r_peak_timing_error_provided": True,
            "q_harmonic_envelope_provided": True,
        }
    )
    config["bounds"]["phase_anchor"]["phase_error_disjoint_from_hrv_term"] = True
    config["bounds"]["window_leakage"].update(
        {
            "alias_factors_computed": True,
            "harmonic_tail_bound_provided": True,
            "front_end_transfer_envelope_provided": True,
        }
    )
    config["bounds"]["sampling"].update(
        {
            "sample_rate_hz": 250.0,
            "sample_rate_provided": True,
            "anti_alias_certified": True,
            "first_derivative_bound_signal_per_s": 5.0,
            "first_derivative_bound_provided": True,
            "implementation": "discrete_wls",
        }
    )
    config["bounds"]["delay"].update(
        {
            "ecg_to_rpeak_channel_alignment_included": True,
            "nominal_transfer_magnitude_lower_bound": 1.0,
            "target_harmonic_abs_upper_bound": 2.0,
        }
    )
    config["bounds"]["measurement_noise"].update(
        {
            "bound_type": "weighted_l2",
            "weighted_l2_bound": 0.0004,
        }
    )
    config["bounds"]["model_residual"].update(
        {
            "residual_location": "output",
            "output_weighted_l2_bound": 0.0004,
        }
    )

    source_radius = 0.0005
    for name in (
        "history",
        "dwell",
        "hrv",
        "phase_anchor",
        "sampling",
        "measurement_noise",
        "model_residual",
    ):
        config["bounds"][name]["computed_w_norm_bound"] = source_radius
    for name in ("window_leakage", "delay", "filter_initialization"):
        config["bounds"][name]["computed_coefficient_radius"] = source_radius
    config["bounds"]["sampling"]["quadrature_radius"] = 0.0

    thresholds = config["decision_thresholds"]
    thresholds.update(
        {
            "deterministic_total_radius_max": 0.1,
            "deterministic_relative_radius_max": 0.1,
            "thresholds_provided": True,
            "target_harmonic_abs_lower_bound": 0.2,
            "target_lower_bound_provided": True,
            "minimum_energy_snr": 20.0,
            "minimum_energy_snr_provided": True,
            "maximum_joint_failure_probability": 0.0,
            "probability_threshold_provided": True,
            "minimum_burn_in_s": 10.0,
            "minimum_dwell_s": 20.0,
            "dwell_thresholds_provided": True,
        }
    )

    probability = config["probability_accounting"]
    probability.update(
        {
            "term_failure_probabilities": [0.0] * len(BOUND_NAMES),
            "joint_failure_probability": 0.0,
            "joint_failure_probability_provided": True,
            "simultaneous_coverage_method": "deterministic_zero",
        }
    )

    gram_min = 0.9
    gram_sqrt_lower = math.nextafter(math.sqrt(gram_min), -math.inf)
    wls_row_norm_upper = _up(1.0 / gram_sqrt_lower)
    amplification = _up(wls_row_norm_upper / 1.0)
    amplified_radius = _up(amplification * source_radius)
    coefficient_radius = _up(source_radius)
    radii = {
        "radius_history": amplified_radius,
        "radius_dwell": amplified_radius,
        "radius_hrv": amplified_radius,
        "radius_phase_anchor": amplified_radius,
        "radius_window_leakage": coefficient_radius,
        "radius_sampling": amplified_radius,
        "radius_quadrature": 0.0,
        "radius_delay": coefficient_radius,
        "radius_filter_initialization": coefficient_radius,
        "radius_measurement_noise": amplified_radius,
        "radius_model_residual": amplified_radius,
    }
    deterministic_total = _up_sum(list(radii.values()))
    config["window_result"] = [
        {
            "window_id": "window-001",
            "rate_id": "rate-001",
            "subject_pseudonym": "subject-001",
            "lead_id": "lead-II",
            "physical_gain_id": "gain-calibration-001",
            "record_manifest_id": "record-manifest-001",
            "record_manifest_sha256": "e" * 64,
            "estimator_manifest_sha256": "f" * 64,
            "wls_execution_sha256": "1" * 64,
            "wls_window_sha256": "2" * 64,
            "wls_window_start_sample": 25_000,
            "wls_window_end_sample_exclusive": 27_000,
            "wls_sample_index_sha256": "3" * 64,
            "wls_digital_sample_sha256": "4" * 64,
            "wls_time_seconds_sha256": "5" * 64,
            "wls_normalized_weights_sha256": "6" * 64,
            "wls_interval_certificate_sha256": "7" * 64,
            "wls_exact_sample_functional_coverage": True,
            "wls_exact_functional_radius_upper": 1.0e-12,
            "protocol_sha256": "",
            "result_sha256": "",
            "left_r_peak_time_s": 100.0,
            "right_r_peak_time_s": 108.0,
            "complete_rr_intervals": 8,
            "window_duration_s": 8.0,
            "nominal_rate_rad_s": 2.0 * 3.141592653589793,
            "sample_count": 2000,
            "gram_min_eigenvalue": gram_min,
            "gram_condition_number": 1.2,
            "gram_bounds_outward_certified": True,
            "gram_certificate_sha256": "d" * 64,
            "wls_row_norm_upper_bound": wls_row_norm_upper,
            "target_harmonic": 1,
            "nominal_transfer_real": 1.0,
            "nominal_transfer_imag": 0.0,
            "z_tilde_real": 0.5,
            "z_tilde_imag": 0.1,
            **radii,
            "radius_total_deterministic": deterministic_total,
            "radius_total_probabilistic": 0.0,
            "joint_failure_probability": 0.0,
            "minimum_burn_in_pass": True,
            "minimum_snr_pass": True,
            "result_status": "PASS_DETERMINISTIC",
            "reason_codes": [],
        }
    ]

    _rehash(config)
    return config


def _rehash(config: dict) -> None:
    protocol_hash = canonical_protocol_sha256(config)
    config["protocol"]["config_sha256"] = protocol_hash
    for result in config.get("window_result", []):
        result["protocol_sha256"] = protocol_hash
        result["result_sha256"] = ""
        result["result_sha256"] = canonical_window_result_sha256(result)


def _valid_probability_config() -> dict:
    config = _valid_deterministic_config()
    config["protocol"]["requested_output_semantics"] = "probabilistic"
    noise = config["bounds"]["measurement_noise"]
    noise.update(
        {
            "semantics": "probabilistic",
            "source_kind": "confidence_interval",
            "failure_probability": 0.02,
        }
    )
    probabilities = [config["bounds"][name]["failure_probability"] for name in BOUND_NAMES]
    joint_probability = _up_sum(probabilities)
    config["evidence_policy"]["joint_probability_method"] = "union_bound"
    config["probability_accounting"].update(
        {
            "term_failure_probabilities": probabilities,
            "joint_failure_probability": joint_probability,
            "simultaneous_coverage_method": "union_bound",
        }
    )
    config["decision_thresholds"]["maximum_joint_failure_probability"] = 0.05
    window = config["window_result"][0]
    all_radius_fields = [
        key for key in window if key.startswith("radius_") and not key.startswith("radius_total_")
    ]
    deterministic_fields = [key for key in all_radius_fields if key != "radius_measurement_noise"]
    window.update(
        {
            "radius_total_deterministic": _up_sum([window[key] for key in deterministic_fields]),
            "radius_total_probabilistic": _up_sum([window[key] for key in all_radius_fields]),
            "joint_failure_probability": joint_probability,
            "result_status": "PASS_PROBABILISTIC",
        }
    )
    _rehash(config)
    return config


def test_unfilled_template_is_not_certifiable() -> None:
    result = validate_r2_protocol_file(TEMPLATE)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert "REJECT_MODEL" not in result.status.value
    assert result.reasons


def test_outward_helpers_do_not_round_positive_subnormal_results_to_zero() -> None:
    tiny = math.nextafter(0.0, math.inf)

    assert _upward_product(tiny, 0.5) >= tiny
    assert _upward_ratio(tiny, 2.0) >= tiny
    assert _upward_sum((tiny, tiny)) >= tiny


def test_deep_or_overflowing_protocol_material_fails_closed() -> None:
    deep = _valid_deterministic_config()
    cursor = deep
    for index in range(100):
        child: dict = {}
        cursor[f"nested-{index}"] = child
        cursor = child

    deep_result = validate_r2_protocol(deep)

    overflowing = _valid_deterministic_config()
    for key in tuple(overflowing["window_result"][0]):
        if key.startswith("radius_"):
            overflowing["window_result"][0][key] = 1.0e308
    _rehash(overflowing)
    overflow_result = validate_r2_protocol(overflowing)

    assert deep_result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert deep_result.reasons == ("root:resource_depth",)
    assert overflow_result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert overflow_result.reasons


def test_complete_deterministic_protocol_passes() -> None:
    result = validate_r2_protocol(_valid_deterministic_config())

    assert result.status is R2ProtocolStatus.PASS_DETERMINISTIC
    assert result.reasons == ()


def test_complete_probability_protocol_passes_with_joint_accounting() -> None:
    config = _valid_probability_config()

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.PASS_PROBABILISTIC
    assert result.reasons == ()


@pytest.mark.parametrize(
    ("field", "reason"),
    (
        (
            "calibration_disjoint_from_validation",
            "data_boundary.calibration_disjoint_from_validation:must_be_true",
        ),
        (
            "identification_disjoint_from_validation",
            "data_boundary.identification_disjoint_from_validation:must_be_true",
        ),
        (
            "subject_specific_identification_precedes_validation",
            "data_boundary.subject_specific_identification_precedes_validation:must_be_true",
        ),
        (
            "record_interval_hashes_required",
            "data_boundary.record_interval_hashes_required:must_be_true",
        ),
        (
            "calibration_subject_disjoint_from_evaluation",
            "data_boundary.calibration_subject_disjoint_from_evaluation:must_be_true",
        ),
        (
            "identification_validation_same_subject",
            "data_boundary.identification_validation_same_subject:must_be_true",
        ),
    ),
)
def test_all_split_and_temporal_leakage_guards_are_mandatory(
    field: str,
    reason: str,
) -> None:
    config = _valid_deterministic_config()
    config["data_boundary"][field] = False
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert reason in result.reasons


def test_record_interval_hash_and_temporal_guard_are_mandatory() -> None:
    missing_hash = _valid_deterministic_config()
    missing_hash["data_boundary"]["record_interval_manifest_sha256"] = ""
    _rehash(missing_hash)
    missing_guard = _valid_deterministic_config()
    missing_guard["data_boundary"]["temporal_guard_provided"] = False
    _rehash(missing_guard)

    hash_result = validate_r2_protocol(missing_hash)
    guard_result = validate_r2_protocol(missing_guard)

    assert hash_result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert "data_boundary.record_interval_manifest_sha256:invalid_sha256" in hash_result.reasons
    assert guard_result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert "data_boundary.temporal_guard_provided:must_be_true" in guard_result.reasons


def test_continuous_l2_noise_cannot_be_sampled_by_discrete_wls_without_a_bridge() -> None:
    config = _valid_deterministic_config()
    config["bounds"]["measurement_noise"].update(
        {
            "bound_type": "continuous_l2",
            "continuous_l2_bound_signal_sqrt_s": 0.001,
        }
    )
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert "bounds.measurement_noise.continuous_l2:requires_continuous_lockin" in result.reasons


@pytest.mark.parametrize(
    "field",
    (
        "record_manifest_sha256",
        "estimator_manifest_sha256",
        "wls_execution_sha256",
        "wls_window_sha256",
        "wls_sample_index_sha256",
        "wls_digital_sample_sha256",
        "wls_time_seconds_sha256",
        "wls_normalized_weights_sha256",
        "wls_interval_certificate_sha256",
    ),
)
def test_each_window_binds_record_and_estimator_manifests(field: str) -> None:
    config = _valid_deterministic_config()
    config["window_result"][0][field] = ""
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert f"window_result[0].{field}:invalid_sha256" in result.reasons


def test_deterministic_pass_requires_exact_wls_coverage_and_numeric_radius_floor() -> None:
    missing_coverage = _valid_deterministic_config()
    missing_coverage["window_result"][0]["wls_exact_sample_functional_coverage"] = False
    _rehash(missing_coverage)
    understated_sampling = _valid_deterministic_config()
    understated_sampling["window_result"][0]["wls_exact_functional_radius_upper"] = 1.0
    _rehash(understated_sampling)
    invalid_geometry = _valid_deterministic_config()
    invalid_geometry["window_result"][0]["wls_window_end_sample_exclusive"] = invalid_geometry[
        "window_result"
    ][0]["wls_window_start_sample"]
    _rehash(invalid_geometry)

    coverage_result = validate_r2_protocol(missing_coverage)
    radius_result = validate_r2_protocol(understated_sampling)
    geometry_result = validate_r2_protocol(invalid_geometry)

    assert coverage_result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert (
        "window_result[0].wls_exact_sample_functional_coverage:must_be_true"
        in coverage_result.reasons
    )
    assert radius_result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert "window_result[0].radius_sampling:below_exact_functional_radius" in radius_result.reasons
    assert geometry_result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert "window_result[0].wls_integer_window_geometry:invalid" in geometry_result.reasons


def test_target_fit_residual_bound_is_not_certifiable() -> None:
    config = _valid_deterministic_config()
    config["bounds"]["model_residual"]["uses_target_fit_residuals"] = True
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("uses_target_fit_residuals" in reason for reason in result.reasons)


def test_surface_complete_protocol_with_unset_window_is_not_certifiable() -> None:
    config = _valid_deterministic_config()
    config["window_result"][0]["window_id"] = "UNSET"
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("window_id" in reason for reason in result.reasons)


def test_missing_bound_specific_evidence_is_not_certifiable() -> None:
    config = _valid_deterministic_config()
    config["bounds"]["window_leakage"]["harmonic_tail_bound_provided"] = False
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("harmonic_tail_bound_provided" in reason for reason in result.reasons)


def test_explicit_acausal_processing_has_priority_over_missing_fields() -> None:
    config = load_r2_protocol(TEMPLATE)
    config["preprocessing"]["all_steps_causal"] = False
    config["bounds"]["model_residual"]["uses_target_fit_residuals"] = True

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.EXCLUDE_PROTOCOL
    assert result.reasons == ("preprocessing.all_steps_causal:explicitly_false",)


def test_nonfinite_numeric_bound_is_not_certifiable() -> None:
    config = _valid_deterministic_config()
    config["bounds"]["history"]["computed_w_norm_bound"] = float("nan")
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any(reason.endswith(":nonfinite") for reason in result.reasons)


def test_canonical_hash_ignores_hash_field_and_mapping_order() -> None:
    config = _valid_deterministic_config()
    baseline = canonical_protocol_sha256(config)

    config["protocol"]["config_sha256"] = "f" * 64
    assert canonical_protocol_sha256(config) == baseline

    reversed_top_level = dict(reversed(list(config.items())))
    assert canonical_protocol_sha256(reversed_top_level) == baseline

    config["units"]["tau_star_s"] = 2.0
    assert canonical_protocol_sha256(config) != baseline


def test_protocol_hash_excludes_results_but_result_hash_detects_mutation() -> None:
    config = _valid_deterministic_config()
    protocol_hash = canonical_protocol_sha256(config)
    result_hash = config["window_result"][0]["result_sha256"]

    config["window_result"][0]["z_tilde_real"] += 0.1

    assert canonical_protocol_sha256(config) == protocol_hash
    invalid = validate_r2_protocol(config)
    assert invalid.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("result_sha256:mismatch" in reason for reason in invalid.reasons)

    _rehash(config)
    repaired = validate_r2_protocol(config)
    assert repaired.status is R2ProtocolStatus.PASS_DETERMINISTIC
    assert config["window_result"][0]["protocol_sha256"] == protocol_hash
    assert config["window_result"][0]["result_sha256"] != result_hash


def test_source_bound_cannot_be_hidden_by_smaller_reported_radius() -> None:
    config = _valid_deterministic_config()
    config["bounds"]["history"]["computed_w_norm_bound"] = 1.0e9
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("radius_history:below_recomputed_floor" in reason for reason in result.reasons)


def test_all_zero_reported_radii_do_not_pass_nonzero_source_bounds() -> None:
    config = _valid_deterministic_config()
    for key in config["window_result"][0]:
        if key.startswith("radius_"):
            config["window_result"][0][key] = 0.0
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("below_recomputed_floor" in reason for reason in result.reasons)


def test_total_radius_may_not_round_below_outward_sum() -> None:
    config = _valid_deterministic_config()
    total = config["window_result"][0]["radius_total_deterministic"]
    config["window_result"][0]["radius_total_deterministic"] = math.nextafter(total, -math.inf)
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("radius_total_deterministic:below_sum" in reason for reason in result.reasons)


def test_wls_row_norm_must_dominate_outward_gram_spectral_floor() -> None:
    config = _valid_deterministic_config()
    gram_min = config["window_result"][0]["gram_min_eigenvalue"]
    config["window_result"][0]["wls_row_norm_upper_bound"] = 1.0 / math.sqrt(gram_min)
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("wls_row_norm_upper_bound:below_gram_floor" in reason for reason in result.reasons)


@pytest.mark.parametrize(
    ("field", "changed_value"),
    [
        ("subject_pseudonym", "different-subject"),
        ("lead_id", "different-lead"),
        ("physical_gain_id", "different-gain"),
    ],
)
def test_joint_rows_must_share_subject_lead_and_gain(
    field: str,
    changed_value: str,
) -> None:
    config = _valid_deterministic_config()
    second = deepcopy(config["window_result"][0])
    second["window_id"] = "window-002"
    second["rate_id"] = "rate-002"
    second[field] = changed_value
    config["window_result"].append(second)
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("joint_identity_mismatch" in reason for reason in result.reasons)


def test_window_harmonic_key_must_be_unique_but_same_window_other_harmonic_is_allowed() -> None:
    config = _valid_deterministic_config()
    duplicate = deepcopy(config["window_result"][0])
    config["window_result"].append(duplicate)
    _rehash(config)

    rejected = validate_r2_protocol(config)
    assert rejected.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("duplicate_window_harmonic" in reason for reason in rejected.reasons)

    config = _valid_deterministic_config()
    other_harmonic = deepcopy(config["window_result"][0])
    other_harmonic["target_harmonic"] = 2
    config["window_result"].append(other_harmonic)
    _rehash(config)

    accepted = validate_r2_protocol(config)
    assert accepted.status is R2ProtocolStatus.PASS_DETERMINISTIC


def test_calibration_bound_may_not_use_identification_or_validation_split() -> None:
    for leaked_split_key in ("identification_split_id", "validation_split_id"):
        config = _valid_deterministic_config()
        bound = config["bounds"]["measurement_noise"]
        bound["source_kind"] = "independent_bounded_calibration"
        bound["calibration_split_id"] = config["data_boundary"][leaked_split_key]
        _rehash(config)

        result = validate_r2_protocol(config)

        assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
        assert any("calibration_split_id:data_leakage" in reason for reason in result.reasons)


def test_configured_source_allowlist_cannot_authorize_target_fit_residuals() -> None:
    config = _valid_deterministic_config()
    config["evidence_policy"]["deterministic_source_kinds"] = ["same_target_fit_residual"]
    for bound in config["bounds"].values():
        bound["source_kind"] = "same_target_fit_residual"
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("not_fixed_allowlist" in reason for reason in result.reasons)


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("model", "q_bounds_unit"), "kg"),
        (("bounds", "phase_anchor", "common_radian_gauge_absorbed_into_q"), False),
        (("bounds", "hrv", "phase_defined_operationally_from_r_peaks"), False),
        (("bounds", "window_leakage", "noninteger_phase_length_retained"), False),
        (("estimator", "kind"), "UNSET"),
        (("estimator", "phase_model"), "UNSET"),
        (("window_design", "left_anchor"), "UNSET"),
        (("window_design", "right_anchor"), "UNSET"),
        (("window_design", "fixed_seconds_window_allowed"), True),
        (("preprocessing", "filter_state_policy"), "UNSET"),
        (("audit_outputs", "write_bound_breakdown"), False),
        (("units", "physical_time_unit"), "UNSET"),
        (("bounds", "sampling", "missing_sample_policy"), "interpolate"),
        (("bounds", "history", "input_and_state_unit"), "kg"),
        (("decision_rules", "missing_required_bound"), "PASS_DETERMINISTIC"),
        (("window_result", 0, "nominal_transfer_real"), 0.5),
        (("window_result", 0, "gram_bounds_outward_certified"), False),
        (("window_result", 0, "gram_certificate_sha256"), "not-a-hash"),
        (("evidence_policy", "bound_numeric_scope"), "single_window_only"),
    ],
)
def test_fixed_theory_and_audit_contracts_fail_closed(
    path: tuple[str | int, ...],
    value: object,
) -> None:
    config = _valid_deterministic_config()
    _set_path(config, path, value)
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("frozen_at_utc", "not-a-date-after-target-access"),
        ("git_commit", "x"),
    ],
)
def test_freeze_metadata_must_be_machine_verifiable(field: str, value: str) -> None:
    config = _valid_deterministic_config()
    config["protocol"][field] = value
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE


def test_decision_dwell_thresholds_must_cover_frozen_window_design() -> None:
    config = _valid_deterministic_config()
    config["window_design"]["burn_in_s_min"] = 1.0e9
    config["window_design"]["dwell_s_min"] = 1.0e9
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("below_window_design" in reason for reason in result.reasons)


def test_probability_method_must_match_preregistered_method() -> None:
    config = _valid_probability_config()
    config["probability_accounting"]["simultaneous_coverage_method"] = (
        "directly_calibrated_joint_bound"
    )
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("not_preregistered" in reason for reason in result.reasons)


def test_direct_joint_method_requires_independent_evidence() -> None:
    config = _valid_probability_config()
    config["evidence_policy"]["joint_probability_method"] = "directly_calibrated_joint_bound"
    config["probability_accounting"]["simultaneous_coverage_method"] = (
        "directly_calibrated_joint_bound"
    )
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("direct_joint" in reason for reason in result.reasons)


def test_direct_joint_method_passes_only_with_bound_independent_evidence() -> None:
    config = _valid_probability_config()
    config["evidence_policy"]["joint_probability_method"] = "directly_calibrated_joint_bound"
    probability = config["probability_accounting"]
    probability.update(
        {
            "simultaneous_coverage_method": "directly_calibrated_joint_bound",
            "joint_failure_probability": 0.01,
            "direct_joint_reference": "independent simultaneous calibration",
            "direct_joint_calibration_split_id": config["data_boundary"]["calibration_split_id"],
            "direct_joint_independent_of_identification_fit": True,
            "direct_joint_evidence_sha256": "e" * 64,
        }
    )
    config["window_result"][0]["joint_failure_probability"] = 0.01
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.PASS_PROBABILISTIC


def test_union_probability_may_not_round_below_outward_term_sum() -> None:
    config = _valid_probability_config()
    config["probability_accounting"]["joint_failure_probability"] = 0.02
    config["window_result"][0]["joint_failure_probability"] = 0.02
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert "probability_accounting:invalid_union_bound" in result.reasons


def test_probability_coverage_scope_must_cover_all_selected_disks() -> None:
    config = _valid_probability_config()
    config["probability_accounting"]["coverage_scope"] = "per_window_only"
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any("coverage_scope" in reason for reason in result.reasons)


@pytest.mark.parametrize(
    ("bound_name", "raw_field", "raw_value", "expected_reason"),
    [
        (
            "measurement_noise",
            "weighted_l2_bound",
            1.0,
            "bounds.measurement_noise.computed_w_norm_bound:below_raw_floor",
        ),
        (
            "model_residual",
            "output_weighted_l2_bound",
            1.0,
            "bounds.model_residual.computed_w_norm_bound:below_raw_floor",
        ),
        (
            "sampling",
            "anti_alias_residual_bound",
            1.0,
            "bounds.sampling.computed_w_norm_bound:below_raw_floor",
        ),
        (
            "delay",
            "relative_group_delay_error_s",
            1.0,
            "bounds.delay.computed_coefficient_radius:below_raw_floor",
        ),
    ],
)
def test_primitive_bounds_cannot_be_hidden_by_zero_or_small_computed_bounds(
    bound_name: str,
    raw_field: str,
    raw_value: float,
    expected_reason: str,
) -> None:
    config = _valid_deterministic_config()
    config["bounds"][bound_name][raw_field] = raw_value
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert expected_reason in result.reasons


@pytest.mark.parametrize(
    ("threshold_field", "threshold_value", "result_flag", "expected_reason"),
    [
        (
            "minimum_burn_in_s",
            13.0,
            "minimum_burn_in_pass",
            "minimum_burn_in_or_dwell_not_met",
        ),
        (
            "minimum_energy_snr",
            1.0e9,
            "minimum_snr_pass",
            "minimum_snr_not_met",
        ),
    ],
)
def test_consistent_false_minimum_gate_cannot_receive_pass_status(
    threshold_field: str,
    threshold_value: float,
    result_flag: str,
    expected_reason: str,
) -> None:
    config = _valid_deterministic_config()
    config["decision_thresholds"][threshold_field] = threshold_value
    config["window_result"][0][result_flag] = False
    _rehash(config)

    result = validate_r2_protocol(config)

    assert result.status is R2ProtocolStatus.NOT_CERTIFIABLE
    assert any(expected_reason in reason for reason in result.reasons)
