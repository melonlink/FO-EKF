"""Self-contained deterministic material for the exact certificate-chain demo."""

from __future__ import annotations

import cmath
import hashlib
import json
import math
import re
from copy import deepcopy
from dataclasses import dataclass, replace
from fractions import Fraction
from pathlib import Path

import numpy as np

from .certified_set import (
    CertifiedHarmonicData,
    CertifiedJointProblem,
    CertifiedParameterBox,
    ClosedInterval,
    RateOffsetInterval,
)
from .r2_bound_certificate import (
    BOUND_NUMERIC_SCOPE,
    DelayBoundInput,
    DwellBoundInput,
    EvidenceReference,
    HistoryBoundInput,
    R2BoundRelation,
    R2BoundRequest,
    SamplingBoundInput,
    certify_r2_bounds,
)
from .r2_protocol import (
    BOUND_NAMES,
    R2ProtocolStatus,
    canonical_protocol_sha256,
    canonical_window_result_sha256,
    load_r2_protocol,
    validate_r2_protocol,
)
from .r2_t1_bundle import R2DiskCertificateLink
from .tube_certificate import TubeRelation, certify_alpha_tube
from .wls_execution import build_wls_execution_record, estimator_identity
from .wls_interval_certificate import build_wls_interval_certificate

ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "config" / "r2_protocol_template.toml"


@dataclass(frozen=True)
class SyntheticCertificateCase:
    """All in-memory inputs needed to issue and replay one bundle."""

    problem: CertifiedJointProblem
    box: CertifiedParameterBox
    baseline_certificate_json: str
    protocol: dict
    links: tuple[R2DiskCertificateLink, ...]


def _up(value: float) -> float:
    return 0.0 if value == 0.0 else math.nextafter(value, math.inf)


def _up_sum(values: list[float]) -> float:
    return _up(math.fsum(values))


def _down_fraction(value: Fraction) -> float:
    candidate = float(value)
    if Fraction.from_float(candidate) > value:
        candidate = math.nextafter(candidate, -math.inf)
    return candidate


def _up_fraction(value: Fraction) -> float:
    candidate = float(value)
    if Fraction.from_float(candidate) < value:
        candidate = math.nextafter(candidate, math.inf)
    return candidate


def _rehash(config: dict) -> None:
    protocol_hash = canonical_protocol_sha256(config)
    config["protocol"]["config_sha256"] = protocol_hash
    for result in config.get("window_result", []):
        result["protocol_sha256"] = protocol_hash
        result["result_sha256"] = ""
        result["result_sha256"] = canonical_window_result_sha256(result)


def _calibration_evidence(label: str) -> EvidenceReference:
    return EvidenceReference(
        source_kind="proved_engineering_envelope",
        source_reference=f"synthetic_generator_contract/{label}",
        evidence_sha256=hashlib.sha256(label.encode()).hexdigest(),
        calibration_split_id="calibration-v1",
        applicability_domain=(
            "directly generated steady-state retained-harmonic synthetic case; not real ECG"
        ),
        signal_unit="mV",
        time_unit="s",
        independent_of_identification_fit=True,
        uses_target_fit_residuals=False,
        numeric_scope=BOUND_NUMERIC_SCOPE,
    )


def _base_request() -> R2BoundRequest:
    sample_times = tuple(20.0 + 0.25 * index for index in range(8))
    count = len(sample_times)
    return R2BoundRequest(
        window_id="window-001",
        protocol_sha256="a" * 64,
        record_manifest_sha256="b" * 64,
        estimator_manifest_sha256="c" * 64,
        calibration_split_id="calibration-v1",
        identification_split_id="identification-v1",
        frozen_before_target_ecg_access=True,
        signal_unit="mV",
        dynamic_front_end_mode="identity",
        tau_star_s=1.0,
        alpha_interval=(0.6, 0.9),
        lambda_interval=(0.4, 1.2),
        left_r_peak_time_s=20.0,
        right_r_peak_time_s=22.0,
        complete_rr_intervals=2,
        sample_times_s=sample_times,
        raw_weights=(1.0,) * count,
        retained_harmonics=(-1, 0, 1),
        target_harmonic=1,
        gram_min_eigenvalue_threshold=0.5,
        gram_condition_number_max=2.0,
        history=HistoryBoundInput(
            record_origin_s=0.0,
            rate_segment_start_s=0.0,
            initial_state_abs_bound=0.0,
            new_input_sup_bound=0.0,
            pre_record_history_coefficient=0.0,
            evidence=_calibration_evidence("history"),
        ),
        dwell=DwellBoundInput(
            previous_input_sup_bound=0.0,
            new_input_sup_bound=0.0,
            evidence=_calibration_evidence("dwell"),
        ),
        sampling=SamplingBoundInput(
            timestamp_error_s=(0.0,) * count,
            interpolation_error_bound=(0.0,) * count,
            anti_alias_error_bound=(0.0,) * count,
            adc_quantization_step=1.0e-9,
            first_derivative_bound_signal_per_s=2.0,
            implementation="discrete_wls",
            quadrature_decomposition_target="not_applicable",
            panel_second_derivative_bound=(),
            evidence=_calibration_evidence("sampling"),
        ),
        delay=DelayBoundInput(
            relative_amplitude_error_bound=0.0,
            relative_group_delay_error_s=0.0,
            residual_phase_calibration_error_rad=0.0,
            target_harmonic_abs_upper_bound=2.0,
            nominal_transfer_magnitude_lower_bound=1.0,
            ecg_to_rpeak_channel_alignment_included=True,
            evidence=_calibration_evidence("delay"),
        ),
        sample_indices=tuple(range(80, 88)),
        window_start_sample=80,
        window_end_sample_exclusive=88,
    )


def _synthetic_problem() -> tuple[CertifiedJointProblem, CertifiedParameterBox]:
    order = 0.7
    damping = 0.5
    frequency = 2.0 * math.pi
    morphology = 1.0 - 0.2j
    fractional_frequency = frequency**order * cmath.exp(0.5j * math.pi * order)
    response = morphology / (damping + fractional_frequency)
    harmonic = CertifiedHarmonicData(
        rate_ids=("rest",),
        base_frequencies=(frequency,),
        harmonic_index=1,
        measured_responses=(response,),
        response_radii=(1.5e-9,),
        morphology_drift_radii=(0.0,),
        morphology_bounds=(0.1, 5.0),
        label="m1",
    )
    return (
        CertifiedJointProblem((harmonic,), "rest"),
        CertifiedParameterBox(
            ClosedInterval(order, order),
            ClosedInterval(damping, damping),
            (RateOffsetInterval("rest", ClosedInterval(0.0, 0.0)),),
        ),
    )


def _protocol_without_results(
    request: R2BoundRequest,
    *,
    baseline_certificate_sha256: str,
    baseline_input_sha256: str,
    git_commit: str,
) -> dict:
    config = deepcopy(load_r2_protocol(TEMPLATE))
    config["protocol"].update(
        {
            "protocol_id": "r2-synthetic-exact-chain",
            "frozen": True,
            "frozen_at_utc": "2026-07-12T00:00:00Z",
            "frozen_before_target_ecg_access": True,
            "git_commit": git_commit,
            "owner": "synthetic-certificate-runner",
            "requested_output_semantics": "deterministic",
            "baseline_t1_certificate_sha256": baseline_certificate_sha256,
            "baseline_t1_input_sha256": baseline_input_sha256,
        }
    )
    config["data_boundary"].update(
        {
            "calibration_split_id": request.calibration_split_id,
            "identification_split_id": request.identification_split_id,
            "validation_split_id": "validation-v1",
            "calibration_disjoint_from_identification": True,
            "calibration_disjoint_from_validation": True,
            "identification_disjoint_from_validation": True,
            "subject_specific_identification_precedes_validation": True,
            "record_interval_hashes_required": True,
            "minimum_temporal_guard_s": 1.0,
            "temporal_guard_provided": True,
            "split_manifest_sha256": "1" * 64,
            "record_interval_manifest_sha256": "2" * 64,
        }
    )
    config["units"].update(
        {
            "signal_unit": request.signal_unit,
            "tau_star_s": request.tau_star_s,
            "tau_star_provided": True,
        }
    )
    config["estimator"].update(
        {
            "signal_representation": "real_two_sided_fourier",
            "retained_harmonics": list(request.retained_harmonics),
            "retained_harmonics_provided": True,
            "include_dc_if_real_signal": True,
            "include_negative_harmonics_if_real_signal": True,
            "weights_sum_to_one": True,
            "weight_rule": "normalized_identity",
            "gram_min_eigenvalue_threshold": request.gram_min_eigenvalue_threshold,
            "gram_condition_number_max": request.gram_condition_number_max,
            "gram_thresholds_provided": True,
        }
    )
    config["model"].update(
        {
            "alpha_interval": list(request.alpha_interval),
            "alpha_interval_provided": True,
            "lambda_interval": list(request.lambda_interval),
            "lambda_interval_provided": True,
            "q_magnitude_bounds_source": "independent synthetic calibration",
            "q_bounds_independent_of_target_fit": True,
            "q_abs_lower_by_harmonic": [0.1, 0.0, 0.1],
            "q_abs_upper_by_harmonic": [5.0, 2.0, 5.0],
            "t1_rate_ids": ["rest"],
            "damping_anchor_rate_id": "rest",
        }
    )
    config["window_design"].update(
        {
            "complete_rr_intervals": request.complete_rr_intervals,
            "complete_rr_intervals_provided": True,
            "burn_in_s_min": 10.0,
            "burn_in_provided": True,
            "dwell_s_min": 14.0,
            "dwell_provided": True,
            "maximum_observed_rr_phase_deviation_rad": 0.2,
            "phase_deviation_threshold_provided": True,
        }
    )
    config["preprocessing"].update(
        {
            "pipeline_manifest": "identity-front-end-v1",
            "pipeline_manifest_sha256": "3" * 64,
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
                "source_kind": "proved_engineering_envelope",
                "source_reference": f"synthetic_generator_contract/{name}",
                "calibration_method": f"direct-construction-{name}-upper-bound",
                "independent_of_identification_fit": True,
                "uses_target_fit_residuals": False,
                "failure_probability": 0.0,
            }
        )
        if "calibration_split_id" in bound:
            bound["calibration_split_id"] = request.calibration_split_id
    config["bounds"]["dwell"].update(
        {
            "burn_in_s": 12.0,
            "window_duration_s": request.right_r_peak_time_s - request.left_r_peak_time_s,
        }
    )
    config["bounds"]["hrv"].update(
        {
            "r_peak_timing_error_s": 0.0,
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
            "sample_rate_hz": 4.0,
            "sample_rate_provided": True,
            "anti_alias_certified": True,
            "anti_alias_residual_bound": 0.0,
            "adc_quantization_step": 1.0e-9,
            "timestamp_jitter_s": 0.0,
            "r_peak_endpoint_interpolation_bound": 0.0,
            "first_derivative_bound_signal_per_s": 2.0,
            "first_derivative_bound_provided": True,
            "implementation": "discrete_wls",
            "quadrature_decomposition_target": "not_applicable",
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
            "weighted_l2_bound": 0.0,
        }
    )
    config["bounds"]["model_residual"].update(
        {
            "residual_location": "output",
            "output_weighted_l2_bound": 0.0,
        }
    )
    source_radius = 0.0
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
    config["decision_thresholds"].update(
        {
            "deterministic_total_radius_max": 100.0,
            "deterministic_relative_radius_max": 1000.0,
            "thresholds_provided": True,
            "target_harmonic_abs_lower_bound": 0.2,
            "target_lower_bound_provided": True,
            "minimum_energy_snr": 20.0,
            "minimum_energy_snr_provided": True,
            "maximum_joint_failure_probability": 0.0,
            "probability_threshold_provided": True,
            "minimum_burn_in_s": 10.0,
            "minimum_dwell_s": 14.0,
            "dwell_thresholds_provided": True,
        }
    )
    config["probability_accounting"].update(
        {
            "term_failure_probabilities": [0.0] * len(BOUND_NAMES),
            "joint_failure_probability": 0.0,
            "joint_failure_probability_provided": True,
            "simultaneous_coverage_method": "deterministic_zero",
        }
    )
    for name, group in (
        ("history", request.history),
        ("dwell", request.dwell),
        ("sampling", request.sampling),
        ("delay", request.delay),
    ):
        assert group is not None
        evidence = group.evidence
        config["bounds"][name].update(
            {
                "source_kind": evidence.source_kind,
                "source_reference": evidence.source_reference,
                "calibration_split_id": evidence.calibration_split_id,
                "independent_of_identification_fit": evidence.independent_of_identification_fit,
                "uses_target_fit_residuals": evidence.uses_target_fit_residuals,
                "evidence_sha256": evidence.evidence_sha256,
            }
        )
    config["window_result"] = []
    return config


def _proof(certificate_json: str) -> dict:
    return json.loads(certificate_json)["proof"]


def build_synthetic_certificate_case(*, git_commit: str) -> SyntheticCertificateCase:
    """Build one small, frozen case using production certificate APIs only."""

    if type(git_commit) is not str or re.fullmatch(r"[0-9a-f]{40}", git_commit) is None:
        raise ValueError("git_commit must be an exact 40-character lowercase hex digest")

    problem, box = _synthetic_problem()
    baseline = certify_alpha_tube(problem, box)
    if baseline.relation is not TubeRelation.ROBUST_INNER:
        raise RuntimeError("synthetic T1 baseline did not produce ROBUST_INNER")
    baseline_document = json.loads(baseline.to_json())

    estimator_hash = estimator_identity()["estimator_sha256"]
    if not isinstance(estimator_hash, str):
        raise RuntimeError("WLS estimator identity is malformed")
    request0 = replace(_base_request(), estimator_manifest_sha256=estimator_hash)
    config = _protocol_without_results(
        request0,
        baseline_certificate_sha256=baseline_document["certificate_sha256"],
        baseline_input_sha256=baseline_document["input_sha256"],
        git_commit=git_commit,
    )

    provisional_hash = canonical_protocol_sha256(config)
    provisional_request = replace(request0, protocol_sha256=provisional_hash)
    provisional = certify_r2_bounds(provisional_request)
    if provisional.relation is not R2BoundRelation.CERTIFIED_BOUND:
        raise RuntimeError(f"provisional R2 component failed: {provisional.reason}")
    provisional_proof = _proof(provisional.to_json())
    dynamic = provisional_proof["dynamic"]
    sampling = provisional_proof["sampling"]
    delay = provisional_proof["delay"]
    config["bounds"]["history"]["computed_w_norm_bound"] = _up_fraction(
        Fraction(dynamic["history"]["weighted_residual_upper"])
    )
    config["bounds"]["dwell"]["computed_w_norm_bound"] = _up_fraction(
        Fraction(dynamic["dwell"]["weighted_residual_upper"])
    )
    config["bounds"]["sampling"]["computed_w_norm_bound"] = _up_fraction(
        Fraction(sampling["weighted_residual_upper"])
    )
    config["bounds"]["sampling"]["computed_w_norm_bound"] = max(
        config["bounds"]["sampling"]["computed_w_norm_bound"],
        1.0e-9,
    )
    config["bounds"]["sampling"]["quadrature_radius"] = _up_fraction(
        Fraction(sampling["quadrature_coefficient_radius_upper"])
    )
    config["bounds"]["delay"]["computed_coefficient_radius"] = _up_fraction(
        Fraction(delay["coefficient_radius_upper"])
    )

    _rehash(config)
    protocol_hash = canonical_protocol_sha256(config)
    request = replace(request0, protocol_sha256=protocol_hash)
    component = certify_r2_bounds(request)
    if component.relation is not R2BoundRelation.CERTIFIED_BOUND:
        raise RuntimeError(f"final R2 component failed: {component.reason}")
    component_proof = _proof(component.to_json())
    gram = component_proof["gram"]

    sample_indices = np.asarray(request.sample_indices, dtype=np.int64)
    phase = (
        2.0
        * np.pi
        * request.complete_rr_intervals
        * (sample_indices - int(request.window_start_sample))
        / (int(request.window_end_sample_exclusive) - int(request.window_start_sample))
    )
    true_response = problem.harmonics[0].measured_responses[0]
    adc_gain = 1.0e9
    physical = 2.0 * np.real(true_response * np.exp(1j * phase))
    digital = np.rint(adc_gain * physical).astype(np.int64)
    wls_window_sha256 = hashlib.sha256(b"window-001-wls-window-manifest").hexdigest()
    execution = build_wls_execution_record(
        digital,
        (80, 84, 88),
        sampling_frequency=4.0,
        adc_gain=adc_gain,
        baseline=0,
        harmonics=(1,),
        protocol_sha256=protocol_hash,
        record_sha256=request.record_manifest_sha256,
        window_sha256=wls_window_sha256,
    )
    interval = build_wls_interval_certificate(execution, digital, precision_bits=192)
    if interval["relation"] != "CERTIFIED_INTERVAL":
        raise RuntimeError(f"exact WLS interval failed: {interval['reason']}")
    exact_numeric = interval["response_disks"][0]["absolute_error_upper"]["float_upper"]

    gram_min = _down_fraction(Fraction(gram["lambda_min_lower"]))
    gram_condition = _up_fraction(Fraction(gram["condition_number_upper"]))
    gram_floor = _up(1.0 / math.nextafter(math.sqrt(gram_min), -math.inf))
    wls_norm = max(_up_fraction(Fraction(gram["wls_row_norm_upper"])), gram_floor)
    amplification = _up(wls_norm)
    radii = {
        "radius_history": _up(amplification * config["bounds"]["history"]["computed_w_norm_bound"]),
        "radius_dwell": _up(amplification * config["bounds"]["dwell"]["computed_w_norm_bound"]),
        "radius_hrv": _up(amplification * config["bounds"]["hrv"]["computed_w_norm_bound"]),
        "radius_phase_anchor": _up(
            amplification * config["bounds"]["phase_anchor"]["computed_w_norm_bound"]
        ),
        "radius_window_leakage": _up(
            config["bounds"]["window_leakage"]["computed_coefficient_radius"]
        ),
        "radius_sampling": _up(
            math.fsum(
                (
                    amplification * config["bounds"]["sampling"]["computed_w_norm_bound"],
                    exact_numeric,
                )
            )
        ),
        "radius_quadrature": _up(config["bounds"]["sampling"]["quadrature_radius"]),
        "radius_delay": _up(config["bounds"]["delay"]["computed_coefficient_radius"]),
        "radius_filter_initialization": _up(
            config["bounds"]["filter_initialization"]["computed_coefficient_radius"]
        ),
        "radius_measurement_noise": _up(
            amplification * config["bounds"]["measurement_noise"]["computed_w_norm_bound"]
        ),
        "radius_model_residual": _up(
            amplification * config["bounds"]["model_residual"]["computed_w_norm_bound"]
        ),
    }
    total = _up(math.fsum(radii.values()))
    target = execution["result"]["z_tilde"][0]
    payload_hashes = execution["payload_hashes"]
    config["window_result"] = [
        {
            "window_id": request.window_id,
            "rate_id": "rest",
            "subject_pseudonym": "synthetic-subject-001",
            "lead_id": "synthetic-lead-II",
            "physical_gain_id": "synthetic-gain-001",
            "record_manifest_id": "synthetic-record-manifest-001",
            "record_manifest_sha256": request.record_manifest_sha256,
            "estimator_manifest_sha256": request.estimator_manifest_sha256,
            "wls_window_sha256": wls_window_sha256,
            "wls_window_start_sample": 80,
            "wls_window_end_sample_exclusive": 88,
            "wls_execution_sha256": execution["execution_sha256"],
            "wls_sample_index_sha256": payload_hashes["sample_index_sha256"],
            "wls_digital_sample_sha256": payload_hashes["digital_sample_sha256"],
            "wls_time_seconds_sha256": payload_hashes["time_seconds_sha256"],
            "wls_normalized_weights_sha256": payload_hashes["normalized_weights_sha256"],
            "wls_interval_certificate_sha256": interval["certificate_sha256"],
            "wls_exact_sample_functional_coverage": True,
            "wls_exact_functional_radius_upper": exact_numeric,
            "protocol_sha256": protocol_hash,
            "result_sha256": "",
            "left_r_peak_time_s": request.left_r_peak_time_s,
            "right_r_peak_time_s": request.right_r_peak_time_s,
            "complete_rr_intervals": request.complete_rr_intervals,
            "window_duration_s": request.right_r_peak_time_s - request.left_r_peak_time_s,
            "nominal_rate_rad_s": problem.base_frequencies[0],
            "sample_count": len(request.sample_times_s),
            "gram_min_eigenvalue": gram_min,
            "gram_condition_number": gram_condition,
            "gram_bounds_outward_certified": True,
            "gram_certificate_sha256": component.certificate_sha256,
            "wls_row_norm_upper_bound": wls_norm,
            "target_harmonic": request.target_harmonic,
            "nominal_transfer_real": 1.0,
            "nominal_transfer_imag": 0.0,
            "z_tilde_real": target["real"],
            "z_tilde_imag": target["imag"],
            **radii,
            "radius_total_deterministic": total,
            "radius_total_probabilistic": 0.0,
            "joint_failure_probability": 0.0,
            "minimum_burn_in_pass": True,
            "minimum_snr_pass": True,
            "result_status": "PASS_DETERMINISTIC",
            "reason_codes": [],
        }
    ]
    _rehash(config)
    validation = validate_r2_protocol(config)
    if validation.status is not R2ProtocolStatus.PASS_DETERMINISTIC:
        raise RuntimeError(f"synthetic protocol did not pass: {validation.reasons}")
    link = R2DiskCertificateLink(
        request.target_harmonic,
        "rest",
        request,
        component,
        execution,
        tuple(int(value) for value in digital),
        interval,
    )
    return SyntheticCertificateCase(
        problem,
        box,
        baseline.to_json(),
        config,
        (link,),
    )


__all__ = ["SyntheticCertificateCase", "build_synthetic_certificate_case"]
