from __future__ import annotations

import cmath
import hashlib
import json
import math
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction

import numpy as np

from fo_ekf import wls_execution as wls_execution_module
from fo_ekf.certified_set import (
    CertifiedHarmonicData,
    CertifiedJointProblem,
    CertifiedParameterBox,
    ClosedInterval,
    RateOffsetInterval,
)
from fo_ekf.r2_bound_certificate import R2BoundRelation, certify_r2_bounds
from fo_ekf.r2_protocol import (
    R2ProtocolStatus,
    canonical_protocol_sha256,
    canonical_window_result_sha256,
    validate_r2_protocol,
)
from fo_ekf.r2_t1_bundle import (
    R2DiskCertificateLink,
    R2T1BundleRelation,
    certify_r2_t1_bundle,
    replay_r2_t1_bundle,
)
from fo_ekf.tube_certificate import TubeRelation, certify_alpha_tube
from fo_ekf.wls_execution import (
    bounded_sample_indices,
    build_wls_execution_record,
    canonical_hash,
    estimator_identity,
)
from fo_ekf.wls_interval_certificate import build_wls_interval_certificate
from tests.test_r2_bound_certificate import _request
from tests.test_r2_protocol import _rehash, _valid_deterministic_config


def _up(value: float) -> float:
    return 0.0 if value == 0.0 else math.nextafter(value, math.inf)


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


def _fractional_frequency(order: float, frequency: float) -> complex:
    return frequency**order * cmath.exp(0.5j * math.pi * order)


def _problem() -> tuple[CertifiedJointProblem, CertifiedParameterBox]:
    order = 0.7
    damping = 0.5
    frequency = 2.0 * math.pi
    morphology = 1.0 - 0.2j
    response = morphology / (damping + _fractional_frequency(order, frequency))
    harmonic = CertifiedHarmonicData(
        rate_ids=("rest",),
        base_frequencies=(frequency,),
        harmonic_index=1,
        measured_responses=(response,),
        response_radii=(0.05,),
        morphology_drift_radii=(0.0,),
        morphology_bounds=(0.1, 5.0),
        label="m1",
    )
    return (
        CertifiedJointProblem((harmonic,), "rest"),
        CertifiedParameterBox(
            ClosedInterval(0.69, 0.71),
            ClosedInterval(damping, damping),
            (RateOffsetInterval("rest", ClosedInterval(0.0, 0.0)),),
        ),
    )


def _proof(certificate_json: str) -> dict[str, object]:
    return json.loads(certificate_json)["proof"]


def _configure_protocol_without_result(
    request=None,
    *,
    sample_rate_hz: float = 4.0,
) -> dict:
    config = _valid_deterministic_config()
    request = _request() if request is None else request
    data = config["data_boundary"]
    data.update(
        {
            "calibration_split_id": request.calibration_split_id,
            "identification_split_id": request.identification_split_id,
        }
    )
    for bound in config["bounds"].values():
        if "calibration_split_id" in bound:
            bound["calibration_split_id"] = request.calibration_split_id
    config["units"]["tau_star_s"] = request.tau_star_s
    config["model"].update(
        {
            "alpha_interval": list(request.alpha_interval),
            "lambda_interval": list(request.lambda_interval),
            "q_abs_lower_by_harmonic": [0.1, 0.0, 0.1],
            "q_abs_upper_by_harmonic": [5.0, 2.0, 5.0],
            "t1_rate_ids": ["rest"],
            "damping_anchor_rate_id": "rest",
        }
    )
    config["estimator"].update(
        {
            "retained_harmonics": list(request.retained_harmonics),
            "gram_min_eigenvalue_threshold": request.gram_min_eigenvalue_threshold,
            "gram_condition_number_max": request.gram_condition_number_max,
        }
    )
    config["window_design"].update(
        {
            "complete_rr_intervals": request.complete_rr_intervals,
            "dwell_s_min": 14.0,
        }
    )
    config["bounds"]["dwell"]["window_duration_s"] = 2.0
    config["bounds"]["sampling"].update(
        {
            "sample_rate_hz": sample_rate_hz,
            "implementation": request.sampling.implementation,
        }
    )
    config["bounds"]["measurement_noise"]["computed_w_norm_bound"] = 1.0
    config["decision_thresholds"].update(
        {
            "deterministic_total_radius_max": 100.0,
            "deterministic_relative_radius_max": 1000.0,
            "minimum_dwell_s": 14.0,
        }
    )
    for name, group in (
        ("history", request.history),
        ("dwell", request.dwell),
        ("sampling", request.sampling),
        ("delay", request.delay),
    ):
        evidence = group.evidence
        config["bounds"][name].update(
            {
                "source_kind": evidence.source_kind,
                "source_reference": evidence.source_reference,
                "calibration_split_id": evidence.calibration_split_id,
                "independent_of_identification_fit": (evidence.independent_of_identification_fit),
                "uses_target_fit_residuals": evidence.uses_target_fit_residuals,
                "evidence_sha256": evidence.evidence_sha256,
            }
        )
    return config


def _material(
    *, large_index_window: bool = False
) -> tuple[
    CertifiedJointProblem,
    CertifiedParameterBox,
    str,
    dict,
    tuple[R2DiskCertificateLink, ...],
]:
    problem, box = _problem()
    baseline = certify_alpha_tube(problem, box)
    assert baseline.relation is TubeRelation.ROBUST_INNER

    estimator_hash = estimator_identity()["estimator_sha256"]
    assert isinstance(estimator_hash, str)
    base_request = replace(_request(), estimator_manifest_sha256=estimator_hash)
    if large_index_window:
        start_sample = 1_234_567
        end_sample = start_sample + 4097
        sample_indices = tuple(
            int(value) for value in bounded_sample_indices(start_sample, end_sample)
        )
        sampling_frequency = (end_sample - start_sample) / 2.0
        sample_times = tuple(value / sampling_frequency for value in sample_indices)
        count = len(sample_indices)
        base_request = replace(
            base_request,
            left_r_peak_time_s=start_sample / sampling_frequency,
            right_r_peak_time_s=end_sample / sampling_frequency,
            sample_times_s=sample_times,
            raw_weights=(1.0,) * count,
            sample_indices=sample_indices,
            window_start_sample=start_sample,
            window_end_sample_exclusive=end_sample,
            sampling=replace(
                base_request.sampling,
                timestamp_error_s=(1.0e-4,) * count,
                interpolation_error_bound=(2.0e-4,) * count,
                anti_alias_error_bound=(1.0e-4,) * count,
            ),
        )
    else:
        start_sample = 8
        end_sample = 16
        sample_indices = tuple(range(start_sample, end_sample))
        sampling_frequency = 4.0
    config = _configure_protocol_without_result(
        base_request,
        sample_rate_hz=sampling_frequency,
    )
    baseline_document = json.loads(baseline.to_json())
    config["protocol"].update(
        {
            "baseline_t1_certificate_sha256": baseline_document["certificate_sha256"],
            "baseline_t1_input_sha256": baseline_document["input_sha256"],
        }
    )
    provisional_hash = canonical_protocol_sha256(config)
    provisional_request = replace(base_request, protocol_sha256=provisional_hash)
    provisional = certify_r2_bounds(provisional_request)
    assert provisional.relation is R2BoundRelation.CERTIFIED_BOUND
    proof = _proof(provisional.to_json())
    dynamic = proof["dynamic"]
    sampling = proof["sampling"]
    delay = proof["delay"]
    assert isinstance(dynamic, dict) and isinstance(sampling, dict) and isinstance(delay, dict)

    config["bounds"]["history"]["computed_w_norm_bound"] = _up_fraction(
        Fraction(dynamic["history"]["weighted_residual_upper"])
    )
    config["bounds"]["dwell"]["computed_w_norm_bound"] = _up_fraction(
        Fraction(dynamic["dwell"]["weighted_residual_upper"])
    )
    config["bounds"]["sampling"]["computed_w_norm_bound"] = _up_fraction(
        Fraction(sampling["weighted_residual_upper"])
    )
    config["bounds"]["sampling"]["quadrature_radius"] = _up_fraction(
        Fraction(sampling["quadrature_coefficient_radius_upper"])
    )
    config["bounds"]["delay"]["computed_coefficient_radius"] = _up_fraction(
        Fraction(delay["coefficient_radius_upper"])
    )

    _rehash(config)
    protocol_hash = canonical_protocol_sha256(config)
    request = replace(base_request, protocol_sha256=protocol_hash)
    component = certify_r2_bounds(request)
    assert component.relation is R2BoundRelation.CERTIFIED_BOUND
    proof = _proof(component.to_json())
    gram = proof["gram"]
    dynamic = proof["dynamic"]
    sampling = proof["sampling"]
    delay = proof["delay"]
    assert all(isinstance(value, dict) for value in (gram, dynamic, sampling, delay))

    gram_min = _down_fraction(Fraction(gram["lambda_min_lower"]))
    gram_condition = _up_fraction(Fraction(gram["condition_number_upper"]))
    gram_floor = _up(1.0 / math.nextafter(math.sqrt(gram_min), -math.inf))
    wls_norm = max(_up_fraction(Fraction(gram["wls_row_norm_upper"])), gram_floor)
    protocol_amplification = _up(wls_norm)
    sample_index_array = np.asarray(sample_indices, dtype=np.int64)
    phase = (
        2.0
        * np.pi
        * request.complete_rr_intervals
        * (sample_index_array - start_sample)
        / (end_sample - start_sample)
    )
    digital = np.rint(1000.0 * (0.4 * np.cos(phase) - 0.2 * np.sin(phase))).astype(np.int64)
    wls_window_sha256 = hashlib.sha256(b"window-001-wls-window-manifest").hexdigest()
    execution = build_wls_execution_record(
        digital,
        (start_sample, start_sample + (end_sample - start_sample) // 2, end_sample),
        sampling_frequency=sampling_frequency,
        adc_gain=1000.0,
        baseline=0,
        harmonics=(1,),
        protocol_sha256=protocol_hash,
        record_sha256=request.record_manifest_sha256,
        window_sha256=wls_window_sha256,
    )
    interval_certificate = (
        None
        if large_index_window
        else build_wls_interval_certificate(execution, digital, precision_bits=192)
    )
    target_numeric = (
        math.nextafter(0.0, math.inf)
        if interval_certificate is None
        else interval_certificate["response_disks"][0]["absolute_error_upper"]["float_upper"]
    )
    interval_certificate_sha256 = (
        "7" * 64 if interval_certificate is None else interval_certificate["certificate_sha256"]
    )
    radii = {
        "radius_history": _up(
            protocol_amplification * config["bounds"]["history"]["computed_w_norm_bound"]
        ),
        "radius_dwell": _up(
            protocol_amplification * config["bounds"]["dwell"]["computed_w_norm_bound"]
        ),
        "radius_hrv": _up(
            protocol_amplification * config["bounds"]["hrv"]["computed_w_norm_bound"]
        ),
        "radius_phase_anchor": _up(
            protocol_amplification * config["bounds"]["phase_anchor"]["computed_w_norm_bound"]
        ),
        "radius_window_leakage": _up(
            config["bounds"]["window_leakage"]["computed_coefficient_radius"]
        ),
        "radius_sampling": _up(
            math.fsum(
                (
                    protocol_amplification * config["bounds"]["sampling"]["computed_w_norm_bound"],
                    target_numeric,
                )
            )
        ),
        "radius_quadrature": _up(config["bounds"]["sampling"]["quadrature_radius"]),
        "radius_delay": _up(config["bounds"]["delay"]["computed_coefficient_radius"]),
        "radius_filter_initialization": _up(
            config["bounds"]["filter_initialization"]["computed_coefficient_radius"]
        ),
        "radius_measurement_noise": _up(
            protocol_amplification * config["bounds"]["measurement_noise"]["computed_w_norm_bound"]
        ),
        "radius_model_residual": _up(
            protocol_amplification * config["bounds"]["model_residual"]["computed_w_norm_bound"]
        ),
    }
    total = _up(math.fsum(radii.values()))
    target = execution["result"]["z_tilde"][0]
    center = complex(target["real"], target["imag"])
    payload_hashes = execution["payload_hashes"]
    config["window_result"] = [
        {
            "window_id": request.window_id,
            "rate_id": "rest",
            "subject_pseudonym": "subject-001",
            "lead_id": "lead-II",
            "physical_gain_id": "gain-calibration-001",
            "record_manifest_id": "record-manifest-001",
            "record_manifest_sha256": request.record_manifest_sha256,
            "estimator_manifest_sha256": request.estimator_manifest_sha256,
            "wls_window_sha256": wls_window_sha256,
            "wls_window_start_sample": start_sample,
            "wls_window_end_sample_exclusive": end_sample,
            "wls_execution_sha256": execution["execution_sha256"],
            "wls_sample_index_sha256": payload_hashes["sample_index_sha256"],
            "wls_digital_sample_sha256": payload_hashes["digital_sample_sha256"],
            "wls_time_seconds_sha256": payload_hashes["time_seconds_sha256"],
            "wls_normalized_weights_sha256": payload_hashes["normalized_weights_sha256"],
            "wls_interval_certificate_sha256": interval_certificate_sha256,
            "wls_exact_sample_functional_coverage": True,
            "wls_exact_functional_radius_upper": target_numeric,
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
            "target_harmonic": 1,
            "nominal_transfer_real": 1.0,
            "nominal_transfer_imag": 0.0,
            "z_tilde_real": center.real,
            "z_tilde_imag": center.imag,
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
    assert canonical_protocol_sha256(config) == protocol_hash
    validation = validate_r2_protocol(config)
    assert validation.status is R2ProtocolStatus.PASS_DETERMINISTIC, validation.reasons
    link = R2DiskCertificateLink(
        1,
        "rest",
        request,
        component,
        execution,
        tuple(int(value) for value in digital),
        interval_certificate,
    )
    return problem, box, baseline.to_json(), config, (link,)


def test_seals_replays_and_reconstructs_exact_sample_response_family() -> None:
    problem, box, baseline, config, links = _material()

    certificate = certify_r2_t1_bundle(problem, box, config, baseline, links)
    replay = replay_r2_t1_bundle(problem, box, config, baseline, links, certificate)

    assert certificate.relation is R2T1BundleRelation.PRESERVED_ROBUST
    assert len(certificate.response_disks) == 1
    assert certificate.updated_problem is not None
    assert replay.valid
    assert replay.relation is R2T1BundleRelation.PRESERVED_ROBUST
    document = json.loads(certificate.to_json())
    assert document["proof"]["high_level_relation"] == "PRESERVED_ROBUST"
    assert document["proof"]["mathematical_exact_sample_wls_coverage"] is True
    assert document["proof"]["uncovered_numeric_term"] is None
    assert document["input_manifest"]["uncovered_numeric_term"] is None
    assert document["proof"]["certified_margin_transfer_relation"] == "PRESERVED_ROBUST"
    assert len(document["proof"]["certified_response_disks"]) == 1
    assert document["proof"]["embedded_evidence"]["baseline_t1_certificate_json"] == baseline
    assert len(document["proof"]["embedded_evidence"]["r2_component_certificates"]) == 1


def test_missing_disk_certificate_is_not_certifiable_never_outer() -> None:
    problem, box, baseline, config, _ = _material()

    certificate = certify_r2_t1_bundle(problem, box, config, baseline, ())
    replay = replay_r2_t1_bundle(problem, box, config, baseline, (), certificate)

    assert certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert certificate.reason == "component_disk_coverage_mismatch"
    assert "OUTER" not in certificate.to_json()
    assert "REJECT" not in certificate.to_json()
    assert replay.valid
    assert replay.relation is R2T1BundleRelation.NOT_CERTIFIABLE


def test_component_hash_and_record_hash_mismatches_fail_closed() -> None:
    problem, box, baseline, config, links = _material()
    wrong_gram = deepcopy(config)
    wrong_gram["window_result"][0]["gram_certificate_sha256"] = "e" * 64
    wrong_gram["window_result"][0]["result_sha256"] = ""
    wrong_gram["window_result"][0]["result_sha256"] = canonical_window_result_sha256(
        wrong_gram["window_result"][0]
    )
    wrong_record = deepcopy(config)
    wrong_record["window_result"][0]["record_manifest_sha256"] = "f" * 64
    wrong_record["window_result"][0]["result_sha256"] = ""
    wrong_record["window_result"][0]["result_sha256"] = canonical_window_result_sha256(
        wrong_record["window_result"][0]
    )

    gram_certificate = certify_r2_t1_bundle(problem, box, wrong_gram, baseline, links)
    record_certificate = certify_r2_t1_bundle(problem, box, wrong_record, baseline, links)

    assert gram_certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert gram_certificate.reason == "window_gram_certificate_hash_mismatch"
    assert record_certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert record_certificate.reason == "component_record_manifest_hash_mismatch"


def test_protocol_center_tamper_is_caught_by_replayed_wls_before_transfer() -> None:
    problem, box, baseline, config, links = _material()
    changed = deepcopy(config)
    changed["window_result"][0]["z_tilde_real"] += 100.0
    changed["window_result"][0]["result_sha256"] = ""
    changed["window_result"][0]["result_sha256"] = canonical_window_result_sha256(
        changed["window_result"][0]
    )
    assert validate_r2_protocol(changed).status is R2ProtocolStatus.PASS_DETERMINISTIC

    certificate = certify_r2_t1_bundle(problem, box, changed, baseline, links)

    assert certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert certificate.reason == "wls_target_center_mismatch"
    assert "OUTER" not in certificate.to_json()


def test_rehashed_tampered_bundle_proof_does_not_replay() -> None:
    problem, box, baseline, config, links = _material()
    certificate = certify_r2_t1_bundle(problem, box, config, baseline, links)
    document = json.loads(certificate.to_json())
    document["proof"]["certified_response_disks"][0]["new_response_radius"] = "0"
    unsigned = dict(document)
    del unsigned["certificate_sha256"]
    payload = json.dumps(
        unsigned,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode()
    document["certificate_sha256"] = hashlib.sha256(payload).hexdigest()
    tampered = json.dumps(document, sort_keys=True, separators=(",", ":"))

    replay = replay_r2_t1_bundle(problem, box, config, baseline, links, tampered)

    assert not replay.valid
    assert replay.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert replay.reason == "bundle_proof_mismatch"


def test_evidence_semantics_mismatch_is_not_certifiable() -> None:
    problem, box, baseline, config, links = _material()
    changed = deepcopy(config)
    changed["bounds"]["history"]["evidence_sha256"] = "0" * 64
    _rehash(changed)

    certificate = certify_r2_t1_bundle(problem, box, changed, baseline, links)

    assert certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert certificate.reason in {
        "component_protocol_hash_mismatch",
        "component_history_evidence_sha256_mismatch",
    }


def test_frozen_protocol_must_bind_baseline_and_matching_q_prior() -> None:
    problem, box, baseline, config, links = _material()
    wrong_baseline = deepcopy(config)
    wrong_baseline["protocol"]["baseline_t1_input_sha256"] = "0" * 64
    _rehash(wrong_baseline)
    wrong_q = deepcopy(config)
    wrong_q["model"]["q_abs_upper_by_harmonic"][-1] = 4.0
    _rehash(wrong_q)

    baseline_certificate = certify_r2_t1_bundle(problem, box, wrong_baseline, baseline, links)
    q_certificate = certify_r2_t1_bundle(problem, box, wrong_q, baseline, links)

    assert baseline_certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert baseline_certificate.reason == "baseline_t1_input_not_bound_by_frozen_protocol"
    assert q_certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert q_certificate.reason == "baseline_t1_q_upper_mismatch"


def test_duplicate_bundle_keys_are_rejected_before_recomputation() -> None:
    problem, box, baseline, config, links = _material()
    certificate = certify_r2_t1_bundle(problem, box, config, baseline, links)
    duplicate = certificate.to_json()[:-1] + f',"schema":"{certificate.relation.value}"}}'

    replay = replay_r2_t1_bundle(problem, box, config, baseline, links, duplicate)

    assert not replay.valid
    assert replay.reason == "invalid_bundle_envelope"


def test_unlinked_protocol_result_cannot_cross_replay_or_issue() -> None:
    problem, box, baseline, config, links = _material()
    certificate = certify_r2_t1_bundle(problem, box, config, baseline, links)
    changed = deepcopy(config)
    extra = deepcopy(changed["window_result"][0])
    extra["window_id"] = "unlinked-window"
    extra["target_harmonic"] = -1
    extra["z_tilde_real"] += 0.01
    changed["window_result"].append(extra)
    _rehash(changed)
    assert validate_r2_protocol(changed).status is R2ProtocolStatus.PASS_DETERMINISTIC

    changed_certificate = certify_r2_t1_bundle(problem, box, changed, baseline, links)
    replay = replay_r2_t1_bundle(
        problem,
        box,
        changed,
        baseline,
        links,
        certificate,
    )

    assert changed_certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert changed_certificate.reason == "protocol_result_set_not_exactly_bound"
    assert not replay.valid
    assert replay.reason == "bundle_proof_mismatch"


def test_deep_protocol_and_hostile_repr_fail_closed_without_recursion() -> None:
    problem, box, baseline, config, links = _material()
    valid_certificate = certify_r2_t1_bundle(problem, box, config, baseline, links)
    deep = deepcopy(config)
    cursor = deep
    for index in range(100):
        child: dict = {}
        cursor[f"nested-{index}"] = child
        cursor = child

    deep_certificate = certify_r2_t1_bundle(problem, box, deep, baseline, links)
    deep_replay = replay_r2_t1_bundle(
        problem,
        box,
        deep,
        baseline,
        links,
        valid_certificate,
    )

    class HostileRepr:
        def __repr__(self) -> str:
            raise RuntimeError("repr must never be called")

    hostile_certificate = certify_r2_t1_bundle(
        problem,
        box,
        config,
        baseline,
        (HostileRepr(),),  # type: ignore[arg-type]
    )

    assert deep_certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert deep_certificate.reason == "protocol_resource_depth"
    assert not deep_replay.valid
    assert hostile_certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert hostile_certificate.reason == "component_link_type_invalid"


def test_legacy_link_without_wls_evidence_is_explicitly_not_certifiable() -> None:
    problem, box, baseline, config, links = _material()
    current = links[0]
    legacy = R2DiskCertificateLink(
        current.harmonic_index,
        current.rate_id,
        current.request,
        current.certificate_json,
    )

    certificate = certify_r2_t1_bundle(problem, box, config, baseline, (legacy,))

    assert certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert certificate.reason == "wls_execution_record_missing"


def test_external_digital_payload_and_resigned_record_identity_tamper_fail_closed() -> None:
    problem, box, baseline, config, links = _material()
    current = links[0]
    changed_samples = list(current.digital_samples or ())
    changed_samples[0] += 1
    payload_tamper = replace(current, digital_samples=tuple(changed_samples))

    record = deepcopy(current.wls_execution_record)
    assert isinstance(record, dict)
    record["binding"]["record_sha256"] = "e" * 64
    unsigned = dict(record)
    unsigned.pop("execution_sha256")
    record["execution_sha256"] = canonical_hash(unsigned)
    identity_tamper = replace(current, wls_execution_record=record)

    payload_certificate = certify_r2_t1_bundle(
        problem,
        box,
        config,
        baseline,
        (payload_tamper,),
    )
    identity_certificate = certify_r2_t1_bundle(
        problem,
        box,
        config,
        baseline,
        (identity_tamper,),
    )

    assert payload_certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert payload_certificate.reason == "wls_execution_record_not_replayable"
    assert identity_certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert identity_certificate.reason == "wls_record_manifest_hash_mismatch"


def test_execution_hash_and_exact_functional_floor_cannot_be_understated() -> None:
    problem, box, baseline, config, links = _material()
    wrong_hash = deepcopy(config)
    wrong_hash["window_result"][0]["wls_execution_sha256"] = "f" * 64
    wrong_hash["window_result"][0]["result_sha256"] = ""
    wrong_hash["window_result"][0]["result_sha256"] = canonical_window_result_sha256(
        wrong_hash["window_result"][0]
    )

    component_document = json.loads(links[0].json_text())
    component_floor = Fraction(component_document["proof"]["sampling"]["coefficient_radius_upper"])
    interval_certificate = links[0].wls_interval_certificate
    assert isinstance(interval_certificate, dict)
    target_bound = interval_certificate["response_disks"][0]["absolute_error_upper"]["dyadic"]
    exact_functional_floor = Fraction(
        int(target_bound["numerator"]),
        int(target_bound["denominator"]),
    )
    understated = deepcopy(config)
    understated["window_result"][0]["radius_sampling"] = _down_fraction(
        component_floor + exact_functional_floor
    )
    radius_fields = [
        key
        for key in understated["window_result"][0]
        if key.startswith("radius_")
        and key not in {"radius_total_deterministic", "radius_total_probabilistic"}
    ]
    understated["window_result"][0]["radius_total_deterministic"] = _up(
        math.fsum(understated["window_result"][0][key] for key in radius_fields)
    )
    understated["window_result"][0]["result_sha256"] = ""
    understated["window_result"][0]["result_sha256"] = canonical_window_result_sha256(
        understated["window_result"][0]
    )
    assert validate_r2_protocol(understated).status is R2ProtocolStatus.PASS_DETERMINISTIC

    hash_certificate = certify_r2_t1_bundle(problem, box, wrong_hash, baseline, links)
    floor_certificate = certify_r2_t1_bundle(problem, box, understated, baseline, links)

    assert hash_certificate.reason == "window_wls_execution_hash_mismatch"
    assert floor_certificate.reason == "window_radius_sampling_below_component_certificate"


def test_missing_unknown_and_tampered_interval_certificates_fail_closed(monkeypatch) -> None:
    problem, box, baseline, config, links = _material()
    current = links[0]
    missing = replace(current, wls_interval_certificate=None)
    missing_certificate = certify_r2_t1_bundle(problem, box, config, baseline, (missing,))

    tampered_document = deepcopy(current.wls_interval_certificate)
    assert isinstance(tampered_document, dict)
    tampered_document["response_disks"][0]["absolute_error_upper"]["float_upper"] *= 2.0
    tampered = replace(current, wls_interval_certificate=tampered_document)
    tampered_certificate = certify_r2_t1_bundle(problem, box, config, baseline, (tampered,))

    monkeypatch.setattr(
        "fo_ekf.r2_t1_bundle.replay_wls_interval_certificate",
        lambda *_: {"status": "MATCH", "relation": "UNKNOWN"},
    )
    unknown_certificate = certify_r2_t1_bundle(problem, box, config, baseline, links)

    assert missing_certificate.reason == "wls_interval_certificate_missing"
    assert tampered_certificate.reason == "wls_interval_certificate_not_replayable"
    assert unknown_certificate.reason == "wls_interval_certificate_not_certified"
    assert all(
        certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
        for certificate in (missing_certificate, tampered_certificate, unknown_certificate)
    )


def test_wls_record_resource_depth_is_bounded_before_replay() -> None:
    problem, box, baseline, config, links = _material()
    current = links[0]
    record = deepcopy(current.wls_execution_record)
    assert isinstance(record, dict)
    cursor = record
    for _ in range(70):
        cursor["nested"] = {}
        cursor = cursor["nested"]
    malformed = replace(current, wls_execution_record=record)

    certificate = certify_r2_t1_bundle(problem, box, config, baseline, (malformed,))

    assert certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert certificate.reason == "wls_execution_resource_depth"


def test_resigned_component_with_different_exact_phase_operator_fails_closed(monkeypatch) -> None:
    monkeypatch.setattr(
        wls_execution_module,
        "_joint_solve_rounding_enclosure",
        lambda *_args, **_kwargs: {
            "status": "UNAVAILABLE",
            "scope": "test_preflight_only",
            "reason": "not needed because exact-index mismatch is rejected before replay",
        },
    )
    problem, box, baseline, config, links = _material(large_index_window=True)
    current = links[0]
    shifted_indices = tuple(value + 1 for value in current.request.sample_indices or ())
    shifted_request = replace(current.request, sample_indices=shifted_indices)
    shifted_component = certify_r2_bounds(shifted_request)
    assert shifted_component.relation is R2BoundRelation.CERTIFIED_BOUND
    shifted_link = replace(
        current,
        request=shifted_request,
        certificate_json=shifted_component,
    )

    start_request = replace(
        current.request,
        window_start_sample=(current.request.window_start_sample or 0) - 1,
    )
    start_component = certify_r2_bounds(start_request)
    assert start_component.relation is R2BoundRelation.CERTIFIED_BOUND
    start_link = replace(current, request=start_request, certificate_json=start_component)

    shifted_result = certify_r2_t1_bundle(
        problem,
        box,
        config,
        baseline,
        (shifted_link,),
    )
    start_result = certify_r2_t1_bundle(
        problem,
        box,
        config,
        baseline,
        (start_link,),
    )

    assert shifted_result.reason == "wls_r2_exact_sample_indices_mismatch"
    assert start_result.reason == "wls_r2_integer_window_geometry_mismatch"
    assert shifted_result.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert start_result.relation is R2T1BundleRelation.NOT_CERTIFIABLE


def test_protocol_integer_window_geometry_must_match_replayed_execution() -> None:
    problem, box, baseline, config, links = _material()
    changed = deepcopy(config)
    changed["window_result"][0]["wls_window_start_sample"] -= 1
    changed["window_result"][0]["result_sha256"] = ""
    changed["window_result"][0]["result_sha256"] = canonical_window_result_sha256(
        changed["window_result"][0]
    )
    assert validate_r2_protocol(changed).status is R2ProtocolStatus.PASS_DETERMINISTIC

    certificate = certify_r2_t1_bundle(problem, box, changed, baseline, links)

    assert certificate.relation is R2T1BundleRelation.NOT_CERTIFIABLE
    assert certificate.reason == "wls_r2_integer_window_geometry_mismatch"
