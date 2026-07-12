from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from fractions import Fraction

from flint import ctx

from fo_ekf.r2_bound_certificate import (
    BOUND_NUMERIC_SCOPE,
    DelayBoundInput,
    DwellBoundInput,
    EvidenceReference,
    HistoryBoundInput,
    R2BoundRelation,
    R2BoundRequest,
    SamplingBoundInput,
    certify_r2_bounds,
    replay_r2_bound_certificate,
)


def _evidence(label: str) -> EvidenceReference:
    return EvidenceReference(
        source_kind="independent_bounded_calibration",
        source_reference=f"calibration/{label}.json",
        evidence_sha256=hashlib.sha256(label.encode()).hexdigest(),
        calibration_split_id="calibration-v1",
        applicability_domain="all frozen synthetic windows",
        signal_unit="mV",
        time_unit="s",
        independent_of_identification_fit=True,
        uses_target_fit_residuals=False,
        numeric_scope=BOUND_NUMERIC_SCOPE,
    )


def _request(
    *,
    trapezoidal: bool = False,
    front_end_mode: str = "identity",
) -> R2BoundRequest:
    if trapezoidal:
        sample_times = (2.0, 2.125, 2.375, 2.75, 3.0, 3.25, 3.625, 3.875, 4.0)
        second_derivatives = (3.0,) * 8
        implementation = "trapezoidal_continuous_lockin"
    else:
        sample_times = tuple(2.0 + 0.25 * index for index in range(8))
        second_derivatives = ()
        implementation = "discrete_wls"
    count = len(sample_times)
    if trapezoidal:
        panels = tuple(
            right - left for left, right in zip(sample_times[:-1], sample_times[1:], strict=True)
        )
        raw_weights = (
            (panels[0],)
            + tuple(panels[index - 1] + panels[index] for index in range(1, count - 1))
            + (panels[-1],)
        )
    else:
        raw_weights = (1.0,) * count
    if front_end_mode == "identity":
        history_output = None
        dwell_output = None
        amplitude_error = 0.0
        group_delay_error = 0.0
        phase_error = 0.0
        transfer_lower = 1.0
    else:
        history_output = (0.03,) * count
        dwell_output = (0.02,) * count
        amplitude_error = 0.01
        group_delay_error = 0.002
        phase_error = 0.003
        transfer_lower = 0.9
    return R2BoundRequest(
        window_id="window-001",
        protocol_sha256="a" * 64,
        record_manifest_sha256="b" * 64,
        estimator_manifest_sha256="c" * 64,
        calibration_split_id="calibration-v1",
        identification_split_id="identification-v1",
        frozen_before_target_ecg_access=True,
        signal_unit="mV",
        dynamic_front_end_mode=front_end_mode,
        tau_star_s=1.0,
        alpha_interval=(0.6, 0.9),
        lambda_interval=(0.4, 1.2),
        left_r_peak_time_s=2.0,
        right_r_peak_time_s=4.0,
        complete_rr_intervals=2,
        sample_times_s=sample_times,
        raw_weights=raw_weights,
        retained_harmonics=(-1, 0, 1),
        target_harmonic=1,
        gram_min_eigenvalue_threshold=0.5,
        gram_condition_number_max=2.0,
        history=HistoryBoundInput(
            record_origin_s=0.0,
            rate_segment_start_s=0.0,
            initial_state_abs_bound=0.5,
            new_input_sup_bound=0.2,
            pre_record_history_coefficient=0.1,
            evidence=_evidence("history"),
            preconvolved_output_envelope=history_output,
        ),
        dwell=DwellBoundInput(
            previous_input_sup_bound=0.3,
            new_input_sup_bound=0.2,
            evidence=_evidence("dwell"),
            preconvolved_output_envelope=dwell_output,
        ),
        sampling=SamplingBoundInput(
            timestamp_error_s=(1.0e-4,) * count,
            interpolation_error_bound=(2.0e-4,) * count,
            anti_alias_error_bound=(1.0e-4,) * count,
            adc_quantization_step=1.0e-3,
            first_derivative_bound_signal_per_s=2.0,
            implementation=implementation,
            panel_second_derivative_bound=second_derivatives,
            evidence=_evidence("sampling"),
        ),
        delay=DelayBoundInput(
            relative_amplitude_error_bound=amplitude_error,
            relative_group_delay_error_s=group_delay_error,
            residual_phase_calibration_error_rad=phase_error,
            target_harmonic_abs_upper_bound=2.0,
            nominal_transfer_magnitude_lower_bound=transfer_lower,
            ecg_to_rpeak_channel_alignment_included=True,
            evidence=_evidence("delay"),
        ),
    )


def _proof(certificate_json: str) -> dict[str, object]:
    return json.loads(certificate_json)["proof"]


def _fraction_at(mapping: dict[str, object], key: str) -> Fraction:
    return Fraction(str(mapping[key]))


def test_certifies_all_six_requested_components_and_replays() -> None:
    request = _request()
    certificate = certify_r2_bounds(request)

    assert certificate.relation is R2BoundRelation.CERTIFIED_BOUND
    assert certificate.precision_bits >= 128
    document = json.loads(certificate.to_json())
    proof = document["proof"]
    assert Fraction(proof["gram"]["lambda_min_lower"]) >= Fraction(1, 2)
    assert Fraction(proof["gram"]["condition_number_upper"]) <= 2
    assert Fraction(proof["dynamic"]["history"]["coefficient_radius_upper"]) > 0
    assert Fraction(proof["dynamic"]["dwell"]["coefficient_radius_upper"]) > 0
    assert Fraction(proof["sampling"]["coefficient_radius_upper"]) > 0
    assert Fraction(proof["sampling"]["quadrature_coefficient_radius_upper"]) == 0
    assert Fraction(proof["delay"]["coefficient_radius_upper"]) == 0
    assert proof["dynamic"]["dynamic_front_end_mode"] == "identity"
    assert Fraction(proof["partial_radius_upper"]) > 0

    replay = replay_r2_bound_certificate(request, certificate)
    assert replay.valid
    assert replay.relation is R2BoundRelation.CERTIFIED_BOUND
    assert replay.reason == "replay_verified"


def test_missing_primitive_group_is_sealed_unknown() -> None:
    request = replace(_request(), history=None)
    certificate = certify_r2_bounds(request)

    assert certificate.relation is R2BoundRelation.UNKNOWN
    assert certificate.reason == "history_bound_missing"
    replay = replay_r2_bound_certificate(request, certificate)
    assert replay.valid
    assert replay.relation is R2BoundRelation.UNKNOWN


def test_singular_aliasing_gram_is_unknown_not_an_exclusion() -> None:
    request = replace(
        _request(),
        retained_harmonics=(0, 4),
        target_harmonic=4,
    )
    certificate = certify_r2_bounds(request)

    assert certificate.relation is R2BoundRelation.UNKNOWN
    assert certificate.reason == "gram_not_positive_by_gershgorin"
    assert "REJECT" not in certificate.to_json()


def test_gram_threshold_and_condition_are_fail_closed() -> None:
    request = _request(trapezoidal=True)
    too_strict_minimum = certify_r2_bounds(replace(request, gram_min_eigenvalue_threshold=1.01))
    too_strict_condition = certify_r2_bounds(replace(_request(), gram_condition_number_max=1.0))

    assert too_strict_minimum.relation is R2BoundRelation.UNKNOWN
    assert too_strict_minimum.reason == "gram_min_below_threshold"
    assert too_strict_condition.relation is R2BoundRelation.UNKNOWN
    assert too_strict_condition.reason == "gram_condition_above_threshold"


def test_trapezoidal_quadrature_has_positive_replayable_radius() -> None:
    request = _request(trapezoidal=True)
    certificate = certify_r2_bounds(request)

    assert certificate.relation is R2BoundRelation.CERTIFIED_BOUND
    sampling = _proof(certificate.to_json())["sampling"]
    gram = _proof(certificate.to_json())["gram"]
    assert isinstance(sampling, dict)
    assert isinstance(gram, dict)
    assert gram["method"] == "continuous_full_cycle_orthogonality_v1"
    assert gram["lambda_min_lower"] == "1"
    assert gram["condition_number_upper"] == "1"
    assert gram["wls_row_norm_upper"] == "1"
    assert Fraction(sampling["quadrature_coefficient_radius_upper"]) > 0
    assert replay_r2_bound_certificate(request, certificate).valid


def test_later_window_reduces_conservative_history_and_dwell_envelopes() -> None:
    early = _request()
    shift = 10.0
    late = replace(
        early,
        left_r_peak_time_s=early.left_r_peak_time_s + shift,
        right_r_peak_time_s=early.right_r_peak_time_s + shift,
        sample_times_s=tuple(value + shift for value in early.sample_times_s),
    )
    early_certificate = certify_r2_bounds(early)
    late_certificate = certify_r2_bounds(late)

    assert early_certificate.relation is R2BoundRelation.CERTIFIED_BOUND
    assert late_certificate.relation is R2BoundRelation.CERTIFIED_BOUND
    early_dynamic = _proof(early_certificate.to_json())["dynamic"]
    late_dynamic = _proof(late_certificate.to_json())["dynamic"]
    assert isinstance(early_dynamic, dict) and isinstance(late_dynamic, dict)
    for name in ("history", "dwell"):
        early_bound = Fraction(early_dynamic[name]["weighted_residual_upper"])
        late_bound = Fraction(late_dynamic[name]["weighted_residual_upper"])
        assert late_bound < early_bound


def test_tampered_proof_fails_even_if_attacker_rehashes_document() -> None:
    request = _request()
    certificate = certify_r2_bounds(request)
    document = json.loads(certificate.to_json())
    document["proof"]["partial_radius_upper"] = "0"
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

    replay = replay_r2_bound_certificate(request, tampered)
    assert not replay.valid
    assert replay.relation is R2BoundRelation.UNKNOWN
    assert replay.reason == "proof_mismatch"


def test_certificate_is_bound_to_request_and_evidence_hashes() -> None:
    request = _request()
    certificate = certify_r2_bounds(request)
    changed_delay = replace(
        request.delay,
        evidence=replace(request.delay.evidence, evidence_sha256="f" * 64),
    )
    changed = replace(request, delay=changed_delay)

    replay = replay_r2_bound_certificate(changed, certificate)
    assert not replay.valid
    assert replay.relation is R2BoundRelation.UNKNOWN
    assert replay.reason == "invalid_certificate_envelope"


def test_target_fit_provenance_and_channel_alignment_fail_closed() -> None:
    request = _request()
    leaked_history = replace(
        request.history,
        evidence=replace(request.history.evidence, uses_target_fit_residuals=True),
    )
    missing_alignment = replace(
        request.delay,
        ecg_to_rpeak_channel_alignment_included=False,
    )

    leaked = certify_r2_bounds(replace(request, history=leaked_history))
    unaligned = certify_r2_bounds(replace(request, delay=missing_alignment))
    assert leaked.relation is R2BoundRelation.UNKNOWN
    assert leaked.reason == "history_uses_target_fit_residuals"
    assert unaligned.relation is R2BoundRelation.UNKNOWN
    assert unaligned.reason == "delay_channel_alignment_missing"


def test_exact_weight_normalization_is_recorded() -> None:
    request = replace(_request(), raw_weights=(0.1,) * 8)
    certificate = certify_r2_bounds(request)

    assert certificate.relation is R2BoundRelation.CERTIFIED_BOUND
    weights = _proof(certificate.to_json())["normalized_weights"]
    assert sum((Fraction(value) for value in weights), Fraction(0)) == 1


def test_issuance_rejects_precision_above_replay_resource_limit() -> None:
    certificate = certify_r2_bounds(
        _request(),
        precision_schedule=(128, 8192),
    )

    assert certificate.relation is R2BoundRelation.UNKNOWN
    assert certificate.reason == "precision_resource_limit"
    assert certificate.attempted_precisions == ()


def test_replay_rejects_duplicate_json_keys() -> None:
    request = _request()
    certificate = certify_r2_bounds(request)
    duplicate = certificate.to_json()[:-1] + ',"schema":"fo-ekf.r2-bound-engine.v2"}'

    replay = replay_r2_bound_certificate(request, duplicate)
    assert not replay.valid
    assert replay.relation is R2BoundRelation.UNKNOWN
    assert replay.reason == "invalid_certificate_envelope"


def test_replay_rejects_deeply_nested_json_without_crashing() -> None:
    request = _request()
    deeply_nested = "[" * 1100 + "0" + "]" * 1100

    replay = replay_r2_bound_certificate(request, deeply_nested)

    assert not replay.valid
    assert replay.relation is R2BoundRelation.UNKNOWN
    assert replay.reason == "invalid_certificate_envelope"


def test_replay_binds_underlying_flint_version() -> None:
    request = _request()
    certificate = certify_r2_bounds(request)
    document = json.loads(certificate.to_json())
    document["environment"]["flint"] = "0.0.0"
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

    replay = replay_r2_bound_certificate(request, tampered)
    assert not replay.valid
    assert replay.relation is R2BoundRelation.UNKNOWN
    assert replay.reason == "invalid_certificate_envelope"


def test_missing_or_insufficient_evidence_is_unknown() -> None:
    request = _request()
    missing = replace(request.history, evidence=None)  # type: ignore[arg-type]
    geometry_only = replace(
        request.history,
        evidence=replace(request.history.evidence, source_kind="observable_window_geometry"),
    )

    missing_certificate = certify_r2_bounds(replace(request, history=missing))
    geometry_certificate = certify_r2_bounds(replace(request, history=geometry_only))
    assert missing_certificate.relation is R2BoundRelation.UNKNOWN
    assert missing_certificate.reason == "history_evidence_missing"
    assert geometry_certificate.relation is R2BoundRelation.UNKNOWN
    assert geometry_certificate.reason == "history_source_kind_not_allowed"


def test_global_arb_precision_is_restored_after_issuance_and_replay() -> None:
    request = _request()
    initial_precision = ctx.prec

    certificate = certify_r2_bounds(request, precision_schedule=(192,))
    replay = replay_r2_bound_certificate(request, certificate)

    assert certificate.relation is R2BoundRelation.CERTIFIED_BOUND
    assert replay.valid
    assert ctx.prec == initial_precision


def test_identity_front_end_contract_rejects_hidden_nonidentity_terms() -> None:
    request = _request()
    nonzero_error = replace(request.delay, relative_amplitude_error_bound=1.0e-3)
    nonunit_transfer = replace(request.delay, nominal_transfer_magnitude_lower_bound=0.99)
    external_envelope = replace(
        request.history,
        preconvolved_output_envelope=(0.1,) * len(request.sample_times_s),
    )

    error_certificate = certify_r2_bounds(replace(request, delay=nonzero_error))
    transfer_certificate = certify_r2_bounds(replace(request, delay=nonunit_transfer))
    envelope_certificate = certify_r2_bounds(replace(request, history=external_envelope))
    assert error_certificate.relation is R2BoundRelation.UNKNOWN
    assert error_certificate.reason == "identity_mode_requires_zero_front_end_error"
    assert transfer_certificate.relation is R2BoundRelation.UNKNOWN
    assert transfer_certificate.reason == "identity_mode_requires_unit_transfer"
    assert envelope_certificate.relation is R2BoundRelation.UNKNOWN
    assert envelope_certificate.reason == "identity_mode_forbids_preconvolved_output_envelopes"


def test_preconvolved_mode_uses_independent_output_envelopes_not_ml_diagnostics() -> None:
    request = _request(front_end_mode="preconvolved_output_envelopes")
    certificate = certify_r2_bounds(request)

    assert certificate.relation is R2BoundRelation.CERTIFIED_BOUND
    proof = _proof(certificate.to_json())
    dynamic = proof["dynamic"]
    assert isinstance(dynamic, dict)
    assert dynamic["dynamic_front_end_mode"] == "preconvolved_output_envelopes"
    assert dynamic["output_bound_source"] == "independent_preconvolved_output_envelope"
    history = dynamic["history"]
    dwell = dynamic["dwell"]
    assert history["front_end_output_envelope_by_sample"] == [str(Fraction.from_float(0.03))] * len(
        request.sample_times_s
    )
    assert dwell["front_end_output_envelope_by_sample"] == [str(Fraction.from_float(0.02))] * len(
        request.sample_times_s
    )
    assert (
        history["input_envelope_diagnostic_by_sample"]
        != history["front_end_output_envelope_by_sample"]
    )
    assert Fraction(history["weighted_residual_upper"]) < Fraction(31, 1000)
    assert Fraction(dwell["weighted_residual_upper"]) < Fraction(21, 1000)
    assert Fraction(proof["delay"]["coefficient_radius_upper"]) > 0
    assert replay_r2_bound_certificate(request, certificate).valid


def test_preconvolved_mode_requires_complete_per_sample_output_envelopes() -> None:
    request = _request(front_end_mode="preconvolved_output_envelopes")
    missing = replace(request.history, preconvolved_output_envelope=None)
    short = replace(request.dwell, preconvolved_output_envelope=(0.1,))

    missing_certificate = certify_r2_bounds(replace(request, history=missing))
    short_certificate = certify_r2_bounds(replace(request, dwell=short))
    assert missing_certificate.relation is R2BoundRelation.UNKNOWN
    assert missing_certificate.reason == "history_preconvolved_output_envelope_invalid"
    assert short_certificate.relation is R2BoundRelation.UNKNOWN
    assert short_certificate.reason == "dwell_preconvolved_output_envelope_invalid"


def test_continuous_lockin_requires_exact_composite_trapezoid_weights() -> None:
    request = _request(trapezoidal=True)
    wrong = replace(request, raw_weights=(1.0,) * len(request.sample_times_s))

    certificate = certify_r2_bounds(wrong)
    assert certificate.relation is R2BoundRelation.UNKNOWN
    assert certificate.reason == "trapezoidal_weights_must_match_composite_rule"


def test_continuous_lockin_does_not_reuse_an_aliased_discrete_gram() -> None:
    request = replace(
        _request(trapezoidal=True),
        retained_harmonics=(0, 8),
        target_harmonic=8,
    )
    certificate = certify_r2_bounds(request)

    assert certificate.relation is R2BoundRelation.CERTIFIED_BOUND
    gram = _proof(certificate.to_json())["gram"]
    assert gram["method"] == "continuous_full_cycle_orthogonality_v1"
    assert gram["wls_row_norm_upper"] == "1"


def test_typed_evidence_units_and_group_source_allowlists_are_enforced() -> None:
    request = _request()
    wrong_signal = replace(
        request.history,
        evidence=replace(request.history.evidence, signal_unit="V"),
    )
    wrong_time = replace(
        request.sampling,
        evidence=replace(request.sampling.evidence, time_unit="ms"),
    )
    wrong_source = replace(
        request.history,
        evidence=replace(request.history.evidence, source_kind="sensor_spec"),
    )
    allowed_sampling_source = replace(
        request.sampling,
        evidence=replace(request.sampling.evidence, source_kind="sensor_spec"),
    )

    signal_certificate = certify_r2_bounds(replace(request, history=wrong_signal))
    time_certificate = certify_r2_bounds(replace(request, sampling=wrong_time))
    source_certificate = certify_r2_bounds(replace(request, history=wrong_source))
    allowed_certificate = certify_r2_bounds(replace(request, sampling=allowed_sampling_source))
    assert signal_certificate.reason == "history_signal_unit_mismatch"
    assert time_certificate.reason == "sampling_time_unit_must_be_s"
    assert source_certificate.reason == "history_source_kind_not_allowed"
    assert allowed_certificate.relation is R2BoundRelation.CERTIFIED_BOUND


def test_non_request_top_level_is_a_replayable_sealed_unknown() -> None:
    malformed = {"window_id": "not-a-typed-request"}
    certificate = certify_r2_bounds(malformed)

    assert certificate.relation is R2BoundRelation.UNKNOWN
    assert certificate.reason == "request_type_invalid"
    replay = replay_r2_bound_certificate(malformed, certificate)  # type: ignore[arg-type]
    assert replay.valid
    assert replay.relation is R2BoundRelation.UNKNOWN
