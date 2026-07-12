import cmath
import hashlib
import json
import math
from dataclasses import replace
from fractions import Fraction

import pytest
from flint import acb, arb, ctx

from fo_ekf import margin_transfer as margin_transfer_module
from fo_ekf.certified_set import (
    CertifiedHarmonicData,
    CertifiedJointProblem,
    CertifiedParameterBox,
    ClosedInterval,
    RateOffsetInterval,
)
from fo_ekf.margin_transfer import (
    DiskMarginPerturbation,
    MarginTransferRelation,
    certify_margin_transfer,
    replay_margin_transfer_certificate,
)
from fo_ekf.tube_certificate import TubeRelation, certify_alpha_tube


def _fractional_frequency(order: float, frequency: float) -> complex:
    return frequency**order * cmath.exp(0.5j * math.pi * order)


def _problem(
    *,
    response_radius: float = 0.05,
    morphology_radius: float = 0.0,
) -> tuple[CertifiedJointProblem, CertifiedParameterBox]:
    order = 0.7
    damping = 0.5
    frequency = 1.2
    morphology = 1.0 - 0.2j
    response = morphology / (damping + _fractional_frequency(order, frequency))
    harmonic = CertifiedHarmonicData(
        rate_ids=("rest",),
        base_frequencies=(frequency,),
        harmonic_index=1,
        measured_responses=(response,),
        response_radii=(response_radius,),
        morphology_drift_radii=(morphology_radius,),
        morphology_bounds=(0.1, 5.0),
        label="m1",
    )
    problem = CertifiedJointProblem((harmonic,), "rest")
    box = CertifiedParameterBox(
        ClosedInterval(0.69, 0.71),
        ClosedInterval(damping, damping),
        (RateOffsetInterval("rest", ClosedInterval(0.0, 0.0)),),
    )
    return problem, box


def _free_damping_problem() -> tuple[CertifiedJointProblem, CertifiedParameterBox]:
    rate_ids = ("rest", "mid", "high")
    frequencies = (0.8, 1.25, 1.9)
    true_order = 0.72
    damping = 0.5
    offsets = (0.0, 0.03, -0.02)
    morphology = 1.0 - 0.25j
    responses = tuple(
        morphology / (damping + offset + _fractional_frequency(true_order, frequency))
        for frequency, offset in zip(frequencies, offsets, strict=True)
    )
    harmonic = CertifiedHarmonicData(
        rate_ids=rate_ids,
        base_frequencies=frequencies,
        harmonic_index=1,
        measured_responses=responses,
        response_radii=(0.05,) * 3,
        morphology_drift_radii=(0.0,) * 3,
        morphology_bounds=(0.5, 1.5),
        label="m1",
    )
    problem = CertifiedJointProblem((harmonic,), "rest")
    box = CertifiedParameterBox(
        ClosedInterval(0.7, 0.74),
        ClosedInterval(0.45, 0.55),
        tuple(
            RateOffsetInterval(rate_id, ClosedInterval(offset, offset))
            for rate_id, offset in zip(rate_ids, offsets, strict=True)
        ),
    )
    return problem, box


def _perturbation(
    problem: CertifiedJointProblem,
    *,
    center_shift: float = 0.0,
    new_response_radius: float | None = None,
    new_morphology_radius: float | None = None,
) -> tuple[DiskMarginPerturbation, ...]:
    harmonic = problem.harmonics[0]
    old_response = harmonic.response_radii[0]
    old_morphology = harmonic.morphology_drift_radii[0]
    return (
        DiskMarginPerturbation(
            harmonic_index=harmonic.harmonic_index,
            rate_id=problem.rate_ids[0],
            center_shift_abs_upper=center_shift,
            old_response_radius=old_response,
            new_response_radius=(
                old_response if new_response_radius is None else new_response_radius
            ),
            old_morphology_radius=old_morphology,
            new_morphology_radius=(
                old_morphology if new_morphology_radius is None else new_morphology_radius
            ),
        ),
    )


def _canonical_json(value: object) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _rat(value: str) -> Fraction:
    return Fraction(value)


def _arb_fraction(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def _denominator_bounds_from_stored_damping(
    damping_image: list[str],
    box: CertifiedParameterBox,
    base_frequency: float,
    harmonic_index: int,
    precision: int,
) -> tuple[float, float]:
    with ctx.workprec(precision):
        damping = _arb_fraction(Fraction(damping_image[0])).union(
            _arb_fraction(Fraction(damping_image[1]))
        )
        order = _arb_fraction(Fraction.from_float(box.order.lower)).union(
            _arb_fraction(Fraction.from_float(box.order.upper))
        )
        frequency = _arb_fraction(Fraction.from_float(base_frequency)) * harmonic_index
        amplitude = (order * frequency.log()).exp()
        angle = order * arb.pi() / 2
        denominator = acb(damping) + acb(amplitude * angle.cos(), amplitude * angle.sin())
        return float(denominator.abs_lower()), float(denominator.abs_upper())


def test_unchanged_disks_preserve_robust_certificate_and_replay() -> None:
    problem, box = _problem()
    baseline = certify_alpha_tube(problem, box)
    perturbations = _perturbation(problem)

    transfer = certify_margin_transfer(problem, box, baseline.to_json(), perturbations)
    replay = replay_margin_transfer_certificate(
        problem,
        box,
        baseline.to_json(),
        perturbations,
        transfer.to_json(),
    )

    assert baseline.relation is TubeRelation.ROBUST_INNER
    assert transfer.relation is MarginTransferRelation.PRESERVED_ROBUST
    assert replay.valid
    assert replay.relation is MarginTransferRelation.PRESERVED_ROBUST
    document = json.loads(transfer.to_json())
    row = document["proof"]["per_disk"][0]
    assert row["adverse_load_abs_upper"] == "0"
    assert _rat(row["remaining_margin_lower"]) > 0
    assert "OUTER" not in transfer.to_json()


def test_guaranteed_radius_increase_offsets_a_larger_center_shift() -> None:
    problem, box = _problem()
    baseline = certify_alpha_tube(problem, box)
    perturbations = _perturbation(
        problem,
        center_shift=0.10,
        new_response_radius=0.16,
    )

    transfer = certify_margin_transfer(problem, box, baseline.to_json(), perturbations)
    row = json.loads(transfer.to_json())["proof"]["per_disk"][0]

    assert transfer.relation is MarginTransferRelation.PRESERVED_ROBUST
    assert _rat(row["adverse_coefficient"]) < 0
    assert _rat(row["adverse_load_abs_upper"]) < 0


def test_morphology_radius_increase_is_monotone_beneficial() -> None:
    problem, box = _problem()
    baseline = certify_alpha_tube(problem, box)
    perturbations = _perturbation(
        problem,
        center_shift=0.05,
        new_morphology_radius=0.2,
    )

    transfer = certify_margin_transfer(problem, box, baseline.to_json(), perturbations)
    row = json.loads(transfer.to_json())["proof"]["per_disk"][0]

    assert transfer.relation is MarginTransferRelation.PRESERVED_ROBUST
    assert _rat(row["morphology_radius_delta"]) > 0
    assert _rat(row["remaining_margin_lower"]) > 0


def test_v2_selector_image_drives_free_damping_bounds_and_both_load_signs() -> None:
    problem, box = _free_damping_problem()
    baseline = certify_alpha_tube(problem, box)
    baseline_document = json.loads(baseline.to_json())
    damping_index = baseline_document["proof"]["variables"].index("damping")
    direct_damping = baseline_document["proof"]["krawczyk_image"][damping_index]
    selector_damping = baseline_document["proof"]["selector_image"][damping_index]
    harmonic = problem.harmonics[0]
    perturbations = tuple(
        DiskMarginPerturbation(
            harmonic_index=harmonic.harmonic_index,
            rate_id=rate_id,
            center_shift_abs_upper=center_shift,
            old_response_radius=harmonic.response_radii[rate_index],
            new_response_radius=new_radius,
            old_morphology_radius=harmonic.morphology_drift_radii[rate_index],
            new_morphology_radius=harmonic.morphology_drift_radii[rate_index],
        )
        for rate_index, (rate_id, center_shift, new_radius) in enumerate(
            zip(
                problem.rate_ids,
                (0.002, 0.001, 0.0),
                (0.051, 0.052, 0.05),
                strict=True,
            )
        )
    )

    transfer = certify_margin_transfer(problem, box, baseline.to_json(), perturbations)
    replay = replay_margin_transfer_certificate(
        problem,
        box,
        baseline.to_json(),
        perturbations,
        transfer.to_json(),
    )
    rows = {row["name"]: row for row in json.loads(transfer.to_json())["proof"]["per_disk"]}
    positive = rows["disk:m1:rest"]
    negative = rows["disk:m1:mid"]

    assert baseline.relation is TubeRelation.ROBUST_INNER
    assert transfer.relation is MarginTransferRelation.PRESERVED_ROBUST
    assert replay.valid
    positive_coefficient = _rat(positive["adverse_coefficient"])
    positive_denominator = tuple(_rat(value) for value in positive["denominator_abs"])
    positive_load = _rat(positive["adverse_load_abs_upper"])
    assert positive_coefficient > 0
    assert positive_load == positive_denominator[1] * positive_coefficient
    negative_coefficient = _rat(negative["adverse_coefficient"])
    negative_denominator = tuple(_rat(value) for value in negative["denominator_abs"])
    negative_load = _rat(negative["adverse_load_abs_upper"])
    assert negative_coefficient < 0
    assert negative_load == negative_denominator[0] * negative_coefficient

    selector_bounds = _denominator_bounds_from_stored_damping(
        selector_damping,
        box,
        harmonic.base_frequencies[0],
        harmonic.harmonic_index,
        baseline.precision_bits,
    )
    direct_bounds = _denominator_bounds_from_stored_damping(
        direct_damping,
        box,
        harmonic.base_frequencies[0],
        harmonic.harmonic_index,
        baseline.precision_bits,
    )
    transfer_bounds = tuple(float(value) for value in positive_denominator)
    assert transfer_bounds[0] <= selector_bounds[0] <= selector_bounds[1] <= transfer_bounds[1]
    assert selector_bounds[0] - transfer_bounds[0] < 1.0e-12
    assert transfer_bounds[1] - selector_bounds[1] < 1.0e-12
    assert transfer_bounds[0] > direct_bounds[0]
    assert transfer_bounds[1] < direct_bounds[1]


def test_exact_claim_consumption_downgrades_to_closed() -> None:
    problem, box = _problem(morphology_radius=1.0)
    baseline = certify_alpha_tube(problem, box)
    baseline_document = json.loads(baseline.to_json())
    disk_entry = next(
        value
        for value in baseline_document["proof"]["constraint_margins"]
        if value["name"] == "disk:m1:rest"
    )
    claimed_margin = _rat(disk_entry["lower"])
    target = Fraction(1) - claimed_margin
    assert target >= 0
    assert Fraction.from_float(float(target)) == target
    perturbations = _perturbation(problem, new_morphology_radius=float(target))

    transfer = certify_margin_transfer(problem, box, baseline.to_json(), perturbations)
    row = json.loads(transfer.to_json())["proof"]["per_disk"][0]

    assert transfer.relation is MarginTransferRelation.PRESERVED_CLOSED
    assert row["remaining_margin_lower"] == "0"
    assert replay_margin_transfer_certificate(
        problem,
        box,
        baseline.to_json(),
        perturbations,
        transfer.to_json(),
    ).valid


def test_exceeded_budget_is_replayable_unknown_never_outer() -> None:
    problem, box = _problem()
    baseline = certify_alpha_tube(problem, box)
    perturbations = _perturbation(problem, center_shift=10.0)

    transfer = certify_margin_transfer(problem, box, baseline.to_json(), perturbations)
    replay = replay_margin_transfer_certificate(
        problem,
        box,
        baseline.to_json(),
        perturbations,
        transfer.to_json(),
    )

    assert transfer.relation is MarginTransferRelation.UNKNOWN
    assert transfer.reason == "margin_budget_exceeded"
    assert replay.valid
    assert replay.relation is MarginTransferRelation.UNKNOWN
    assert "OUTER" not in transfer.to_json()


def test_closed_baseline_cannot_be_promoted_by_transfer() -> None:
    problem, box = _problem(response_radius=0.0, morphology_radius=0.0)
    baseline = certify_alpha_tube(problem, box)

    transfer = certify_margin_transfer(
        problem,
        box,
        baseline.to_json(),
        _perturbation(problem, new_response_radius=1.0),
    )

    assert baseline.relation is TubeRelation.CLOSED_INNER
    assert transfer.relation is MarginTransferRelation.UNKNOWN
    assert transfer.reason == "baseline_not_robust"


def test_old_radii_and_complete_canonical_layout_are_mandatory() -> None:
    problem, box = _problem()
    baseline = certify_alpha_tube(problem, box)
    valid = _perturbation(problem)[0]

    with pytest.raises(ValueError, match="one perturbation"):
        certify_margin_transfer(problem, box, baseline.to_json(), ())
    with pytest.raises(ValueError, match="old_response_radius"):
        certify_margin_transfer(
            problem,
            box,
            baseline.to_json(),
            (replace(valid, old_response_radius=valid.old_response_radius + 0.001),),
        )
    with pytest.raises(ValueError, match="nonnegative"):
        replace(valid, center_shift_abs_upper=-0.1)


def test_replay_binds_the_perturbation_manifest() -> None:
    problem, box = _problem()
    baseline = certify_alpha_tube(problem, box)
    original = _perturbation(problem)
    transfer = certify_margin_transfer(problem, box, baseline.to_json(), original)
    changed = _perturbation(problem, center_shift=1.0e-4)

    replay = replay_margin_transfer_certificate(
        problem,
        box,
        baseline.to_json(),
        changed,
        transfer.to_json(),
    )

    assert not replay.valid
    assert replay.relation is MarginTransferRelation.UNKNOWN
    assert replay.reason == "invalid_certificate_envelope"


def test_tampered_margin_is_rejected_even_after_digest_recomputation() -> None:
    problem, box = _problem()
    baseline = certify_alpha_tube(problem, box)
    perturbations = _perturbation(problem)
    transfer = certify_margin_transfer(problem, box, baseline.to_json(), perturbations)
    document = json.loads(transfer.to_json())
    document["proof"]["per_disk"][0]["remaining_margin_lower"] = "999"
    unsigned = dict(document)
    del unsigned["certificate_sha256"]
    document["certificate_sha256"] = hashlib.sha256(
        _canonical_json(unsigned).encode("utf-8")
    ).hexdigest()

    replay = replay_margin_transfer_certificate(
        problem,
        box,
        baseline.to_json(),
        perturbations,
        _canonical_json(document),
    )

    assert not replay.valid
    assert replay.relation is MarginTransferRelation.UNKNOWN
    assert replay.reason == "transfer_proof_mismatch"


def test_duplicate_json_key_is_rejected_before_replay() -> None:
    problem, box = _problem()
    baseline = certify_alpha_tube(problem, box)
    perturbations = _perturbation(problem)
    transfer = certify_margin_transfer(problem, box, baseline.to_json(), perturbations)
    duplicated = transfer.to_json().replace(
        '"algorithm":"frozen-selector-monotone-transfer-v1",',
        '"algorithm":"frozen-selector-monotone-transfer-v1",'
        '"algorithm":"frozen-selector-monotone-transfer-v1",',
        1,
    )

    replay = replay_margin_transfer_certificate(
        problem,
        box,
        baseline.to_json(),
        perturbations,
        duplicated,
    )

    assert duplicated != transfer.to_json()
    assert not replay.valid
    assert replay.relation is MarginTransferRelation.UNKNOWN
    assert replay.reason == "invalid_certificate_envelope"


def test_deeply_nested_certificate_is_rejected_before_replay() -> None:
    problem, box = _problem()
    baseline = certify_alpha_tube(problem, box)
    perturbations = _perturbation(problem)
    deeply_nested = "[" * 1100 + "0" + "]" * 1100

    replay = replay_margin_transfer_certificate(
        problem,
        box,
        baseline.to_json(),
        perturbations,
        deeply_nested,
    )

    assert not replay.valid
    assert replay.relation is MarginTransferRelation.UNKNOWN
    assert replay.reason == "invalid_certificate_envelope"


def test_baseline_material_prefers_selector_image_with_field_fallback() -> None:
    base_proof = {
        "variables": ["damping"],
        "krawczyk_image": [["3", "4"]],
        "constraint_margins": [{"name": "disk:m1:rest", "lower": "1/2"}],
    }
    fallback_document = {"precision_bits": 128, "proof": base_proof}
    _, fallback_image, _ = margin_transfer_module._baseline_material(
        fallback_document,
        ("disk:m1:rest",),
    )
    v2_document = {
        "precision_bits": 128,
        "proof": {**base_proof, "selector_image": [["1", "2"]]},
    }
    _, v2_image, _ = margin_transfer_module._baseline_material(
        v2_document,
        ("disk:m1:rest",),
    )

    assert fallback_image["damping"] == (Fraction(3), Fraction(4))
    assert v2_image["damping"] == (Fraction(1), Fraction(2))


def test_invalid_baseline_is_only_unknown_and_replays_as_unknown() -> None:
    problem, box = _problem()
    baseline = certify_alpha_tube(problem, box)
    baseline_document = json.loads(baseline.to_json())
    baseline_document["certificate_sha256"] = "0" * 64
    invalid_baseline = _canonical_json(baseline_document)
    perturbations = _perturbation(problem)

    transfer = certify_margin_transfer(problem, box, invalid_baseline, perturbations)
    replay = replay_margin_transfer_certificate(
        problem,
        box,
        invalid_baseline,
        perturbations,
        transfer.to_json(),
    )

    assert transfer.relation is MarginTransferRelation.UNKNOWN
    assert transfer.reason == "baseline_certificate_invalid"
    assert replay.valid
    assert replay.relation is MarginTransferRelation.UNKNOWN
