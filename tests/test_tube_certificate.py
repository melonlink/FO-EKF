import cmath
import hashlib
import json
import math
from fractions import Fraction

import pytest

from fo_ekf.certified_set import (
    CertifiedHarmonicData,
    CertifiedJointProblem,
    CertifiedParameterBox,
    CertifiedRelation,
    ClosedInterval,
    RateOffsetInterval,
    certify_parameter_box,
)
from fo_ekf.tube_certificate import (
    TubeRelation,
    branch_alpha_tube,
    certify_alpha_tube,
    replay_alpha_tube_branch,
    replay_tube_certificate,
)


def _fractional_frequency(order: float, frequency: float) -> complex:
    return frequency**order * cmath.exp(0.5j * math.pi * order)


def _single_rate_problem(
    *,
    scale: float = 1.0,
    response_radius: float = 0.05,
    order_interval: tuple[float, float] = (0.6, 0.8),
    morphology_bounds: tuple[float, float] = (0.1, 5.0),
) -> tuple[CertifiedJointProblem, CertifiedParameterBox]:
    true_order = 0.7
    damping = 0.5
    frequency = 1.2
    morphology = scale * (1.0 - 0.2j)
    response = morphology / (damping + _fractional_frequency(true_order, frequency))
    harmonic = CertifiedHarmonicData(
        rate_ids=("rest",),
        base_frequencies=(frequency,),
        harmonic_index=1,
        measured_responses=(response,),
        response_radii=(scale * response_radius,),
        morphology_drift_radii=(0.0,),
        morphology_bounds=(
            scale * morphology_bounds[0],
            scale * morphology_bounds[1],
        ),
        label="m1",
    )
    problem = CertifiedJointProblem((harmonic,), "rest")
    box = CertifiedParameterBox(
        ClosedInterval(*order_interval),
        ClosedInterval(damping, damping),
        (RateOffsetInterval("rest", ClosedInterval(0.0, 0.0)),),
    )
    return problem, box


def _multirate_free_damping_problem(
    *,
    scale: float = 1.0,
    order_interval: tuple[float, float] = (0.7199, 0.7201),
) -> tuple[CertifiedJointProblem, CertifiedParameterBox]:
    rate_ids = ("rest", "mid", "high")
    frequencies = (0.8, 1.25, 1.9)
    true_order = 0.72
    damping = 0.5
    offsets = (0.0, 0.03, -0.02)
    morphology = scale * (1.0 - 0.25j)
    responses = tuple(
        morphology / (damping + offset + _fractional_frequency(true_order, frequency))
        for frequency, offset in zip(frequencies, offsets, strict=True)
    )
    harmonic = CertifiedHarmonicData(
        rate_ids=rate_ids,
        base_frequencies=frequencies,
        harmonic_index=1,
        measured_responses=responses,
        response_radii=(scale * 0.05,) * 3,
        morphology_drift_radii=(0.0,) * 3,
        morphology_bounds=(scale * 0.5, scale * 1.5),
        label="m1",
    )
    problem = CertifiedJointProblem((harmonic,), "rest")
    box = CertifiedParameterBox(
        ClosedInterval(*order_interval),
        ClosedInterval(0.49, 0.51),
        tuple(
            RateOffsetInterval(rate_id, ClosedInterval(offset, offset))
            for rate_id, offset in zip(rate_ids, offsets, strict=True)
        ),
    )
    return problem, box


def _canonical_json(value: object) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def test_robust_tube_replays_and_binds_exact_input_hash() -> None:
    problem, box = _single_rate_problem()

    certificate = certify_alpha_tube(problem, box)
    replay = replay_tube_certificate(problem, box, certificate.to_json())

    assert certificate.relation is TubeRelation.ROBUST_INNER
    assert replay.valid
    assert replay.relation is TubeRelation.ROBUST_INNER
    document = json.loads(certificate.to_json())
    assert document["input_sha256"] == certificate.input_sha256
    assert document["certificate_sha256"] == certificate.certificate_sha256
    assert document["proof"]["contraction_inf_norm"] == "0"
    assert all(isinstance(value, str) for row in document["proof"]["variable_box"] for value in row)

    changed_problem, changed_box = _single_rate_problem(response_radius=0.051)
    changed_replay = replay_tube_certificate(
        changed_problem,
        changed_box,
        certificate.to_json(),
    )
    assert not changed_replay.valid
    assert changed_replay.relation is TubeRelation.UNKNOWN


def test_tamper_is_rejected_even_if_attacker_recomputes_document_digest() -> None:
    problem, box = _single_rate_problem()
    certificate = certify_alpha_tube(problem, box)
    document = json.loads(certificate.to_json())
    document["proof"]["preconditioner_c"][0][0] = "2"
    unsigned = dict(document)
    del unsigned["certificate_sha256"]
    document["certificate_sha256"] = hashlib.sha256(
        _canonical_json(unsigned).encode("utf-8")
    ).hexdigest()

    replay = replay_tube_certificate(problem, box, _canonical_json(document))

    assert not replay.valid
    assert replay.relation is TubeRelation.UNKNOWN
    assert replay.reason == "malformed_proof"


@pytest.mark.parametrize("scale", (1.0e-9, 1.0, 1.0e9))
def test_common_measurement_scale_preserves_robust_tube(scale: float) -> None:
    problem, box = _multirate_free_damping_problem(scale=scale)

    certificate = certify_alpha_tube(problem, box)

    assert certificate.relation is TubeRelation.ROBUST_INNER
    assert replay_tube_certificate(problem, box, certificate.to_json()).valid


def test_variable_q_tube_closes_fixed_q_counterexample() -> None:
    """No constant q works over the box, but q(alpha) does for every alpha."""

    problem, box = _single_rate_problem(response_radius=1.0e-4)

    old_uniform_q = certify_parameter_box(problem, box)
    tube = certify_alpha_tube(problem, box)

    assert old_uniform_q.relation is CertifiedRelation.UNKNOWN
    assert tube.relation is TubeRelation.ROBUST_INNER
    replay = replay_tube_certificate(problem, box, tube.to_json())
    assert replay.valid
    assert replay.relation is TubeRelation.ROBUST_INNER

    harmonic = problem.harmonics[0]
    response = harmonic.measured_responses[0]
    endpoint_denominators = tuple(
        box.damping.lower + _fractional_frequency(order, harmonic.base_frequencies[0])
        for order in (box.order.lower, box.order.upper)
    )
    endpoint_centers = tuple(value * response for value in endpoint_denominators)
    endpoint_radii = tuple(
        abs(value) * harmonic.response_radii[0] for value in endpoint_denominators
    )
    assert abs(endpoint_centers[0] - endpoint_centers[1]) > sum(endpoint_radii)


def test_zero_radius_point_selector_is_closed_inner() -> None:
    problem, broad_box = _single_rate_problem(
        response_radius=0.0,
        order_interval=(0.7, 0.7),
    )

    certificate = certify_alpha_tube(problem, broad_box)

    assert certificate.relation is TubeRelation.CLOSED_INNER
    replay = replay_tube_certificate(problem, broad_box, certificate.to_json())
    assert replay.valid
    assert replay.relation is TubeRelation.CLOSED_INNER


def test_exact_physical_boundary_is_unknown_not_false_robust() -> None:
    """A transcendental boundary that Arb cannot sign stays UNKNOWN."""

    harmonic = CertifiedHarmonicData(
        rate_ids=("rest",),
        base_frequencies=(1.0,),
        harmonic_index=1,
        measured_responses=(0.5 - 0.5j,),
        response_radii=(0.0,),
        morphology_drift_radii=(0.0,),
        morphology_bounds=(1.0, 2.0),
        label="boundary",
    )
    problem = CertifiedJointProblem((harmonic,), "rest")
    box = CertifiedParameterBox(
        ClosedInterval(1.0, 1.0),
        ClosedInterval(1.0, 1.0),
        (RateOffsetInterval("rest", ClosedInterval(0.0, 0.0)),),
    )

    certificate = certify_alpha_tube(problem, box)

    assert certificate.relation is TubeRelation.UNKNOWN
    assert replay_tube_certificate(problem, box, certificate.to_json()).valid


def test_shared_free_damping_multirate_tube_replays() -> None:
    problem, box = _multirate_free_damping_problem()

    certificate = certify_alpha_tube(problem, box)

    assert certificate.relation is TubeRelation.ROBUST_INNER
    assert "damping" in json.loads(certificate.to_json())["proof"]["variables"]
    assert replay_tube_certificate(problem, box, certificate.to_json()).valid


def test_real_separable_multirate_selector_is_exactly_rank_deficient() -> None:
    rate_ids = ("a", "b")
    first = CertifiedHarmonicData(
        rate_ids=rate_ids,
        base_frequencies=(1.0, 2.0),
        harmonic_index=1,
        measured_responses=(1.0 + 1.0j, 2.0 + 2.0j),
        response_radii=(0.5, 0.5),
        morphology_drift_radii=(0.0, 0.0),
        morphology_bounds=(0.1, 10.0),
        label="m1",
    )
    second = CertifiedHarmonicData(
        rate_ids=rate_ids,
        base_frequencies=(1.0, 2.0),
        harmonic_index=2,
        measured_responses=(1.0 - 2.0j, 2.0 - 4.0j),
        response_radii=(0.5, 0.5),
        morphology_drift_radii=(0.0, 0.0),
        morphology_bounds=(0.1, 10.0),
        label="m2",
    )
    problem = CertifiedJointProblem((first, second), "a")
    box = CertifiedParameterBox(
        ClosedInterval(0.7, 0.7),
        ClosedInterval(0.4, 0.6),
        (
            RateOffsetInterval("a", ClosedInterval(0.0, 0.0)),
            RateOffsetInterval("b", ClosedInterval(-0.1, 0.1)),
        ),
    )

    certificate = certify_alpha_tube(problem, box)

    assert certificate.relation is TubeRelation.UNKNOWN
    assert certificate.reason == "selector_rank_deficient"
    assert replay_tube_certificate(problem, box, certificate.to_json()).valid
    branch = branch_alpha_tube(problem, box, max_depth=8, max_leaves=256)
    assert len(branch.leaves) == 1
    assert branch.leaves[0].terminal_reason == "structural_unknown"
    assert not branch.budget_exhausted
    assert replay_alpha_tube_branch(problem, box, branch).valid


def test_lens_intersection_uses_gram_point_not_observation_center() -> None:
    """Neither disk center is feasible, while their LS midpoint is strict."""

    rate_ids = ("left", "right")
    harmonic = CertifiedHarmonicData(
        rate_ids=rate_ids,
        base_frequencies=(1.0, 1.0),
        harmonic_index=1,
        measured_responses=(0.5 - 0.5j, 1.5 - 1.5j),
        response_radii=(0.0, 0.0),
        morphology_drift_radii=(1.1, 1.1),
        morphology_bounds=(0.5, 3.5),
        label="lens",
    )
    problem = CertifiedJointProblem((harmonic,), "left")
    box = CertifiedParameterBox(
        ClosedInterval(1.0, 1.0),
        ClosedInterval(1.0, 1.0),
        (
            RateOffsetInterval("left", ClosedInterval(0.0, 0.0)),
            RateOffsetInterval("right", ClosedInterval(0.0, 0.0)),
        ),
    )
    centers = (1.0 + 0.0j, 3.0 + 0.0j)
    radius = 1.1
    assert abs(centers[0] - centers[1]) > radius
    assert abs(2.0 - centers[0]) < radius
    assert abs(2.0 - centers[1]) < radius

    certificate = certify_alpha_tube(problem, box)

    assert certificate.relation is TubeRelation.ROBUST_INNER
    document = json.loads(certificate.to_json())
    q_real_index = document["proof"]["variables"].index("q_re:m1")
    q_real_image = document["proof"]["krawczyk_image"][q_real_index]
    assert float(Fraction(q_real_image[0])) < 2.0 < float(Fraction(q_real_image[1]))
    assert replay_tube_certificate(problem, box, certificate.to_json()).valid


def test_alpha_only_subdivision_covers_inner_and_unknown_and_replays() -> None:
    problem, initial = _multirate_free_damping_problem(order_interval=(0.7, 0.74))

    result = branch_alpha_tube(problem, initial, max_depth=5, max_leaves=64)
    replay = replay_alpha_tube_branch(problem, initial, result)

    assert result.coverage_certified
    assert result.budget_exhausted
    assert result.certified_inner
    assert result.unknown
    assert any(leaf.certificate.relation is TubeRelation.ROBUST_INNER for leaf in result.leaves)
    assert any(leaf.certificate.relation is TubeRelation.UNKNOWN for leaf in result.leaves)
    assert replay.valid
    assert len(replay.leaf_replays) == len(result.leaves)
    assert all(value.valid for value in replay.leaf_replays)
    assert sum(Fraction(1, 2 ** len(leaf.path)) for leaf in result.leaves) == 1

    repeated = branch_alpha_tube(problem, initial, max_depth=5, max_leaves=64)
    assert tuple(leaf.path for leaf in repeated.leaves) == tuple(
        leaf.path for leaf in result.leaves
    )
    assert tuple(leaf.certificate.certificate_sha256 for leaf in repeated.leaves) == tuple(
        leaf.certificate.certificate_sha256 for leaf in result.leaves
    )


def test_one_dyadic_split_closes_parent_and_leaf_proofs_are_domain_bound() -> None:
    rate_ids = ("slow", "fast")
    frequencies = (0.7, 2.0)
    true_order = 0.7
    damping = 0.5
    morphology = 1.0 - 0.2j
    responses = tuple(
        morphology / (damping + _fractional_frequency(true_order, frequency))
        for frequency in frequencies
    )
    harmonic = CertifiedHarmonicData(
        rate_ids=rate_ids,
        base_frequencies=frequencies,
        harmonic_index=1,
        measured_responses=responses,
        response_radii=(0.0, 0.0),
        morphology_drift_radii=(0.08, 0.08),
        morphology_bounds=(0.2, 2.0),
        label="split",
    )
    problem = CertifiedJointProblem((harmonic,), "slow")
    initial = CertifiedParameterBox(
        ClosedInterval(0.65, 0.75),
        ClosedInterval(damping, damping),
        (
            RateOffsetInterval("slow", ClosedInterval(0.0, 0.0)),
            RateOffsetInterval("fast", ClosedInterval(0.0, 0.0)),
        ),
    )
    assert certify_alpha_tube(problem, initial).relation is TubeRelation.UNKNOWN

    result = branch_alpha_tube(problem, initial, max_depth=1, max_leaves=2)

    assert tuple(leaf.path for leaf in result.leaves) == ((0,), (1,))
    assert all(leaf.certificate.relation is TubeRelation.ROBUST_INNER for leaf in result.leaves)
    assert result.leaves[0].box.order.upper == result.leaves[1].box.order.lower
    assert result.certified_inner == (initial.order,)
    assert result.unknown == ()
    assert replay_alpha_tube_branch(problem, initial, result).valid

    left, right = result.leaves
    assert not replay_tube_certificate(
        problem,
        right.box,
        left.certificate.to_json(),
    ).valid
    assert not replay_tube_certificate(
        problem,
        initial,
        left.certificate.to_json(),
    ).valid


def test_alpha_subdivision_budgets_keep_root_unknown() -> None:
    problem, initial = _multirate_free_damping_problem(order_interval=(0.7, 0.74))

    depth_limited = branch_alpha_tube(problem, initial, max_depth=0, max_leaves=4)
    leaf_limited = branch_alpha_tube(problem, initial, max_depth=4, max_leaves=1)

    assert depth_limited.leaves[0].terminal_reason == "max_depth_unknown"
    assert leaf_limited.leaves[0].terminal_reason == "max_leaves_unknown"
    for result in (depth_limited, leaf_limited):
        assert result.budget_exhausted
        assert result.certified_inner == ()
        assert result.unknown == (initial.order,)
        assert replay_alpha_tube_branch(problem, initial, result).valid


def test_rehashed_excessive_precision_is_rejected_before_arb_replay() -> None:
    problem, box = _single_rate_problem()
    certificate = certify_alpha_tube(problem, box)
    document = json.loads(certificate.to_json())
    document["precision_bits"] = 1_000_000
    document["attempted_precisions"] = [1_000_000]
    unsigned = dict(document)
    del unsigned["certificate_sha256"]
    document["certificate_sha256"] = hashlib.sha256(
        _canonical_json(unsigned).encode("utf-8")
    ).hexdigest()

    replay = replay_tube_certificate(problem, box, _canonical_json(document))

    assert not replay.valid
    assert replay.relation is TubeRelation.UNKNOWN
    assert replay.reason == "invalid_certificate_envelope"


def test_rehashed_closed_to_robust_upgrade_is_rejected() -> None:
    problem, box = _single_rate_problem(
        response_radius=0.0,
        order_interval=(0.7, 0.7),
    )
    certificate = certify_alpha_tube(problem, box)
    assert certificate.relation is TubeRelation.CLOSED_INNER
    document = json.loads(certificate.to_json())
    document["relation"] = TubeRelation.ROBUST_INNER.value
    unsigned = dict(document)
    del unsigned["certificate_sha256"]
    document["certificate_sha256"] = hashlib.sha256(
        _canonical_json(unsigned).encode("utf-8")
    ).hexdigest()

    replay = replay_tube_certificate(problem, box, _canonical_json(document))

    assert not replay.valid
    assert replay.relation is TubeRelation.UNKNOWN
    assert replay.reason == "robust_claim_without_margin"


def test_underdetermined_selector_is_unknown_and_never_outer() -> None:
    problem, broad_box = _single_rate_problem()
    free_damping_box = CertifiedParameterBox(
        broad_box.order,
        ClosedInterval(0.4, 0.6),
        broad_box.damping_offsets,
    )

    certificate = certify_alpha_tube(problem, free_damping_box)

    assert certificate.relation is TubeRelation.UNKNOWN
    assert certificate.reason == "selector_underdetermined"
    assert "OUTER" not in certificate.to_json()
    replay = replay_tube_certificate(problem, free_damping_box, certificate.to_json())
    assert replay.valid
    assert replay.relation is TubeRelation.UNKNOWN
