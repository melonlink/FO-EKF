import cmath
import math

import flint
import pytest

import fo_ekf.certified_set as certified_set
from fo_ekf.certified_set import (
    CertifiedBoxResult,
    CertifiedHarmonicData,
    CertifiedJointProblem,
    CertifiedLeaf,
    CertifiedParameterBox,
    CertifiedRelation,
    ClosedInterval,
    RateOffsetInterval,
    branch_parameter_box,
    certify_parameter_box,
)

RATE_IDS = ("rest", "mid", "high")
BASE_FREQUENCIES = (0.8, 1.25, 1.9)
ANCHOR = "rest"


def test_python_flint_090_negation_regression() -> None:
    assert flint.__version__ == "0.9.0"
    assert -flint.arb(2) == -2


def _fractional_frequency(order: float, frequency: float) -> complex:
    amplitude = math.exp(order * math.log(frequency))
    angle = 0.5 * math.pi * order
    return amplitude * complex(math.cos(angle), math.sin(angle))


def _responses(
    order: float,
    damping: float,
    offsets: tuple[float, ...],
    morphology: complex,
    harmonic_index: int,
    base_frequencies: tuple[float, ...] = BASE_FREQUENCIES,
) -> tuple[complex, ...]:
    return tuple(
        morphology
        / (damping + offset + _fractional_frequency(order, harmonic_index * base_frequency))
        for base_frequency, offset in zip(base_frequencies, offsets, strict=True)
    )


def _harmonic(
    order: float,
    damping: float,
    offsets: tuple[float, ...],
    morphology: complex,
    harmonic_index: int,
    *,
    response_radius: float = 0.05,
    drift_radius: float = 0.0,
    q_margin: float = 0.3,
    rate_ids: tuple[str, ...] = RATE_IDS,
    base_frequencies: tuple[float, ...] = BASE_FREQUENCIES,
    scale: float = 1.0,
) -> CertifiedHarmonicData:
    magnitude = abs(morphology)
    return CertifiedHarmonicData(
        rate_ids=rate_ids,
        base_frequencies=base_frequencies,
        harmonic_index=harmonic_index,
        measured_responses=tuple(
            scale * value
            for value in _responses(
                order,
                damping,
                offsets,
                morphology,
                harmonic_index,
                base_frequencies,
            )
        ),
        response_radii=(scale * response_radius,) * len(rate_ids),
        morphology_drift_radii=(scale * drift_radius,) * len(rate_ids),
        morphology_bounds=(
            scale * (1.0 - q_margin) * magnitude,
            scale * (1.0 + q_margin) * magnitude,
        ),
        label=f"m{harmonic_index}",
    )


def _problem(
    order: float = 0.72,
    damping: float = 0.5,
    offsets: tuple[float, ...] = (0.0, 0.03, -0.02),
    *,
    response_radius: float = 0.05,
    scale: float = 1.0,
) -> CertifiedJointProblem:
    return CertifiedJointProblem(
        (
            _harmonic(
                order,
                damping,
                offsets,
                1.0 - 0.25j,
                1,
                response_radius=response_radius,
                scale=scale,
            ),
            _harmonic(
                order,
                damping,
                offsets,
                0.55 + 0.65j,
                2,
                response_radius=response_radius,
                scale=scale,
            ),
        ),
        ANCHOR,
    )


def _box(
    order: ClosedInterval,
    damping: ClosedInterval,
    offsets: tuple[ClosedInterval, ...],
) -> CertifiedParameterBox:
    return CertifiedParameterBox(
        order,
        damping,
        tuple(
            RateOffsetInterval(rate_id, interval)
            for rate_id, interval in zip(RATE_IDS, offsets, strict=True)
        ),
    )


def _point_box(
    order: float,
    damping: float,
    offsets: tuple[float, ...],
) -> CertifiedParameterBox:
    return _box(
        ClosedInterval(order, order),
        ClosedInterval(damping, damping),
        tuple(ClosedInterval(value, value) for value in offsets),
    )


def test_rate_permutation_and_wrong_rate_id_are_rejected() -> None:
    order = 0.72
    damping = 0.5
    offsets = (0.0, 0.03, -0.02)
    first = _harmonic(order, damping, offsets, 1.0 - 0.25j, 1)
    permutation = (2, 0, 1)
    permuted_ids = tuple(RATE_IDS[index] for index in permutation)
    permuted_frequencies = tuple(BASE_FREQUENCIES[index] for index in permutation)
    permuted_offsets = tuple(offsets[index] for index in permutation)
    second_permuted = _harmonic(
        order,
        damping,
        permuted_offsets,
        0.55 + 0.65j,
        2,
        rate_ids=permuted_ids,
        base_frequencies=permuted_frequencies,
    )
    with pytest.raises(ValueError, match="permutation"):
        CertifiedJointProblem((first, second_permuted), ANCHOR)

    wrong_ids = ("rest", "mid", "other")
    second_wrong = _harmonic(
        order,
        damping,
        offsets,
        0.55 + 0.65j,
        2,
        rate_ids=wrong_ids,
    )
    with pytest.raises(ValueError, match="same rate ids"):
        CertifiedJointProblem((first, second_wrong), ANCHOR)


def test_degenerate_morphology_annulus_is_rejected() -> None:
    with pytest.raises(ValueError, match="q_min < q_max"):
        CertifiedHarmonicData(
            rate_ids=(ANCHOR,),
            base_frequencies=(1.0,),
            harmonic_index=1,
            measured_responses=(1.0 + 0.0j,),
            response_radii=(0.0,),
            morphology_drift_radii=(0.0,),
            morphology_bounds=(1.0, 1.0),
        )


def test_wrong_anchor_and_nonzero_anchor_offset_are_rejected() -> None:
    harmonic = _harmonic(0.72, 0.5, (0.0, 0.03, -0.02), 1.0 - 0.25j, 1)
    with pytest.raises(ValueError, match="anchor"):
        CertifiedJointProblem((harmonic,), "missing")

    problem = CertifiedJointProblem((harmonic,), ANCHOR)
    invalid = _box(
        ClosedInterval(0.7, 0.8),
        ClosedInterval(0.4, 0.6),
        (
            ClosedInterval(1.0e-12, 1.0e-12),
            ClosedInterval(0.02, 0.04),
            ClosedInterval(-0.03, -0.01),
        ),
    )
    with pytest.raises(ValueError, match="exactly"):
        certify_parameter_box(problem, invalid)


def test_offset_rate_permutation_is_rejected() -> None:
    problem = _problem()
    offsets = (
        RateOffsetInterval("mid", ClosedInterval(0.02, 0.04)),
        RateOffsetInterval("rest", ClosedInterval(0.0, 0.0)),
        RateOffsetInterval("high", ClosedInterval(-0.03, -0.01)),
    )
    box = CertifiedParameterBox(
        ClosedInterval(0.7, 0.8),
        ClosedInterval(0.4, 0.6),
        offsets,
    )

    with pytest.raises(ValueError, match="permutation"):
        certify_parameter_box(problem, box)


def test_alpha_one_uses_principal_branch_and_certifies_true_point() -> None:
    order = 1.0
    damping = 0.62
    offsets = (0.0, 0.04, -0.03)
    problem = _problem(
        order,
        damping,
        offsets,
        response_radius=1.0e-7,
    )

    result = certify_parameter_box(problem, _point_box(order, damping, offsets))

    assert result.relation is CertifiedRelation.INNER_INCLUDED
    assert result.precision_bits in (128, 256, 512)
    assert all(witness is not None for witness in result.morphology_witnesses)


def test_tolerance_false_inner_counterexample_is_not_inner_certified() -> None:
    tolerance = 1.0e-10
    order = 0.7
    damping = 0.5
    rate_ids = ("rest",)
    base_frequencies = (0.8,)
    offsets = (0.0,)
    morphology = complex(1.0 - 0.5 * tolerance, 0.0)
    harmonic = CertifiedHarmonicData(
        rate_ids=rate_ids,
        base_frequencies=base_frequencies,
        harmonic_index=1,
        measured_responses=_responses(
            order,
            damping,
            offsets,
            morphology,
            1,
            base_frequencies,
        ),
        response_radii=(0.0,),
        morphology_drift_radii=(0.0,),
        morphology_bounds=(1.0, 2.0),
    )
    problem = CertifiedJointProblem((harmonic,), ANCHOR)
    box = CertifiedParameterBox(
        ClosedInterval(order, order),
        ClosedInterval(damping, damping),
        (RateOffsetInterval(ANCHOR, ClosedInterval(0.0, 0.0)),),
    )

    result = certify_parameter_box(problem, box)

    assert result.relation is CertifiedRelation.OUTER_EXCLUDED
    assert result.reason.startswith("outer_disk_inside_annulus_hole")


def test_outer_pair_separation_is_arb_certified() -> None:
    order = 0.68
    damping = 0.45
    offsets = (0.0, 0.0, 0.0)
    rate_morphologies = (1.0 + 0.0j, 2.0 + 0.0j, 1.5 + 0.0j)
    responses = tuple(
        morphology / (damping + _fractional_frequency(order, base_frequency))
        for morphology, base_frequency in zip(
            rate_morphologies,
            BASE_FREQUENCIES,
            strict=True,
        )
    )
    harmonic = CertifiedHarmonicData(
        rate_ids=RATE_IDS,
        base_frequencies=BASE_FREQUENCIES,
        harmonic_index=1,
        measured_responses=responses,
        response_radii=(0.0, 0.0, 0.0),
        morphology_drift_radii=(0.0, 0.0, 0.0),
        morphology_bounds=(0.5, 3.0),
    )
    problem = CertifiedJointProblem((harmonic,), ANCHOR)

    result = certify_parameter_box(problem, _point_box(order, damping, offsets))

    assert result.relation is CertifiedRelation.OUTER_EXCLUDED
    assert result.reason.startswith("outer_pair_separated")


def test_outer_single_disk_beyond_qmax_is_arb_certified() -> None:
    order = 0.68
    damping = 0.45
    rate_ids = ("rest",)
    base_frequencies = (0.8,)
    offsets = (0.0,)
    harmonic = CertifiedHarmonicData(
        rate_ids=rate_ids,
        base_frequencies=base_frequencies,
        harmonic_index=1,
        measured_responses=_responses(
            order,
            damping,
            offsets,
            2.0 + 0.0j,
            1,
            base_frequencies,
        ),
        response_radii=(0.0,),
        morphology_drift_radii=(0.0,),
        morphology_bounds=(0.2, 1.0),
    )
    problem = CertifiedJointProblem((harmonic,), ANCHOR)
    box = CertifiedParameterBox(
        ClosedInterval(order, order),
        ClosedInterval(damping, damping),
        (RateOffsetInterval(ANCHOR, ClosedInterval(0.0, 0.0)),),
    )

    result = certify_parameter_box(problem, box)

    assert result.relation is CertifiedRelation.OUTER_EXCLUDED
    assert result.reason.startswith("outer_disk_outside_qmax")


def test_qbox_certifies_pairwise_overlap_with_empty_triple_intersection() -> None:
    """The q-box cover handles a classical failure of pairwise disk tests."""

    order = 0.68
    damping = 0.45
    offsets = (0.0, 0.0, 0.0)
    side = 1.9
    height = 0.5 * math.sqrt(3.0) * side
    centers = (
        complex(-0.5 * side, -height / 3.0),
        complex(0.5 * side, -height / 3.0),
        complex(0.0, 2.0 * height / 3.0),
    )
    denominators = tuple(
        damping + _fractional_frequency(order, base_frequency)
        for base_frequency in BASE_FREQUENCIES
    )
    harmonic = CertifiedHarmonicData(
        rate_ids=RATE_IDS,
        base_frequencies=BASE_FREQUENCIES,
        harmonic_index=1,
        measured_responses=tuple(
            center / denominator for center, denominator in zip(centers, denominators, strict=True)
        ),
        response_radii=tuple(1.0 / abs(value) for value in denominators),
        morphology_drift_radii=(0.0, 0.0, 0.0),
        morphology_bounds=(0.01, 3.0),
    )
    problem = CertifiedJointProblem((harmonic,), ANCHOR)

    result = certify_parameter_box(problem, _point_box(order, damping, offsets))

    assert result.relation is CertifiedRelation.OUTER_EXCLUDED
    assert result.reason.startswith("outer_qbox_cover")


def test_qbox_budget_and_closed_boundaries_remain_unknown() -> None:
    feasible = certified_set._HarmonicEnclosure(
        outer_disks=(certified_set._BallDisk(center=flint.acb(1), radius=flint.arb("0.25")),),
        inner_disks=None,
        morphology_bounds=(flint.arb("0.5"), flint.arb(2)),
    )
    assert (
        certified_set._qbox_outer_exclusion_reason(
            feasible,
            max_depth=0,
            max_nodes=1,
        )
        is None
    )

    tangent_pair = certified_set._HarmonicEnclosure(
        outer_disks=(
            certified_set._BallDisk(center=flint.acb(0), radius=flint.arb(1)),
            certified_set._BallDisk(center=flint.acb(2), radius=flint.arb(1)),
        ),
        inner_disks=None,
        morphology_bounds=(flint.arb("0.5"), flint.arb(3)),
    )
    assert certified_set._fast_outer_exclusion_reason(tangent_pair) is None

    qmax_tangent = certified_set._HarmonicEnclosure(
        outer_disks=(certified_set._BallDisk(center=flint.acb(2), radius=flint.arb(1)),),
        inner_disks=None,
        morphology_bounds=(flint.arb("0.5"), flint.arb(1)),
    )
    assert certified_set._fast_outer_exclusion_reason(qmax_tangent) is None


@pytest.mark.parametrize("scale", (1.0e-9, 1.0, 1.0e9))
def test_common_measurement_scaling_preserves_inner_certificate(scale: float) -> None:
    order = 0.72
    damping = 0.5
    offsets = (0.0, 0.03, -0.02)
    problem = _problem(
        order,
        damping,
        offsets,
        response_radius=0.08,
        scale=scale,
    )
    box = _box(
        ClosedInterval(order - 1.0e-5, order + 1.0e-5),
        ClosedInterval(damping - 1.0e-5, damping + 1.0e-5),
        (
            ClosedInterval(0.0, 0.0),
            ClosedInterval(offsets[1] - 1.0e-5, offsets[1] + 1.0e-5),
            ClosedInterval(offsets[2] - 1.0e-5, offsets[2] + 1.0e-5),
        ),
    )

    result = certify_parameter_box(problem, box)

    assert result.relation is CertifiedRelation.INNER_INCLUDED


def test_true_parameter_in_box_is_never_outer_excluded() -> None:
    order = 0.77
    damping = 0.56
    offsets = (0.0, 0.025, -0.015)
    problem = _problem(
        order,
        damping,
        offsets,
        response_radius=0.02,
    )
    box = _box(
        ClosedInterval(order - 0.03, order + 0.03),
        ClosedInterval(damping - 0.05, damping + 0.05),
        (
            ClosedInterval(0.0, 0.0),
            ClosedInterval(offsets[1] - 0.02, offsets[1] + 0.02),
            ClosedInterval(offsets[2] - 0.02, offsets[2] + 0.02),
        ),
    )

    result = certify_parameter_box(problem, box)

    assert result.relation is not CertifiedRelation.OUTER_EXCLUDED


def test_branch_leaves_cover_initial_box_and_form_alpha_sandwich() -> None:
    order = 0.72
    damping = 0.5
    offsets = (0.0, 0.03, -0.02)
    problem = _problem(
        order,
        damping,
        offsets,
        response_radius=0.015,
    )
    initial = _box(
        ClosedInterval(0.45, 0.95),
        ClosedInterval(0.25, 0.9),
        (
            ClosedInterval(0.0, 0.0),
            ClosedInterval(-0.08, 0.12),
            ClosedInterval(-0.1, 0.08),
        ),
    )

    result = branch_parameter_box(problem, initial, max_depth=2, max_leaves=8)

    assert result.coverage_certified
    assert 1 < len(result.leaves) <= 8
    for value in (0.45, 0.5, 0.72, 0.9, 0.95):
        assert any(leaf.box.order.contains(value) for leaf in result.leaves)
    parameter_samples = (
        (0.45, 0.25, 0.0, -0.08, -0.1),
        (0.72, 0.575, 0.0, 0.02, -0.01),
        (0.95, 0.9, 0.0, 0.12, 0.08),
    )
    for order_value, damping_value, anchor_value, mid_value, high_value in parameter_samples:
        assert any(
            leaf.box.order.contains(order_value)
            and leaf.box.damping.contains(damping_value)
            and all(
                offset.interval.contains(value)
                for offset, value in zip(
                    leaf.box.damping_offsets,
                    (anchor_value, mid_value, high_value),
                    strict=True,
                )
            )
            for leaf in result.leaves
        )
    assert result.alpha.outer_contains(order)
    for interval in result.alpha.inner:
        assert any(
            outer.lower <= interval.lower and interval.upper <= outer.upper
            for outer in result.alpha.outer
        )


def test_budget_exhaustion_keeps_unknown_leaf_and_outer_alpha() -> None:
    problem = _problem(response_radius=0.01)
    initial = _box(
        ClosedInterval(0.4, 1.0),
        ClosedInterval(0.2, 1.0),
        (
            ClosedInterval(0.0, 0.0),
            ClosedInterval(-0.1, 0.15),
            ClosedInterval(-0.1, 0.1),
        ),
    )

    result = branch_parameter_box(problem, initial, max_depth=0, max_leaves=1)

    assert result.budget_exhausted
    assert len(result.leaves) == 1
    assert result.leaves[0].certificate.relation is CertifiedRelation.UNKNOWN
    assert result.leaves[0].certificate.attempted_precisions == (128, 256, 512)
    assert result.leaves[0].certificate.precision_bits == 512
    assert result.leaves[0].terminal_reason == "max_depth_unknown"
    assert result.alpha.inner == ()
    assert result.alpha.outer == (initial.order,)


def test_internal_contradiction_is_sticky_across_precision_schedule(monkeypatch) -> None:
    problem = _problem()
    box = _point_box(0.72, 0.5, (0.0, 0.03, -0.02))
    attempted: list[int] = []

    def contradictory_result(
        _problem_value: CertifiedJointProblem,
        _box_value: CertifiedParameterBox,
        precision: int,
    ) -> CertifiedBoxResult:
        attempted.append(precision)
        relation = (
            CertifiedRelation.UNKNOWN if precision == 128 else CertifiedRelation.OUTER_EXCLUDED
        )
        reason = "arb_outer_inner_contradiction" if precision == 128 else "outer_pair_separated:0:1"
        return CertifiedBoxResult(relation, reason, precision, (precision,), (None, None))

    monkeypatch.setattr(
        certified_set,
        "_certify_at_current_precision",
        contradictory_result,
    )

    result = certify_parameter_box(problem, box)

    assert result.relation is CertifiedRelation.UNKNOWN
    assert result.reason == "arb_outer_inner_contradiction"
    assert result.attempted_precisions == (128,)
    assert attempted == [128]

    branched = branch_parameter_box(problem, box, max_depth=4, max_leaves=16)
    assert len(branched.leaves) == 1
    assert branched.leaves[0].terminal_reason == "integrity_failure_unknown"
    assert branched.leaves[0].certificate.relation is CertifiedRelation.UNKNOWN
    assert attempted == [128, 128]


def test_coverage_audit_rejects_equal_volume_overlap_and_gap() -> None:
    problem = _problem()
    initial = _box(
        ClosedInterval(0.6, 0.8),
        ClosedInterval(0.4, 0.6),
        (
            ClosedInterval(0.0, 0.0),
            ClosedInterval(0.02, 0.04),
            ClosedInterval(-0.03, -0.01),
        ),
    )
    certificate = CertifiedBoxResult(
        CertifiedRelation.UNKNOWN,
        "test",
        128,
        (128,),
        (None, None),
    )
    malformed_left = _box(
        ClosedInterval(0.6, 0.7),
        ClosedInterval(0.4, 0.55),
        tuple(offset.interval for offset in initial.damping_offsets),
    )
    malformed_right = _box(
        ClosedInterval(0.7, 0.8),
        ClosedInterval(0.45, 0.6),
        tuple(offset.interval for offset in initial.damping_offsets),
    )
    leaves = [
        CertifiedLeaf(malformed_left, 1, certificate, "test", (0,)),
        CertifiedLeaf(malformed_right, 1, certificate, "test", (1,)),
    ]

    assert not certified_set._leaves_cover_full_box(problem, initial, leaves)


def test_certified_inner_leaf_has_matching_alpha_inner_and_outer() -> None:
    order = 0.72
    damping = 0.5
    offsets = (0.0, 0.03, -0.02)
    problem = _problem(response_radius=0.08)
    initial = _box(
        ClosedInterval(order - 1.0e-5, order + 1.0e-5),
        ClosedInterval(damping - 1.0e-5, damping + 1.0e-5),
        (
            ClosedInterval(0.0, 0.0),
            ClosedInterval(offsets[1] - 1.0e-5, offsets[1] + 1.0e-5),
            ClosedInterval(offsets[2] - 1.0e-5, offsets[2] + 1.0e-5),
        ),
    )

    result = branch_parameter_box(problem, initial, max_depth=2, max_leaves=8)

    assert len(result.leaves) == 1
    assert result.leaves[0].certificate.relation is CertifiedRelation.INNER_INCLUDED
    assert result.alpha.inner == (initial.order,)
    assert result.alpha.outer == (initial.order,)


def test_main_branch_formula_matches_cmath_at_point_values() -> None:
    for order in (0.4, 0.73, 1.0):
        for frequency in (0.6, 1.0, 3.2):
            expected = (1j * frequency) ** order
            actual = _fractional_frequency(order, frequency)
            assert actual == pytest.approx(expected, abs=2.0e-15)
            assert cmath.phase(actual) == pytest.approx(0.5 * math.pi * order)
