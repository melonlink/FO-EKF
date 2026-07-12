import math

import pytest

from fo_ekf.joint_set import (
    BoxRelation,
    ComplexDisk,
    FeasibilityStatus,
    HarmonicData,
    JointDiskProblem,
    ParameterBox,
    RealInterval,
    branch_parameter_box,
    certify_parameter_box,
    disk_intersection_annulus,
    evaluate_fixed_parameters,
    fractional_frequency,
)


def _responses(
    order: float,
    damping: float,
    morphology: complex,
    frequencies: tuple[float, ...],
    offsets: tuple[float, ...] | None = None,
) -> tuple[complex, ...]:
    if offsets is None:
        offsets = (0.0,) * len(frequencies)
    return tuple(
        morphology / (damping + offset + fractional_frequency(order, frequency))
        for frequency, offset in zip(frequencies, offsets, strict=True)
    )


def _harmonic(
    order: float,
    damping: float,
    morphology: complex,
    frequencies: tuple[float, ...],
    *,
    offsets: tuple[float, ...] | None = None,
    response_radius: float = 0.0,
    drift_radii: tuple[float, ...] | None = None,
    q_margin: float = 0.25,
    label: str = "",
) -> HarmonicData:
    if drift_radii is None:
        drift_radii = (0.0,) * len(frequencies)
    magnitude = abs(morphology)
    return HarmonicData(
        frequencies=frequencies,
        measured_responses=_responses(
            order,
            damping,
            morphology,
            frequencies,
            offsets,
        ),
        response_radii=(response_radius,) * len(frequencies),
        morphology_drift_radii=drift_radii,
        morphology_bounds=((1.0 - q_margin) * magnitude, (1.0 + q_margin) * magnitude),
        label=label,
    )


def test_fixed_parameters_jointly_share_order_damping_and_rate_offsets() -> None:
    order = 0.71
    damping = 0.58
    offsets = (0.06, -0.04, 0.02)
    first = _harmonic(
        order,
        damping,
        1.2 - 0.4j,
        (0.8, 1.3, 1.9),
        offsets=offsets,
        label="m1",
    )
    second = _harmonic(
        order,
        damping,
        -0.45 + 0.8j,
        (1.6, 2.6, 3.8),
        offsets=offsets,
        label="m2",
    )
    problem = JointDiskProblem((first, second))

    result = evaluate_fixed_parameters(
        problem,
        order,
        damping,
        damping_offsets=offsets,
        tolerance=1.0e-9,
    )
    wrong_offsets = evaluate_fixed_parameters(
        problem,
        order,
        damping,
        damping_offsets=(0.0, 0.0, 0.0),
        tolerance=1.0e-9,
    )

    assert result.status is FeasibilityStatus.FEASIBLE
    assert result.morphology_witnesses[0] == pytest.approx(1.2 - 0.4j, abs=1.0e-9)
    assert result.morphology_witnesses[1] == pytest.approx(-0.45 + 0.8j, abs=1.0e-9)
    assert wrong_offsets.status is FeasibilityStatus.INFEASIBLE


def test_independent_complex_morphology_drift_is_an_exact_radius_sum() -> None:
    order = 0.66
    damping = 0.43
    frequencies = (0.75, 1.1, 1.7, 2.2)
    nominal = 0.9 + 0.35j
    drifts = (0.04 - 0.01j, -0.02 + 0.03j, 0.01 + 0.02j, -0.03 - 0.02j)
    responses = tuple(
        (nominal + drift) / (damping + fractional_frequency(order, frequency))
        for frequency, drift in zip(frequencies, drifts, strict=True)
    )
    harmonic = HarmonicData(
        frequencies=frequencies,
        measured_responses=responses,
        response_radii=(0.0,) * len(frequencies),
        morphology_drift_radii=tuple(abs(drift) for drift in drifts),
        morphology_bounds=(0.8, 1.1),
    )

    result = evaluate_fixed_parameters(JointDiskProblem((harmonic,)), order, damping)

    assert result.status is FeasibilityStatus.FEASIBLE
    assert result.morphology_witnesses[0] == pytest.approx(nominal, abs=1.0e-9)


def test_annulus_hole_rejects_a_nonempty_disk_intersection() -> None:
    result = disk_intersection_annulus(
        (ComplexDisk(0.0j, 0.1), ComplexDisk(0.02j, 0.08)),
        (0.5, 1.0),
    )

    assert result.status is FeasibilityStatus.INFEASIBLE
    assert result.reason == "disk_intersection_inside_annulus_hole"
    assert result.maximum_modulus is not None
    assert result.maximum_modulus < 0.5


def test_connected_disk_intersection_gets_a_constructive_annulus_witness() -> None:
    result = disk_intersection_annulus((ComplexDisk(0.0j, 2.0),), (0.9, 1.1))

    assert result.status is FeasibilityStatus.FEASIBLE
    assert result.witness is not None
    assert abs(result.witness) == pytest.approx(1.0, abs=2.0e-10)


def test_pairwise_overlapping_disks_can_have_empty_triple_intersection() -> None:
    side = 1.9
    height = math.sqrt(3.0) * side / 2.0
    disks = (
        ComplexDisk(0.0j, 1.0),
        ComplexDisk(complex(side, 0.0), 1.0),
        ComplexDisk(complex(side / 2.0, height), 1.0),
    )
    assert all(
        abs(left.center - right.center) < left.radius + right.radius
        for index, left in enumerate(disks)
        for right in disks[index + 1 :]
    )

    result = disk_intersection_annulus(disks, (0.1, 10.0))

    assert result.status is FeasibilityStatus.INFEASIBLE
    assert result.reason == "empty_disk_intersection"


def test_small_parameter_box_has_a_uniform_inner_witness() -> None:
    order = 0.74
    damping = 0.52
    offsets = (0.03, -0.02, 0.01)
    problem = JointDiskProblem(
        (
            _harmonic(
                order,
                damping,
                1.0 - 0.2j,
                (0.8, 1.2, 1.8),
                offsets=offsets,
                response_radius=0.08,
                q_margin=0.3,
            ),
            _harmonic(
                order,
                damping,
                0.55 + 0.65j,
                (1.6, 2.4, 3.6),
                offsets=offsets,
                response_radius=0.08,
                q_margin=0.3,
            ),
        )
    )
    box = ParameterBox(
        RealInterval(order - 1.0e-4, order + 1.0e-4),
        RealInterval(damping - 1.0e-4, damping + 1.0e-4),
        tuple(RealInterval(value - 1.0e-4, value + 1.0e-4) for value in offsets),
    )

    certificate = certify_parameter_box(problem, box)

    assert certificate.relation is BoxRelation.INNER_INCLUDED
    assert all(result.status is FeasibilityStatus.FEASIBLE for result in certificate.inner_results)


def test_tolerance_expanded_witness_is_not_promoted_to_inner_inclusion() -> None:
    tolerance = 1.0e-10
    order = 0.7
    damping = 0.5
    frequency = 1.2
    morphology = complex(1.0 - 0.5 * tolerance, 0.0)
    denominator = damping + fractional_frequency(order, frequency)
    harmonic = HarmonicData(
        frequencies=(frequency,),
        measured_responses=(morphology / denominator,),
        response_radii=(0.0,),
        morphology_drift_radii=(0.0,),
        morphology_bounds=(1.0, 2.0),
    )
    problem = JointDiskProblem((harmonic,))
    point_box = ParameterBox(
        RealInterval(order, order),
        RealInterval(damping, damping),
    )

    certificate = certify_parameter_box(problem, point_box, tolerance=tolerance)

    assert certificate.inner_results[0].status is FeasibilityStatus.UNKNOWN
    assert certificate.relation is BoxRelation.UNKNOWN


def test_wrong_zero_width_parameter_box_is_certifiably_excluded() -> None:
    order = 0.69
    damping = 0.47
    harmonic = _harmonic(order, damping, 1.1 + 0.3j, (0.75, 1.25, 2.0))
    problem = JointDiskProblem((harmonic,))
    wrong_box = ParameterBox(
        RealInterval(0.95, 0.95),
        RealInterval(0.8, 0.8),
    )

    certificate = certify_parameter_box(problem, wrong_box, tolerance=1.0e-11)

    assert certificate.relation is BoxRelation.OUTER_EXCLUDED


def test_broad_box_remains_unknown_instead_of_being_called_exact() -> None:
    order = 0.69
    damping = 0.47
    problem = JointDiskProblem((_harmonic(order, damping, 1.1 + 0.3j, (0.75, 1.25, 2.0)),))
    broad_box = ParameterBox(
        RealInterval(0.4, 1.0),
        RealInterval(0.2, 1.0),
    )

    certificate = certify_parameter_box(problem, broad_box)

    assert certificate.relation is BoxRelation.UNKNOWN
    assert certificate.inner_results[0].reason == "empty_uniform_inner_core"


def test_bounded_branching_preserves_unknown_leaves() -> None:
    order = 0.69
    damping = 0.47
    problem = JointDiskProblem((_harmonic(order, damping, 1.1 + 0.3j, (0.75, 1.25, 2.0)),))
    broad_box = ParameterBox(
        RealInterval(0.4, 1.0),
        RealInterval(0.2, 1.0),
    )

    leaves = branch_parameter_box(problem, broad_box, max_depth=2, max_leaves=8)

    assert 1 < len(leaves) <= 8
    assert all(leaf.depth <= 2 for leaf in leaves)
    assert any(leaf.certificate.relation is BoxRelation.UNKNOWN for leaf in leaves)


def test_invalid_drift_box_cannot_cross_nonpositive_rate_damping() -> None:
    problem = JointDiskProblem((_harmonic(0.7, 0.5, 1.0 + 0.2j, (0.8, 1.3)),))
    invalid_box = ParameterBox(
        RealInterval(0.6, 0.8),
        RealInterval(0.1, 0.2),
        (RealInterval(-0.15, 0.0), RealInterval(0.0, 0.1)),
    )

    with pytest.raises(ValueError, match="positive"):
        certify_parameter_box(problem, invalid_box)


def test_parameter_box_enforces_declared_order_and_nominal_damping_domain() -> None:
    with pytest.raises(ValueError, match="exceed one"):
        ParameterBox(RealInterval(0.8, 1.1), RealInterval(0.2, 0.4))

    with pytest.raises(ValueError, match="nominal damping"):
        ParameterBox(RealInterval(0.6, 0.8), RealInterval(-0.1, 0.4))
