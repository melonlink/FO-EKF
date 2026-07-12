"""Certified building blocks for the joint forward-disk order set.

The declared scalar harmonic surrogate is

``Z[r, m] = (q_bar[m] + delta_q[r, m]) /
             (lambda_bar + delta_lambda[r] + (1j * nu[r, m]) ** alpha)``.

For fixed ``alpha``, ``lambda_bar``, and real ``delta_lambda[r]``, bounded
response error and *independent complex-disk* morphology drift give one disk
for the shared nominal ``q_bar[m]`` per observation.  Feasibility is exactly
the intersection of those disks with a calibrated annulus, up to the stated
floating-point geometry tolerance.

Unknown real damping drift is deliberately retained as one parameter per
heart-rate index and shared by every harmonic.  It is never absorbed into an
independent complex radius.  In exact arithmetic, the parameter-box routines
below give sufficient outer and inner enclosures.  The implementation uses
ordinary floating-point arithmetic with an explicit geometry tolerance, not
directed-rounding interval arithmetic; its statuses are therefore an auditable
algorithm prototype, not machine-certified propositions.  An UNKNOWN result
is intentional: the enclosures are not claimed to be a complete interval
solver.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum


class FeasibilityStatus(str, Enum):
    """Outcome of fixed-parameter disk/annulus feasibility."""

    FEASIBLE = "FEASIBLE"
    INFEASIBLE = "INFEASIBLE"
    UNKNOWN = "UNKNOWN"


class BoxRelation(str, Enum):
    """Certified relation between a parameter box and the feasible set.

    ``OUTER_EXCLUDED`` means the entire box is outside the feasible set.
    ``INNER_INCLUDED`` means every point in the box is feasible (possibly with
    a different morphology coefficient for each harmonic, although the proof
    implemented here actually finds a constant coefficient over the box).
    """

    OUTER_EXCLUDED = "OUTER_EXCLUDED"
    INNER_INCLUDED = "INNER_INCLUDED"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class ComplexDisk:
    """Closed disk in the complex plane."""

    center: complex
    radius: float

    def __post_init__(self) -> None:
        center = complex(self.center)
        radius = float(self.radius)
        if not (math.isfinite(center.real) and math.isfinite(center.imag)):
            raise ValueError("disk center must be finite")
        if not math.isfinite(radius) or radius < 0.0:
            raise ValueError("disk radius must be finite and nonnegative")
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "radius", radius)


@dataclass(frozen=True)
class HarmonicData:
    """All rate-indexed observations for one retained harmonic.

    ``morphology_drift_radii[r]`` bounds an independent complex variable
    ``delta_q[r, m]``.  The annulus applies to the shared nominal coefficient
    ``q_bar[m]``, not to every drifted coefficient.
    """

    frequencies: tuple[float, ...]
    measured_responses: tuple[complex, ...]
    response_radii: tuple[float, ...]
    morphology_drift_radii: tuple[float, ...]
    morphology_bounds: tuple[float, float]
    label: str = ""

    def __post_init__(self) -> None:
        frequencies = tuple(float(value) for value in self.frequencies)
        responses = tuple(complex(value) for value in self.measured_responses)
        response_radii = tuple(float(value) for value in self.response_radii)
        drift_radii = tuple(float(value) for value in self.morphology_drift_radii)
        q_min, q_max = (float(value) for value in self.morphology_bounds)
        count = len(frequencies)
        if count == 0:
            raise ValueError("at least one heart-rate observation is required")
        if not (len(responses) == len(response_radii) == len(drift_radii) == count):
            raise ValueError("all harmonic observation arrays must have equal length")
        if any(not math.isfinite(value) or value <= 0.0 for value in frequencies):
            raise ValueError("frequencies must be finite and positive")
        if any(
            not (math.isfinite(value.real) and math.isfinite(value.imag)) for value in responses
        ):
            raise ValueError("measured responses must be finite")
        if any(not math.isfinite(value) or value < 0.0 for value in response_radii):
            raise ValueError("response radii must be finite and nonnegative")
        if any(not math.isfinite(value) or value < 0.0 for value in drift_radii):
            raise ValueError("morphology drift radii must be finite and nonnegative")
        if not (0.0 < q_min <= q_max < math.inf):
            raise ValueError("morphology bounds must satisfy 0 < q_min <= q_max")
        object.__setattr__(self, "frequencies", frequencies)
        object.__setattr__(self, "measured_responses", responses)
        object.__setattr__(self, "response_radii", response_radii)
        object.__setattr__(self, "morphology_drift_radii", drift_radii)
        object.__setattr__(self, "morphology_bounds", (q_min, q_max))


@dataclass(frozen=True)
class JointDiskProblem:
    """Harmonics sharing one order and the same rate-indexed damping drift."""

    harmonics: tuple[HarmonicData, ...]

    def __post_init__(self) -> None:
        harmonics = tuple(self.harmonics)
        if not harmonics:
            raise ValueError("at least one harmonic is required")
        rate_count = len(harmonics[0].frequencies)
        if any(len(harmonic.frequencies) != rate_count for harmonic in harmonics):
            raise ValueError("all harmonics must use the same number of rate indices")
        object.__setattr__(self, "harmonics", harmonics)

    @property
    def rate_count(self) -> int:
        return len(self.harmonics[0].frequencies)


@dataclass(frozen=True)
class DiskAnnulusResult:
    """Geometry result for ``intersection(disks) intersect annulus``."""

    status: FeasibilityStatus
    witness: complex | None
    disk_witness: complex | None
    minimum_modulus: float | None
    maximum_modulus: float | None
    reason: str


@dataclass(frozen=True)
class JointFeasibilityResult:
    """Fixed-parameter result across all harmonics."""

    status: FeasibilityStatus
    harmonic_results: tuple[DiskAnnulusResult, ...]
    morphology_witnesses: tuple[complex | None, ...]
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class RealInterval:
    """Closed finite real interval."""

    lower: float
    upper: float

    def __post_init__(self) -> None:
        lower = float(self.lower)
        upper = float(self.upper)
        if not (math.isfinite(lower) and math.isfinite(upper) and lower <= upper):
            raise ValueError("interval endpoints must be finite and ordered")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @property
    def midpoint(self) -> float:
        return 0.5 * (self.lower + self.upper)

    @property
    def radius(self) -> float:
        return 0.5 * (self.upper - self.lower)


@dataclass(frozen=True)
class ParameterBox:
    """Box in order, nominal damping, and shared real rate-drift variables."""

    order: RealInterval
    damping: RealInterval
    damping_offsets: tuple[RealInterval, ...] = ()

    def __post_init__(self) -> None:
        offsets = tuple(self.damping_offsets)
        if self.order.lower <= 0.0:
            raise ValueError("order box must be strictly positive")
        if self.order.upper > 1.0:
            raise ValueError("the declared fractional-order box must not exceed one")
        if self.damping.lower <= 0.0:
            raise ValueError("nominal damping box must be strictly positive")
        object.__setattr__(self, "damping_offsets", offsets)


@dataclass(frozen=True)
class BoxCertificate:
    """Analytic enclosure result from the floating-point prototype."""

    relation: BoxRelation
    outer_results: tuple[DiskAnnulusResult, ...]
    inner_results: tuple[DiskAnnulusResult, ...]
    reason: str


@dataclass(frozen=True)
class BranchLeaf:
    """One leaf produced by bounded-depth branch and bound."""

    box: ParameterBox
    depth: int
    certificate: BoxCertificate


def fractional_frequency(order: float, frequency: float) -> complex:
    """Principal-branch value of ``(1j * frequency) ** order``."""

    if not math.isfinite(order) or not 0.0 < order <= 1.0:
        raise ValueError("order must be finite and in (0, 1]")
    if not math.isfinite(frequency) or frequency <= 0.0:
        raise ValueError("frequency must be finite and positive")
    return frequency**order * cmath.exp(0.5j * math.pi * order)


def forward_morphology_disks(
    harmonic: HarmonicData,
    order: float,
    damping: float,
    *,
    damping_offsets: Sequence[float] | None = None,
) -> tuple[ComplexDisk, ...]:
    """Return the exact nominal-morphology disks at fixed parameters.

    Multiplying ``|Z_hat - Z| <= epsilon`` by the fixed denominator gives
    center ``d * Z_hat`` and radius ``|d| * epsilon``.  The Minkowski sum
    with an independent complex disk ``|delta_q| <= rho_q`` adds ``rho_q``
    exactly.  A real damping offset is shared by all harmonics through its
    common rate index; callers must pass its concrete value here.
    """

    order = float(order)
    damping = float(damping)
    if not math.isfinite(damping) or damping <= 0.0:
        raise ValueError("nominal damping must be finite and positive")
    if damping_offsets is None:
        offsets = (0.0,) * len(harmonic.frequencies)
    else:
        offsets = tuple(float(value) for value in damping_offsets)
    if len(offsets) != len(harmonic.frequencies):
        raise ValueError("one damping offset is required per rate index")
    if any(not math.isfinite(value) for value in offsets):
        raise ValueError("damping offsets must be finite")

    disks: list[ComplexDisk] = []
    for frequency, response, error, q_drift, offset in zip(
        harmonic.frequencies,
        harmonic.measured_responses,
        harmonic.response_radii,
        harmonic.morphology_drift_radii,
        offsets,
        strict=True,
    ):
        rate_damping = damping + offset
        if rate_damping <= 0.0:
            raise ValueError("every rate-specific damping value must be positive")
        denominator = rate_damping + fractional_frequency(order, frequency)
        disks.append(
            ComplexDisk(
                center=denominator * response,
                radius=abs(denominator) * error + q_drift,
            )
        )
    return tuple(disks)


def _disk_residual(point: complex, disk: ComplexDisk) -> float:
    return abs(point - disk.center) - disk.radius


def _inside_all(point: complex, disks: Sequence[ComplexDisk], tolerance: float) -> bool:
    return all(_disk_residual(point, disk) <= tolerance for disk in disks)


def _circle_intersections(
    left: ComplexDisk,
    right: ComplexDisk,
    tolerance: float,
) -> tuple[complex, ...]:
    delta = right.center - left.center
    distance = abs(delta)
    if distance <= tolerance:
        return ()
    if distance > left.radius + right.radius + tolerance:
        return ()
    if distance < abs(left.radius - right.radius) - tolerance:
        return ()

    along = (left.radius**2 - right.radius**2 + distance**2) / (2.0 * distance)
    height_squared = left.radius**2 - along**2
    scale = max(1.0, left.radius**2, right.radius**2, distance**2)
    if height_squared < -tolerance * scale:
        return ()
    height = math.sqrt(max(0.0, height_squared))
    direction = delta / distance
    base = left.center + along * direction
    perpendicular = 1j * direction
    if height <= tolerance:
        return (base,)
    return (base + height * perpendicular, base - height * perpendicular)


def _deduplicate(points: Sequence[complex], tolerance: float) -> tuple[complex, ...]:
    unique: list[complex] = []
    for point in points:
        if not any(abs(point - previous) <= tolerance for previous in unique):
            unique.append(point)
    return tuple(unique)


def _structural_candidates(
    disks: Sequence[ComplexDisk],
    tolerance: float,
) -> tuple[complex, ...]:
    points: list[complex] = [disk.center for disk in disks]
    for left_index, left in enumerate(disks):
        for right in disks[left_index + 1 :]:
            points.extend(_circle_intersections(left, right, tolerance))
    return _deduplicate(points, tolerance)


def _extremum_candidates(
    disks: Sequence[ComplexDisk],
    tolerance: float,
) -> tuple[complex, ...]:
    points = list(_structural_candidates(disks, tolerance))
    points.append(0.0j)
    for disk in disks:
        center_magnitude = abs(disk.center)
        if center_magnitude <= tolerance:
            points.extend((complex(disk.radius, 0.0), complex(-disk.radius, 0.0)))
        else:
            direction = disk.center / center_magnitude
            points.extend(
                (
                    disk.center - disk.radius * direction,
                    disk.center + disk.radius * direction,
                )
            )
    return _deduplicate(points, tolerance)


def _segment_modulus_witness(
    lower_point: complex,
    upper_point: complex,
    target: float,
    tolerance: float,
) -> complex | None:
    lower_value = abs(lower_point) - target
    upper_value = abs(upper_point) - target
    if abs(lower_value) <= tolerance:
        return lower_point
    if abs(upper_value) <= tolerance:
        return upper_point
    if lower_value * upper_value > 0.0:
        return None
    left = lower_point
    right = upper_point
    for _ in range(100):
        middle = 0.5 * (left + right)
        value = abs(middle) - target
        if abs(value) <= tolerance:
            return middle
        if lower_value * value <= 0.0:
            right = middle
            upper_value = value
        else:
            left = middle
            lower_value = value
    return 0.5 * (left + right)


def disk_intersection_annulus(
    disks: Sequence[ComplexDisk],
    morphology_bounds: tuple[float, float],
    *,
    tolerance: float = 1.0e-10,
) -> DiskAnnulusResult:
    """Decide finite disk-intersection/annulus feasibility to a tolerance.

    In exact arithmetic, a nonempty finite disk intersection contains either
    a disk center or a pairwise circle-boundary intersection.  Conditional on
    nonemptiness, extrema of ``|q|`` occur at the origin, a radial extremum of
    one active circle, or a pairwise boundary intersection.  These finite
    candidate sets therefore decide whether the connected convex disk
    intersection meets the closed annulus.

    Candidate harvesting and disk-residual checks use the absolute tolerance;
    the morphology annulus itself is not expanded.  Strict separations larger
    than the tolerance produce ``INFEASIBLE``.  Numerically unresolved
    near-contact configurations return ``UNKNOWN``.  Because these are
    mixed floating guards rather than directed rounding, callers must not
    promote a fixed-parameter status into a machine certificate.
    """

    values = tuple(disks)
    if not values:
        raise ValueError("at least one disk is required")
    q_min, q_max = (float(value) for value in morphology_bounds)
    if not (0.0 < q_min <= q_max < math.inf):
        raise ValueError("morphology bounds must satisfy 0 < q_min <= q_max")
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be finite and positive")

    for left_index, left in enumerate(values):
        for right in values[left_index + 1 :]:
            if abs(left.center - right.center) > left.radius + right.radius + tolerance:
                return DiskAnnulusResult(
                    status=FeasibilityStatus.INFEASIBLE,
                    witness=None,
                    disk_witness=None,
                    minimum_modulus=None,
                    maximum_modulus=None,
                    reason="strictly_disjoint_disk_pair",
                )

    structural = _structural_candidates(values, tolerance)
    disk_points = tuple(point for point in structural if _inside_all(point, values, tolerance))
    if not disk_points:
        best_violation = min(
            max(_disk_residual(point, disk) for disk in values) for point in structural
        )
        status = (
            FeasibilityStatus.INFEASIBLE
            if best_violation > 10.0 * tolerance
            else FeasibilityStatus.UNKNOWN
        )
        return DiskAnnulusResult(
            status=status,
            witness=None,
            disk_witness=None,
            minimum_modulus=None,
            maximum_modulus=None,
            reason=(
                "empty_disk_intersection"
                if status is FeasibilityStatus.INFEASIBLE
                else "near_contact_disk_geometry"
            ),
        )

    extrema = tuple(
        point
        for point in _extremum_candidates(values, tolerance)
        if _inside_all(point, values, tolerance)
    )
    if not extrema:
        return DiskAnnulusResult(
            status=FeasibilityStatus.UNKNOWN,
            witness=None,
            disk_witness=disk_points[0],
            minimum_modulus=None,
            maximum_modulus=None,
            reason="unresolved_modulus_extrema",
        )

    minimum_point = min(extrema, key=abs)
    maximum_point = max(extrema, key=abs)
    minimum_modulus = abs(minimum_point)
    maximum_modulus = abs(maximum_point)
    disk_witness = disk_points[0]

    for point in extrema:
        magnitude = abs(point)
        if _inside_all(point, values, tolerance) and q_min <= magnitude <= q_max:
            return DiskAnnulusResult(
                status=FeasibilityStatus.FEASIBLE,
                witness=point,
                disk_witness=disk_witness,
                minimum_modulus=minimum_modulus,
                maximum_modulus=maximum_modulus,
                reason="checked_geometry_witness",
            )

    if minimum_modulus > q_max + tolerance:
        return DiskAnnulusResult(
            status=FeasibilityStatus.INFEASIBLE,
            witness=None,
            disk_witness=disk_witness,
            minimum_modulus=minimum_modulus,
            maximum_modulus=maximum_modulus,
            reason="disk_intersection_outside_outer_morphology_radius",
        )
    if maximum_modulus < q_min - tolerance:
        return DiskAnnulusResult(
            status=FeasibilityStatus.INFEASIBLE,
            witness=None,
            disk_witness=disk_witness,
            minimum_modulus=minimum_modulus,
            maximum_modulus=maximum_modulus,
            reason="disk_intersection_inside_annulus_hole",
        )

    overlap_lower = max(minimum_modulus, q_min)
    overlap_upper = min(maximum_modulus, q_max)
    target = 0.5 * (overlap_lower + overlap_upper)
    witness = _segment_modulus_witness(
        minimum_point,
        maximum_point,
        target,
        tolerance,
    )
    if witness is not None and _inside_all(witness, values, 2.0 * tolerance):
        if q_min <= abs(witness) <= q_max:
            return DiskAnnulusResult(
                status=FeasibilityStatus.FEASIBLE,
                witness=witness,
                disk_witness=disk_witness,
                minimum_modulus=minimum_modulus,
                maximum_modulus=maximum_modulus,
                reason="connected_intersection_modulus_witness",
            )
    return DiskAnnulusResult(
        status=FeasibilityStatus.UNKNOWN,
        witness=None,
        disk_witness=disk_witness,
        minimum_modulus=minimum_modulus,
        maximum_modulus=maximum_modulus,
        reason="annulus_contact_within_geometry_tolerance",
    )


def evaluate_fixed_parameters(
    problem: JointDiskProblem,
    order: float,
    damping: float,
    *,
    damping_offsets: Sequence[float] | None = None,
    tolerance: float = 1.0e-10,
) -> JointFeasibilityResult:
    """Evaluate the fixed-parameter forward-disk predicate numerically.

    A concrete real offset is shared by all harmonics at each rate index.
    Unknown offset intervals belong in :func:`certify_parameter_box` instead.
    """

    if damping_offsets is None:
        offsets = (0.0,) * problem.rate_count
    else:
        offsets = tuple(float(value) for value in damping_offsets)
    if len(offsets) != problem.rate_count:
        raise ValueError("one damping offset is required per rate index")

    results = tuple(
        disk_intersection_annulus(
            forward_morphology_disks(
                harmonic,
                order,
                damping,
                damping_offsets=offsets,
            ),
            harmonic.morphology_bounds,
            tolerance=tolerance,
        )
        for harmonic in problem.harmonics
    )
    if any(result.status is FeasibilityStatus.INFEASIBLE for result in results):
        status = FeasibilityStatus.INFEASIBLE
    elif all(result.status is FeasibilityStatus.FEASIBLE for result in results):
        status = FeasibilityStatus.FEASIBLE
    else:
        status = FeasibilityStatus.UNKNOWN
    return JointFeasibilityResult(
        status=status,
        harmonic_results=results,
        morphology_witnesses=tuple(result.witness for result in results),
        reasons=tuple(result.reason for result in results),
    )


def _validate_box(problem: JointDiskProblem, box: ParameterBox) -> tuple[RealInterval, ...]:
    if box.damping_offsets:
        offsets = box.damping_offsets
    else:
        offsets = (RealInterval(0.0, 0.0),) * problem.rate_count
    if len(offsets) != problem.rate_count:
        raise ValueError("one damping-offset interval is required per rate index")
    if any(box.damping.lower + offset.lower <= 0.0 for offset in offsets):
        raise ValueError("the full box must keep every rate-specific damping positive")
    return offsets


def _fractional_frequency_variation(order: RealInterval, frequency: float) -> float:
    """Mean-value enclosure around the order midpoint."""

    derivative_factor = math.hypot(math.log(frequency), 0.5 * math.pi)
    endpoint_magnitude = max(frequency**order.lower, frequency**order.upper)
    return order.radius * derivative_factor * endpoint_magnitude


def _has_guarded_inner_witness(
    result: DiskAnnulusResult,
    disks: Sequence[ComplexDisk],
    morphology_bounds: tuple[float, float],
    tolerance: float,
) -> bool:
    """Require an inward numerical margin before promoting inner inclusion.

    The geometry routine uses tolerance-expanded candidate tests to avoid
    losing near-boundary configurations.  Such a witness is unsuitable for
    the stronger claim that one core lies inside every member of a parameter
    box.  This guard prevents promotion unless the witness lies inside every
    original constraint by at least tolerance.  It is a defensive check, not
    a substitute for directed rounding.
    """

    witness = result.witness
    if result.status is not FeasibilityStatus.FEASIBLE or witness is None:
        return False
    q_min, q_max = morphology_bounds
    magnitude = abs(witness)
    return q_min + tolerance <= magnitude <= q_max - tolerance and all(
        _disk_residual(witness, disk) <= -tolerance for disk in disks
    )


def certify_parameter_box(
    problem: JointDiskProblem,
    box: ParameterBox,
    *,
    tolerance: float = 1.0e-10,
) -> BoxCertificate:
    """Return an outer/inner/unknown prototype classification for a box.

    The outer disks contain the union of each forward disk over the box.  If
    even their intersection misses the annulus, no point of the parameter box
    can be feasible.  The inner cores are contained in every corresponding
    forward disk.  A core-intersection witness therefore proves every point in
    the box feasible.  Failure of either sufficient test is only ``UNKNOWN``.

    The analytic enclosures use a global mean-value bound in alpha.  The
    floating implementation adds a strict guard before promoting inner
    inclusion, but outer exclusion still lacks directed-rounding error bounds.
    Consequently this is not a machine certificate, completeness theorem, or
    converged projection algorithm.
    """

    offsets = _validate_box(problem, box)
    outer_results: list[DiskAnnulusResult] = []
    inner_results: list[DiskAnnulusResult] = []
    guarded_inner_witnesses: list[bool] = []
    inner_core_empty = False
    for harmonic in problem.harmonics:
        # Construct inline so an empty core can be represented without an
        # invalid numerical disk center.
        outer_disks: list[ComplexDisk] = []
        inner_disks: list[ComplexDisk] = []
        harmonic_core_empty = False
        for frequency, response, error, q_drift, offset in zip(
            harmonic.frequencies,
            harmonic.measured_responses,
            harmonic.response_radii,
            harmonic.morphology_drift_radii,
            offsets,
            strict=True,
        ):
            denominator_midpoint = (
                box.damping.midpoint
                + offset.midpoint
                + fractional_frequency(box.order.midpoint, frequency)
            )
            denominator_variation = (
                box.damping.radius
                + offset.radius
                + _fractional_frequency_variation(box.order, frequency)
            )
            center_midpoint = denominator_midpoint * response
            center_variation = abs(response) * denominator_variation
            denominator_lower = max(0.0, abs(denominator_midpoint) - denominator_variation)
            denominator_upper = abs(denominator_midpoint) + denominator_variation
            radius_lower = denominator_lower * error + q_drift
            radius_upper = denominator_upper * error + q_drift
            outer_disks.append(ComplexDisk(center_midpoint, center_variation + radius_upper))
            core_radius = radius_lower - center_variation
            if core_radius < 0.0:
                harmonic_core_empty = True
            else:
                inner_disks.append(ComplexDisk(center_midpoint, core_radius))

        outer_result = disk_intersection_annulus(
            outer_disks,
            harmonic.morphology_bounds,
            tolerance=tolerance,
        )
        outer_results.append(outer_result)
        if harmonic_core_empty:
            inner_core_empty = True
            guarded_inner_witnesses.append(False)
            inner_results.append(
                DiskAnnulusResult(
                    status=FeasibilityStatus.UNKNOWN,
                    witness=None,
                    disk_witness=None,
                    minimum_modulus=None,
                    maximum_modulus=None,
                    reason="empty_uniform_inner_core",
                )
            )
        else:
            inner_result = disk_intersection_annulus(
                inner_disks,
                harmonic.morphology_bounds,
                tolerance=tolerance,
            )
            inner_results.append(inner_result)
            guarded_inner_witnesses.append(
                _has_guarded_inner_witness(
                    inner_result,
                    inner_disks,
                    harmonic.morphology_bounds,
                    tolerance,
                )
            )

    outer_excluded = any(result.status is FeasibilityStatus.INFEASIBLE for result in outer_results)
    inner_included = (
        not inner_core_empty
        and all(guarded_inner_witnesses)
        and all(result.status is FeasibilityStatus.FEASIBLE for result in inner_results)
    )
    if outer_excluded and inner_included:
        relation = BoxRelation.UNKNOWN
        reason = "contradictory_float_outer_inner_results"
    elif outer_excluded:
        relation = BoxRelation.OUTER_EXCLUDED
        reason = "an_outer_union_disk_set_is_infeasible"
    elif inner_included:
        relation = BoxRelation.INNER_INCLUDED
        reason = "constant_morphology_witnesses_cover_the_entire_box"
    else:
        relation = BoxRelation.UNKNOWN
        reason = "outer_and_inner_sufficient_tests_do_not_close_the_box"
    return BoxCertificate(
        relation=relation,
        outer_results=tuple(outer_results),
        inner_results=tuple(inner_results),
        reason=reason,
    )


def bisect_parameter_box(box: ParameterBox) -> tuple[ParameterBox, ParameterBox] | None:
    """Bisect the widest absolute-width coordinate of a parameter box."""

    coordinates = [box.order, box.damping, *box.damping_offsets]
    widths = [coordinate.upper - coordinate.lower for coordinate in coordinates]
    index = max(range(len(widths)), key=widths.__getitem__)
    if widths[index] == 0.0:
        return None
    selected = coordinates[index]
    midpoint = selected.midpoint
    left_interval = RealInterval(selected.lower, midpoint)
    right_interval = RealInterval(midpoint, selected.upper)

    def replace(interval: RealInterval) -> ParameterBox:
        if index == 0:
            return ParameterBox(interval, box.damping, box.damping_offsets)
        if index == 1:
            return ParameterBox(box.order, interval, box.damping_offsets)
        offsets = list(box.damping_offsets)
        offsets[index - 2] = interval
        return ParameterBox(box.order, box.damping, tuple(offsets))

    return replace(left_interval), replace(right_interval)


def branch_parameter_box(
    problem: JointDiskProblem,
    initial_box: ParameterBox,
    *,
    max_depth: int = 8,
    max_leaves: int = 4096,
    tolerance: float = 1.0e-10,
) -> tuple[BranchLeaf, ...]:
    """Run a bounded-depth conservative branch-and-bound skeleton.

    Unknown leaves remain unknown.  In particular, exhausting ``max_depth``
    or ``max_leaves`` is not evidence for feasibility or infeasibility.
    """

    if max_depth < 0:
        raise ValueError("max_depth must be nonnegative")
    if max_leaves < 1:
        raise ValueError("max_leaves must be positive")
    pending: list[tuple[ParameterBox, int]] = [(initial_box, 0)]
    leaves: list[BranchLeaf] = []
    while pending:
        box, depth = pending.pop()
        certificate = certify_parameter_box(problem, box, tolerance=tolerance)
        split = bisect_parameter_box(box)
        can_add_children = len(leaves) + len(pending) + 2 <= max_leaves
        if (
            certificate.relation is BoxRelation.UNKNOWN
            and depth < max_depth
            and split is not None
            and can_add_children
        ):
            pending.append((split[1], depth + 1))
            pending.append((split[0], depth + 1))
        else:
            leaves.append(BranchLeaf(box=box, depth=depth, certificate=certificate))
    return tuple(leaves)
