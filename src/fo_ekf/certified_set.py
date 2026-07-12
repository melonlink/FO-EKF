"""Arb-certified outer/inner boxes for the joint forward-disk model.

This module is the rigorous counterpart of :mod:`fo_ekf.joint_set`.  Ordinary
floating-point geometry is used only to *propose* a nominal morphology
coefficient.  Every issued ``INNER_INCLUDED`` or ``OUTER_EXCLUDED`` result is
proved again with python-flint ``arb``/``acb`` balls.

``OUTER_EXCLUDED`` first uses three constant-time sufficient certificate
classes: pairwise outer-disk separation, one outer disk lying strictly inside
the annulus hole, or one outer disk lying strictly beyond ``q_max``.  If none
applies, a rigorous Cartesian q-box cover can prove that every point in the
annulus square violates at least one constraint.  Failure to finish that cover
within its fixed resource budget is ``UNKNOWN``, never evidence of feasibility.

The declared model for rate ``r`` and harmonic ``m`` is

``Z[r,m] = (q_bar[m] + delta_q[r,m]) /
            (lambda_bar + delta_lambda[r] + (1j*m*omega[r])**alpha)``.

``delta_q`` is an independent complex disk.  ``delta_lambda`` is real, shared
by all harmonics with the same explicit rate id, and fixed to exact zero at a
pre-registered anchor rate.  The anchor removes the translation redundancy
between nominal damping and the rate offsets.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from enum import Enum
from fractions import Fraction
from threading import RLock
from typing import TypeAlias

import flint
from flint import acb, arb, ctx

from .joint_set import ComplexDisk as FloatDisk
from .joint_set import disk_intersection_annulus as propose_disk_intersection

DEFAULT_PRECISION_SCHEDULE = (128, 256, 512)
REQUIRED_PYTHON_FLINT_VERSION = "0.9.0"
DEFAULT_QBOX_MAX_DEPTH = 18
DEFAULT_QBOX_MAX_NODES = 100_000

if flint.__version__ != REQUIRED_PYTHON_FLINT_VERSION:
    raise RuntimeError(
        "R1-CERT was audited against python-flint "
        f"{REQUIRED_PYTHON_FLINT_VERSION}, but {flint.__version__} is installed"
    )

# python-flint's precision context is process-global.  All public certificate
# issuance through this module is therefore serialized across the full
# precision-escalation schedule.
_CERTIFICATION_LOCK = RLock()


class CertifiedRelation(str, Enum):
    """Rigorous relation between a closed parameter box and feasibility."""

    OUTER_EXCLUDED = "OUTER_EXCLUDED"
    INNER_INCLUDED = "INNER_INCLUDED"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class ClosedInterval:
    """Closed interval whose endpoints are finite Python binary floats."""

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
        return self.lower + 0.5 * (self.upper - self.lower)

    @property
    def width(self) -> float:
        return self.upper - self.lower

    def contains(self, value: float) -> bool:
        return self.lower <= value <= self.upper


@dataclass(frozen=True)
class CertifiedHarmonicData:
    """Rate-indexed observations for one explicit harmonic.

    Frequencies are not accepted independently.  They are derived as
    ``harmonic_index * base_frequencies[r]`` so that a rate id cannot silently
    receive a frequency or damping offset belonging to another rate.
    """

    rate_ids: tuple[str, ...]
    base_frequencies: tuple[float, ...]
    harmonic_index: int
    measured_responses: tuple[complex, ...]
    response_radii: tuple[float, ...]
    morphology_drift_radii: tuple[float, ...]
    morphology_bounds: tuple[float, float]
    label: str = ""

    def __post_init__(self) -> None:
        rate_ids = tuple(str(value) for value in self.rate_ids)
        base_frequencies = tuple(float(value) for value in self.base_frequencies)
        responses = tuple(complex(value) for value in self.measured_responses)
        response_radii = tuple(float(value) for value in self.response_radii)
        drift_radii = tuple(float(value) for value in self.morphology_drift_radii)
        q_min, q_max = (float(value) for value in self.morphology_bounds)
        count = len(rate_ids)
        if count == 0:
            raise ValueError("at least one rate is required")
        if any(not value for value in rate_ids) or len(set(rate_ids)) != count:
            raise ValueError("rate ids must be nonempty and unique")
        if isinstance(self.harmonic_index, bool) or not isinstance(self.harmonic_index, int):
            raise ValueError("harmonic index must be a positive integer")
        if self.harmonic_index < 1:
            raise ValueError("harmonic index must be a positive integer")
        if not (
            len(base_frequencies)
            == len(responses)
            == len(response_radii)
            == len(drift_radii)
            == count
        ):
            raise ValueError("all harmonic arrays must have one entry per rate id")
        if any(not math.isfinite(value) or value <= 0.0 for value in base_frequencies):
            raise ValueError("base frequencies must be finite and positive")
        if any(
            not (math.isfinite(value.real) and math.isfinite(value.imag)) for value in responses
        ):
            raise ValueError("measured responses must be finite")
        if any(not math.isfinite(value) or value < 0.0 for value in response_radii):
            raise ValueError("response radii must be finite and nonnegative")
        if any(not math.isfinite(value) or value < 0.0 for value in drift_radii):
            raise ValueError("morphology drift radii must be finite and nonnegative")
        if not (0.0 < q_min < q_max < math.inf):
            raise ValueError("morphology bounds must satisfy 0 < q_min < q_max")
        object.__setattr__(self, "rate_ids", rate_ids)
        object.__setattr__(self, "base_frequencies", base_frequencies)
        object.__setattr__(self, "measured_responses", responses)
        object.__setattr__(self, "response_radii", response_radii)
        object.__setattr__(self, "morphology_drift_radii", drift_radii)
        object.__setattr__(self, "morphology_bounds", (q_min, q_max))

    @property
    def frequencies(self) -> tuple[float, ...]:
        return tuple(self.harmonic_index * value for value in self.base_frequencies)


@dataclass(frozen=True)
class CertifiedJointProblem:
    """Harmonics with a common canonical rate layout and damping anchor."""

    harmonics: tuple[CertifiedHarmonicData, ...]
    anchor_rate_id: str

    def __post_init__(self) -> None:
        harmonics = tuple(self.harmonics)
        anchor = str(self.anchor_rate_id)
        if not harmonics:
            raise ValueError("at least one harmonic is required")
        canonical_ids = harmonics[0].rate_ids
        canonical_frequencies = harmonics[0].base_frequencies
        if anchor not in canonical_ids:
            raise ValueError("anchor rate id must be present in every harmonic")
        indices: set[int] = set()
        for harmonic in harmonics:
            if harmonic.rate_ids != canonical_ids:
                if set(harmonic.rate_ids) == set(canonical_ids):
                    raise ValueError("rate permutation is forbidden; use the canonical order")
                raise ValueError("all harmonics must use exactly the same rate ids")
            if harmonic.base_frequencies != canonical_frequencies:
                raise ValueError("base frequencies must align exactly across harmonics")
            if harmonic.harmonic_index in indices:
                raise ValueError("harmonic indices must be unique")
            indices.add(harmonic.harmonic_index)
        object.__setattr__(self, "harmonics", harmonics)
        object.__setattr__(self, "anchor_rate_id", anchor)

    @property
    def rate_ids(self) -> tuple[str, ...]:
        return self.harmonics[0].rate_ids

    @property
    def base_frequencies(self) -> tuple[float, ...]:
        return self.harmonics[0].base_frequencies

    @property
    def anchor_index(self) -> int:
        return self.rate_ids.index(self.anchor_rate_id)


@dataclass(frozen=True)
class RateOffsetInterval:
    """Real damping-offset interval associated with an explicit rate id."""

    rate_id: str
    interval: ClosedInterval

    def __post_init__(self) -> None:
        rate_id = str(self.rate_id)
        if not rate_id:
            raise ValueError("rate id must be nonempty")
        object.__setattr__(self, "rate_id", rate_id)


@dataclass(frozen=True)
class CertifiedParameterBox:
    """Closed box in order, nominal damping, and anchored rate offsets."""

    order: ClosedInterval
    damping: ClosedInterval
    damping_offsets: tuple[RateOffsetInterval, ...]

    def __post_init__(self) -> None:
        offsets = tuple(self.damping_offsets)
        if not 0.0 < self.order.lower <= self.order.upper <= 1.0:
            raise ValueError("order box must lie in (0, 1]")
        if self.damping.lower <= 0.0:
            raise ValueError("nominal damping box must be strictly positive")
        ids = tuple(value.rate_id for value in offsets)
        if len(set(ids)) != len(ids):
            raise ValueError("damping-offset rate ids must be unique")
        object.__setattr__(self, "damping_offsets", offsets)


@dataclass(frozen=True)
class CertifiedBoxResult:
    """Final Arb-issued relation for one parameter box."""

    relation: CertifiedRelation
    reason: str
    precision_bits: int
    attempted_precisions: tuple[int, ...]
    morphology_witnesses: tuple[complex | None, ...]


@dataclass(frozen=True)
class CertifiedLeaf:
    """Terminal branch-and-bound leaf."""

    box: CertifiedParameterBox
    depth: int
    certificate: CertifiedBoxResult
    terminal_reason: str
    path: tuple[int, ...]


@dataclass(frozen=True)
class AlphaSandwich:
    """Inner and outer enclosures of the feasible-order projection."""

    inner: tuple[ClosedInterval, ...]
    outer: tuple[ClosedInterval, ...]

    def inner_contains(self, value: float) -> bool:
        return any(interval.contains(value) for interval in self.inner)

    def outer_contains(self, value: float) -> bool:
        return any(interval.contains(value) for interval in self.outer)


@dataclass(frozen=True)
class CertifiedBranchResult:
    """Conservative branch result with leaf coverage and alpha sandwich."""

    leaves: tuple[CertifiedLeaf, ...]
    alpha: AlphaSandwich
    coverage_certified: bool
    budget_exhausted: bool


Ball: TypeAlias = arb
ComplexBall: TypeAlias = acb


@dataclass(frozen=True)
class _BallDisk:
    """Complex center ball with a certified nonnegative scalar radius."""

    center: ComplexBall
    radius: Ball


@dataclass(frozen=True)
class _HarmonicEnclosure:
    outer_disks: tuple[_BallDisk, ...]
    inner_disks: tuple[_BallDisk, ...] | None
    morphology_bounds: tuple[Ball, Ball]


@dataclass(frozen=True)
class _UnitQBox:
    """Exact rational coordinates scaled by the harmonic's q_max.

    Keeping the subdivision coordinates as :class:`fractions.Fraction`
    makes child coverage independent of binary floating-point rounding.  Arb
    performs the final outward-rounded multiplication by q_max.
    """

    real: tuple[Fraction, Fraction]
    imag: tuple[Fraction, Fraction]


def _as_fraction(value: float) -> Fraction:
    return Fraction.from_float(float(value))


def _arb_interval(interval: ClosedInterval) -> Ball:
    """Return an Arb ball containing both Python-float endpoint values."""

    return arb(interval.lower).union(arb(interval.upper))


def _containing_acb(value: complex) -> ComplexBall:
    """Return an ACB ball containing both components of a Python complex."""

    return acb(arb(value.real), arb(value.imag))


def _validate_precision_schedule(schedule: tuple[int, ...]) -> tuple[int, ...]:
    values = tuple(int(value) for value in schedule)
    if not values or any(value < 64 for value in values):
        raise ValueError("at least one precision of 64 bits or more is required")
    if any(right <= left for left, right in zip(values, values[1:], strict=False)):
        raise ValueError("precision schedule must be strictly increasing")
    return values


def _validate_box(problem: CertifiedJointProblem, box: CertifiedParameterBox) -> None:
    offset_ids = tuple(value.rate_id for value in box.damping_offsets)
    if offset_ids != problem.rate_ids:
        if set(offset_ids) == set(problem.rate_ids):
            raise ValueError("offset rate permutation is forbidden; use the canonical order")
        raise ValueError("one damping-offset interval is required for every canonical rate id")
    anchor = box.damping_offsets[problem.anchor_index].interval
    if anchor.lower != 0.0 or anchor.upper != 0.0:
        raise ValueError("anchor-rate damping offset must be exactly [0, 0]")
    damping_lower = _as_fraction(box.damping.lower)
    for offset in box.damping_offsets:
        if damping_lower + _as_fraction(offset.interval.lower) <= 0:
            raise ValueError("every rate-specific damping value must remain positive")


def _principal_fractional_frequency(
    order: Ball,
    base_frequency: float,
    harmonic: int,
) -> ComplexBall:
    """Evaluate the declared principal branch using real Arb formulas.

    For positive ``nu = harmonic * base_frequency`` and ``0 < alpha <= 1``,

    ``(i*nu)**alpha = exp(alpha*log(nu)) *
                       (cos(pi*alpha/2) + i*sin(pi*alpha/2))``.

    No generic complex power or implicit branch selection is used.
    """

    frequency = arb(base_frequency) * harmonic
    amplitude = (order * frequency.log()).exp()
    angle = order * arb.pi() / 2
    return acb(amplitude * angle.cos(), amplitude * angle.sin())


def _rate_denominator(
    order: Ball,
    damping: Ball,
    offset: Ball,
    base_frequency: float,
    harmonic: int,
) -> ComplexBall:
    fractional = _principal_fractional_frequency(order, base_frequency, harmonic)
    return acb(damping + offset) + fractional


def _build_enclosures(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
) -> tuple[_HarmonicEnclosure, ...]:
    order_box = _arb_interval(box.order)
    damping_box = _arb_interval(box.damping)
    order_point = arb(box.order.midpoint)
    damping_point = arb(box.damping.midpoint)
    offset_boxes = tuple(_arb_interval(value.interval) for value in box.damping_offsets)
    offset_points = tuple(arb(value.interval.midpoint) for value in box.damping_offsets)
    enclosures: list[_HarmonicEnclosure] = []

    for harmonic in problem.harmonics:
        outer_disks: list[_BallDisk] = []
        inner_disks: list[_BallDisk] = []
        inner_exists = True
        observations = zip(
            harmonic.base_frequencies,
            harmonic.measured_responses,
            harmonic.response_radii,
            harmonic.morphology_drift_radii,
            offset_boxes,
            offset_points,
            strict=True,
        )
        for (
            base_frequency,
            response_value,
            response_error_value,
            q_drift_value,
            offset_box,
            offset_point,
        ) in observations:
            denominator_box = _rate_denominator(
                order_box,
                damping_box,
                offset_box,
                base_frequency,
                harmonic.harmonic_index,
            )
            denominator_point = _rate_denominator(
                order_point,
                damping_point,
                offset_point,
                base_frequency,
                harmonic.harmonic_index,
            )
            response = _containing_acb(response_value)
            response_error = arb(response_error_value)
            q_drift = arb(q_drift_value)
            center = denominator_point * response

            denominator_deviation = (denominator_box - denominator_point).abs_upper()
            center_deviation_upper = (denominator_deviation * response.abs_upper()).upper()
            forward_radius_upper = (denominator_box.abs_upper() * response_error + q_drift).upper()
            forward_radius_lower = (denominator_box.abs_lower() * response_error + q_drift).lower()
            outer_radius = (center_deviation_upper + forward_radius_upper).upper()
            outer_disks.append(_BallDisk(center=center, radius=outer_radius))

            core_radius = (forward_radius_lower - center_deviation_upper).lower()
            if core_radius < 0:
                inner_exists = False
            else:
                inner_disks.append(_BallDisk(center=center, radius=core_radius))

        q_min, q_max = harmonic.morphology_bounds
        enclosures.append(
            _HarmonicEnclosure(
                outer_disks=tuple(outer_disks),
                inner_disks=tuple(inner_disks) if inner_exists else None,
                morphology_bounds=(arb(q_min), arb(q_max)),
            )
        )
    return tuple(enclosures)


def _fast_outer_exclusion_reason(enclosure: _HarmonicEnclosure) -> str | None:
    """Return one of three constant-time OUTER certificates, if proved.

    The implemented classes are pair separation, strict containment of one
    disk in the annulus hole, and strict exclusion of one disk beyond qmax.
    ``None`` means only that these sufficient tests did not decide the box; it
    is not a complete multidisk-feasibility result.
    """

    disks = enclosure.outer_disks
    q_min, q_max = enclosure.morphology_bounds
    for left_index, left in enumerate(disks):
        for right_index, right in enumerate(disks[left_index + 1 :], start=left_index + 1):
            distance_lower = (left.center - right.center).abs_lower()
            radius_sum_upper = (left.radius + right.radius).upper()
            if distance_lower > radius_sum_upper:
                return f"outer_pair_separated:{left_index}:{right_index}"
    for index, disk in enumerate(disks):
        maximum_modulus = (disk.center.abs_upper() + disk.radius).upper()
        if maximum_modulus < q_min:
            return f"outer_disk_inside_annulus_hole:{index}"
        minimum_modulus = (disk.center.abs_lower() - disk.radius).lower()
        if minimum_modulus > q_max:
            return f"outer_disk_outside_qmax:{index}"
    return None


def _fraction_ball(value: Fraction) -> Ball:
    """Embed one exact rational in an outward-rounded Arb ball."""

    return arb(value.numerator) / value.denominator


def _scaled_q_interval(
    endpoints: tuple[Fraction, Fraction],
    q_max: Ball,
) -> Ball:
    """Enclose a rational unit interval after multiplication by q_max."""

    lower, upper = endpoints
    return (_fraction_ball(lower) * q_max).union(_fraction_ball(upper) * q_max)


def _qbox_is_excluded(
    box: _UnitQBox,
    enclosure: _HarmonicEnclosure,
) -> bool:
    """Prove that every point of one q-box violates a declared constraint."""

    q_min, q_max = enclosure.morphology_bounds
    q = acb(
        _scaled_q_interval(box.real, q_max),
        _scaled_q_interval(box.imag, q_max),
    )
    # Strict tests preserve the closed annulus and disk boundaries.  Each
    # compared value is itself an outward-rounded bound; an indeterminate Arb
    # comparison therefore falls through rather than issuing a certificate.
    if q.abs_lower() > q_max:
        return True
    if q.abs_upper() < q_min:
        return True
    return any((q - disk.center).abs_lower() > disk.radius for disk in enclosure.outer_disks)


def _bisect_qbox(box: _UnitQBox) -> tuple[_UnitQBox, _UnitQBox]:
    """Bisect one exact rational q-box without gaps."""

    real_width = box.real[1] - box.real[0]
    imag_width = box.imag[1] - box.imag[0]
    if real_width >= imag_width:
        midpoint = (box.real[0] + box.real[1]) / 2
        return (
            _UnitQBox((box.real[0], midpoint), box.imag),
            _UnitQBox((midpoint, box.real[1]), box.imag),
        )
    midpoint = (box.imag[0] + box.imag[1]) / 2
    return (
        _UnitQBox(box.real, (box.imag[0], midpoint)),
        _UnitQBox(box.real, (midpoint, box.imag[1])),
    )


def _qbox_outer_exclusion_reason(
    enclosure: _HarmonicEnclosure,
    *,
    max_depth: int = DEFAULT_QBOX_MAX_DEPTH,
    max_nodes: int = DEFAULT_QBOX_MAX_NODES,
) -> str | None:
    """Certify empty outer-disk/annulus intersection by a finite q-box cover.

    The initial square ``[-q_max,q_max]^2`` contains the complete annulus.
    Every removed box is rigorously outside the annulus or at least one outer
    disk.  Consequently, exhausting the exact binary partition proves empty
    intersection.  A depth or node limit returns ``None`` and hence UNKNOWN.
    """

    unit = Fraction(1)
    pending: list[tuple[_UnitQBox, int]] = [(_UnitQBox((-unit, unit), (-unit, unit)), 0)]
    visited = 0
    maximum_depth_seen = 0
    while pending:
        box, depth = pending.pop()
        visited += 1
        maximum_depth_seen = max(maximum_depth_seen, depth)
        if visited > max_nodes:
            return None
        if _qbox_is_excluded(box, enclosure):
            continue
        if depth >= max_depth:
            return None
        left, right = _bisect_qbox(box)
        pending.append((right, depth + 1))
        pending.append((left, depth + 1))
    return f"outer_qbox_cover:{visited}:{maximum_depth_seen}"


def _outer_exclusion_reason(enclosure: _HarmonicEnclosure) -> str | None:
    """Return a rigorous fast-geometry or finite q-box OUTER certificate."""

    fast_reason = _fast_outer_exclusion_reason(enclosure)
    if fast_reason is not None:
        return fast_reason
    return _qbox_outer_exclusion_reason(enclosure)


def _float_value(value: Ball) -> float:
    return float(value.mid())


def _float_candidates(enclosure: _HarmonicEnclosure) -> tuple[complex, ...]:
    """Generate finite float candidates without certificate authority."""

    if enclosure.inner_disks is None:
        return ()
    q_min = _float_value(enclosure.morphology_bounds[0])
    q_max = _float_value(enclosure.morphology_bounds[1])
    disks: list[FloatDisk] = []
    for disk in enclosure.inner_disks:
        center = complex(_float_value(disk.center.real), _float_value(disk.center.imag))
        radius = max(0.0, _float_value(disk.radius))
        disks.append(FloatDisk(center, radius))
    scale = max(
        1.0e-300,
        q_max,
        *(abs(disk.center) + disk.radius for disk in disks),
    )
    tolerance = max(64.0 * math.ulp(scale), scale * 1.0e-14)
    candidates: list[complex] = []
    try:
        proposed = propose_disk_intersection(
            disks,
            (q_min, q_max),
            tolerance=tolerance,
        )
        if proposed.witness is not None:
            candidates.append(complex(proposed.witness))
    except (ArithmeticError, ValueError):
        pass

    raw_points = [disk.center for disk in disks]
    raw_points.append(sum(raw_points, 0.0j) / len(raw_points))
    target = q_min + 0.5 * (q_max - q_min)
    for point in raw_points:
        candidates.append(point)
        magnitude = abs(point)
        if magnitude == 0.0:
            candidates.append(complex(target, 0.0))
        else:
            candidates.append(point * (target / magnitude))

    unique: list[complex] = []
    for candidate in candidates:
        if not (math.isfinite(candidate.real) and math.isfinite(candidate.imag)):
            continue
        if candidate not in unique:
            unique.append(candidate)
    return tuple(unique)


def _verify_inner_candidate(
    candidate: complex,
    enclosure: _HarmonicEnclosure,
) -> bool:
    disks = enclosure.inner_disks
    if disks is None:
        return False
    q = _containing_acb(candidate)
    q_min, q_max = enclosure.morphology_bounds
    if not (q.abs_lower() >= q_min and q.abs_upper() <= q_max):
        return False
    return all((q - disk.center).abs_upper() <= disk.radius for disk in disks)


def _certify_at_current_precision(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    precision: int,
) -> CertifiedBoxResult:
    enclosures = _build_enclosures(problem, box)
    outer_reasons = tuple(_outer_exclusion_reason(enclosure) for enclosure in enclosures)

    witnesses: list[complex | None] = []
    all_inner = True
    for enclosure in enclosures:
        witness = next(
            (
                candidate
                for candidate in _float_candidates(enclosure)
                if _verify_inner_candidate(candidate, enclosure)
            ),
            None,
        )
        witnesses.append(witness)
        if witness is None:
            all_inner = False

    outer_reason = next((value for value in outer_reasons if value is not None), None)
    if outer_reason is not None and all_inner:
        return CertifiedBoxResult(
            relation=CertifiedRelation.UNKNOWN,
            reason="arb_outer_inner_contradiction",
            precision_bits=precision,
            attempted_precisions=(precision,),
            morphology_witnesses=tuple(witnesses),
        )
    if outer_reason is not None:
        return CertifiedBoxResult(
            relation=CertifiedRelation.OUTER_EXCLUDED,
            reason=outer_reason,
            precision_bits=precision,
            attempted_precisions=(precision,),
            morphology_witnesses=tuple(witnesses),
        )
    if all_inner:
        return CertifiedBoxResult(
            relation=CertifiedRelation.INNER_INCLUDED,
            reason="arb_uniform_inner_witnesses",
            precision_bits=precision,
            attempted_precisions=(precision,),
            morphology_witnesses=tuple(witnesses),
        )
    return CertifiedBoxResult(
        relation=CertifiedRelation.UNKNOWN,
        reason="arb_bounds_do_not_decide",
        precision_bits=precision,
        attempted_precisions=(precision,),
        morphology_witnesses=tuple(witnesses),
    )


def certify_parameter_box(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    *,
    precision_schedule: tuple[int, ...] = DEFAULT_PRECISION_SCHEDULE,
) -> CertifiedBoxResult:
    """Certify a box using 128/256/512-bit Arb escalation by default.

    Issuance is serialized because python-flint's precision context is global.
    OUTER first uses the three fast sufficient certificate classes documented
    by :func:`_fast_outer_exclusion_reason`, then a bounded rigorous q-box
    cover.  Undecided multidisk geometry stays UNKNOWN.
    """

    _validate_box(problem, box)
    schedule = _validate_precision_schedule(precision_schedule)
    with _CERTIFICATION_LOCK:
        attempted: list[int] = []
        final: CertifiedBoxResult | None = None
        for precision in schedule:
            attempted.append(precision)
            with ctx.workprec(precision):
                result = _certify_at_current_precision(problem, box, precision)
            final = replace(result, attempted_precisions=tuple(attempted))
            if final.reason == "arb_outer_inner_contradiction":
                # A contradiction means an implementation or backend invariant
                # has failed.  Precision escalation must never turn that same
                # issuance attempt into a trusted certificate.
                return final
            if final.relation is not CertifiedRelation.UNKNOWN:
                return final
        if final is None:  # pragma: no cover - rejected by schedule validation
            raise RuntimeError("precision schedule was unexpectedly empty")
        return final


def _replace_box_coordinate(
    box: CertifiedParameterBox,
    index: int,
    interval: ClosedInterval,
) -> CertifiedParameterBox:
    if index == 0:
        return CertifiedParameterBox(interval, box.damping, box.damping_offsets)
    if index == 1:
        return CertifiedParameterBox(box.order, interval, box.damping_offsets)
    offsets = list(box.damping_offsets)
    selected = offsets[index - 2]
    offsets[index - 2] = RateOffsetInterval(selected.rate_id, interval)
    return CertifiedParameterBox(box.order, box.damping, tuple(offsets))


def bisect_parameter_box(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
) -> tuple[CertifiedParameterBox, CertifiedParameterBox] | None:
    """Bisect the widest normalized non-anchor coordinate."""

    _validate_box(problem, box)
    coordinates = [box.order, box.damping]
    coordinates.extend(value.interval for value in box.damping_offsets)
    scores: list[float] = []
    for index, interval in enumerate(coordinates):
        if index >= 2 and box.damping_offsets[index - 2].rate_id == problem.anchor_rate_id:
            scores.append(-1.0)
            continue
        scale = max(1.0, abs(interval.midpoint))
        scores.append(interval.width / scale)
    index = max(range(len(scores)), key=scores.__getitem__)
    selected = coordinates[index]
    if scores[index] <= 0.0:
        return None
    midpoint = selected.midpoint
    if midpoint in (selected.lower, selected.upper):
        return None
    left = ClosedInterval(selected.lower, midpoint)
    right = ClosedInterval(midpoint, selected.upper)
    return (
        _replace_box_coordinate(box, index, left),
        _replace_box_coordinate(box, index, right),
    )


def _merge_intervals(intervals: list[ClosedInterval]) -> tuple[ClosedInterval, ...]:
    if not intervals:
        return ()
    ordered = sorted(intervals, key=lambda value: (value.lower, value.upper))
    merged: list[ClosedInterval] = [ordered[0]]
    for interval in ordered[1:]:
        previous = merged[-1]
        if interval.lower <= previous.upper:
            merged[-1] = ClosedInterval(previous.lower, max(previous.upper, interval.upper))
        else:
            merged.append(interval)
    return tuple(merged)


def _intervals_are_subset(
    inner: tuple[ClosedInterval, ...],
    outer: tuple[ClosedInterval, ...],
) -> bool:
    return all(
        any(
            container.lower <= value.lower and value.upper <= container.upper for container in outer
        )
        for value in inner
    )


def _leaves_cover_full_box(
    problem: CertifiedJointProblem,
    initial: CertifiedParameterBox,
    leaves: list[CertifiedLeaf],
) -> bool:
    """Replay and audit the exact complete binary partition tree.

    A volume sum alone cannot distinguish a gap from an equal-volume overlap.
    Here every leaf carries its left/right path.  Replaying each path verifies
    its geometry, prefix-freeness rules out ancestor/descendant duplication,
    and exact Kraft equality proves that both branches are represented wherever
    a split occurred.
    """

    paths = [leaf.path for leaf in leaves]
    if len(paths) != len(set(paths)):
        return False
    leaves_by_path = {leaf.path: leaf for leaf in leaves}
    path_set = set(paths)
    for path in paths:
        if any(step not in (0, 1) for step in path):
            return False
        if any(path[:length] in path_set for length in range(len(path))):
            return False

        replayed = initial
        for step in path:
            children = bisect_parameter_box(problem, replayed)
            if children is None:
                return False
            replayed = children[step]
        matching_geometry = replayed == leaves_by_path[path].box
        if not matching_geometry:
            return False

    kraft_sum = sum((Fraction(1, 2 ** len(path)) for path in paths), Fraction(0))
    return kraft_sum == 1


def branch_parameter_box(
    problem: CertifiedJointProblem,
    initial_box: CertifiedParameterBox,
    *,
    max_depth: int = 8,
    max_leaves: int = 4096,
    precision_schedule: tuple[int, ...] = DEFAULT_PRECISION_SCHEDULE,
) -> CertifiedBranchResult:
    """Run conservative branching and return an alpha inner/outer sandwich.

    Every terminal unknown box remains in the outer alpha projection.  Hitting
    a depth or leaf budget never changes its relation from ``UNKNOWN``.
    """

    _validate_box(problem, initial_box)
    schedule = _validate_precision_schedule(precision_schedule)
    if max_depth < 0:
        raise ValueError("max_depth must be nonnegative")
    if max_leaves < 1:
        raise ValueError("max_leaves must be positive")

    pending: list[tuple[CertifiedParameterBox, int, tuple[int, ...]]] = [(initial_box, 0, ())]
    leaves: list[CertifiedLeaf] = []
    budget_exhausted = False
    while pending:
        box, depth, path = pending.pop()
        certificate = certify_parameter_box(problem, box, precision_schedule=schedule)
        split = bisect_parameter_box(problem, box)
        can_add_children = len(leaves) + len(pending) + 2 <= max_leaves
        integrity_failure = certificate.reason == "arb_outer_inner_contradiction"
        if (
            certificate.relation is CertifiedRelation.UNKNOWN
            and split is not None
            and not integrity_failure
        ):
            if depth >= max_depth:
                budget_exhausted = True
                leaves.append(CertifiedLeaf(box, depth, certificate, "max_depth_unknown", path))
            elif not can_add_children:
                budget_exhausted = True
                leaves.append(CertifiedLeaf(box, depth, certificate, "max_leaves_unknown", path))
            else:
                pending.append((split[1], depth + 1, (*path, 1)))
                pending.append((split[0], depth + 1, (*path, 0)))
        else:
            if integrity_failure:
                terminal = "integrity_failure_unknown"
            elif certificate.relation is CertifiedRelation.UNKNOWN:
                terminal = "atomic_unknown"
            else:
                terminal = "certified"
            leaves.append(CertifiedLeaf(box, depth, certificate, terminal, path))

    leaf_alpha_cover = _merge_intervals([leaf.box.order for leaf in leaves])
    alpha_coverage = (
        len(leaf_alpha_cover) == 1
        and leaf_alpha_cover[0].lower == initial_box.order.lower
        and leaf_alpha_cover[0].upper == initial_box.order.upper
    )
    coverage_certified = alpha_coverage and _leaves_cover_full_box(
        problem,
        initial_box,
        leaves,
    )
    if not coverage_certified:
        raise RuntimeError("branch leaves lost coverage of the initial alpha interval")

    alpha_inner = _merge_intervals(
        [
            leaf.box.order
            for leaf in leaves
            if leaf.certificate.relation is CertifiedRelation.INNER_INCLUDED
        ]
    )
    alpha_outer = _merge_intervals(
        [
            leaf.box.order
            for leaf in leaves
            if leaf.certificate.relation is not CertifiedRelation.OUTER_EXCLUDED
        ]
    )
    if not _intervals_are_subset(alpha_inner, alpha_outer):
        raise RuntimeError("certified alpha inner set is not contained in the outer set")
    return CertifiedBranchResult(
        leaves=tuple(leaves),
        alpha=AlphaSandwich(inner=alpha_inner, outer=alpha_outer),
        coverage_certified=True,
        budget_exhausted=budget_exhausted,
    )


__all__ = [
    "AlphaSandwich",
    "CertifiedBoxResult",
    "CertifiedBranchResult",
    "CertifiedHarmonicData",
    "CertifiedJointProblem",
    "CertifiedLeaf",
    "CertifiedParameterBox",
    "CertifiedRelation",
    "ClosedInterval",
    "DEFAULT_PRECISION_SCHEDULE",
    "DEFAULT_QBOX_MAX_DEPTH",
    "DEFAULT_QBOX_MAX_NODES",
    "RateOffsetInterval",
    "bisect_parameter_box",
    "branch_parameter_box",
    "certify_parameter_box",
]
