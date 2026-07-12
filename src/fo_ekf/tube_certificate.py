"""Replayable parametric Gram/Krawczyk inner certificates.

This module proves an *inner* statement only.  For every fractional order in
one closed interval, a continuously selected damping/morphology tuple exists
and satisfies all declared forward disks and morphology annuli.  Failure of
the selector, ball arithmetic, or a resource limit is ``UNKNOWN`` and is
never an outer exclusion.

For fixed order, stack the real and imaginary parts of

``q[m] - Z[r,m] * (lambda_bar + delta_lambda[r] + (1j*nu[r,m])**alpha)``.

The result is ``g(alpha, x) = B*x - c(alpha)`` in the free damping
coordinates and the real/imaginary morphology coordinates.  The stationary
equations of the unweighted central least-squares selector are

``F(alpha, x) = B.T * g = H*x - B.T*c(alpha) = 0``.

``H = B.T*B`` is an exact rational matrix independent of ``alpha``.  When it
is nonsingular, the certificate stores the exact rational inverse ``C``.
Consequently ``I-C*H`` is exactly zero: the parametric Krawczyk map is the
constant map ``C*B.T*c(alpha)``.  Arb encloses that tube, while the original
disk/annulus inequalities -- not the least-squares equations -- carry the
physical feasibility claim.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from typing import Any

import flint
from flint import acb, arb, ctx, fmpq, fmpq_mat

from .certified_set import (
    _CERTIFICATION_LOCK,
    DEFAULT_PRECISION_SCHEDULE,
    CertifiedJointProblem,
    CertifiedParameterBox,
    ClosedInterval,
    _principal_fractional_frequency,
    _validate_box,
    _validate_precision_schedule,
)

SCHEMA = "fo-ekf.parametric-gram-krawczyk.v1"
ALGORITHM = "central-unweighted-gram-selector-v1"
MAX_SELECTOR_DIMENSION = 32
MAX_RESIDUAL_ROWS = 512
MAX_RATIONAL_BITS = 16_384
MAX_RATIONAL_CHARACTERS = 20_000
MAX_CERTIFICATE_BYTES = 8_000_000
MAX_REPLAY_PRECISION = 4096
MAX_ALPHA_DEPTH = 64
MAX_ALPHA_LEAVES = 4096
_SUBDIVIDABLE_UNKNOWN_REASON = "arb_parametric_tube_undecided"
_RESOURCE_UNKNOWN_REASONS = frozenset(
    {
        "precision_resource_limit",
        "certificate_resource_bytes",
        "selector_inverse_resource_bits",
        "selector_resource_dimension",
    }
)


class _CertificateResourceExceeded(RuntimeError):
    pass


class TubeRelation(str, Enum):
    """Only sound inner relations are expressible by this module."""

    ROBUST_INNER = "ROBUST_INNER"
    CLOSED_INNER = "CLOSED_INNER"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class TubeCertificate:
    """A sealed JSON certificate or an auditable UNKNOWN result."""

    relation: TubeRelation
    reason: str
    precision_bits: int
    attempted_precisions: tuple[int, ...]
    input_sha256: str
    certificate_sha256: str
    certificate_json: str

    def to_json(self) -> str:
        return self.certificate_json


@dataclass(frozen=True)
class TubeReplayResult:
    """Independent replay result; invalid material never becomes an inner."""

    valid: bool
    relation: TubeRelation
    reason: str


@dataclass(frozen=True)
class AlphaTubeLeaf:
    """One terminal alpha-only subdivision leaf with its own replay proof."""

    box: CertifiedParameterBox
    depth: int
    path: tuple[int, ...]
    certificate: TubeCertificate
    terminal_reason: str


@dataclass(frozen=True)
class AlphaTubeBranchResult:
    """Gap-free inner/UNKNOWN cover of the initial order interval."""

    leaves: tuple[AlphaTubeLeaf, ...]
    certified_inner: tuple[ClosedInterval, ...]
    robust_inner: tuple[ClosedInterval, ...]
    closed_inner: tuple[ClosedInterval, ...]
    unknown: tuple[ClosedInterval, ...]
    coverage_certified: bool
    budget_exhausted: bool
    max_depth: int
    max_leaves: int


@dataclass(frozen=True)
class AlphaTubeBranchReplayResult:
    """Audit result for the partition and every leaf JSON proof."""

    valid: bool
    reason: str
    leaf_replays: tuple[TubeReplayResult, ...]


@dataclass(frozen=True)
class _Selector:
    variable_names: tuple[str, ...]
    variable_box: tuple[tuple[Fraction, Fraction], ...]
    free_damping_columns: tuple[tuple[str, int], ...]
    fixed_rate_damping: tuple[Fraction, ...]
    q_columns: tuple[tuple[int, int], ...]
    b: tuple[tuple[Fraction, ...], ...]
    c_inverse: tuple[tuple[Fraction, ...], ...]
    selector_map: tuple[tuple[Fraction, ...], ...]
    residual_map: tuple[tuple[Fraction, ...], ...]


def _fraction(value: float) -> Fraction:
    return Fraction.from_float(float(value))


def _rat(value: Fraction) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _parse_rat(value: Any) -> Fraction:
    if not isinstance(value, str) or not value:
        raise ValueError("rational values must be nonempty strings")
    if len(value) > MAX_RATIONAL_CHARACTERS:
        raise ValueError("rational string exceeds the replay resource limit")
    if value.strip() != value or value.startswith("+"):
        raise ValueError("rational strings must be canonical")
    try:
        parsed = Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("invalid rational string") from exc
    if max(abs(parsed.numerator).bit_length(), parsed.denominator.bit_length()) > MAX_RATIONAL_BITS:
        raise ValueError("rational value exceeds the replay bit limit")
    if _rat(parsed) != value:
        raise ValueError("rational strings must be reduced and canonical")
    return parsed


def _arb_fraction(value: Fraction) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def _arb_fraction_interval(endpoints: tuple[Fraction, Fraction]) -> arb:
    return _arb_fraction(endpoints[0]).union(_arb_fraction(endpoints[1]))


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _problem_manifest(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
) -> dict[str, Any]:
    return {
        "model": "anchored-forward-disks-principal-branch-v1",
        "anchor_rate_id": problem.anchor_rate_id,
        "rate_ids": list(problem.rate_ids),
        "harmonics": [
            {
                "label": harmonic.label,
                "harmonic_index": harmonic.harmonic_index,
                "base_frequencies": [_rat(_fraction(value)) for value in harmonic.base_frequencies],
                "measured_responses": [
                    [_rat(_fraction(value.real)), _rat(_fraction(value.imag))]
                    for value in harmonic.measured_responses
                ],
                "response_radii": [_rat(_fraction(value)) for value in harmonic.response_radii],
                "morphology_drift_radii": [
                    _rat(_fraction(value)) for value in harmonic.morphology_drift_radii
                ],
                "morphology_bounds": [
                    _rat(_fraction(harmonic.morphology_bounds[0])),
                    _rat(_fraction(harmonic.morphology_bounds[1])),
                ],
            }
            for harmonic in problem.harmonics
        ],
        "parameter_box": {
            "order": [_rat(_fraction(box.order.lower)), _rat(_fraction(box.order.upper))],
            "damping": [
                _rat(_fraction(box.damping.lower)),
                _rat(_fraction(box.damping.upper)),
            ],
            "damping_offsets": [
                {
                    "rate_id": value.rate_id,
                    "interval": [
                        _rat(_fraction(value.interval.lower)),
                        _rat(_fraction(value.interval.upper)),
                    ],
                }
                for value in box.damping_offsets
            ],
        },
    }


def _matrix_transpose(matrix: tuple[tuple[Fraction, ...], ...]) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(row[index] for row in matrix) for index in range(len(matrix[0])))


def _matrix_product(
    left: tuple[tuple[Fraction, ...], ...],
    right: tuple[tuple[Fraction, ...], ...],
) -> tuple[tuple[Fraction, ...], ...]:
    right_t = _matrix_transpose(right)
    return tuple(
        tuple(
            sum((a * b for a, b in zip(row, column, strict=True)), Fraction(0))
            for column in right_t
        )
        for row in left
    )


def _matrix_inverse(
    matrix: tuple[tuple[Fraction, ...], ...],
) -> tuple[tuple[Fraction, ...], ...] | None:
    values = fmpq_mat(
        [[fmpq(value.numerator, value.denominator) for value in row] for row in matrix]
    )
    if values.det() == 0:
        return None
    inverse = values.inv().tolist()
    return tuple(
        tuple(Fraction(int(value.numerator), int(value.denominator)) for value in row)
        for row in inverse
    )


def _identity(size: int) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Fraction(int(row == column)) for column in range(size)) for row in range(size)
    )


def _build_selector(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
) -> tuple[_Selector | None, str]:
    _validate_box(problem, box)
    variable_names: list[str] = []
    variable_box: list[tuple[Fraction, Fraction]] = []
    free_damping: list[tuple[str, int]] = []

    if box.damping.width > 0.0:
        free_damping.append(("damping", len(variable_names)))
        variable_names.append("damping")
        variable_box.append((_fraction(box.damping.lower), _fraction(box.damping.upper)))

    for value in box.damping_offsets:
        if value.rate_id == problem.anchor_rate_id:
            continue
        if value.interval.width > 0.0:
            free_damping.append((f"offset:{value.rate_id}", len(variable_names)))
            variable_names.append(f"offset:{value.rate_id}")
            variable_box.append((_fraction(value.interval.lower), _fraction(value.interval.upper)))

    fixed_rate_damping: list[Fraction] = []
    damping_is_free = box.damping.width > 0.0
    free_offset_ids = {
        name.removeprefix("offset:") for name, _ in free_damping if name.startswith("offset:")
    }
    for offset in box.damping_offsets:
        fixed = Fraction(0) if damping_is_free else _fraction(box.damping.lower)
        if offset.rate_id not in free_offset_ids:
            fixed += _fraction(offset.interval.lower)
        fixed_rate_damping.append(fixed)

    q_columns: list[tuple[int, int]] = []
    for harmonic in problem.harmonics:
        q_max = _fraction(harmonic.morphology_bounds[1])
        real_index = len(variable_names)
        imag_index = real_index + 1
        q_columns.append((real_index, imag_index))
        variable_names.extend(
            (
                f"q_re:m{harmonic.harmonic_index}",
                f"q_im:m{harmonic.harmonic_index}",
            )
        )
        variable_box.extend(((-q_max, q_max), (-q_max, q_max)))

    dimension = len(variable_names)
    row_count = 2 * len(problem.rate_ids) * len(problem.harmonics)
    if dimension > MAX_SELECTOR_DIMENSION or row_count > MAX_RESIDUAL_ROWS:
        return None, "selector_resource_dimension"
    if row_count < dimension:
        return None, "selector_underdetermined"

    free_column_by_name = dict(free_damping)
    rows: list[tuple[Fraction, ...]] = []
    for harmonic_index, harmonic in enumerate(problem.harmonics):
        q_real, q_imag = q_columns[harmonic_index]
        for _rate_index, (rate_id, response) in enumerate(
            zip(problem.rate_ids, harmonic.measured_responses, strict=True)
        ):
            real = [Fraction(0) for _ in range(dimension)]
            imag = [Fraction(0) for _ in range(dimension)]
            z_real = _fraction(response.real)
            z_imag = _fraction(response.imag)
            if "damping" in free_column_by_name:
                column = free_column_by_name["damping"]
                real[column] -= z_real
                imag[column] -= z_imag
            offset_name = f"offset:{rate_id}"
            if offset_name in free_column_by_name:
                column = free_column_by_name[offset_name]
                real[column] -= z_real
                imag[column] -= z_imag
            real[q_real] = Fraction(1)
            imag[q_imag] = Fraction(1)
            rows.extend((tuple(real), tuple(imag)))

    b = tuple(rows)
    b_t = _matrix_transpose(b)
    hessian = _matrix_product(b_t, b)
    c_inverse = _matrix_inverse(hessian)
    if c_inverse is None:
        return None, "selector_rank_deficient"
    if any(
        max(abs(value.numerator).bit_length(), value.denominator.bit_length()) > MAX_RATIONAL_BITS
        for row in c_inverse
        for value in row
    ):
        return None, "selector_inverse_resource_bits"
    if _matrix_product(c_inverse, hessian) != _identity(dimension):
        return None, "selector_inverse_integrity"

    selector_map = _matrix_product(c_inverse, b_t)
    projection = _matrix_product(b, selector_map)
    residual_map = tuple(
        tuple(value - Fraction(int(row == column)) for column, value in enumerate(values))
        for row, values in enumerate(projection)
    )
    return (
        _Selector(
            variable_names=tuple(variable_names),
            variable_box=tuple(variable_box),
            free_damping_columns=tuple(free_damping),
            fixed_rate_damping=tuple(fixed_rate_damping),
            q_columns=tuple(q_columns),
            b=b,
            c_inverse=c_inverse,
            selector_map=selector_map,
            residual_map=residual_map,
        ),
        "ok",
    )


def _linear_ball_map(
    matrix: tuple[tuple[Fraction, ...], ...],
    vector: tuple[arb, ...],
) -> tuple[arb, ...]:
    results: list[arb] = []
    for row in matrix:
        total = arb(0)
        for coefficient, value in zip(row, vector, strict=True):
            if coefficient:
                total += _arb_fraction(coefficient) * value
        results.append(total)
    return tuple(results)


def _central_vector(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    selector: _Selector,
) -> tuple[arb, ...]:
    order = _arb_fraction_interval((_fraction(box.order.lower), _fraction(box.order.upper)))
    values: list[arb] = []
    for harmonic in problem.harmonics:
        for rate_index, (frequency, response) in enumerate(
            zip(harmonic.base_frequencies, harmonic.measured_responses, strict=True)
        ):
            fractional = _principal_fractional_frequency(
                order,
                frequency,
                harmonic.harmonic_index,
            )
            fixed = _arb_fraction(selector.fixed_rate_damping[rate_index])
            z = acb(
                _arb_fraction(_fraction(response.real)),
                _arb_fraction(_fraction(response.imag)),
            )
            center = z * (acb(fixed) + fractional)
            values.extend((center.real, center.imag))
    return tuple(values)


def _physical_rate_dampings(
    problem: CertifiedJointProblem,
    selector: _Selector,
    selected: tuple[arb, ...],
) -> tuple[arb, ...]:
    free = dict(selector.free_damping_columns)
    values: list[arb] = []
    for rate_index, rate_id in enumerate(problem.rate_ids):
        value = _arb_fraction(selector.fixed_rate_damping[rate_index])
        if "damping" in free:
            value += selected[free["damping"]]
        offset_name = f"offset:{rate_id}"
        if offset_name in free:
            value += selected[free[offset_name]]
        values.append(value)
    return tuple(values)


def _constraint_balls(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    selector: _Selector,
    central: tuple[arb, ...],
    selected: tuple[arb, ...],
) -> tuple[tuple[str, arb], ...]:
    residual = _linear_ball_map(selector.residual_map, central)
    order = _arb_fraction_interval((_fraction(box.order.lower), _fraction(box.order.upper)))
    rate_dampings = _physical_rate_dampings(problem, selector, selected)
    margins: list[tuple[str, arb]] = []
    row = 0
    for harmonic_index, harmonic in enumerate(problem.harmonics):
        q_real_index, q_imag_index = selector.q_columns[harmonic_index]
        q = acb(selected[q_real_index], selected[q_imag_index])
        q_min = _arb_fraction(_fraction(harmonic.morphology_bounds[0]))
        q_max = _arb_fraction(_fraction(harmonic.morphology_bounds[1]))
        margins.append((f"annulus_lower:m{harmonic.harmonic_index}", q.abs_lower() - q_min))
        margins.append((f"annulus_upper:m{harmonic.harmonic_index}", q_max - q.abs_upper()))
        for rate_index, rate_id in enumerate(problem.rate_ids):
            fractional = _principal_fractional_frequency(
                order,
                harmonic.base_frequencies[rate_index],
                harmonic.harmonic_index,
            )
            denominator = acb(rate_dampings[rate_index]) + fractional
            disk_radius_lower = (
                denominator.abs_lower() * arb(harmonic.response_radii[rate_index])
                + arb(harmonic.morphology_drift_radii[rate_index])
            ).lower()
            residual_value = acb(residual[row], residual[row + 1])
            margins.append(
                (
                    f"disk:m{harmonic.harmonic_index}:{rate_id}",
                    disk_radius_lower - residual_value.abs_upper(),
                )
            )
            row += 2
    return tuple(margins)


def _outer_dyadic(value: arb) -> tuple[Fraction, Fraction] | None:
    if not value.is_finite():
        return None
    try:
        lower_float = math.nextafter(float(value.lower()), -math.inf)
        upper_float = math.nextafter(float(value.upper()), math.inf)
    except (OverflowError, ValueError):
        return None
    if not (math.isfinite(lower_float) and math.isfinite(upper_float)):
        return None
    for _ in range(64):
        lower = Fraction.from_float(lower_float)
        if value >= _arb_fraction(lower):
            break
        lower_float = math.nextafter(lower_float, -math.inf)
    else:
        return None
    for _ in range(64):
        upper = Fraction.from_float(upper_float)
        if value <= _arb_fraction(upper):
            break
        upper_float = math.nextafter(upper_float, math.inf)
    else:
        return None
    return lower, upper


def _positive_margin_claim(value: arb) -> Fraction | None:
    if not value > 0:
        return None
    try:
        lower = float(value.lower())
    except (OverflowError, ValueError):
        return None
    if not math.isfinite(lower) or lower <= 0.0:
        return None
    candidate = Fraction.from_float(lower) / 2
    if candidate > 0 and value > _arb_fraction(candidate):
        return candidate
    return None


def _matrix_payload(matrix: tuple[tuple[Fraction, ...], ...]) -> list[list[str]]:
    return [[_rat(value) for value in row] for row in matrix]


def _seal(payload: dict[str, Any]) -> tuple[str, str]:
    digest = _sha256_json(payload)
    document = dict(payload)
    document["certificate_sha256"] = digest
    return _canonical_json(document), digest


def _unknown_certificate(
    manifest: dict[str, Any],
    reason: str,
    attempted: tuple[int, ...],
    precision: int,
) -> TubeCertificate:
    input_hash = _sha256_json(manifest)
    payload = {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "input_manifest": manifest,
        "input_sha256": input_hash,
        "relation": TubeRelation.UNKNOWN.value,
        "reason": reason,
        "precision_bits": precision,
        "attempted_precisions": list(attempted),
        "backend": {
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
    }
    certificate_json, digest = _seal(payload)
    return TubeCertificate(
        TubeRelation.UNKNOWN,
        reason,
        precision,
        attempted,
        input_hash,
        digest,
        certificate_json,
    )


def _proof_at_precision(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    selector: _Selector,
    manifest: dict[str, Any],
    attempted: tuple[int, ...],
    precision: int,
) -> TubeCertificate | None:
    central = _central_vector(problem, box, selector)
    selected = _linear_ball_map(selector.selector_map, central)
    krawczyk = tuple(_outer_dyadic(value) for value in selected)
    if any(value is None for value in krawczyk):
        return None
    image = tuple(value for value in krawczyk if value is not None)

    closed_containment = all(
        domain[0] <= enclosure[0] and enclosure[1] <= domain[1]
        for domain, enclosure in zip(selector.variable_box, image, strict=True)
    )
    if not closed_containment:
        return None
    gaps = tuple(
        gap
        for domain, enclosure in zip(selector.variable_box, image, strict=True)
        for gap in (enclosure[0] - domain[0], domain[1] - enclosure[1])
    )
    krawczyk_margin = min(gaps)
    strict_containment = krawczyk_margin > 0

    constraint_balls = _constraint_balls(problem, box, selector, central, selected)
    if not all(value >= 0 for _, value in constraint_balls):
        return None
    claimed_margins: list[dict[str, str]] = []
    all_strict = True
    for name, value in constraint_balls:
        claim = _positive_margin_claim(value)
        if claim is None:
            all_strict = False
            claim = Fraction(0)
        claimed_margins.append({"name": name, "lower": _rat(claim)})

    relation = (
        TubeRelation.ROBUST_INNER
        if strict_containment and all_strict
        else TubeRelation.CLOSED_INNER
    )
    input_hash = _sha256_json(manifest)
    x_midpoints = tuple((lower + upper) / 2 for lower, upper in selector.variable_box)
    payload = {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "input_manifest": manifest,
        "input_sha256": input_hash,
        "relation": relation.value,
        "reason": "arb_parametric_gram_tube",
        "precision_bits": precision,
        "attempted_precisions": list(attempted),
        "backend": {
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
        "proof": {
            "variables": list(selector.variable_names),
            "variable_box": [[_rat(a), _rat(b)] for a, b in selector.variable_box],
            "center": [_rat(value) for value in x_midpoints],
            "preconditioner_c": _matrix_payload(selector.c_inverse),
            "contraction_inf_norm": "0",
            "krawczyk_image": [[_rat(a), _rat(b)] for a, b in image],
            "krawczyk_margin": _rat(krawczyk_margin),
            "constraint_margins": claimed_margins,
        },
    }
    certificate_json, digest = _seal(payload)
    if len(certificate_json.encode("utf-8")) > MAX_CERTIFICATE_BYTES:
        raise _CertificateResourceExceeded
    return TubeCertificate(
        relation,
        "arb_parametric_gram_tube",
        precision,
        attempted,
        input_hash,
        digest,
        certificate_json,
    )


def certify_alpha_tube(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    *,
    precision_schedule: tuple[int, ...] = DEFAULT_PRECISION_SCHEDULE,
) -> TubeCertificate:
    """Prove ``forall alpha in box.order, exists nuisance parameters``.

    The central Gram selector is only a deterministic candidate rule.  Every
    physical claim is rechecked against all original disks and annuli.  No
    failure path in this function has outer-exclusion semantics.
    """

    schedule = _validate_precision_schedule(precision_schedule)
    manifest = _problem_manifest(problem, box)
    if any(value > MAX_REPLAY_PRECISION for value in schedule):
        return _unknown_certificate(manifest, "precision_resource_limit", (), 0)
    try:
        selector, reason = _build_selector(problem, box)
    except ValueError:
        raise
    if selector is None:
        return _unknown_certificate(manifest, reason, (), 0)

    attempted: list[int] = []
    with _CERTIFICATION_LOCK:
        for precision in schedule:
            attempted.append(precision)
            try:
                with ctx.workprec(precision):
                    certificate = _proof_at_precision(
                        problem,
                        box,
                        selector,
                        manifest,
                        tuple(attempted),
                        precision,
                    )
            except _CertificateResourceExceeded:
                return _unknown_certificate(
                    manifest,
                    "certificate_resource_bytes",
                    tuple(attempted),
                    precision,
                )
            if certificate is not None:
                return certificate
    return _unknown_certificate(
        manifest,
        "arb_parametric_tube_undecided",
        tuple(attempted),
        attempted[-1],
    )


def _no_duplicate_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key")
        result[key] = value
    return result


def _load_document(certificate_json: str) -> dict[str, Any]:
    if not isinstance(certificate_json, str):
        raise ValueError("certificate must be JSON text")
    if len(certificate_json.encode("utf-8")) > MAX_CERTIFICATE_BYTES:
        raise ValueError("certificate exceeds the replay byte limit")
    value = json.loads(certificate_json, object_pairs_hook=_no_duplicate_object)
    if not isinstance(value, dict):
        raise ValueError("certificate root must be an object")
    return value


def _parse_matrix(value: Any, rows: int, columns: int) -> tuple[tuple[Fraction, ...], ...]:
    if not isinstance(value, list) or len(value) != rows:
        raise ValueError("matrix row count mismatch")
    parsed: list[tuple[Fraction, ...]] = []
    for row in value:
        if not isinstance(row, list) or len(row) != columns:
            raise ValueError("matrix column count mismatch")
        parsed.append(tuple(_parse_rat(item) for item in row))
    return tuple(parsed)


def _verify_proof(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    selector: _Selector,
    document: dict[str, Any],
) -> TubeReplayResult:
    proof = document.get("proof")
    if not isinstance(proof, dict):
        return TubeReplayResult(False, TubeRelation.UNKNOWN, "missing_proof")
    try:
        relation = TubeRelation(document["relation"])
        if document.get("reason") != "arb_parametric_gram_tube":
            raise ValueError("proof reason mismatch")
        precision = document["precision_bits"]
        if not isinstance(precision, int) or isinstance(precision, bool) or precision < 64:
            raise ValueError("invalid precision")
        variables = proof["variables"]
        if variables != list(selector.variable_names):
            raise ValueError("variable order mismatch")
        variable_box = tuple(
            (_parse_rat(value[0]), _parse_rat(value[1])) for value in proof["variable_box"]
        )
        if variable_box != selector.variable_box:
            raise ValueError("variable box mismatch")
        expected_center = tuple((a + b) / 2 for a, b in selector.variable_box)
        center = tuple(_parse_rat(value) for value in proof["center"])
        if center != expected_center:
            raise ValueError("center mismatch")
        dimension = len(selector.variable_names)
        preconditioner = _parse_matrix(proof["preconditioner_c"], dimension, dimension)
        if preconditioner != selector.c_inverse:
            raise ValueError("preconditioner mismatch")
        if _parse_rat(proof["contraction_inf_norm"]) != 0:
            raise ValueError("contraction mismatch")
        image = tuple(
            (_parse_rat(value[0]), _parse_rat(value[1])) for value in proof["krawczyk_image"]
        )
        if len(image) != dimension or any(a > b for a, b in image):
            raise ValueError("invalid Krawczyk image")
        stored_k_margin = _parse_rat(proof["krawczyk_margin"])
        margins = proof["constraint_margins"]
        if not isinstance(margins, list):
            raise ValueError("invalid margin table")
    except (KeyError, TypeError, ValueError, IndexError):
        return TubeReplayResult(False, TubeRelation.UNKNOWN, "malformed_proof")

    with _CERTIFICATION_LOCK:
        with ctx.workprec(precision):
            central = _central_vector(problem, box, selector)
            selected = _linear_ball_map(selector.selector_map, central)
            constraints = _constraint_balls(problem, box, selector, central, selected)
            for ball, enclosure in zip(selected, image, strict=True):
                if not (
                    ball >= _arb_fraction(enclosure[0]) and ball <= _arb_fraction(enclosure[1])
                ):
                    return TubeReplayResult(False, TubeRelation.UNKNOWN, "krawczyk_image_mismatch")

            closed_containment = all(
                domain[0] <= enclosure[0] and enclosure[1] <= domain[1]
                for domain, enclosure in zip(selector.variable_box, image, strict=True)
            )
            if not closed_containment:
                return TubeReplayResult(False, TubeRelation.UNKNOWN, "krawczyk_not_contained")
            gaps = tuple(
                gap
                for domain, enclosure in zip(selector.variable_box, image, strict=True)
                for gap in (enclosure[0] - domain[0], domain[1] - enclosure[1])
            )
            actual_k_margin = min(gaps)
            if stored_k_margin != actual_k_margin:
                return TubeReplayResult(False, TubeRelation.UNKNOWN, "krawczyk_margin_mismatch")

            if len(margins) != len(constraints):
                return TubeReplayResult(False, TubeRelation.UNKNOWN, "constraint_count_mismatch")
            all_claimed_strict = True
            for stored, (expected_name, ball) in zip(margins, constraints, strict=True):
                if not isinstance(stored, dict) or stored.get("name") != expected_name:
                    return TubeReplayResult(False, TubeRelation.UNKNOWN, "constraint_name_mismatch")
                try:
                    claimed = _parse_rat(stored.get("lower"))
                except ValueError:
                    return TubeReplayResult(
                        False, TubeRelation.UNKNOWN, "constraint_margin_malformed"
                    )
                if claimed < 0:
                    return TubeReplayResult(
                        False, TubeRelation.UNKNOWN, "constraint_margin_negative"
                    )
                if claimed > 0:
                    if not ball > _arb_fraction(claimed):
                        return TubeReplayResult(
                            False, TubeRelation.UNKNOWN, "constraint_margin_false"
                        )
                else:
                    all_claimed_strict = False
                    if not ball >= 0:
                        return TubeReplayResult(
                            False, TubeRelation.UNKNOWN, "closed_constraint_false"
                        )

    strict_k = actual_k_margin > 0
    if relation is TubeRelation.ROBUST_INNER:
        if not (strict_k and all_claimed_strict):
            return TubeReplayResult(False, TubeRelation.UNKNOWN, "robust_claim_without_margin")
    elif relation is not TubeRelation.CLOSED_INNER:
        return TubeReplayResult(False, TubeRelation.UNKNOWN, "proof_relation_invalid")
    return TubeReplayResult(True, relation, "replay_verified")


def replay_tube_certificate(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    certificate_json: str,
) -> TubeReplayResult:
    """Replay a sealed certificate without trusting any derived field."""

    try:
        document = _load_document(certificate_json)
        digest = document.get("certificate_sha256")
        if not isinstance(digest, str) or len(digest) != 64:
            raise ValueError("invalid digest")
        unsigned = dict(document)
        del unsigned["certificate_sha256"]
        if _sha256_json(unsigned) != digest:
            raise ValueError("digest mismatch")
        if document.get("schema") != SCHEMA or document.get("algorithm") != ALGORITHM:
            raise ValueError("schema mismatch")
        if document.get("backend") != {
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        }:
            raise ValueError("backend mismatch")
        precision = document.get("precision_bits")
        attempted = document.get("attempted_precisions")
        if (
            not isinstance(precision, int)
            or isinstance(precision, bool)
            or not 0 <= precision <= MAX_REPLAY_PRECISION
            or not isinstance(attempted, list)
            or any(
                not isinstance(value, int)
                or isinstance(value, bool)
                or not 64 <= value <= MAX_REPLAY_PRECISION
                for value in attempted
            )
            or any(right <= left for left, right in zip(attempted, attempted[1:], strict=False))
            or (attempted and attempted[-1] != precision)
            or (not attempted and precision != 0)
        ):
            raise ValueError("invalid precision evidence")
        if not isinstance(document.get("reason"), str) or not document["reason"]:
            raise ValueError("invalid reason")
        manifest = _problem_manifest(problem, box)
        input_hash = _sha256_json(manifest)
        if document.get("input_manifest") != manifest or document.get("input_sha256") != input_hash:
            raise ValueError("input hash mismatch")
        relation = TubeRelation(document["relation"])
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return TubeReplayResult(False, TubeRelation.UNKNOWN, "invalid_certificate_envelope")

    if relation is TubeRelation.UNKNOWN:
        return TubeReplayResult(True, TubeRelation.UNKNOWN, "replay_unknown")
    try:
        selector, reason = _build_selector(problem, box)
    except ValueError:
        return TubeReplayResult(False, TubeRelation.UNKNOWN, "invalid_problem")
    if selector is None:
        return TubeReplayResult(False, TubeRelation.UNKNOWN, reason)
    return _verify_proof(problem, box, selector, document)


def _replace_order(
    box: CertifiedParameterBox,
    order: ClosedInterval,
) -> CertifiedParameterBox:
    return CertifiedParameterBox(order, box.damping, box.damping_offsets)


def _bisect_alpha(
    box: CertifiedParameterBox,
) -> tuple[CertifiedParameterBox, CertifiedParameterBox] | None:
    """Deterministically bisect only the binary-float order coordinate."""

    order = box.order
    if order.lower == order.upper:
        return None
    midpoint = order.midpoint
    if midpoint in (order.lower, order.upper):
        return None
    return (
        _replace_order(box, ClosedInterval(order.lower, midpoint)),
        _replace_order(box, ClosedInterval(midpoint, order.upper)),
    )


def _merge_order_intervals(intervals: list[ClosedInterval]) -> tuple[ClosedInterval, ...]:
    if not intervals:
        return ()
    ordered = sorted(intervals, key=lambda value: (value.lower, value.upper))
    merged = [ordered[0]]
    for interval in ordered[1:]:
        previous = merged[-1]
        if interval.lower <= previous.upper:
            merged[-1] = ClosedInterval(previous.lower, max(previous.upper, interval.upper))
        else:
            merged.append(interval)
    return tuple(merged)


def _alpha_partition_covers(
    initial: CertifiedParameterBox,
    leaves: tuple[AlphaTubeLeaf, ...],
) -> bool:
    """Replay every binary path and prove a complete prefix-free tree."""

    paths = tuple(leaf.path for leaf in leaves)
    if len(paths) != len(set(paths)):
        return False
    path_set = set(paths)
    for leaf in leaves:
        path = leaf.path
        if leaf.depth != len(path) or any(step not in (0, 1) for step in path):
            return False
        if any(path[:length] in path_set for length in range(len(path))):
            return False
        replayed = initial
        for step in path:
            children = _bisect_alpha(replayed)
            if children is None:
                return False
            replayed = children[step]
        if replayed != leaf.box:
            return False
    return sum((Fraction(1, 2 ** len(path)) for path in paths), Fraction(0)) == 1


def _branch_intervals(
    leaves: tuple[AlphaTubeLeaf, ...],
    relation: TubeRelation | None,
) -> tuple[ClosedInterval, ...]:
    if relation is None:
        selected = [
            leaf.box.order
            for leaf in leaves
            if leaf.certificate.relation in (TubeRelation.ROBUST_INNER, TubeRelation.CLOSED_INNER)
        ]
    else:
        selected = [leaf.box.order for leaf in leaves if leaf.certificate.relation is relation]
    return _merge_order_intervals(selected)


def branch_alpha_tube(
    problem: CertifiedJointProblem,
    initial_box: CertifiedParameterBox,
    *,
    max_depth: int = 8,
    max_leaves: int = 256,
    precision_schedule: tuple[int, ...] = DEFAULT_PRECISION_SCHEDULE,
) -> AlphaTubeBranchResult:
    """Subdivide alpha only and retain every undecided interval as UNKNOWN."""

    _validate_box(problem, initial_box)
    schedule = _validate_precision_schedule(precision_schedule)
    if not 0 <= max_depth <= MAX_ALPHA_DEPTH:
        raise ValueError(f"max_depth must lie in [0, {MAX_ALPHA_DEPTH}]")
    if not 1 <= max_leaves <= MAX_ALPHA_LEAVES:
        raise ValueError(f"max_leaves must lie in [1, {MAX_ALPHA_LEAVES}]")

    pending: list[tuple[CertifiedParameterBox, int, tuple[int, ...]]] = [(initial_box, 0, ())]
    leaves: list[AlphaTubeLeaf] = []
    budget_exhausted = False
    while pending:
        box, depth, path = pending.pop()
        certificate = certify_alpha_tube(problem, box, precision_schedule=schedule)
        children = _bisect_alpha(box)
        subdividable = (
            certificate.relation is TubeRelation.UNKNOWN
            and certificate.reason == _SUBDIVIDABLE_UNKNOWN_REASON
            and children is not None
        )
        if subdividable:
            can_add_children = len(leaves) + len(pending) + 2 <= max_leaves
            if depth >= max_depth:
                budget_exhausted = True
                leaves.append(AlphaTubeLeaf(box, depth, path, certificate, "max_depth_unknown"))
            elif not can_add_children:
                budget_exhausted = True
                leaves.append(AlphaTubeLeaf(box, depth, path, certificate, "max_leaves_unknown"))
            else:
                pending.append((children[1], depth + 1, (*path, 1)))
                pending.append((children[0], depth + 1, (*path, 0)))
        else:
            if certificate.relation is not TubeRelation.UNKNOWN:
                terminal = "certified"
            elif certificate.reason in _RESOURCE_UNKNOWN_REASONS:
                terminal = "resource_unknown"
            elif certificate.reason != _SUBDIVIDABLE_UNKNOWN_REASON:
                terminal = "structural_unknown"
            else:
                terminal = "atomic_unknown"
            leaves.append(AlphaTubeLeaf(box, depth, path, certificate, terminal))

    frozen_leaves = tuple(leaves)
    coverage = _alpha_partition_covers(initial_box, frozen_leaves)
    if not coverage:
        raise RuntimeError("alpha tube leaves lost binary partition coverage")
    return AlphaTubeBranchResult(
        leaves=frozen_leaves,
        certified_inner=_branch_intervals(frozen_leaves, None),
        robust_inner=_branch_intervals(frozen_leaves, TubeRelation.ROBUST_INNER),
        closed_inner=_branch_intervals(frozen_leaves, TubeRelation.CLOSED_INNER),
        unknown=_branch_intervals(frozen_leaves, TubeRelation.UNKNOWN),
        coverage_certified=True,
        budget_exhausted=budget_exhausted,
        max_depth=max_depth,
        max_leaves=max_leaves,
    )


def replay_alpha_tube_branch(
    problem: CertifiedJointProblem,
    initial_box: CertifiedParameterBox,
    result: AlphaTubeBranchResult,
) -> AlphaTubeBranchReplayResult:
    """Replay the partition geometry and every terminal leaf certificate."""

    if (
        not 0 <= result.max_depth <= MAX_ALPHA_DEPTH
        or not 1 <= result.max_leaves <= MAX_ALPHA_LEAVES
    ):
        return AlphaTubeBranchReplayResult(False, "invalid_branch_budget", ())
    if len(result.leaves) > result.max_leaves or not _alpha_partition_covers(
        initial_box, result.leaves
    ):
        return AlphaTubeBranchReplayResult(False, "invalid_branch_partition", ())

    leaf_replays: list[TubeReplayResult] = []
    for leaf in result.leaves:
        replay = replay_tube_certificate(
            problem,
            leaf.box,
            leaf.certificate.to_json(),
        )
        leaf_replays.append(replay)
        if not replay.valid or replay.relation is not leaf.certificate.relation:
            return AlphaTubeBranchReplayResult(
                False,
                "invalid_leaf_certificate",
                tuple(leaf_replays),
            )
        try:
            sealed_reason = _load_document(leaf.certificate.to_json())["reason"]
        except (KeyError, TypeError, ValueError, json.JSONDecodeError):
            return AlphaTubeBranchReplayResult(
                False,
                "invalid_leaf_envelope",
                tuple(leaf_replays),
            )
        if sealed_reason != leaf.certificate.reason:
            return AlphaTubeBranchReplayResult(
                False,
                "leaf_metadata_mismatch",
                tuple(leaf_replays),
            )
        children = _bisect_alpha(leaf.box)
        if replay.relation is TubeRelation.UNKNOWN:
            expected_terminal = {
                "atomic_unknown"
                if sealed_reason == _SUBDIVIDABLE_UNKNOWN_REASON and children is None
                else "",
                "resource_unknown" if sealed_reason in _RESOURCE_UNKNOWN_REASONS else "",
                "structural_unknown"
                if sealed_reason != _SUBDIVIDABLE_UNKNOWN_REASON
                and sealed_reason not in _RESOURCE_UNKNOWN_REASONS
                else "",
                "max_depth_unknown"
                if sealed_reason == _SUBDIVIDABLE_UNKNOWN_REASON
                and children is not None
                and leaf.depth >= result.max_depth
                else "",
                "max_leaves_unknown"
                if sealed_reason == _SUBDIVIDABLE_UNKNOWN_REASON
                and children is not None
                and len(result.leaves) >= result.max_leaves
                else "",
            }
            expected_terminal.discard("")
            if leaf.terminal_reason not in expected_terminal:
                return AlphaTubeBranchReplayResult(
                    False,
                    "invalid_unknown_terminal",
                    tuple(leaf_replays),
                )
        elif leaf.terminal_reason != "certified":
            return AlphaTubeBranchReplayResult(
                False,
                "invalid_inner_terminal",
                tuple(leaf_replays),
            )

    leaves = result.leaves
    expected = (
        _branch_intervals(leaves, None),
        _branch_intervals(leaves, TubeRelation.ROBUST_INNER),
        _branch_intervals(leaves, TubeRelation.CLOSED_INNER),
        _branch_intervals(leaves, TubeRelation.UNKNOWN),
    )
    actual = (
        result.certified_inner,
        result.robust_inner,
        result.closed_inner,
        result.unknown,
    )
    if expected != actual:
        return AlphaTubeBranchReplayResult(
            False,
            "branch_interval_summary_mismatch",
            tuple(leaf_replays),
        )
    expected_budget = any(
        leaf.terminal_reason in ("max_depth_unknown", "max_leaves_unknown") for leaf in leaves
    )
    if not result.coverage_certified or result.budget_exhausted != expected_budget:
        return AlphaTubeBranchReplayResult(
            False,
            "branch_status_mismatch",
            tuple(leaf_replays),
        )
    return AlphaTubeBranchReplayResult(True, "branch_replay_verified", tuple(leaf_replays))


__all__ = [
    "ALGORITHM",
    "AlphaTubeBranchReplayResult",
    "AlphaTubeBranchResult",
    "AlphaTubeLeaf",
    "SCHEMA",
    "TubeCertificate",
    "TubeRelation",
    "TubeReplayResult",
    "branch_alpha_tube",
    "certify_alpha_tube",
    "replay_alpha_tube_branch",
    "replay_tube_certificate",
]
