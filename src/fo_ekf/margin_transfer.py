"""Replayable one-way transfer of a certified T1 selector margin.

The transfer freezes the selector path proved by an existing
``ROBUST_INNER`` tube certificate.  It never rebuilds the selector from a
perturbed response centre.  For every rate--harmonic disk it checks the
absolute-margin implication

``new_margin >= old_margin
                 - |d| * center_shift
                 + |d| * (new_epsilon - old_epsilon)
                 + (new_rho - old_rho)``.

Only preservation statements are expressible.  Failure of this sufficient
test is ``UNKNOWN`` and has no outer-exclusion meaning.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from typing import Any

import flint
from flint import acb, arb, ctx

from .certified_set import (
    _CERTIFICATION_LOCK,
    CertifiedJointProblem,
    CertifiedParameterBox,
    _validate_box,
)
from .tube_certificate import TubeRelation, replay_tube_certificate

SCHEMA = "fo-ekf.r2-to-t1-margin-transfer.v1"
ALGORITHM = "frozen-selector-monotone-transfer-v1"
MAX_CERTIFICATE_BYTES = 8_000_000
MAX_RATIONAL_BITS = 16_384
MAX_RATIONAL_CHARACTERS = 20_000
MAX_REPLAY_PRECISION = 4096


class MarginTransferRelation(str, Enum):
    """Sound relations supported by the one-way transfer."""

    PRESERVED_ROBUST = "PRESERVED_ROBUST"
    PRESERVED_CLOSED = "PRESERVED_CLOSED"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class DiskMarginPerturbation:
    """Marginal perturbation bounds for one canonical response disk.

    ``center_shift_abs_upper`` has response units and is universally
    quantified: the new centre may be any complex value within that distance
    of the old centre.  Radius values are exact declared inputs, not sampled
    estimates.
    """

    harmonic_index: int
    rate_id: str
    center_shift_abs_upper: float
    old_response_radius: float
    new_response_radius: float
    old_morphology_radius: float
    new_morphology_radius: float

    def __post_init__(self) -> None:
        if isinstance(self.harmonic_index, bool) or not isinstance(self.harmonic_index, int):
            raise ValueError("harmonic_index must be a positive integer")
        if self.harmonic_index < 1:
            raise ValueError("harmonic_index must be a positive integer")
        rate_id = str(self.rate_id)
        if not rate_id:
            raise ValueError("rate_id must be nonempty")
        values = (
            self.center_shift_abs_upper,
            self.old_response_radius,
            self.new_response_radius,
            self.old_morphology_radius,
            self.new_morphology_radius,
        )
        converted = tuple(float(value) for value in values)
        if any(not math.isfinite(value) or value < 0.0 for value in converted):
            raise ValueError("all centre-shift and radius values must be finite and nonnegative")
        object.__setattr__(self, "rate_id", rate_id)
        (
            center_shift,
            old_response,
            new_response,
            old_morphology,
            new_morphology,
        ) = converted
        object.__setattr__(self, "center_shift_abs_upper", center_shift)
        object.__setattr__(self, "old_response_radius", old_response)
        object.__setattr__(self, "new_response_radius", new_response)
        object.__setattr__(self, "old_morphology_radius", old_morphology)
        object.__setattr__(self, "new_morphology_radius", new_morphology)


@dataclass(frozen=True)
class MarginTransferCertificate:
    """A sealed preservation certificate or an auditable UNKNOWN result."""

    relation: MarginTransferRelation
    reason: str
    precision_bits: int
    input_sha256: str
    certificate_sha256: str
    certificate_json: str

    def to_json(self) -> str:
        return self.certificate_json


@dataclass(frozen=True)
class MarginTransferReplayResult:
    """Result of rebuilding every transfer quantity from bound inputs."""

    valid: bool
    relation: MarginTransferRelation
    reason: str


@dataclass(frozen=True)
class _Evaluation:
    relation: MarginTransferRelation
    reason: str
    precision_bits: int
    proof: dict[str, Any] | None


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


def _arb_interval(lower: Fraction, upper: Fraction) -> arb:
    return _arb_fraction(lower).union(_arb_fraction(upper))


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


def _seal(payload: dict[str, Any]) -> tuple[str, str]:
    digest = _sha256_json(payload)
    document = dict(payload)
    document["certificate_sha256"] = digest
    return _canonical_json(document), digest


def _load_document(value: str) -> dict[str, Any]:
    if not isinstance(value, str) or len(value.encode("utf-8")) > MAX_CERTIFICATE_BYTES:
        raise ValueError("certificate exceeds the replay resource limit")

    def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        document_object: dict[str, Any] = {}
        for key, item in pairs:
            if key in document_object:
                raise ValueError("duplicate JSON object key")
            document_object[key] = item
        return document_object

    try:
        document = json.loads(value, object_pairs_hook=reject_duplicate_keys)
    except RecursionError as exc:
        raise ValueError("certificate nesting exceeds the replay resource limit") from exc
    if not isinstance(document, dict):
        raise ValueError("certificate must be a JSON object")
    return document


def _perturbation_payload(value: DiskMarginPerturbation) -> dict[str, Any]:
    return {
        "harmonic_index": value.harmonic_index,
        "rate_id": value.rate_id,
        "center_shift_abs_upper": _rat(_fraction(value.center_shift_abs_upper)),
        "old_response_radius": _rat(_fraction(value.old_response_radius)),
        "new_response_radius": _rat(_fraction(value.new_response_radius)),
        "old_morphology_radius": _rat(_fraction(value.old_morphology_radius)),
        "new_morphology_radius": _rat(_fraction(value.new_morphology_radius)),
    }


def _input_manifest(
    baseline_certificate_json: str,
    perturbations: tuple[DiskMarginPerturbation, ...],
) -> dict[str, Any]:
    return {
        "semantics": "frozen-old-selector-forall-bounded-center-shifts",
        "baseline_certificate_json_sha256": hashlib.sha256(
            baseline_certificate_json.encode("utf-8")
        ).hexdigest(),
        "perturbations": [_perturbation_payload(value) for value in perturbations],
    }


def _expected_layout(
    problem: CertifiedJointProblem,
) -> tuple[tuple[int, str, float, float], ...]:
    return tuple(
        (
            harmonic.harmonic_index,
            rate_id,
            harmonic.response_radii[rate_index],
            harmonic.morphology_drift_radii[rate_index],
        )
        for harmonic in problem.harmonics
        for rate_index, rate_id in enumerate(problem.rate_ids)
    )


def _validate_perturbations(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    perturbations: Sequence[DiskMarginPerturbation],
) -> tuple[DiskMarginPerturbation, ...]:
    _validate_box(problem, box)
    frozen = tuple(perturbations)
    expected = _expected_layout(problem)
    actual_ids = tuple((value.harmonic_index, value.rate_id) for value in frozen)
    expected_ids = tuple((harmonic, rate_id) for harmonic, rate_id, _, _ in expected)
    if actual_ids != expected_ids:
        raise ValueError("one perturbation is required per disk in canonical harmonic/rate order")
    for value, (_, _, response_radius, morphology_radius) in zip(frozen, expected, strict=True):
        if value.old_response_radius != response_radius:
            raise ValueError("old_response_radius must match the baseline T1 problem exactly")
        if value.old_morphology_radius != morphology_radius:
            raise ValueError("old_morphology_radius must match the baseline T1 problem exactly")
    return frozen


def _outer_dyadic(value: arb) -> tuple[Fraction, Fraction] | None:
    """Return deterministic binary-float rational bounds containing an Arb value."""

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


def _baseline_material(
    baseline_document: dict[str, Any],
    expected_names: tuple[str, ...],
) -> tuple[
    int,
    dict[str, tuple[Fraction, Fraction]],
    dict[str, Fraction],
]:
    precision = baseline_document["precision_bits"]
    if (
        not isinstance(precision, int)
        or isinstance(precision, bool)
        or not 64 <= precision <= MAX_REPLAY_PRECISION
    ):
        raise ValueError("invalid baseline precision")
    proof = baseline_document["proof"]
    if not isinstance(proof, dict):
        raise ValueError("missing baseline proof")
    variables = proof["variables"]
    image_values = proof["selector_image"] if "selector_image" in proof else proof["krawczyk_image"]
    if not isinstance(variables, list) or not isinstance(image_values, list):
        raise ValueError("invalid baseline selector image")
    if len(variables) != len(image_values) or any(not isinstance(name, str) for name in variables):
        raise ValueError("invalid baseline selector dimensions")
    if len(set(variables)) != len(variables):
        raise ValueError("duplicate baseline selector variable")
    image: dict[str, tuple[Fraction, Fraction]] = {}
    for name, endpoints in zip(variables, image_values, strict=True):
        if not isinstance(endpoints, list) or len(endpoints) != 2:
            raise ValueError("invalid baseline selector interval")
        lower, upper = (_parse_rat(item) for item in endpoints)
        if lower > upper:
            raise ValueError("reversed baseline selector interval")
        image[name] = (lower, upper)

    constraints = proof["constraint_margins"]
    if not isinstance(constraints, list):
        raise ValueError("invalid baseline constraint table")
    claims: dict[str, Fraction] = {}
    for entry in constraints:
        if not isinstance(entry, dict) or not isinstance(entry.get("name"), str):
            raise ValueError("invalid baseline constraint entry")
        name = entry["name"]
        if name in claims:
            raise ValueError("duplicate baseline constraint name")
        claims[name] = _parse_rat(entry.get("lower"))
    disk_claims: dict[str, Fraction] = {}
    for name in expected_names:
        claim = claims.get(name)
        if claim is None or claim <= 0:
            raise ValueError("baseline robust disk margin is missing")
        disk_claims[name] = claim
    return precision, image, disk_claims


def _physical_damping_interval(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    image: dict[str, tuple[Fraction, Fraction]],
    rate_index: int,
) -> tuple[Fraction, Fraction]:
    if "damping" in image:
        lower, upper = image["damping"]
    else:
        lower = upper = _fraction(box.damping.lower)
    rate_id = problem.rate_ids[rate_index]
    offset_name = f"offset:{rate_id}"
    if offset_name in image:
        offset_lower, offset_upper = image[offset_name]
    else:
        offset = box.damping_offsets[rate_index].interval
        if offset.lower != offset.upper:
            raise ValueError("free damping offset is absent from the baseline selector image")
        offset_lower = offset_upper = _fraction(offset.lower)
    return lower + offset_lower, upper + offset_upper


def _denominator_abs_bounds(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    image: dict[str, tuple[Fraction, Fraction]],
    harmonic_position: int,
    rate_index: int,
) -> tuple[Fraction, Fraction] | None:
    damping_lower, damping_upper = _physical_damping_interval(problem, box, image, rate_index)
    order = _arb_interval(_fraction(box.order.lower), _fraction(box.order.upper))
    harmonic = problem.harmonics[harmonic_position]
    frequency = _arb_fraction(_fraction(harmonic.base_frequencies[rate_index]))
    frequency *= harmonic.harmonic_index
    amplitude = (order * frequency.log()).exp()
    angle = order * arb.pi() / 2
    fractional = acb(amplitude * angle.cos(), amplitude * angle.sin())
    damping = _arb_interval(damping_lower, damping_upper)
    denominator = acb(damping) + fractional
    lower_outer = _outer_dyadic(denominator.abs_lower())
    upper_outer = _outer_dyadic(denominator.abs_upper())
    if lower_outer is None or upper_outer is None:
        return None
    lower = max(Fraction(0), lower_outer[0])
    upper = upper_outer[1]
    if upper < lower:
        return None
    return lower, upper


def _evaluate_transfer(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    baseline_certificate_json: str,
    perturbations: tuple[DiskMarginPerturbation, ...],
) -> _Evaluation:
    baseline_replay = replay_tube_certificate(problem, box, baseline_certificate_json)
    if not baseline_replay.valid:
        return _Evaluation(
            MarginTransferRelation.UNKNOWN,
            "baseline_certificate_invalid",
            0,
            None,
        )
    if baseline_replay.relation is not TubeRelation.ROBUST_INNER:
        return _Evaluation(
            MarginTransferRelation.UNKNOWN,
            "baseline_not_robust",
            0,
            None,
        )
    try:
        baseline_document = _load_document(baseline_certificate_json)
        expected_names = tuple(
            f"disk:m{harmonic.harmonic_index}:{rate_id}"
            for harmonic in problem.harmonics
            for rate_id in problem.rate_ids
        )
        precision, image, claims = _baseline_material(baseline_document, expected_names)
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return _Evaluation(
            MarginTransferRelation.UNKNOWN,
            "baseline_margin_material_invalid",
            0,
            None,
        )

    rows: list[dict[str, Any]] = []
    all_positive = True
    any_negative = False
    with _CERTIFICATION_LOCK:
        with ctx.workprec(precision):
            position = 0
            for harmonic_position, harmonic in enumerate(problem.harmonics):
                for rate_index, rate_id in enumerate(problem.rate_ids):
                    perturbation = perturbations[position]
                    position += 1
                    bounds = _denominator_abs_bounds(
                        problem,
                        box,
                        image,
                        harmonic_position,
                        rate_index,
                    )
                    if bounds is None:
                        return _Evaluation(
                            MarginTransferRelation.UNKNOWN,
                            "denominator_bound_undecided",
                            precision,
                            None,
                        )
                    denominator_lower, denominator_upper = bounds
                    name = f"disk:m{harmonic.harmonic_index}:{rate_id}"
                    baseline_margin = claims[name]
                    center_shift = _fraction(perturbation.center_shift_abs_upper)
                    response_delta = _fraction(perturbation.new_response_radius) - _fraction(
                        perturbation.old_response_radius
                    )
                    morphology_delta = _fraction(perturbation.new_morphology_radius) - _fraction(
                        perturbation.old_morphology_radius
                    )
                    coefficient = center_shift - response_delta
                    selected_denominator = (
                        denominator_upper if coefficient >= 0 else denominator_lower
                    )
                    adverse_load = selected_denominator * coefficient - morphology_delta
                    remaining = baseline_margin - adverse_load
                    if remaining <= 0:
                        all_positive = False
                    if remaining < 0:
                        any_negative = True
                    rows.append(
                        {
                            "name": name,
                            "denominator_abs": [
                                _rat(denominator_lower),
                                _rat(denominator_upper),
                            ],
                            "baseline_margin_abs_lower": _rat(baseline_margin),
                            "center_shift_abs_upper": _rat(center_shift),
                            "response_radius_delta": _rat(response_delta),
                            "morphology_radius_delta": _rat(morphology_delta),
                            "adverse_coefficient": _rat(coefficient),
                            "adverse_load_abs_upper": _rat(adverse_load),
                            "remaining_margin_lower": _rat(remaining),
                        }
                    )

    proof = {
        "semantics": "frozen-old-selector-forall-bounded-center-shifts",
        "baseline_relation": TubeRelation.ROBUST_INNER.value,
        "per_disk": rows,
    }
    if any_negative:
        return _Evaluation(
            MarginTransferRelation.UNKNOWN,
            "margin_budget_exceeded",
            precision,
            proof,
        )
    if all_positive:
        return _Evaluation(
            MarginTransferRelation.PRESERVED_ROBUST,
            "strict_margin_preserved",
            precision,
            proof,
        )
    return _Evaluation(
        MarginTransferRelation.PRESERVED_CLOSED,
        "closed_margin_preserved",
        precision,
        proof,
    )


def certify_margin_transfer(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    baseline_certificate_json: str,
    perturbations: Sequence[DiskMarginPerturbation],
) -> MarginTransferCertificate:
    """Transfer one robust T1 witness to bounded R2 disk perturbations.

    The function proves a universal statement over every complex centre shift
    satisfying each supplied magnitude bound.  A failed budget is only
    ``UNKNOWN``; this module cannot issue an outer relation.
    """

    if not isinstance(baseline_certificate_json, str):
        raise TypeError("baseline_certificate_json must be a string")
    frozen = _validate_perturbations(problem, box, perturbations)
    manifest = _input_manifest(baseline_certificate_json, frozen)
    evaluation = _evaluate_transfer(problem, box, baseline_certificate_json, frozen)
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "input_manifest": manifest,
        "input_sha256": _sha256_json(manifest),
        "relation": evaluation.relation.value,
        "reason": evaluation.reason,
        "precision_bits": evaluation.precision_bits,
        "backend": {
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
    }
    if evaluation.proof is not None:
        payload["proof"] = evaluation.proof
    certificate_json, digest = _seal(payload)
    if len(certificate_json.encode("utf-8")) > MAX_CERTIFICATE_BYTES:
        evaluation = _Evaluation(
            MarginTransferRelation.UNKNOWN,
            "certificate_resource_bytes",
            0,
            None,
        )
        payload.update(
            {
                "relation": evaluation.relation.value,
                "reason": evaluation.reason,
                "precision_bits": 0,
            }
        )
        payload.pop("proof", None)
        certificate_json, digest = _seal(payload)
    return MarginTransferCertificate(
        evaluation.relation,
        evaluation.reason,
        evaluation.precision_bits,
        payload["input_sha256"],
        digest,
        certificate_json,
    )


def replay_margin_transfer_certificate(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    baseline_certificate_json: str,
    perturbations: Sequence[DiskMarginPerturbation],
    certificate_json: str,
) -> MarginTransferReplayResult:
    """Replay a sealed transfer without trusting any stored derived field."""

    try:
        frozen = _validate_perturbations(problem, box, perturbations)
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
        manifest = _input_manifest(baseline_certificate_json, frozen)
        if document.get("input_manifest") != manifest:
            raise ValueError("input manifest mismatch")
        if document.get("input_sha256") != _sha256_json(manifest):
            raise ValueError("input hash mismatch")
        stored_relation = MarginTransferRelation(document["relation"])
        stored_reason = document["reason"]
        stored_precision = document["precision_bits"]
        if not isinstance(stored_reason, str) or not stored_reason:
            raise ValueError("invalid reason")
        if (
            not isinstance(stored_precision, int)
            or isinstance(stored_precision, bool)
            or not 0 <= stored_precision <= MAX_REPLAY_PRECISION
        ):
            raise ValueError("invalid precision")
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return MarginTransferReplayResult(
            False,
            MarginTransferRelation.UNKNOWN,
            "invalid_certificate_envelope",
        )

    evaluation = _evaluate_transfer(problem, box, baseline_certificate_json, frozen)
    if (
        stored_relation is not evaluation.relation
        or stored_reason != evaluation.reason
        or stored_precision != evaluation.precision_bits
        or document.get("proof") != evaluation.proof
    ):
        return MarginTransferReplayResult(
            False,
            MarginTransferRelation.UNKNOWN,
            "transfer_proof_mismatch",
        )
    return MarginTransferReplayResult(True, stored_relation, "replay_verified")
