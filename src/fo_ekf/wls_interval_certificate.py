"""Rigorous integer-payload certificates for the frozen joint complex WLS center.

The ordinary WLS execution record is intentionally a binary64 execution trace.
Its embedded enclosure starts from the *recorded* binary64 ``G`` and ``b`` and
therefore does not cover phase, sample conversion, or accumulation rounding.
This module closes that narrow numerical gap.  It first requires the execution
record to replay exactly, then rebuilds the mathematical sample functional from
the selected int64 payload with Arb/ACB ball arithmetic:

``y_n = (d_n - baseline) / gain``, ``w_n = 1/n_s`` for ``n_s`` samples,
``phi_n = 2*pi*N*(j_n-a)/(b-a)``, ``G = Psi^* W Psi``, and
``b = Psi^* W y``.

The complete contiguous two-sided basis from the execution record is used.
The ADC gain and front-end responses are interpreted as the exact dyadic values
of their frozen binary64 encodings.  Each positive-harmonic result is a proved
complex disk, centered at the execution record's binary64 ``z_tilde=C/H_hat``,
whose outward dyadic radius contains the exact-functional response.

This is only a numerical certificate for the frozen mathematical sample
functional.  It does not cover acquisition, beat annotation, anti-aliasing,
front-end calibration, model discrepancy, or any physiological claim.
"""

from __future__ import annotations

import hashlib
import math
import platform
import re
from fractions import Fraction
from pathlib import Path
from typing import Any

import flint
import numpy as np
from flint import acb, acb_mat, arb, ctx

from . import wls_execution as _wls_execution
from .wls_execution import (
    MAX_ARB_PRECISION_BITS,
    MAX_RETAINED_HARMONICS,
    MAX_SAMPLES,
    bounded_sample_indices,
    canonical_hash,
    replay_wls_execution_record,
    retained_harmonic_set,
)

SCHEMA_VERSION = "fo-ekf-exact-integer-joint-wls-interval-certificate-v1"
ALGORITHM = "arb-direct-exact-payload-joint-complex-wls-v1"
MIN_PRECISION_BITS = 80
DEFAULT_PRECISION_BITS = 256
MAX_DOCUMENT_DEPTH = 32
MAX_DOCUMENT_NODES = 20_000
MAX_DOCUMENT_TEXT_CHARS = 1_000_000
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _plain_tree_error(value: Any) -> str | None:
    """Reject non-JSON material, cycles, aliases, and resource bombs."""

    stack = [(value, 0)]
    seen_containers: set[int] = set()
    nodes = 0
    text_chars = 0
    while stack:
        current, depth = stack.pop()
        nodes += 1
        if nodes > MAX_DOCUMENT_NODES:
            return "document_node_budget_exceeded"
        if depth > MAX_DOCUMENT_DEPTH:
            return "document_depth_budget_exceeded"
        if type(current) is str:
            text_chars += len(current)
            if text_chars > MAX_DOCUMENT_TEXT_CHARS:
                return "document_text_budget_exceeded"
        elif current is None or type(current) in (bool, int, float):
            continue
        elif type(current) is dict:
            identity = id(current)
            if identity in seen_containers:
                return "document_cycle_or_alias_detected"
            seen_containers.add(identity)
            for key, child in current.items():
                if type(key) is not str:
                    return "document_nonstring_key"
                text_chars += len(key)
                if text_chars > MAX_DOCUMENT_TEXT_CHARS:
                    return "document_text_budget_exceeded"
                stack.append((child, depth + 1))
        elif type(current) is list:
            identity = id(current)
            if identity in seen_containers:
                return "document_cycle_or_alias_detected"
            seen_containers.add(identity)
            stack.extend((child, depth + 1) for child in current)
        else:
            return "document_non_plain_json_value"
    return None


def _validate_sha256(name: str, value: Any) -> str:
    if type(value) is not str or _SHA256_RE.fullmatch(value) is None:
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def _validate_precision(precision_bits: Any) -> int:
    if (
        type(precision_bits) is not int
        or not MIN_PRECISION_BITS <= precision_bits <= MAX_ARB_PRECISION_BITS
    ):
        raise ValueError("Arb precision lies outside the certificate resource contract")
    return precision_bits


def _validated_int64_payload(
    digital_samples: np.ndarray,
    *,
    expected_count: int,
) -> np.ndarray:
    if type(digital_samples) is not np.ndarray:
        raise ValueError("digital payload must be an exact numpy ndarray")
    if digital_samples.ndim != 1 or digital_samples.size != expected_count:
        raise ValueError("digital payload does not match the selected sample count")
    if not 1 <= digital_samples.size <= MAX_SAMPLES:
        raise ValueError("digital payload exceeds the sample resource contract")
    if digital_samples.dtype == np.bool_ or digital_samples.dtype.kind not in "iu":
        raise ValueError("digital payload must have an exact integer dtype, not bool")
    converted = np.asarray(digital_samples, dtype=np.int64)
    if not np.array_equal(digital_samples, converted):
        raise ValueError("digital payload values must lie in the int64 domain")
    return np.array(converted, dtype="<i8", order="C", copy=True)


def _int64_payload_sha256(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values, dtype="<i8").tobytes(order="C")).hexdigest()


def _dyadic_payload(value: float) -> dict[str, str]:
    if type(value) is not float or not math.isfinite(value):
        raise ValueError("a frozen binary64 value must be finite and non-boolean")
    exact = Fraction.from_float(value)
    return {"numerator": str(exact.numerator), "denominator": str(exact.denominator)}


def _arb_binary64(value: float) -> arb:
    if type(value) is not float or not math.isfinite(value):
        raise ValueError("a frozen binary64 value must be finite and exact-type float")
    exact = Fraction.from_float(value)
    return arb(exact.numerator) / exact.denominator


def _acb_binary64(real: float, imag: float) -> acb:
    return acb(_arb_binary64(real), _arb_binary64(imag))


def _outward_dyadic_upper(value: arb) -> dict[str, Any] | None:
    """Extract and re-prove a finite binary64/dyadic upper bound."""

    if not value.is_finite():
        return None
    try:
        candidate = math.nextafter(float(value.upper()), math.inf)
    except (OverflowError, ValueError):
        return None
    if not math.isfinite(candidate) or candidate < 0.0:
        return None
    for _ in range(64):
        exact = Fraction.from_float(candidate)
        claim = arb(exact.numerator) / exact.denominator
        if value <= claim:
            return {
                "float_upper": candidate,
                "dyadic": {
                    "numerator": str(exact.numerator),
                    "denominator": str(exact.denominator),
                },
            }
        candidate = math.nextafter(candidate, math.inf)
        if not math.isfinite(candidate):
            return None
    return None


def _implementation_identity() -> dict[str, str]:
    source_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    payload = {
        "algorithm": ALGORITHM,
        "implementation_source_sha256": source_sha256,
    }
    return payload | {"implementation_sha256": canonical_hash(payload)}


def _runtime_environment() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "python_flint": flint.__version__,
        "flint": flint.__FLINT_VERSION__,
    }


def _complex_row(row: Any, *, name: str) -> tuple[float, float]:
    if type(row) is not dict:
        raise ValueError(f"{name} must be a plain mapping")
    real = row.get("real")
    imag = row.get("imag")
    if (
        type(real) is not float
        or type(imag) is not float
        or not math.isfinite(real)
        or not math.isfinite(imag)
    ):
        raise ValueError(f"{name} components must be finite binary64 values")
    return real, imag


def _extract_execution_contract(execution: dict[str, Any]) -> dict[str, Any]:
    """Extract only fields already proved canonical by execution replay."""

    window = execution["window"]
    positive = tuple(window["positive_harmonics"])
    retained = retained_harmonic_set(positive)
    if tuple(window["retained_harmonics"]) != retained:
        raise ValueError("execution retained basis is not the canonical complete two-sided set")
    if len(retained) > MAX_RETAINED_HARMONICS:
        raise ValueError("retained basis exceeds the resource contract")

    start = window["start_sample_inclusive"]
    end = window["end_sample_exclusive"]
    rr_count = window["rr_interval_count"]
    baseline = window["baseline"]
    if any(type(value) is not int for value in (start, end, rr_count, baseline)):
        raise ValueError("execution integer fields must be exact integers")
    if not 0 <= start < end <= np.iinfo(np.int64).max or rr_count <= 0:
        raise ValueError("execution window lies outside the integer resource contract")
    indices = bounded_sample_indices(start, end)
    if int(window["sample_count"]) != indices.size:
        raise ValueError("execution sample count conflicts with deterministic selection")

    adc_gain = window["adc_gain"]
    if type(adc_gain) is not float or not math.isfinite(adc_gain) or adc_gain <= 0.0:
        raise ValueError("execution ADC gain must be finite positive binary64")

    front_end_rows = execution["front_end"]["H_hat_by_harmonic"]
    center_rows = execution["result"]["z_tilde"]
    if type(front_end_rows) is not list or type(center_rows) is not list:
        raise ValueError("front-end and center rows must be lists")
    front_end: dict[int, tuple[float, float]] = {}
    centers: dict[int, tuple[float, float]] = {}
    for row in front_end_rows:
        if type(row) is not dict or type(row.get("harmonic")) is not int:
            raise ValueError("front-end harmonic rows are malformed")
        harmonic = row["harmonic"]
        if harmonic in front_end:
            raise ValueError("front-end harmonic rows must be unique")
        front_end[harmonic] = _complex_row(row, name="front-end response")
    for row in center_rows:
        if type(row) is not dict or type(row.get("harmonic")) is not int:
            raise ValueError("center harmonic rows are malformed")
        harmonic = row["harmonic"]
        if harmonic in centers:
            raise ValueError("center harmonic rows must be unique")
        centers[harmonic] = _complex_row(row, name="stored response center")
    if set(front_end) != set(positive) or set(centers) != set(positive):
        raise ValueError("front-end and center rows must cover exactly the positive harmonics")
    if any(real == 0.0 and imag == 0.0 for real, imag in front_end.values()):
        raise ValueError("front-end responses must be nonzero")

    return {
        "positive": positive,
        "retained": retained,
        "start": start,
        "end": end,
        "rr_count": rr_count,
        "baseline": baseline,
        "indices": indices,
        "adc_gain": adc_gain,
        "front_end": front_end,
        "centers": centers,
    }


def _direct_interval_solve(
    contract: dict[str, Any],
    digital: np.ndarray,
) -> tuple[acb_mat, acb, dict[str, acb]]:
    """Directly accumulate exact-functional G/b and solve with error bounds."""

    positive: tuple[int, ...] = contract["positive"]
    retained: tuple[int, ...] = contract["retained"]
    indices: np.ndarray = contract["indices"]
    start: int = contract["start"]
    span = contract["end"] - start
    rr_count: int = contract["rr_count"]
    baseline: int = contract["baseline"]
    gain = _arb_binary64(contract["adc_gain"])
    maximum_harmonic = positive[-1]
    maximum_difference = 2 * maximum_harmonic
    moment_sums = {
        difference: acb(0) for difference in range(-maximum_difference, maximum_difference + 1)
    }
    rhs_sums = {harmonic: acb(0) for harmonic in retained}

    for sample_index, digital_value in zip(indices, digital, strict=True):
        phase = arb.pi() * (2 * rr_count * (int(sample_index) - start)) / span
        unit = acb(0, phase).exp()
        powers: dict[int, acb] = {0: acb(1)}
        for harmonic in range(1, maximum_difference + 1):
            powers[harmonic] = powers[harmonic - 1] * unit
            powers[-harmonic] = powers[harmonic].conjugate()
        for difference in moment_sums:
            moment_sums[difference] += powers[difference]
        physical = arb(int(digital_value) - baseline) / gain
        for harmonic in retained:
            rhs_sums[harmonic] += powers[-harmonic] * physical

    count = digital.size
    moments = {difference: total / count for difference, total in moment_sums.items()}
    gram = acb_mat([[moments[right - left] for right in retained] for left in retained])
    rhs = acb_mat([[rhs_sums[harmonic] / count] for harmonic in retained])
    determinant = gram.det()
    if not determinant.is_finite() or determinant.contains(0):
        raise ArithmeticError("interval Gram determinant does not exclude zero")
    solution = gram.solve(rhs)
    if any(not solution[index, 0].is_finite() for index in range(len(retained))):
        raise ArithmeticError("interval WLS solution is nonfinite")

    responses: dict[str, acb] = {}
    for harmonic in positive:
        real, imag = contract["front_end"][harmonic]
        response = _acb_binary64(real, imag)
        if response.contains(0):
            raise ArithmeticError("front-end response interval contains zero")
        retained_index = retained.index(harmonic)
        responses[str(harmonic)] = solution[retained_index, 0] / response
        if not responses[str(harmonic)].is_finite():
            raise ArithmeticError("front-end-corrected response is nonfinite")
    return solution, determinant, responses


def _certified_result(
    contract: dict[str, Any],
    digital: np.ndarray,
    *,
    precision_bits: int,
) -> tuple[str, str, list[dict[str, Any]], dict[str, Any]]:
    try:
        solution, determinant, responses = _direct_interval_solve(contract, digital)
        disks: list[dict[str, Any]] = []
        for harmonic in contract["positive"]:
            center_real, center_imag = contract["centers"][harmonic]
            center = _acb_binary64(center_real, center_imag)
            error = (responses[str(harmonic)] - center).abs_upper()
            radius = _outward_dyadic_upper(error)
            if radius is None:
                raise ArithmeticError("could not extract an outward dyadic response radius")
            disks.append(
                {
                    "harmonic": harmonic,
                    "center_binary64": {"real": center_real, "imag": center_imag},
                    "absolute_error_upper": radius,
                    "radius_floor_semantics": (
                        "when this certificate is used as numerical evidence, the declared "
                        "frozen-sample-functional radius centered at this stored z_tilde must "
                        "be no smaller than this published certified upper bound; the bound is "
                        "sufficient and is not asserted to be the minimal true error"
                    ),
                }
            )
        proof = {
            "precision_bits": precision_bits,
            "gram_shape": [len(contract["retained"]), len(contract["retained"])],
            "determinant_zero_excluded": not determinant.contains(0),
            "all_solution_components_finite": all(
                solution[index, 0].is_finite() for index in range(len(contract["retained"]))
            ),
        }
        return "CERTIFIED_INTERVAL", "exact_integer_payload_joint_wls_interval", disks, proof
    except (ArithmeticError, OverflowError, ZeroDivisionError) as error:
        return (
            "UNKNOWN",
            f"interval_computation_undecided:{type(error).__name__}",
            [],
            {
                "precision_bits": precision_bits,
                "gram_shape": [len(contract["retained"]), len(contract["retained"])],
                "determinant_zero_excluded": False,
                "all_solution_components_finite": False,
            },
        )


def build_wls_interval_certificate(
    execution_record: dict[str, Any],
    digital_samples: np.ndarray,
    *,
    precision_bits: int = DEFAULT_PRECISION_BITS,
) -> dict[str, Any]:
    """Issue a sealed exact-functional interval certificate.

    The supplied execution record must first replay as ``MATCH`` against the
    same selected digital payload.  Validation failures raise ``ValueError``;
    an unresolved interval solve returns a sealed ``UNKNOWN`` document.
    """

    precision = _validate_precision(precision_bits)
    if type(execution_record) is not dict:
        raise ValueError("execution record must be a plain mapping")
    execution_tree_error = _plain_tree_error(execution_record)
    if execution_tree_error is not None:
        raise ValueError(f"execution record is not an exact plain tree: {execution_tree_error}")
    if type(digital_samples) is not np.ndarray:
        raise ValueError("digital payload must be an exact numpy ndarray")
    digital = _validated_int64_payload(
        digital_samples,
        expected_count=int(digital_samples.size),
    )
    execution_replay = replay_wls_execution_record(execution_record, digital)
    if execution_replay.get("status") != "MATCH":
        raise ValueError("execution replay must be MATCH before interval issuance")
    execution_sha256 = _validate_sha256(
        "execution_sha256", execution_record.get("execution_sha256")
    )
    contract = _extract_execution_contract(execution_record)
    digital = _validated_int64_payload(
        digital,
        expected_count=int(execution_record["window"]["sample_count"]),
    )
    digital_sha256 = _int64_payload_sha256(digital)
    if digital_sha256 != execution_record["payload_hashes"]["digital_sample_sha256"]:
        raise ValueError("digital payload hash conflicts with the replayed execution")

    with _wls_execution._ARB_LOCK:
        previous_precision = ctx.prec
        try:
            ctx.prec = precision
            relation, reason, disks, numeric_proof = _certified_result(
                contract,
                digital,
                precision_bits=precision,
            )
        finally:
            ctx.prec = previous_precision

    front_end_dyadics = []
    for harmonic in contract["positive"]:
        real, imag = contract["front_end"][harmonic]
        front_end_dyadics.append(
            {
                "harmonic": harmonic,
                "real": _dyadic_payload(real),
                "imag": _dyadic_payload(imag),
            }
        )
    document: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "relation": relation,
        "reason": reason,
        "binding": {
            "execution_sha256": execution_sha256,
            "protocol_sha256": execution_record["binding"]["protocol_sha256"],
            "record_sha256": execution_record["binding"]["record_sha256"],
            "estimator_sha256": execution_record["binding"]["estimator_sha256"],
            "window_sha256": execution_record["binding"]["window_sha256"],
            "digital_sample_sha256": digital_sha256,
        },
        "implementation": _implementation_identity(),
        "environment": _runtime_environment(),
        "mathematical_contract": {
            "sample_payload": "selected exact signed int64 digital values in deterministic order",
            "physical_samples": "y_n=(digital_n-baseline)/exact_binary64(adc_gain)",
            "phase": "phi_n=2*pi*N*(sample_index_n-a)/(b-a)",
            "weights": "uniform exact rational w_n=1/sample_count",
            "basis": "complete contiguous two-sided Psi_(n,k)=exp(i*k*phi_n)",
            "statistics": "G=Psi^* W Psi and b=Psi^* W y accumulated in ACB",
            "solve": "ACB verified linear solve followed by exact_binary64(H_hat_m) division",
            "positive_harmonics": list(contract["positive"]),
            "retained_harmonics": list(contract["retained"]),
            "sample_count": int(digital.size),
            "start_sample_inclusive": contract["start"],
            "end_sample_exclusive": contract["end"],
            "rr_interval_count": contract["rr_count"],
            "baseline": contract["baseline"],
            "adc_gain_exact_dyadic": _dyadic_payload(contract["adc_gain"]),
            "front_end_exact_dyadics": front_end_dyadics,
        },
        "numeric_proof": numeric_proof,
        "response_disks": disks,
        "coverage": [
            "digital_to_physical_binary64_conversion_rounding",
            "phase_and_complex_exponential_rounding",
            "G_and_b_direct_accumulation_rounding",
            "verified_joint_linear_solve_rounding",
            "front_end_binary64_division_rounding",
        ],
        "exclusions": [
            "acquisition_and_sensor_error",
            "beat_annotation_and_window_timing_error",
            "sampling_interpolation_and_quadrature_error",
            "anti_alias_error",
            "front_end_calibration_uncertainty",
            "model_discrepancy_and_physiological_interpretation",
        ],
        "claim_scope": "frozen_mathematical_sample_functional_only",
    }
    document["certificate_sha256"] = canonical_hash(document)
    return document


def replay_wls_interval_certificate(
    document: dict[str, Any],
    execution_record: dict[str, Any],
    digital_samples: np.ndarray,
) -> dict[str, Any]:
    """Independently rebuild a sealed interval certificate and compare it."""

    tree_error = _plain_tree_error(document)
    if tree_error is not None:
        return {"status": "UNKNOWN", "reason": tree_error}
    try:
        if type(document) is not dict:
            raise ValueError("certificate must be a plain mapping")
        unsigned = dict(document)
        certificate_sha256 = _validate_sha256(
            "certificate_sha256", unsigned.pop("certificate_sha256")
        )
        if canonical_hash(unsigned) != certificate_sha256:
            return {"status": "MISMATCH", "reason": "stored_certificate_hash_mismatch"}
        if unsigned.get("schema_version") != SCHEMA_VERSION:
            return {"status": "UNKNOWN", "reason": "unsupported_schema_version"}
        precision = _validate_precision(unsigned["numeric_proof"]["precision_bits"])
        rebuilt = build_wls_interval_certificate(
            execution_record,
            digital_samples,
            precision_bits=precision,
        )
    except (
        IndexError,
        KeyError,
        OverflowError,
        RecursionError,
        TypeError,
        ValueError,
    ) as error:
        return {
            "status": "UNKNOWN",
            "reason": f"malformed_or_unreplayable_certificate:{type(error).__name__}",
        }
    if rebuilt["certificate_sha256"] != certificate_sha256:
        mismatched_fields = [
            field
            for field in (
                "relation",
                "reason",
                "binding",
                "implementation",
                "environment",
                "mathematical_contract",
                "numeric_proof",
                "response_disks",
                "coverage",
                "exclusions",
                "claim_scope",
            )
            if rebuilt.get(field) != document.get(field)
        ]
        return {
            "status": "MISMATCH",
            "reason": "recomputed_certificate_differs",
            "mismatched_fields": mismatched_fields,
            "recomputed_certificate_sha256": rebuilt["certificate_sha256"],
        }
    return {
        "status": "MATCH",
        "relation": document["relation"],
        "reason": document["reason"],
        "certificate_sha256": certificate_sha256,
    }


def certified_radius_floors(
    document: dict[str, Any],
    execution_record: dict[str, Any],
    digital_samples: np.ndarray,
) -> dict[int, Fraction]:
    """Replay first, then return all proved dyadic radius floors."""

    replay = replay_wls_interval_certificate(document, execution_record, digital_samples)
    if replay.get("status") != "MATCH" or replay.get("relation") != "CERTIFIED_INTERVAL":
        raise ValueError("a replayed CERTIFIED_INTERVAL document is required")
    floors: dict[int, Fraction] = {}
    for row in document["response_disks"]:
        harmonic = row["harmonic"]
        if type(harmonic) is not int or harmonic <= 0 or harmonic in floors:
            raise ValueError("certified harmonic rows must be exact, positive, and unique")
        dyadic = row["absolute_error_upper"]["dyadic"]
        numerator = int(dyadic["numerator"])
        denominator = int(dyadic["denominator"])
        value = Fraction(numerator, denominator)
        if value < 0 or denominator <= 0:
            raise ValueError("stored dyadic radius is invalid")
        floors[harmonic] = value
    if not floors:
        raise ValueError("certified response disks are empty")
    return floors
