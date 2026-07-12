"""Replayable joint complex-WLS records for external ECG payloads.

The estimator follows the frozen R2 convention: one whole-window nominal
phase, a complete retained two-sided Fourier basis (including DC), fixed raw
weights normalized once, the joint solve ``G C = b``, and coefficient-level
front-end correction ``z_tilde_m = C_m / H_hat_m``.  Records bind execution
without copying raw waveforms.  A match proves execution consistency only,
never an R2 error certificate.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import threading
from fractions import Fraction
from functools import lru_cache
from numbers import Complex
from pathlib import Path
from typing import Any

import numpy as np

SCHEMA_VERSION = "fo-ekf-whole-window-joint-complex-wls-execution-v3"
MAX_SAMPLES = 4096
MAX_RETAINED_HARMONICS = 64
MAX_ARB_PRECISION_BITS = 4096
MAX_DOCUMENT_DEPTH = 32
MAX_DOCUMENT_NODES = 20_000
MAX_DOCUMENT_TEXT_CHARS = 1_000_000
ESTIMATOR_SPECIFICATION = {
    "algorithm": "whole_window_two_sided_joint_complex_wls",
    "schema_version": SCHEMA_VERSION,
    "numeric_dtype": "IEEE-754-binary64-little-endian",
    "phase": "phi_n=2*pi*N*(sample_n-a)/(b-a), N=RR_count, n in [a,b)",
    "sample_interval": "half-open [first WFDB beat, last WFDB beat)",
    "retained_basis": "K={-h_max,...,-1,0,1,...,h_max}; Psi_(n,k)=exp(i*k*phi_n)",
    "weights": "fixed raw identity weights, normalized once so sum_n w_n=1",
    "statistics": "G=Psi^* W Psi; b=Psi^* W y",
    "solver": "C_tilde=numpy.linalg.solve(G,b)",
    "front_end_correction": "z_tilde_m=C_tilde_m/H_hat_m",
    "coefficient_convention": (
        "y=A*cos(m*phi)+B*sin(m*phi) implies positive-frequency C_m=(A-i*B)/2 and C_-m=conj(C_m)"
    ),
}
_ARB_LOCK = threading.RLock()


def _canonical_bytes(payload: Any) -> bytes:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")


def canonical_hash(payload: Any) -> str:
    """Return a canonical SHA-256 JSON digest."""

    return hashlib.sha256(_canonical_bytes(payload)).hexdigest()


def _array_hash(values: np.ndarray, dtype: str) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype=dtype))
    return hashlib.sha256(array.tobytes(order="C")).hexdigest()


def _complex_payload(value: complex) -> dict[str, float]:
    return {"real": float(value.real), "imag": float(value.imag)}


def _runtime_environment() -> dict[str, str | None]:
    try:
        import flint
    except ImportError:
        python_flint = None
        flint_version = None
    else:
        python_flint = flint.__version__
        flint_version = flint.__FLINT_VERSION__
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "python_flint": python_flint,
        "flint": flint_version,
    }


def _exact_int_tuple(name: str, values: tuple[int, ...]) -> tuple[int, ...]:
    if not isinstance(values, tuple) or not values:
        raise ValueError(f"{name} must be a nonempty tuple")
    if any(type(value) is not int for value in values):
        raise ValueError(f"{name} values must be exact integers, not bools")
    return values


def _plain_tree_error(value: Any) -> str | None:
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
        elif current is None or type(current) in {bool, int, float}:
            continue
        elif isinstance(current, dict):
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
        elif isinstance(current, list):
            identity = id(current)
            if identity in seen_containers:
                return "document_cycle_or_alias_detected"
            seen_containers.add(identity)
            stack.extend((child, depth + 1) for child in current)
        else:
            return "document_non_plain_json_value"
    return None


def bounded_sample_indices(
    start_sample: int,
    end_sample_exclusive: int,
    *,
    maximum_samples: int = MAX_SAMPLES,
) -> np.ndarray:
    """Return the preregistered value-independent bounded sample set."""

    if type(start_sample) is not int or type(end_sample_exclusive) is not int:
        raise ValueError("sample endpoints must be exact integers")
    if type(maximum_samples) is not int or not 1 <= maximum_samples <= MAX_SAMPLES:
        raise ValueError("maximum_samples lies outside the resource contract")
    int64_max = np.iinfo(np.int64).max
    if not 0 <= start_sample < end_sample_exclusive <= int64_max:
        raise ValueError("sample endpoints lie outside the nonnegative int64 domain")
    span = end_sample_exclusive - start_sample
    if span <= 0:
        raise ValueError("sample endpoints must define a nonempty half-open span")
    if span <= maximum_samples:
        return np.arange(start_sample, end_sample_exclusive, dtype=np.int64)
    indices = np.asarray(
        [start_sample + (index * span // maximum_samples) for index in range(maximum_samples)],
        dtype=np.int64,
    )
    if np.any(np.diff(indices) <= 0):
        raise ValueError("bounded sample policy did not produce unique indices")
    return indices


@lru_cache(maxsize=1)
def estimator_identity() -> dict[str, Any]:
    """Bind the declared algorithm and the source implementing it."""

    source_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    payload = {
        "specification": ESTIMATOR_SPECIFICATION,
        "implementation_source_sha256": source_sha256,
    }
    return payload | {"estimator_sha256": canonical_hash(payload)}


def _validate_sha256(name: str, value: str) -> None:
    if type(value) is not str or len(value) != 64:
        raise ValueError(f"{name} must be a SHA-256 hex digest")
    try:
        int(value, 16)
    except ValueError as error:
        raise ValueError(f"{name} must be a SHA-256 hex digest") from error


def retained_harmonic_set(positive_harmonics: tuple[int, ...]) -> tuple[int, ...]:
    """Return the complete contiguous two-sided set required for real ECG."""

    exact = _exact_int_tuple("positive_harmonics", positive_harmonics)
    if any(value <= 0 for value in exact):
        raise ValueError("positive harmonics must be positive")
    ordered = tuple(sorted(set(exact)))
    if ordered != exact:
        raise ValueError("positive harmonics must be ordered and unique")
    if ordered != tuple(range(1, ordered[-1] + 1)):
        raise ValueError("positive harmonics must be the contiguous set 1..h_max")
    retained = tuple(range(-ordered[-1], ordered[-1] + 1))
    if len(retained) > MAX_RETAINED_HARMONICS:
        raise ValueError("retained harmonic count exceeds the resource contract")
    return retained


def whole_window_phase_and_basis(
    beat_samples: np.ndarray,
    sampling_frequency: float,
    positive_harmonics: tuple[int, ...],
    *,
    sample_indices: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, tuple[int, ...], np.ndarray]:
    """Build frozen time, nominal phase, retained indices, and complex basis."""

    raw_beats = np.asarray(beat_samples)
    if raw_beats.ndim != 1 or raw_beats.size < 2:
        raise ValueError("at least two beat samples are required")
    if raw_beats.dtype == np.bool_ or raw_beats.dtype.kind not in "iu":
        raise ValueError("beat samples must have an exact integer dtype, not bool")
    beats = np.asarray(raw_beats, dtype=np.int64)
    if not np.array_equal(raw_beats, beats):
        raise ValueError("beat samples must be exact integers")
    if np.any(np.diff(beats) <= 0):
        raise ValueError("beat samples must be strictly increasing")
    if sampling_frequency <= 0.0 or not math.isfinite(sampling_frequency):
        raise ValueError("sampling frequency must be finite and positive")
    retained = retained_harmonic_set(positive_harmonics)

    start = int(beats[0])
    end = int(beats[-1])
    rr_count = beats.size - 1
    selected = (
        bounded_sample_indices(start, end) if sample_indices is None else np.asarray(sample_indices)
    )
    if selected.ndim != 1 or selected.size == 0 or selected.size > MAX_SAMPLES:
        raise ValueError("selected sample count lies outside the resource contract")
    if selected.dtype == np.bool_ or selected.dtype.kind not in "iu":
        raise ValueError("sample indices must have an exact integer dtype, not bool")
    integer_selected = np.asarray(selected, dtype=np.int64)
    if not np.array_equal(selected, integer_selected):
        raise ValueError("selected sample indices must be exact integers")
    if (
        integer_selected[0] < start
        or integer_selected[-1] >= end
        or np.any(np.diff(integer_selected) <= 0)
    ):
        raise ValueError("selected sample indices violate the window contract")
    phase = 2.0 * np.pi * rr_count * (integer_selected - start) / (end - start)
    basis = np.column_stack([np.exp(1j * harmonic * phase) for harmonic in retained])
    time_seconds = integer_selected.astype(np.float64) / sampling_frequency
    return time_seconds, phase, retained, basis


def joint_complex_wls_sufficient_statistics(
    physical_samples: np.ndarray,
    basis: np.ndarray,
    raw_weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return normalized weights, full Gram G, RHS b, and joint solution C."""

    values = np.asarray(physical_samples, dtype=np.float64)
    matrix = np.asarray(basis, dtype=np.complex128)
    declared_weights = (
        np.ones(values.size, dtype=np.float64)
        if raw_weights is None
        else np.asarray(raw_weights, dtype=np.float64)
    )
    if matrix.ndim != 2 or values.shape != (matrix.shape[0],):
        raise ValueError("physical samples must align with the basis rows")
    if not 1 <= matrix.shape[0] <= MAX_SAMPLES:
        raise ValueError("sample count exceeds the WLS resource contract")
    if not 1 <= matrix.shape[1] <= MAX_RETAINED_HARMONICS:
        raise ValueError("retained basis exceeds the WLS resource contract")
    if declared_weights.shape != values.shape:
        raise ValueError("weights must align with samples")
    if matrix.shape[0] < matrix.shape[1]:
        raise ValueError("joint WLS design must be overdetermined")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(matrix)):
        raise ValueError("samples and basis must be finite")
    if not np.all(np.isfinite(declared_weights)) or np.any(declared_weights <= 0.0):
        raise ValueError("raw WLS weights must be finite and strictly positive")

    raw_weight_sum = float(np.sum(declared_weights))
    if not math.isfinite(raw_weight_sum) or raw_weight_sum <= 0.0:
        raise ValueError("raw WLS weight sum must be finite and positive")
    normalized_weights = declared_weights / raw_weight_sum
    weighted_basis = normalized_weights[:, None] * matrix
    gram = np.conjugate(matrix).T @ weighted_basis
    right_hand_side = np.conjugate(matrix).T @ (normalized_weights * values)
    coefficients = np.linalg.solve(gram, right_hand_side)
    if not (
        np.all(np.isfinite(gram))
        and np.all(np.isfinite(right_hand_side))
        and np.all(np.isfinite(coefficients))
    ):
        raise ValueError("joint WLS statistics or solution are nonfinite")
    return normalized_weights, gram, right_hand_side, coefficients


def _verified_arb_upper(value: Any, arb_type: Any) -> dict[str, Any] | None:
    """Extract an outward binary64/dyadic upper claim and prove containment."""

    if not value.is_finite():
        return None
    try:
        candidate = math.nextafter(float(value.upper()), math.inf)
    except (OverflowError, ValueError):
        return None
    if not math.isfinite(candidate):
        return None
    for _ in range(64):
        exact = Fraction.from_float(candidate)
        exact_ball = arb_type(exact.numerator) / exact.denominator
        if value <= exact_ball:
            return {
                "float_upper": candidate,
                "dyadic_numerator": str(exact.numerator),
                "dyadic_denominator": str(exact.denominator),
            }
        candidate = math.nextafter(candidate, math.inf)
        if not math.isfinite(candidate):
            return None
    return None


def _joint_solve_rounding_enclosure(
    gram: np.ndarray,
    right_hand_side: np.ndarray,
    coefficients: np.ndarray,
    retained: tuple[int, ...],
    positive_harmonics: tuple[int, ...],
    front_end_response: dict[int, complex],
    z_tilde: np.ndarray,
    *,
    precision_bits: int,
) -> dict[str, Any]:
    """Enclose solve/division rounding relative to recorded binary64 G/b/H."""

    try:
        import flint
        from flint import acb, acb_mat, arb, ctx
    except ImportError:
        return {
            "status": "UNAVAILABLE",
            "scope": "joint_solve_and_H_division_relative_to_recorded_binary64_G_b_H",
            "reason": "python-flint is not installed",
        }
    if type(precision_bits) is not int or not 80 <= precision_bits <= MAX_ARB_PRECISION_BITS:
        raise ValueError("Arb precision lies outside the resource contract")

    with _ARB_LOCK:
        previous = ctx.prec
        try:
            ctx.prec = precision_bits
            gram_ball = acb_mat(
                [
                    [acb(arb(float(value.real)), arb(float(value.imag))) for value in row]
                    for row in gram
                ]
            )
            rhs_ball = acb_mat(
                [[acb(arb(float(value.real)), arb(float(value.imag)))] for value in right_hand_side]
            )
            solution = gram_ball.solve(rhs_ball)
            coefficient_bounds: list[dict[str, Any]] = []
            for index, value in enumerate(coefficients):
                stored = acb(arb(float(value.real)), arb(float(value.imag)))
                claim = _verified_arb_upper((solution[index, 0] - stored).abs_upper(), arb)
                if claim is None:
                    raise ArithmeticError("could not extract a verified C_tilde upper bound")
                coefficient_bounds.append(claim)

            z_bounds: list[dict[str, Any]] = []
            for output_index, harmonic in enumerate(positive_harmonics):
                retained_index = retained.index(harmonic)
                response = front_end_response[harmonic]
                response_ball = acb(arb(float(response.real)), arb(float(response.imag)))
                exact_z = solution[retained_index, 0] / response_ball
                stored_z = z_tilde[output_index]
                stored_z_ball = acb(arb(float(stored_z.real)), arb(float(stored_z.imag)))
                claim = _verified_arb_upper((exact_z - stored_z_ball).abs_upper(), arb)
                if claim is None:
                    raise ArithmeticError("could not extract a verified z_tilde upper bound")
                z_bounds.append(claim)
        except Exception as error:  # pragma: no cover - defensive fail-closed path
            return {
                "status": "UNKNOWN",
                "scope": "joint_solve_and_H_division_relative_to_recorded_binary64_G_b_H",
                "reason": f"{type(error).__name__}: {error}",
            }
        finally:
            ctx.prec = previous
    all_bounds = coefficient_bounds + z_bounds
    maximum = max(
        all_bounds,
        key=lambda item: Fraction(int(item["dyadic_numerator"]), int(item["dyadic_denominator"])),
    )
    return {
        "status": "ENCLOSED",
        "scope": "joint_solve_and_H_division_relative_to_recorded_binary64_G_b_H",
        "precision_bits": precision_bits,
        "C_tilde_component_absolute_error_upper_bounds": coefficient_bounds,
        "z_tilde_component_absolute_error_upper_bounds": z_bounds,
        "maximum_absolute_error_upper_bound": maximum["float_upper"],
        "maximum_absolute_error_upper_bound_dyadic": {
            "numerator": maximum["dyadic_numerator"],
            "denominator": maximum["dyadic_denominator"],
        },
        "environment": {
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
        "exclusions": [
            "G_and_b_accumulation_rounding",
            "annotation_timing_error",
            "sampling_and_quadrature_error",
            "measurement_error",
            "front_end_calibration_error",
        ],
    }


def build_wls_execution_record(
    digital_samples: np.ndarray,
    beat_samples: tuple[int, ...],
    *,
    sampling_frequency: float,
    adc_gain: float,
    baseline: int,
    harmonics: tuple[int, ...],
    protocol_sha256: str,
    record_sha256: str,
    window_sha256: str,
    front_end_response_by_harmonic: dict[int, complex] | None = None,
    arb_precision_bits: int = 128,
) -> dict[str, Any]:
    """Build one frozen joint-WLS record without embedding raw samples."""

    for name, value in (
        ("protocol_sha256", protocol_sha256),
        ("record_sha256", record_sha256),
        ("window_sha256", window_sha256),
    ):
        _validate_sha256(name, value)
    if isinstance(adc_gain, bool) or adc_gain <= 0.0 or not math.isfinite(adc_gain):
        raise ValueError("ADC gain must be finite and positive")
    if isinstance(sampling_frequency, bool):
        raise ValueError("sampling frequency cannot be boolean")
    if type(baseline) is not int:
        raise ValueError("baseline must be an exact integer, not bool")
    if (
        type(arb_precision_bits) is not int
        or not 80 <= arb_precision_bits <= MAX_ARB_PRECISION_BITS
    ):
        raise ValueError("Arb precision lies outside the resource contract")

    positive_harmonics = _exact_int_tuple("harmonics", harmonics)
    retained_harmonic_set(positive_harmonics)
    exact_beats = _exact_int_tuple("beat_samples", beat_samples)
    if len(exact_beats) < 2:
        raise ValueError("at least two beat samples are required")
    if any(value < 0 or value > np.iinfo(np.int64).max for value in exact_beats):
        raise ValueError("beat samples lie outside the nonnegative int64 domain")
    beats = np.asarray(exact_beats, dtype=np.int64)
    sample_indices = bounded_sample_indices(int(beats[0]), int(beats[-1]))
    raw_digital = np.asarray(digital_samples)
    if raw_digital.dtype == np.bool_ or raw_digital.dtype.kind not in "iu":
        raise ValueError("digital samples must have an exact integer dtype, not bool")
    digital = np.asarray(raw_digital, dtype=np.int64)
    if not np.array_equal(raw_digital, digital):
        raise ValueError("digital samples must be exact integers")
    if digital.shape != sample_indices.shape:
        raise ValueError("digital payload must match the preregistered bounded sample set")
    time_seconds, phase, retained, basis = whole_window_phase_and_basis(
        beats,
        sampling_frequency,
        positive_harmonics,
        sample_indices=sample_indices,
    )
    physical = (digital.astype(np.float64) - int(baseline)) / float(adc_gain)
    raw_weights = np.ones(digital.size, dtype=np.float64)
    normalized_weights, gram, rhs, coefficients = joint_complex_wls_sufficient_statistics(
        physical, basis, raw_weights
    )
    front_end_response = (
        {harmonic: 1.0 + 0.0j for harmonic in positive_harmonics}
        if front_end_response_by_harmonic is None
        else dict(front_end_response_by_harmonic)
    )
    if any(type(harmonic) is not int for harmonic in front_end_response):
        raise ValueError("front-end harmonic keys must be exact integers, not bools")
    if set(front_end_response) != set(positive_harmonics):
        raise ValueError("front-end responses must cover exactly the positive harmonics")
    if any(
        isinstance(value, bool) or not isinstance(value, Complex)
        for value in front_end_response.values()
    ):
        raise ValueError("front-end H values must be numeric and not boolean")
    front_end_response = {
        harmonic: complex(value) for harmonic, value in front_end_response.items()
    }
    if any(
        value == 0.0 or not (math.isfinite(value.real) and math.isfinite(value.imag))
        for value in front_end_response.values()
    ):
        raise ValueError("front-end H values must be finite and nonzero")
    z_tilde = np.asarray(
        [
            coefficients[retained.index(harmonic)] / front_end_response[harmonic]
            for harmonic in positive_harmonics
        ],
        dtype=np.complex128,
    )
    if not np.all(np.isfinite(z_tilde)):
        raise ValueError("front-end-corrected targets must be finite")
    estimator = estimator_identity()

    basis_by_harmonic = [
        {
            "harmonic": harmonic,
            "payload_sha256": _array_hash(basis[:, index], "<c16"),
        }
        for index, harmonic in enumerate(retained)
    ]
    payload_hashes = {
        "sample_index_sha256": _array_hash(sample_indices, "<i8"),
        "digital_sample_sha256": _array_hash(digital, "<i8"),
        "physical_sample_sha256": _array_hash(physical, "<f8"),
        "time_seconds_sha256": _array_hash(time_seconds, "<f8"),
        "phase_sha256": _array_hash(phase, "<f8"),
        "raw_weights_sha256": _array_hash(raw_weights, "<f8"),
        "normalized_weights_sha256": _array_hash(normalized_weights, "<f8"),
        "basis_matrix_sha256": _array_hash(basis, "<c16"),
    }
    payload_hashes["combined_payload_sha256"] = canonical_hash(payload_hashes)
    conjugate_errors = [
        abs(
            coefficients[retained.index(-harmonic)]
            - np.conjugate(coefficients[retained.index(harmonic)])
        )
        for harmonic in positive_harmonics
    ]
    result = {
        "schema_version": SCHEMA_VERSION,
        "binding": {
            "protocol_sha256": protocol_sha256,
            "record_sha256": record_sha256,
            "estimator_sha256": estimator["estimator_sha256"],
            "window_sha256": window_sha256,
        },
        "environment": _runtime_environment(),
        "estimator": estimator,
        "window": {
            "start_sample_inclusive": int(beats[0]),
            "end_sample_exclusive": int(beats[-1]),
            "sample_count": int(digital.size),
            "full_span_sample_count": int(beats[-1] - beats[0]),
            "sample_selection": (
                "all half-open samples when span<=4096; otherwise "
                "a+floor(j*(b-a)/4096), j=0,...,4095"
            ),
            "rr_interval_count": int(beats.size - 1),
            "beat_samples": [int(value) for value in beats],
            "sampling_frequency_hz": float(sampling_frequency),
            "adc_gain": float(adc_gain),
            "baseline": int(baseline),
            "positive_harmonics": list(positive_harmonics),
            "retained_harmonics": list(retained),
        },
        "payload_hashes": payload_hashes,
        "weights": {
            "raw_definition": "identity raw weights w_raw_n=1",
            "normalization": "w_n=w_raw_n/sum_j(w_raw_j)",
            "count": int(raw_weights.size),
            "raw_sum": float(np.sum(raw_weights)),
            "normalized_sum": float(np.sum(normalized_weights)),
            "raw_payload_sha256": payload_hashes["raw_weights_sha256"],
            "normalized_payload_sha256": payload_hashes["normalized_weights_sha256"],
        },
        "basis": {
            "phase_definition": ESTIMATOR_SPECIFICATION["phase"],
            "retained_basis_definition": ESTIMATOR_SPECIFICATION["retained_basis"],
            "shape": [int(value) for value in basis.shape],
            "matrix_payload_sha256": payload_hashes["basis_matrix_sha256"],
            "by_harmonic": basis_by_harmonic,
        },
        "sufficient_statistics": {
            "definition": ESTIMATOR_SPECIFICATION["statistics"],
            "G": [[_complex_payload(value) for value in row] for row in gram],
            "b": [_complex_payload(value) for value in rhs],
            "G_sha256": _array_hash(gram, "<c16"),
            "b_sha256": _array_hash(rhs, "<c16"),
        },
        "front_end": {
            "mode": (
                "identity"
                if all(value == 1.0 + 0.0j for value in front_end_response.values())
                else "caller_supplied"
            ),
            "H_hat_by_harmonic": [
                {"harmonic": harmonic, **_complex_payload(front_end_response[harmonic])}
                for harmonic in positive_harmonics
            ],
        },
        "result": {
            "joint_formula": ESTIMATOR_SPECIFICATION["solver"],
            "front_end_formula": ESTIMATOR_SPECIFICATION["front_end_correction"],
            "coefficient_convention": ESTIMATOR_SPECIFICATION["coefficient_convention"],
            "C_tilde": [
                {"harmonic": harmonic, **_complex_payload(value)}
                for harmonic, value in zip(retained, coefficients, strict=True)
            ],
            "C_tilde_sha256": _array_hash(coefficients, "<c16"),
            "z_tilde": [
                {"harmonic": harmonic, **_complex_payload(value)}
                for harmonic, value in zip(positive_harmonics, z_tilde, strict=True)
            ],
            "z_tilde_sha256": _array_hash(z_tilde, "<c16"),
            "maximum_conjugate_symmetry_error": float(max(conjugate_errors)),
            "joint_solve_rounding_enclosure": _joint_solve_rounding_enclosure(
                gram,
                rhs,
                coefficients,
                retained,
                positive_harmonics,
                front_end_response,
                z_tilde,
                precision_bits=arb_precision_bits,
            ),
            "certification_status": "NOT_CERTIFIABLE",
        },
    }
    result["execution_sha256"] = canonical_hash(result)
    return result


def replay_wls_execution_record(
    document: dict[str, Any],
    digital_samples: np.ndarray,
) -> dict[str, Any]:
    """Rebuild one record from external samples and compare every bound field."""

    tree_error = _plain_tree_error(document)
    if tree_error is not None:
        return {"status": "UNKNOWN", "reason": tree_error}
    try:
        stored = dict(document)
        execution_sha256 = str(stored.pop("execution_sha256"))
        _validate_sha256("execution_sha256", execution_sha256)
        if canonical_hash(stored) != execution_sha256:
            return {
                "status": "MISMATCH",
                "reason": "stored_execution_hash_mismatch",
            }
        if stored.get("schema_version") != SCHEMA_VERSION:
            return {"status": "UNKNOWN", "reason": "unsupported_schema_version"}
        binding = stored["binding"]
        window = stored["window"]
        front_end_rows = stored["front_end"]["H_hat_by_harmonic"]
        if not isinstance(front_end_rows, list):
            raise ValueError("front-end response rows must be a list")
        front_end: dict[int, complex] = {}
        for row in front_end_rows:
            harmonic = row["harmonic"]
            if type(harmonic) is not int or harmonic in front_end:
                raise ValueError("front-end harmonic keys must be exact and unique")
            front_end[harmonic] = complex(row["real"], row["imag"])
        precision = stored["result"]["joint_solve_rounding_enclosure"].get("precision_bits", 128)
        rebuilt = build_wls_execution_record(
            digital_samples,
            tuple(window["beat_samples"]),
            sampling_frequency=float(window["sampling_frequency_hz"]),
            adc_gain=float(window["adc_gain"]),
            baseline=window["baseline"],
            harmonics=tuple(window["positive_harmonics"]),
            protocol_sha256=binding["protocol_sha256"],
            record_sha256=binding["record_sha256"],
            window_sha256=binding["window_sha256"],
            front_end_response_by_harmonic=front_end,
            arb_precision_bits=precision,
        )
    except (
        IndexError,
        KeyError,
        OverflowError,
        RecursionError,
        TypeError,
        ValueError,
        np.linalg.LinAlgError,
    ) as error:
        return {
            "status": "UNKNOWN",
            "reason": f"malformed_or_unreplayable_record:{type(error).__name__}",
        }
    if rebuilt["execution_sha256"] != execution_sha256:
        mismatches = []
        for field in (
            "binding",
            "environment",
            "payload_hashes",
            "weights",
            "basis",
            "sufficient_statistics",
            "front_end",
            "result",
        ):
            if rebuilt.get(field) != document.get(field):
                mismatches.append(field)
        return {
            "status": "MISMATCH",
            "reason": "recomputed_execution_differs",
            "mismatched_fields": mismatches,
            "recomputed_execution_sha256": rebuilt["execution_sha256"],
        }
    return {
        "status": "MATCH",
        "execution_sha256": execution_sha256,
        "certification_status": "NOT_CERTIFIABLE",
    }


def validate_execution_manifest_binding(
    wrapper: dict[str, Any],
    manifest_row: dict[str, Any],
    *,
    external_record_sha256: str,
    protocol_sha256: str,
) -> list[str]:
    """Cross-check one serialized execution against source and manifest bindings."""

    failures: list[str] = []
    try:
        execution = wrapper["execution"]
        binding = execution["binding"]
        window = execution["window"]
        comparisons = {
            "wrapper_record": (wrapper["record"], manifest_row["record"]),
            "wrapper_role": (wrapper["role"], manifest_row["role"]),
            "wrapper_split": (wrapper["split"], manifest_row["split"]),
            "wrapper_window_index": (
                int(wrapper["window_index"]),
                int(manifest_row["window_index"]),
            ),
            "binding_record_sha256": (
                binding["record_sha256"],
                external_record_sha256,
            ),
            "binding_protocol_sha256": (
                binding["protocol_sha256"],
                protocol_sha256,
            ),
            "binding_window_sha256": (
                binding["window_sha256"],
                manifest_row["window_sha256"],
            ),
            "execution_sha256": (
                execution["execution_sha256"],
                manifest_row["wls_execution_sha256"],
            ),
            "start_sample": (
                int(window["start_sample_inclusive"]),
                int(manifest_row["start_sample"]),
            ),
            "end_sample": (
                int(window["end_sample_exclusive"]),
                int(manifest_row["end_sample_exclusive"]),
            ),
            "selected_sample_count": (
                int(window["sample_count"]),
                int(manifest_row["sample_count"]),
            ),
            "full_span_sample_count": (
                int(window["full_span_sample_count"]),
                int(manifest_row["full_span_sample_count"]),
            ),
        }
    except (KeyError, TypeError, ValueError) as error:
        return [f"malformed_cross_binding:{type(error).__name__}"]
    for name, (actual, expected) in comparisons.items():
        if actual != expected:
            failures.append(name)
    return failures
