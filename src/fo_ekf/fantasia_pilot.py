"""Locked Fantasia feasibility pipeline for subject-specific rate response.

The module deliberately produces *empirical* held-out evidence only.  It does
not turn the Fantasia annotations or residuals into certified R2 primitives;
the corresponding output is therefore fail-closed as ``NOT_CERTIFIABLE``.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import platform
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from fo_ekf import wls_execution as _wls_execution
from fo_ekf.wls_execution import (
    bounded_sample_indices,
    build_wls_execution_record,
    joint_complex_wls_sufficient_statistics,
    replay_wls_execution_record,
    whole_window_phase_and_basis,
)

RECORD_ROLES = {
    "f1y01": "engineering_calibration",
    "f1o01": "engineering_calibration",
    "f2y01": "locked_engineering_validation",
    "f2o01": "locked_engineering_validation",
}

# Checked against the official Fantasia v1.0.0 SHA256SUMS.txt.  The pilot
# refuses a same-name file from another database revision.
FANTASIA_V1_SHA256 = {
    "f1y01.hea": "519eff62d1ad082ea1f3ed68e4e90dc474dc35f7c99757b1d00cdaf1738bc9c2",
    "f1y01.dat": "c992e4e69b1869e2da1f94e00a8cca070bc099cfc58911dd5a65aa436003b9bd",
    "f1y01.ecg": "7f1da797a9292ca296ac93affec2bf16f43f0c88cf2a2f5850a3336c3b2d59f0",
    "f1o01.hea": "c0ad734a39b077e5fe5bb3a0f910aa5128f2655f57247219bf02b734a2f0c2ea",
    "f1o01.dat": "5758cb184d4f529383fdafd4dca68478d654e511f15621ff00d04a17c01f80aa",
    "f1o01.ecg": "480a47fd00aba1cf11129e67160b2349afb63b573214a63b89641d14f7a3c808",
    "f2y01.hea": "4074b8fe12804229305d6e84c92453205dab5cce704856dadfc16e94ad6435a0",
    "f2y01.dat": "aec272b391cb2f43785272901233ce87d10ae72486c53108f78a365cecaa7a63",
    "f2y01.ecg": "e8632c4adf93a08819d126781a489287c6687d894ca4321bc6adc3b3ac5cc75d",
    "f2o01.hea": "3d57d614f597f324e5798dca285fc4029c2afb036b65925a60f233d5512acb6d",
    "f2o01.dat": "ee9f7a8f92a4e8d9aa8f209ba3e4d112a7b4090ea08f0de4be02b1656646199a",
    "f2o01.ecg": "2eb9fd71ef5fc197782ce510358c42233db8192f73d697f26e01c57fd276dbe4",
}


@dataclass(frozen=True)
class FantasiaPilotConfig:
    """Immutable choices fixed before held-out engineering evaluation."""

    records: tuple[str, ...] = ("f1y01", "f1o01", "f2y01", "f2o01")
    identification_end_seconds: float = 30.0 * 60.0
    guard_end_seconds: float = 35.0 * 60.0
    rr_intervals_per_window: int = 16
    harmonics: tuple[int, ...] = (1, 2, 3)
    beat_symbols: tuple[str, ...] = ("N", "S", "V")
    ignored_nonbeat_symbols: tuple[str, ...] = ("+",)
    alpha_bounds: tuple[float, float] = (0.05, 1.0)
    lambda_bounds: tuple[float, float] = (1.0e-3, 1.0e3)
    tau_star_seconds: float = 1.0
    shuffle_seed: int = 20260712
    front_end: str = "identity"
    phase_anchor: str = "whole_window_linear_WFDB_endpoint_R_annotations"
    normalization: str = "none"
    wls_weights: str = "raw_identity_normalized_once"
    sample_selection: str = "all_if_span_le_4096_else_floor_j_span_over_4096_for_j_0_to_4095"

    def validate(self) -> None:
        if tuple(self.records) != tuple(RECORD_ROLES):
            raise ValueError("the pilot record order is locked")
        if not 0.0 < self.identification_end_seconds < self.guard_end_seconds:
            raise ValueError("the guard must follow the identification interval")
        if self.guard_end_seconds - self.identification_end_seconds < 300.0:
            raise ValueError("the guard interval must be at least five minutes")
        if self.rr_intervals_per_window <= 0:
            raise ValueError("RR intervals per window must be positive")
        if not self.harmonics or any(value <= 0 for value in self.harmonics):
            raise ValueError("harmonic indices must be positive")
        if len(set(self.harmonics)) != len(self.harmonics):
            raise ValueError("harmonic indices must be unique")
        if not 0.0 < self.alpha_bounds[0] < self.alpha_bounds[1] <= 1.0:
            raise ValueError("alpha bounds must lie in (0, 1]")
        if not 0.0 < self.lambda_bounds[0] < self.lambda_bounds[1]:
            raise ValueError("lambda bounds must be positive and ordered")
        if not math.isfinite(self.tau_star_seconds) or self.tau_star_seconds <= 0.0:
            raise ValueError("tau_star_seconds must be finite and positive")
        if self.front_end != "identity":
            raise ValueError("this pilot permits only the identity front end")
        if self.normalization != "none" or self.wls_weights != "raw_identity_normalized_once":
            raise ValueError("window normalization or non-identity weights are not permitted")

    def canonical_payload(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["record_roles"] = RECORD_ROLES
        payload["blindness_scope"] = (
            "The code/configuration hash is common to all records and no outcome-dependent "
            "adaptation is made on f2 records; source metadata may be known during preflight."
        )
        return payload

    def lock_hash(self) -> str:
        return canonical_json_hash(self.canonical_payload())


@dataclass(frozen=True)
class RRWindow:
    """One non-overlapping block of complete annotated RR intervals."""

    split: str
    index: int
    beat_samples: tuple[int, ...]
    beat_symbols: tuple[str, ...]

    @property
    def start_sample(self) -> int:
        return self.beat_samples[0]

    @property
    def end_sample(self) -> int:
        return self.beat_samples[-1]


@dataclass(frozen=True)
class WindowEstimate:
    """Complex harmonic WLS coefficients and numerical diagnostics."""

    dc: float
    coefficients: tuple[complex, ...]
    residual_rmse: float
    gram_min_eigenvalue: float
    gram_max_eigenvalue: float
    gram_condition_number: float
    gram_rank: int
    sample_count: int


def canonical_json_hash(payload: Any) -> str:
    """Hash JSON with stable ordering and separators."""

    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def sha256_file(path: Path, *, block_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of one explicitly selected source file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def fantasia_implementation_manifest() -> dict[str, Any]:
    """Bind the production sources that determine the derived pilot bundle."""

    repository_root = Path(__file__).resolve().parents[2]
    sources = {
        "src/fo_ekf/fantasia_pilot.py": Path(__file__).resolve(),
        "src/fo_ekf/wls_execution.py": Path(_wls_execution.__file__).resolve(),
        "scripts/run_fantasia_pilot.py": repository_root / "scripts" / "run_fantasia_pilot.py",
    }
    if any(not path.is_file() for path in sources.values()):
        raise RuntimeError("Fantasia production source manifest is incomplete")
    file_hashes = {name: sha256_file(path) for name, path in sorted(sources.items())}
    payload: dict[str, Any] = {
        "algorithm": "locked-fantasia-rate-response-pilot-v2",
        "source_sha256": file_hashes,
    }
    return payload | {"implementation_sha256": canonical_json_hash(payload)}


def fantasia_runtime_environment() -> dict[str, str]:
    """Record the numerical runtime without using it as scientific evidence."""

    import flint
    import scipy
    import wfdb

    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "wfdb": wfdb.__version__,
        "python_flint": flint.__version__,
        "flint": flint.__FLINT_VERSION__,
    }


def build_rr_windows(
    annotation_samples: Sequence[int],
    annotation_symbols: Sequence[str],
    *,
    sampling_frequency: float,
    segment_start_seconds: float,
    segment_end_seconds: float,
    split: str,
    rr_intervals_per_window: int = 16,
    beat_symbols: Iterable[str] = ("N", "S", "V"),
    ignored_nonbeat_symbols: Iterable[str] = ("+",),
) -> tuple[RRWindow, ...]:
    """Partition valid annotation runs into complete, non-overlapping RR windows.

    Unknown/event annotations break a run instead of silently bridging an RR
    interval.  Rhythm markers in ``ignored_nonbeat_symbols`` carry no beat time
    and do not break a run.  Adjacent windows may share their boundary R peak,
    but their half-open signal sample ranges never overlap.
    """

    if len(annotation_samples) != len(annotation_symbols):
        raise ValueError("annotation samples and symbols must have the same length")
    if sampling_frequency <= 0.0:
        raise ValueError("sampling frequency must be positive")
    if segment_start_seconds < 0.0 or segment_end_seconds <= segment_start_seconds:
        raise ValueError("segment bounds must be ordered and nonnegative")
    if rr_intervals_per_window <= 0:
        raise ValueError("RR interval count must be positive")

    allowed = frozenset(beat_symbols)
    ignored = frozenset(ignored_nonbeat_symbols)
    start_sample = math.ceil(segment_start_seconds * sampling_frequency)
    end_sample = math.floor(segment_end_seconds * sampling_frequency)
    runs: list[list[tuple[int, str]]] = []
    current: list[tuple[int, str]] = []

    for raw_sample, raw_symbol in zip(annotation_samples, annotation_symbols, strict=True):
        sample = int(raw_sample)
        symbol = str(raw_symbol)
        if sample < start_sample or sample > end_sample:
            continue
        if symbol in allowed:
            if current and sample <= current[-1][0]:
                raise ValueError("beat annotations must be strictly increasing")
            current.append((sample, symbol))
        elif symbol in ignored:
            continue
        elif current:
            runs.append(current)
            current = []
    if current:
        runs.append(current)

    result: list[RRWindow] = []
    width = rr_intervals_per_window
    for run in runs:
        for offset in range(0, len(run) - width, width):
            block = run[offset : offset + width + 1]
            result.append(
                RRWindow(
                    split=split,
                    index=len(result),
                    beat_samples=tuple(item[0] for item in block),
                    beat_symbols=tuple(item[1] for item in block),
                )
            )
    return tuple(result)


def cardiac_phase(beat_samples: Sequence[int]) -> np.ndarray:
    """Return the frozen whole-window linear phase on the half-open sample set."""

    boundaries = np.asarray(beat_samples, dtype=np.int64)
    if boundaries.ndim != 1 or boundaries.size < 2:
        raise ValueError("at least two beat boundaries are required")
    if np.any(np.diff(boundaries) <= 0):
        raise ValueError("beat boundaries must be strictly increasing")
    samples = np.arange(boundaries[0], boundaries[-1], dtype=np.int64)
    rr_count = boundaries.size - 1
    return 2.0 * np.pi * rr_count * (samples - boundaries[0]) / (boundaries[-1] - boundaries[0])


def estimate_window_coefficients(
    signal: np.ndarray,
    beat_samples: Sequence[int],
    harmonics: Sequence[int],
) -> WindowEstimate:
    """Estimate fixed-phase Fourier coefficients with identity-weight WLS."""

    values = np.asarray(signal, dtype=np.float64)
    boundaries = tuple(int(value) for value in beat_samples)
    if not boundaries or boundaries[0] < 0 or boundaries[-1] > values.size:
        raise ValueError("beat boundaries lie outside the signal")
    harmonic_indices = tuple(int(value) for value in harmonics)
    if not harmonic_indices or any(value <= 0 for value in harmonic_indices):
        raise ValueError("harmonics must be positive")

    selected_indices = bounded_sample_indices(boundaries[0], boundaries[-1])
    response = values[selected_indices]
    _, _phase, retained, basis = whole_window_phase_and_basis(
        np.asarray(boundaries, dtype=np.int64),
        1.0,
        harmonic_indices,
        sample_indices=selected_indices,
    )
    raw_weights = np.ones(response.size, dtype=np.float64)
    _weights, gram, _right_hand_side, retained_coefficients = (
        joint_complex_wls_sufficient_statistics(response, basis, raw_weights)
    )
    dc = float(retained_coefficients[retained.index(0)].real)
    centers = np.asarray(
        [retained_coefficients[retained.index(value)] for value in harmonic_indices]
    )
    fitted = np.real(basis @ retained_coefficients)
    residual_rmse = float(np.sqrt(np.mean(np.square(response - fitted))))
    eigenvalues = np.linalg.eigvalsh(gram)
    minimum = float(eigenvalues[0])
    maximum = float(eigenvalues[-1])
    condition = float(maximum / minimum) if minimum > 0.0 else math.inf
    return WindowEstimate(
        dc=dc,
        coefficients=tuple(complex(value) for value in centers),
        residual_rmse=residual_rmse,
        gram_min_eigenvalue=minimum,
        gram_max_eigenvalue=maximum,
        gram_condition_number=condition,
        gram_rank=int(np.linalg.matrix_rank(gram)),
        sample_count=int(basis.shape[0]),
    )


def _responses(
    base_rates_hz: np.ndarray,
    harmonics: Sequence[int],
    alpha: float,
    damping: float,
    tau_star_seconds: float,
) -> np.ndarray:
    if not math.isfinite(tau_star_seconds) or tau_star_seconds <= 0.0:
        raise ValueError("tau_star_seconds must be finite and positive")
    dimensionless_frequency = (
        2.0 * np.pi * tau_star_seconds * base_rates_hz[:, None] * np.asarray(harmonics)[None, :]
    )
    fractional = np.power(1j * dimensionless_frequency, alpha)
    return 1.0 / (damping + fractional)


def _profile_q(response_basis: np.ndarray, observations: np.ndarray) -> np.ndarray:
    numerator = np.sum(np.conjugate(response_basis) * observations, axis=0)
    denominator = np.sum(np.abs(response_basis) ** 2, axis=0)
    if np.any(denominator <= 0.0):
        raise ValueError("fractional response basis is degenerate")
    return numerator / denominator


def _fit_at_parameters(
    base_rates_hz: np.ndarray,
    observations: np.ndarray,
    harmonics: Sequence[int],
    alpha: float,
    damping: float,
    tau_star_seconds: float,
) -> tuple[float, np.ndarray, np.ndarray]:
    basis = _responses(base_rates_hz, harmonics, alpha, damping, tau_star_seconds)
    morphology = _profile_q(basis, observations)
    predicted = basis * morphology[None, :]
    loss = float(np.sum(np.abs(observations - predicted) ** 2))
    return loss, morphology, predicted


def fit_fractional_model(
    base_rates_hz: Sequence[float],
    observations: np.ndarray,
    harmonics: Sequence[int],
    *,
    alpha_bounds: tuple[float, float] = (0.05, 1.0),
    lambda_bounds: tuple[float, float] = (1.0e-3, 1.0e3),
    tau_star_seconds: float = 1.0,
) -> dict[str, Any]:
    """Fit common alpha/lambda with profiled complex morphology coefficients."""

    rates = np.asarray(base_rates_hz, dtype=np.float64)
    values = np.asarray(observations, dtype=np.complex128)
    if rates.ndim != 1 or values.shape != (rates.size, len(tuple(harmonics))):
        raise ValueError("observation matrix does not match rates and harmonics")
    if rates.size < 3 or np.any(rates <= 0.0) or not np.all(np.isfinite(values)):
        raise ValueError("at least three finite, positive-rate windows are required")

    alpha_lower, alpha_upper = alpha_bounds
    lambda_lower, lambda_upper = lambda_bounds
    log_bounds = (math.log(lambda_lower), math.log(lambda_upper))

    def objective(point: np.ndarray) -> float:
        loss, _, _ = _fit_at_parameters(
            rates,
            values,
            harmonics,
            float(point[0]),
            math.exp(float(point[1])),
            tau_star_seconds,
        )
        return loss

    alpha_grid = np.linspace(alpha_lower, alpha_upper, 17)
    log_grid = np.linspace(log_bounds[0], log_bounds[1], 25)
    candidates = sorted(
        (
            (objective(np.asarray((alpha, log_damping))), alpha, log_damping)
            for alpha in alpha_grid
            for log_damping in log_grid
        ),
        key=lambda item: item[0],
    )
    starts = candidates[:6]
    solutions = []
    for _, alpha, log_damping in starts:
        solutions.append(
            minimize(
                objective,
                x0=np.asarray((alpha, log_damping)),
                method="L-BFGS-B",
                bounds=(alpha_bounds, log_bounds),
                options={"ftol": 1.0e-14, "gtol": 1.0e-10, "maxiter": 1000},
            )
        )
    best = min(solutions, key=lambda item: float(item.fun))
    alpha = float(best.x[0])
    damping = math.exp(float(best.x[1]))
    loss, morphology, _ = _fit_at_parameters(
        rates, values, harmonics, alpha, damping, tau_star_seconds
    )
    tolerance = 1.0e-5
    return {
        "alpha": alpha,
        "lambda": damping,
        "tau_star_seconds": tau_star_seconds,
        "q": tuple(complex(value) for value in morphology),
        "identification_sse": loss,
        "optimizer_success": bool(best.success),
        "optimizer_message": str(best.message),
        "alpha_boundary_hit": (
            alpha - alpha_lower <= tolerance * max(1.0, abs(alpha_lower))
            or alpha_upper - alpha <= tolerance * max(1.0, abs(alpha_upper))
        ),
        "lambda_boundary_hit": (
            damping / lambda_lower <= 1.0 + tolerance or lambda_upper / damping <= 1.0 + tolerance
        ),
    }


def fit_order_one_model(
    base_rates_hz: Sequence[float],
    observations: np.ndarray,
    harmonics: Sequence[int],
    *,
    lambda_bounds: tuple[float, float] = (1.0e-3, 1.0e3),
    tau_star_seconds: float = 1.0,
) -> dict[str, Any]:
    """Fit the exactly nested alpha=1 model."""

    rates = np.asarray(base_rates_hz, dtype=np.float64)
    values = np.asarray(observations, dtype=np.complex128)
    log_bounds = tuple(math.log(value) for value in lambda_bounds)

    def objective(log_damping: float) -> float:
        return _fit_at_parameters(
            rates,
            values,
            harmonics,
            1.0,
            math.exp(log_damping),
            tau_star_seconds,
        )[0]

    result = minimize_scalar(
        objective,
        bounds=log_bounds,
        method="bounded",
        options={"xatol": 1.0e-12, "maxiter": 1000},
    )
    damping = math.exp(float(result.x))
    loss, morphology, _ = _fit_at_parameters(
        rates, values, harmonics, 1.0, damping, tau_star_seconds
    )
    tolerance = 1.0e-5
    return {
        "alpha": 1.0,
        "lambda": damping,
        "tau_star_seconds": tau_star_seconds,
        "q": tuple(complex(value) for value in morphology),
        "identification_sse": loss,
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "lambda_boundary_hit": (
            damping / lambda_bounds[0] <= 1.0 + tolerance
            or lambda_bounds[1] / damping <= 1.0 + tolerance
        ),
    }


def fit_constant_fourier_model(observations: np.ndarray) -> dict[str, Any]:
    """Fit one rate-independent complex coefficient per harmonic."""

    values = np.asarray(observations, dtype=np.complex128)
    if values.ndim != 2 or values.shape[0] < 1:
        raise ValueError("a nonempty window-by-harmonic matrix is required")
    morphology = np.mean(values, axis=0)
    predicted = np.broadcast_to(morphology, values.shape)
    return {
        "q": tuple(complex(value) for value in morphology),
        "identification_sse": float(np.sum(np.abs(values - predicted) ** 2)),
    }


def predict_model(
    model: dict[str, Any],
    base_rates_hz: Sequence[float],
    harmonics: Sequence[int],
    *,
    model_name: str,
) -> np.ndarray:
    """Predict complex coefficients for a fitted model."""

    rates = np.asarray(base_rates_hz, dtype=np.float64)
    morphology = np.asarray(model["q"], dtype=np.complex128)
    if model_name == "constant_fourier":
        return np.broadcast_to(morphology, (rates.size, morphology.size)).copy()
    basis = _responses(
        rates,
        harmonics,
        float(model["alpha"]),
        float(model["lambda"]),
        float(model["tau_star_seconds"]),
    )
    return basis * morphology[None, :]


def coefficient_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    """Return transparent complex-coefficient error metrics."""

    values = np.asarray(observed, dtype=np.complex128)
    estimates = np.asarray(predicted, dtype=np.complex128)
    if values.shape != estimates.shape or values.size == 0:
        raise ValueError("observed and predicted arrays must be nonempty and aligned")
    residual = values - estimates
    sse = float(np.sum(np.abs(residual) ** 2))
    energy = float(np.sum(np.abs(values) ** 2))
    centered = values - np.mean(values, axis=0, keepdims=True)
    centered_energy = float(np.sum(np.abs(centered) ** 2))
    return {
        "count_complex": int(values.size),
        "sse_mV2": sse,
        "rmse_mV": float(np.sqrt(sse / values.size)),
        "mae_mV": float(np.mean(np.abs(residual))),
        "nrmse_energy": float(np.sqrt(sse / energy)) if energy > 0.0 else math.inf,
        "r2_centered": 1.0 - sse / centered_energy if centered_energy > 0.0 else math.nan,
    }


def _complex_payload(values: Sequence[complex]) -> list[dict[str, float]]:
    return [{"real": float(value.real), "imag": float(value.imag)} for value in values]


def _serializable_model(model: dict[str, Any]) -> dict[str, Any]:
    payload = dict(model)
    payload["q"] = _complex_payload(payload["q"])
    return payload


def _window_hash(
    record_name: str,
    window: RRWindow,
    digital_signal: np.ndarray,
    config_hash: str,
) -> str:
    metadata = {
        "record": record_name,
        "split": window.split,
        "index": window.index,
        "beat_samples": window.beat_samples,
        "beat_symbols": window.beat_symbols,
        "config_sha256": config_hash,
        "sample_encoding": "little_endian_int64_half_open",
    }
    digest = hashlib.sha256()
    digest.update(json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode("utf-8"))
    raw = np.asarray(digital_signal[window.start_sample : window.end_sample], dtype="<i8")
    digest.update(raw.tobytes(order="C"))
    return digest.hexdigest()


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: Sequence[str]) -> None:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    _write_text_lf(path, buffer.getvalue())


def _write_csv_gzip(path: Path, rows: list[dict[str, Any]], fields: Sequence[str]) -> None:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    encoded = buffer.getvalue().encode("utf-8")
    with path.open("wb") as raw_stream:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            compresslevel=9,
            fileobj=raw_stream,
            mtime=0,
        ) as stream:
            stream.write(encoded)


def _write_text_lf(path: Path, text: str, *, encoding: str = "utf-8") -> None:
    """Write a sealed text artifact without platform newline translation."""

    if "\r" in text:
        raise ValueError("sealed text artifacts must use LF line endings")
    path.write_bytes(text.encode(encoding))


def _write_json(path: Path, payload: Any) -> None:
    _write_text_lf(
        path,
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n",
    )


def _write_jsonl_gzip(path: Path, payloads: Sequence[dict[str, Any]]) -> None:
    with path.open("wb") as raw_stream:
        stream = gzip.GzipFile(
            filename="",
            mode="wb",
            compresslevel=9,
            fileobj=raw_stream,
            mtime=0,
        )
        for payload in payloads:
            encoded = (
                json.dumps(
                    payload,
                    sort_keys=True,
                    separators=(",", ":"),
                    ensure_ascii=True,
                    allow_nan=False,
                )
                + "\n"
            ).encode("utf-8")
            stream.write(encoded)
        stream.close()


def _finite_json_metrics(metrics: dict[str, float]) -> dict[str, float | int | None]:
    return {
        key: (value if isinstance(value, int) or math.isfinite(value) else None)
        for key, value in metrics.items()
    }


def _evaluate_record_models(
    record_name: str,
    harmonics: tuple[int, ...],
    identification_rows: list[dict[str, Any]],
    validation_rows: list[dict[str, Any]],
    config: FantasiaPilotConfig,
    record_index: int,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    identification_rates = np.asarray(
        [row["base_rate_hz"] for row in identification_rows], dtype=np.float64
    )
    validation_rates = np.asarray(
        [row["base_rate_hz"] for row in validation_rows], dtype=np.float64
    )
    identification_values = np.asarray(
        [row["coefficients"] for row in identification_rows], dtype=np.complex128
    )
    validation_values = np.asarray(
        [row["coefficients"] for row in validation_rows], dtype=np.complex128
    )

    fractional = fit_fractional_model(
        identification_rates,
        identification_values,
        harmonics,
        alpha_bounds=config.alpha_bounds,
        lambda_bounds=config.lambda_bounds,
        tau_star_seconds=config.tau_star_seconds,
    )
    order_one = fit_order_one_model(
        identification_rates,
        identification_values,
        harmonics,
        lambda_bounds=config.lambda_bounds,
        tau_star_seconds=config.tau_star_seconds,
    )
    constant = fit_constant_fourier_model(identification_values)

    rng = np.random.default_rng(config.shuffle_seed + record_index)
    permutation = rng.permutation(identification_rates.size)
    shuffled_rates = identification_rates[permutation]
    shuffled = fit_fractional_model(
        shuffled_rates,
        identification_values,
        harmonics,
        alpha_bounds=config.alpha_bounds,
        lambda_bounds=config.lambda_bounds,
        tau_star_seconds=config.tau_star_seconds,
    )

    models = {
        "fractional_common": fractional,
        "alpha_one_nested": order_one,
        "constant_fourier": constant,
        "rate_shuffle_negative_control": shuffled,
    }
    metric_rows: list[dict[str, Any]] = []
    serializable: dict[str, Any] = {}
    for model_name, model in models.items():
        fit_rates = (
            shuffled_rates
            if model_name == "rate_shuffle_negative_control"
            else identification_rates
        )
        identification_prediction = predict_model(
            model,
            fit_rates,
            harmonics,
            model_name=(
                "fractional_common" if model_name == "rate_shuffle_negative_control" else model_name
            ),
        )
        validation_prediction = predict_model(
            model,
            validation_rates,
            harmonics,
            model_name=(
                "fractional_common" if model_name == "rate_shuffle_negative_control" else model_name
            ),
        )
        split_metrics = {
            "identification": coefficient_metrics(identification_values, identification_prediction),
            "validation": coefficient_metrics(validation_values, validation_prediction),
        }
        serializable[model_name] = _serializable_model(model) | {
            "metrics": {key: _finite_json_metrics(value) for key, value in split_metrics.items()}
        }
        for split, observed, predicted in (
            ("identification", identification_values, identification_prediction),
            ("validation", validation_values, validation_prediction),
        ):
            aggregate = coefficient_metrics(observed, predicted)
            metric_rows.append(
                {
                    "record": record_name,
                    "role": RECORD_ROLES[record_name],
                    "model": model_name,
                    "split": split,
                    "harmonic": "all",
                    **_finite_json_metrics(aggregate),
                }
            )
            for harmonic_index, harmonic in enumerate(harmonics):
                one = coefficient_metrics(
                    observed[:, harmonic_index : harmonic_index + 1],
                    predicted[:, harmonic_index : harmonic_index + 1],
                )
                metric_rows.append(
                    {
                        "record": record_name,
                        "role": RECORD_ROLES[record_name],
                        "model": model_name,
                        "split": split,
                        "harmonic": harmonic,
                        **_finite_json_metrics(one),
                    }
                )

    serializable["rate_shuffle_negative_control"]["permutation_sha256"] = hashlib.sha256(
        np.asarray(permutation, dtype="<i8").tobytes()
    ).hexdigest()
    serializable["rate_shuffle_negative_control"]["permutation_definition"] = (
        f"numpy.default_rng({config.shuffle_seed + record_index}).permutation(n_identification)"
    )
    result = {
        "record": record_name,
        "role": RECORD_ROLES[record_name],
        "harmonics": list(harmonics),
        "identification_window_count": len(identification_rows),
        "validation_window_count": len(validation_rows),
        "identification_rate_hz": {
            "minimum": float(np.min(identification_rates)),
            "median": float(np.median(identification_rates)),
            "maximum": float(np.max(identification_rates)),
        },
        "validation_rate_hz": {
            "minimum": float(np.min(validation_rates)),
            "median": float(np.median(validation_rates)),
            "maximum": float(np.max(validation_rates)),
        },
        "models": serializable,
        "r2_certification": {
            "status": "NOT_CERTIFIABLE",
            "reason": (
                "Point estimates and empirical held-out residuals do not instantiate "
                "the required interval-valued R2 primitive bounds."
            ),
            "missing_primitives": [
                "annotation_timing_error_bound",
                "hrv_and_whole_window_phase_warp_bound",
                "measurement_noise_envelope",
                "model_discrepancy_envelope",
                "finite_dwell_transient_bound",
                "prehistory_or_memory_tail_bound",
                "sampling_and_quadrature_error_bound",
                "subsampling_anti_alias_error_bound",
                "omitted_harmonic_tail_and_window_leakage_bound",
                "adc_gain_and_analog_front_end_calibration_bound",
                "G_and_b_accumulation_rounding_enclosure",
                "end_to_end_WLS_interval_enclosure",
            ],
        },
    }
    return result, metric_rows


def pilot_decision_payload(
    *,
    record_count: int,
    locked_record_count: int,
    alpha_boundary_hits: int,
    materially_better_than_both: int,
    locked_better_than_constant: int,
    locked_better_than_shuffle: int,
) -> dict[str, str]:
    """Derive empirical labels from counts instead of hard-coding an outcome."""

    counts = (
        record_count,
        locked_record_count,
        alpha_boundary_hits,
        materially_better_than_both,
        locked_better_than_constant,
        locked_better_than_shuffle,
    )
    if any(type(value) is not int or value < 0 for value in counts):
        raise ValueError("pilot decision counts must be nonnegative exact integers")
    if (
        record_count == 0
        or locked_record_count > record_count
        or alpha_boundary_hits > record_count
        or materially_better_than_both > record_count
        or locked_better_than_constant > locked_record_count
        or locked_better_than_shuffle > locked_record_count
    ):
        raise ValueError("pilot decision counts are internally inconsistent")

    distinct_not_supported = (
        alpha_boundary_hits == record_count and materially_better_than_both == 0
    )
    confirmatory_not_supported = (
        locked_better_than_constant < locked_record_count
        or locked_better_than_shuffle < locked_record_count
    )
    return {
        "pipeline_feasibility": "PASS",
        "distinct_fractional_order_evidence": (
            "NOT_SUPPORTED" if distinct_not_supported else "INCONCLUSIVE"
        ),
        "confirmatory_rate_response_evidence": (
            "NOT_SUPPORTED" if confirmatory_not_supported else "INCONCLUSIVE"
        ),
        "reason": (
            f"{alpha_boundary_hits}/{record_count} profiled fractional fits reached the "
            f"alpha=1 boundary; {materially_better_than_both}/{record_count} materially "
            "improved on both the nested alpha=1 and constant baselines; among "
            f"{locked_record_count} locked records, {locked_better_than_constant} beat the "
            f"constant baseline and {locked_better_than_shuffle} beat the rate-shuffle "
            "control."
        ),
        "recommended_use": (
            "Use this bundle to validate segmentation and leakage controls only; use a "
            "controlled multi-rate protocol with certified R2 primitives for theory testing."
        ),
    }


def run_fantasia_pilot(
    data_dir: Path,
    output_dir: Path,
    *,
    config: FantasiaPilotConfig | None = None,
) -> dict[str, Any]:
    """Run the locked four-record pilot and write a small reproducibility bundle."""

    import wfdb

    selected = config or FantasiaPilotConfig()
    selected.validate()
    source_root = data_dir.resolve()
    destination = output_dir.resolve()
    if (
        destination == source_root
        or source_root in destination.parents
        or destination in source_root.parents
    ):
        raise ValueError("output and read-only ECG_DATA_DIR must be disjoint")
    destination.mkdir(parents=True, exist_ok=True)
    config_hash = selected.lock_hash()
    implementation_manifest = fantasia_implementation_manifest()
    runtime_environment = fantasia_runtime_environment()

    design_payload = selected.canonical_payload() | {
        "config_sha256": config_hash,
        "implementation_manifest": implementation_manifest,
        "runtime_environment": runtime_environment,
        "evidence_scope": (
            "Subject-specific temporal held-out feasibility; not a certified R2 result, "
            "not prospective validation, and not a cardiac-digital-twin claim."
        ),
        "coefficient_convention": (
            "ECG(t)=dc+2*Re(sum_h Z_h exp(i*h*phase(t))); therefore a real term "
            "A*cos(h*phase)+B*sin(h*phase) has positive-frequency Z_h=(A-i*B)/2."
        ),
        "fractional_model": (
            "Z_(w,h)=q_h/[lambda+(i*2*pi*h*f_w*tau_star)^alpha] on the principal "
            "branch; tau_star=1 second is frozen, alpha/lambda are dimensionless and "
            "shared within one subject, and q_h is harmonic-specific."
        ),
        "parameter_units": {
            "alpha": "dimensionless",
            "lambda": "dimensionless",
            "q_h": "mV",
            "f_w": "Hz",
            "tau_star": "seconds",
        },
        "real_ecg_reconstruction": ("y=Z_0+2*Re(sum_{m>0} Z_m*exp(i*m*phi))"),
    }
    _write_json(destination / "design_lock.json", design_payload)

    manifest_rows: list[dict[str, Any]] = []
    coefficient_rows: list[dict[str, Any]] = []
    source_inventory: list[dict[str, Any]] = []
    record_hashes: dict[str, str] = {}
    model_results: list[dict[str, Any]] = []
    all_metric_rows: list[dict[str, Any]] = []
    execution_records: list[dict[str, Any]] = []
    execution_replay_matches = 0

    for record_index, record_name in enumerate(selected.records):
        paths = {
            extension: source_root / f"{record_name}.{extension}"
            for extension in ("hea", "dat", "ecg")
        }
        missing = [str(path) for path in paths.values() if not path.is_file()]
        if missing:
            raise FileNotFoundError(f"missing locked Fantasia inputs: {missing}")
        record_file_hashes: dict[str, str] = {}
        for extension, path in paths.items():
            file_sha256 = sha256_file(path)
            if file_sha256 != FANTASIA_V1_SHA256[path.name]:
                raise ValueError(
                    f"{path.name} does not match the official Fantasia v1.0.0 checksum"
                )
            record_file_hashes[extension] = file_sha256
            source_inventory.append(
                {
                    "record": record_name,
                    "extension": extension,
                    "path_relative_to_ECG_DATA_DIR": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": file_sha256,
                }
            )
        record_hash = canonical_json_hash({"record": record_name, "files": record_file_hashes})
        record_hashes[record_name] = record_hash

        header = wfdb.rdheader(str(source_root / record_name))
        try:
            channel_index = tuple(header.sig_name).index("ECG")
        except ValueError as error:
            raise ValueError(f"record {record_name} has no channel named ECG") from error
        record = wfdb.rdrecord(
            str(source_root / record_name), channels=[channel_index], physical=False
        )
        annotation = wfdb.rdann(str(source_root / record_name), "ecg")
        if record.d_signal is None or record.d_signal.shape[1] != 1:
            raise ValueError(f"record {record_name} did not return one digital ECG channel")
        digital = np.asarray(record.d_signal[:, 0], dtype=np.int64)
        gain = float(record.adc_gain[0])
        baseline = int(record.baseline[0])
        if not math.isfinite(gain) or gain <= 0.0:
            raise ValueError(f"record {record_name} has invalid ECG gain")
        physical = (digital.astype(np.float64) - baseline) / gain
        sampling_frequency = float(record.fs)
        duration_seconds = digital.size / sampling_frequency

        split_windows = {
            "identification": build_rr_windows(
                annotation.sample,
                annotation.symbol,
                sampling_frequency=sampling_frequency,
                segment_start_seconds=0.0,
                segment_end_seconds=selected.identification_end_seconds,
                split="identification",
                rr_intervals_per_window=selected.rr_intervals_per_window,
                beat_symbols=selected.beat_symbols,
                ignored_nonbeat_symbols=selected.ignored_nonbeat_symbols,
            ),
            "validation": build_rr_windows(
                annotation.sample,
                annotation.symbol,
                sampling_frequency=sampling_frequency,
                segment_start_seconds=selected.guard_end_seconds,
                segment_end_seconds=duration_seconds,
                split="validation",
                rr_intervals_per_window=selected.rr_intervals_per_window,
                beat_symbols=selected.beat_symbols,
                ignored_nonbeat_symbols=selected.ignored_nonbeat_symbols,
            ),
        }
        if any(len(windows) < 3 for windows in split_windows.values()):
            raise ValueError(f"record {record_name} has too few complete windows")

        fit_rows: dict[str, list[dict[str, Any]]] = {
            "identification": [],
            "validation": [],
        }
        for split, windows in split_windows.items():
            previous_end = None
            for window in windows:
                if previous_end is not None and window.start_sample < previous_end:
                    raise AssertionError("RR window sample ranges overlap")
                previous_end = window.end_sample
                estimate = estimate_window_coefficients(
                    physical, window.beat_samples, selected.harmonics
                )
                duration = (window.end_sample - window.start_sample) / sampling_frequency
                base_rate = selected.rr_intervals_per_window / duration
                window_digest = _window_hash(record_name, window, digital, config_hash)
                selected_indices = bounded_sample_indices(window.start_sample, window.end_sample)
                digital_window = digital[selected_indices]
                execution = build_wls_execution_record(
                    digital_window,
                    window.beat_samples,
                    sampling_frequency=sampling_frequency,
                    adc_gain=gain,
                    baseline=baseline,
                    harmonics=selected.harmonics,
                    protocol_sha256=config_hash,
                    record_sha256=record_hash,
                    window_sha256=window_digest,
                )
                execution_z = execution["result"]["z_tilde"]
                if any(
                    coefficient != complex(float(stored["real"]), float(stored["imag"]))
                    for coefficient, stored in zip(estimate.coefficients, execution_z, strict=True)
                ):
                    raise AssertionError("reported coefficients do not match WLS execution")
                replay = replay_wls_execution_record(execution, digital_window)
                if replay["status"] != "MATCH":
                    raise AssertionError(
                        f"WLS replay failed for {record_name}/{split}/{window.index}: {replay}"
                    )
                execution_replay_matches += 1
                execution_records.append(
                    {
                        "record": record_name,
                        "role": RECORD_ROLES[record_name],
                        "split": split,
                        "window_index": window.index,
                        "execution": execution,
                    }
                )
                manifest_rows.append(
                    {
                        "record": record_name,
                        "role": RECORD_ROLES[record_name],
                        "split": split,
                        "window_index": window.index,
                        "start_sample": window.start_sample,
                        "end_sample_exclusive": window.end_sample,
                        "start_seconds": window.start_sample / sampling_frequency,
                        "end_seconds": window.end_sample / sampling_frequency,
                        "duration_seconds": duration,
                        "rr_intervals": selected.rr_intervals_per_window,
                        "sample_count": estimate.sample_count,
                        "full_span_sample_count": window.end_sample - window.start_sample,
                        "base_rate_hz": base_rate,
                        "beat_symbols": "".join(window.beat_symbols),
                        "window_sha256": window_digest,
                        "wls_execution_sha256": execution["execution_sha256"],
                    }
                )
                fit_rows[split].append(
                    {
                        "base_rate_hz": base_rate,
                        "coefficients": estimate.coefficients,
                    }
                )
                for harmonic, coefficient in zip(
                    selected.harmonics, estimate.coefficients, strict=True
                ):
                    coefficient_rows.append(
                        {
                            "record": record_name,
                            "role": RECORD_ROLES[record_name],
                            "split": split,
                            "window_index": window.index,
                            "harmonic": harmonic,
                            "base_rate_hz": base_rate,
                            "angular_frequency_rad_s": 2.0 * math.pi * harmonic * base_rate,
                            "dimensionless_frequency": 2.0
                            * math.pi
                            * harmonic
                            * base_rate
                            * selected.tau_star_seconds,
                            "coefficient_real_mV": coefficient.real,
                            "coefficient_imag_mV": coefficient.imag,
                            "dc_mV": estimate.dc,
                            "wls_residual_rmse_mV": estimate.residual_rmse,
                            "gram_min_eigenvalue": estimate.gram_min_eigenvalue,
                            "gram_max_eigenvalue": estimate.gram_max_eigenvalue,
                            "gram_condition_number": estimate.gram_condition_number,
                            "gram_rank": estimate.gram_rank,
                            "design_columns": 1 + 2 * len(selected.harmonics),
                            "window_sha256": window_digest,
                            "wls_execution_sha256": execution["execution_sha256"],
                        }
                    )

        record_result, metric_rows = _evaluate_record_models(
            record_name,
            selected.harmonics,
            fit_rows["identification"],
            fit_rows["validation"],
            selected,
            record_index,
        )
        record_result["source_channel"] = {
            "name": "ECG",
            "original_channel_index": channel_index,
            "units": str(record.units[0]),
            "adc_gain": gain,
            "baseline": baseline,
            "sampling_frequency_hz": sampling_frequency,
            "signal_samples": int(digital.size),
            "duration_seconds": duration_seconds,
        }
        record_result["window_gram_diagnostics"] = {
            "minimum_eigenvalue_over_all_windows": min(
                row["gram_min_eigenvalue"]
                for row in coefficient_rows
                if row["record"] == record_name
            ),
            "maximum_condition_number_over_all_windows": max(
                row["gram_condition_number"]
                for row in coefficient_rows
                if row["record"] == record_name
            ),
            "all_full_column_rank": all(
                row["gram_rank"] == row["design_columns"]
                for row in coefficient_rows
                if row["record"] == record_name
            ),
        }
        model_results.append(record_result)
        all_metric_rows.extend(metric_rows)

    manifest_fields = (
        "record",
        "role",
        "split",
        "window_index",
        "start_sample",
        "end_sample_exclusive",
        "start_seconds",
        "end_seconds",
        "duration_seconds",
        "rr_intervals",
        "sample_count",
        "full_span_sample_count",
        "base_rate_hz",
        "beat_symbols",
        "window_sha256",
        "wls_execution_sha256",
    )
    coefficient_fields = (
        "record",
        "role",
        "split",
        "window_index",
        "harmonic",
        "base_rate_hz",
        "angular_frequency_rad_s",
        "dimensionless_frequency",
        "coefficient_real_mV",
        "coefficient_imag_mV",
        "dc_mV",
        "wls_residual_rmse_mV",
        "gram_min_eigenvalue",
        "gram_max_eigenvalue",
        "gram_condition_number",
        "gram_rank",
        "design_columns",
        "window_sha256",
        "wls_execution_sha256",
    )
    metric_fields = (
        "record",
        "role",
        "model",
        "split",
        "harmonic",
        "count_complex",
        "sse_mV2",
        "rmse_mV",
        "mae_mV",
        "nrmse_energy",
        "r2_centered",
    )
    _write_csv_gzip(destination / "window_manifest.csv.gz", manifest_rows, manifest_fields)
    _write_csv_gzip(destination / "wls_coefficients.csv.gz", coefficient_rows, coefficient_fields)
    _write_csv(destination / "heldout_metrics.csv", all_metric_rows, metric_fields)
    _write_jsonl_gzip(destination / "wls_execution_records.jsonl.gz", execution_records)
    _write_json(
        destination / "model_fits.json",
        {
            "model_definition": design_payload["fractional_model"],
            "coefficient_convention": design_payload["coefficient_convention"],
            "records": model_results,
        },
    )
    source_payload = {
        "source_root_policy": "resolved externally through ECG_DATA_DIR",
        "raw_data_copied_into_repository": False,
        "selected_files_only": True,
        "dataset": {
            "name": "Fantasia Database",
            "version": "1.0.0",
            "doi": "10.13026/C2RG61",
            "landing_url": "https://physionet.org/content/fantasia/1.0.0/",
            "upstream_sha256_manifest_url": (
                "https://physionet.org/files/fantasia/1.0.0/SHA256SUMS.txt"
            ),
            "selected_files_match_upstream_sha256_manifest": True,
            "file_license": "Open Data Commons Attribution License v1.0",
            "required_citations": [
                "Iyengar et al., Am J Physiol 271 (1996), R1078-R1084",
                "Goldberger et al., Circulation 101 (2000), e215-e220",
            ],
            "derived_content_notice": (
                "The repository bundle contains no waveform samples, but serialized WLS "
                "execution records retain selected beat-sample indices derived from the "
                "database annotations; the dataset attribution and license boundary remain."
            ),
        },
        "files": source_inventory,
        "record_sha256": record_hashes,
        "unique_source_bytes": sum(row["bytes"] for row in source_inventory),
    }
    _write_json(destination / "source_inventory.json", source_payload)
    r2_payload = {
        "overall_status": "NOT_CERTIFIABLE",
        "record_statuses": {
            result["record"]: result["r2_certification"] for result in model_results
        },
        "prohibition": (
            "Empirical WLS/model-fit intervals must not be relabeled as a certified R2 PASS."
        ),
    }
    _write_json(destination / "r2_status.json", r2_payload)

    material_tolerance = 1.0e-3
    materially_better_than_both = 0
    alpha_boundary_hits = 0
    locked_better_than_constant = 0
    locked_better_than_shuffle = 0
    validation_comparisons: dict[str, Any] = {}
    for result in model_results:
        models = result["models"]
        alpha_boundary_hits += int(models["fractional_common"]["alpha_boundary_hit"])
        fractional_rmse = models["fractional_common"]["metrics"]["validation"]["rmse_mV"]
        order_one_rmse = models["alpha_one_nested"]["metrics"]["validation"]["rmse_mV"]
        constant_rmse = models["constant_fourier"]["metrics"]["validation"]["rmse_mV"]
        shuffled_rmse = models["rate_shuffle_negative_control"]["metrics"]["validation"]["rmse_mV"]
        if fractional_rmse < (1.0 - material_tolerance) * min(order_one_rmse, constant_rmse):
            materially_better_than_both += 1
        if result["role"] == "locked_engineering_validation":
            locked_better_than_constant += int(
                fractional_rmse < (1.0 - material_tolerance) * constant_rmse
            )
            locked_better_than_shuffle += int(
                fractional_rmse < (1.0 - material_tolerance) * shuffled_rmse
            )
        validation_comparisons[result["record"]] = {
            "fractional_rmse_mV": fractional_rmse,
            "alpha_one_rmse_mV": order_one_rmse,
            "constant_rmse_mV": constant_rmse,
            "rate_shuffle_rmse_mV": shuffled_rmse,
            "fractional_vs_constant_ratio": fractional_rmse / constant_rmse,
            "fractional_vs_alpha_one_ratio": fractional_rmse / order_one_rmse,
            "fractional_vs_rate_shuffle_ratio": fractional_rmse / shuffled_rmse,
        }
    decision = pilot_decision_payload(
        record_count=len(model_results),
        locked_record_count=2,
        alpha_boundary_hits=alpha_boundary_hits,
        materially_better_than_both=materially_better_than_both,
        locked_better_than_constant=locked_better_than_constant,
        locked_better_than_shuffle=locked_better_than_shuffle,
    )
    rounding_enclosures = [
        row["execution"]["result"]["joint_solve_rounding_enclosure"] for row in execution_records
    ]
    enclosed_rounding = [row for row in rounding_enclosures if row["status"] == "ENCLOSED"]
    summary = {
        "config_sha256": config_hash,
        "implementation_sha256": implementation_manifest["implementation_sha256"],
        "runtime_environment": runtime_environment,
        "record_count": len(model_results),
        "calibration_record_count": 2,
        "locked_engineering_validation_record_count": 2,
        "window_count": len(manifest_rows),
        "coefficient_count": len(coefficient_rows),
        "wls_execution_record_count": len(execution_records),
        "wls_execution_replay_match_count": execution_replay_matches,
        "wls_solve_rounding_enclosed_count": len(enclosed_rounding),
        "wls_solve_rounding_maximum_absolute_error_upper_bound": max(
            (row["maximum_absolute_error_upper_bound"] for row in enclosed_rounding),
            default=None,
        ),
        "wls_execution_rounding_scope": (
            "Arb enclosure covers the joint solve and identity-H division relative to "
            "recorded binary64 G and b; their accumulation and all R2 primitives remain open."
        ),
        "material_improvement_threshold_fraction": material_tolerance,
        "fractional_materially_better_than_alpha_one_and_constant_count": (
            materially_better_than_both
        ),
        "fractional_alpha_boundary_hit_count": alpha_boundary_hits,
        "locked_fractional_better_than_constant_count": locked_better_than_constant,
        "locked_fractional_better_than_rate_shuffle_count": locked_better_than_shuffle,
        "validation_comparisons": validation_comparisons,
        "r2_status": "NOT_CERTIFIABLE",
        "manifest_sha256": sha256_file(destination / "window_manifest.csv.gz"),
        "evidence_scope": design_payload["evidence_scope"],
        "decision": decision,
    }
    _write_json(destination / "summary.json", summary)

    checksum_targets = (
        "design_lock.json",
        "window_manifest.csv.gz",
        "wls_coefficients.csv.gz",
        "wls_execution_records.jsonl.gz",
        "heldout_metrics.csv",
        "model_fits.json",
        "source_inventory.json",
        "r2_status.json",
        "summary.json",
    )
    checksums = {name: sha256_file(destination / name) for name in checksum_targets}
    checksum_lines = [f"{digest}  {name}" for name, digest in checksums.items()]
    _write_text_lf(
        destination / "SHA256SUMS",
        "\n".join(checksum_lines) + "\n",
        encoding="ascii",
    )
    return summary
