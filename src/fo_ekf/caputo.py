"""Independent full-history Caputo truth utilities for S1 validation.

The solver is a Diethelm-Ford-Freed predictor-corrector product-integration
scheme.  It retains every past right-hand-side evaluation and never restarts at
a heart-rate switch.  It is intentionally separate from any future observer or
finite-memory implementation.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import quad, quad_vec

ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class HarmonicSchedule:
    """Uniform-grid piecewise-constant frequency schedule with continuous phase."""

    time: FloatArray
    phase: FloatArray
    forcing: ComplexArray
    segment: NDArray[np.int64]
    frequencies: tuple[float, ...]
    dwell: float


def solve_linear_caputo_pece(
    order: float,
    damping: float,
    time: Sequence[float] | FloatArray,
    forcing: Sequence[complex] | ComplexArray,
    initial: complex,
) -> ComplexArray:
    """Solve ``D^alpha z = -damping*z + forcing`` with full history.

    The grid must be uniform and start at zero.  For ``order == 1`` the same
    formula reduces to a second-order predictor-corrector ODE scheme.
    """

    if not 0.0 < order <= 1.0:
        raise ValueError("order must satisfy 0 < order <= 1")
    if damping <= 0.0:
        raise ValueError("damping must be positive")

    grid = np.asarray(time, dtype=float)
    input_values = np.asarray(forcing, dtype=complex)
    if grid.ndim != 1 or len(grid) < 2:
        raise ValueError("time must be a one-dimensional grid with at least two points")
    if input_values.shape != grid.shape:
        raise ValueError("forcing must have the same shape as time")
    if not math.isclose(grid[0], 0.0, abs_tol=1.0e-14):
        raise ValueError("time grid must start at zero")
    steps = np.diff(grid)
    step = float(steps[0])
    if step <= 0.0 or not np.allclose(steps, step, rtol=1.0e-11, atol=1.0e-13):
        raise ValueError("time grid must be strictly increasing and uniform")

    count = len(grid) - 1
    state = np.empty(count + 1, dtype=complex)
    rhs = np.empty(count + 1, dtype=complex)
    state[0] = complex(initial)
    rhs[0] = -damping * state[0] + input_values[0]

    lag = np.arange(count + 1, dtype=float)
    predictor_weights = (lag + 1.0) ** order - lag**order
    corrector_weights = (
        (lag + 2.0) ** (order + 1.0) + lag ** (order + 1.0) - 2.0 * (lag + 1.0) ** (order + 1.0)
    )
    predictor_factor = step**order / math.gamma(order + 1.0)
    corrector_factor = step**order / math.gamma(order + 2.0)

    for index in range(count):
        predicted = initial + predictor_factor * np.dot(
            predictor_weights[: index + 1], rhs[index::-1]
        )
        first_weight = index ** (order + 1.0) - (index - order) * (index + 1.0) ** order
        history = first_weight * rhs[0]
        if index >= 1:
            history += np.dot(corrector_weights[:index], rhs[index:0:-1])
        predicted_rhs = -damping * predicted + input_values[index + 1]
        state[index + 1] = initial + corrector_factor * (predicted_rhs + history)
        rhs[index + 1] = -damping * state[index + 1] + input_values[index + 1]

    return state


def solve_linear_caputo_product_integration(
    order: float,
    damping: float,
    time: Sequence[float] | FloatArray,
    forcing: Sequence[complex] | ComplexArray,
    initial: complex,
) -> ComplexArray:
    """Solve the linear equation with a full-history implicit product integral."""

    if not 0.0 < order <= 1.0:
        raise ValueError("order must satisfy 0 < order <= 1")
    if damping <= 0.0:
        raise ValueError("damping must be positive")
    grid = np.asarray(time, dtype=float)
    input_values = np.asarray(forcing, dtype=complex)
    if grid.ndim != 1 or len(grid) < 2 or input_values.shape != grid.shape:
        raise ValueError("time and forcing must share a one-dimensional nontrivial shape")
    if not math.isclose(grid[0], 0.0, abs_tol=1.0e-14):
        raise ValueError("time grid must start at zero")
    steps = np.diff(grid)
    step = float(steps[0])
    if step <= 0.0 or not np.allclose(steps, step, rtol=1.0e-11, atol=1.0e-13):
        raise ValueError("time grid must be strictly increasing and uniform")

    count = len(grid) - 1
    state = np.empty(count + 1, dtype=complex)
    rhs = np.empty(count + 1, dtype=complex)
    state[0] = complex(initial)
    rhs[0] = -damping * state[0] + input_values[0]
    lag = np.arange(count + 1, dtype=float)
    corrector_weights = (
        (lag + 2.0) ** (order + 1.0) + lag ** (order + 1.0) - 2.0 * (lag + 1.0) ** (order + 1.0)
    )
    factor = step**order / math.gamma(order + 2.0)

    for index in range(count):
        first_weight = index ** (order + 1.0) - (index - order) * (index + 1.0) ** order
        history = first_weight * rhs[0]
        if index >= 1:
            history += np.dot(corrector_weights[:index], rhs[index:0:-1])
        state[index + 1] = (initial + factor * (history + input_values[index + 1])) / (
            1.0 + factor * damping
        )
        rhs[index + 1] = -damping * state[index + 1] + input_values[index + 1]

    return state


def make_harmonic_schedule(
    frequencies: Sequence[float],
    dwell: float,
    step: float,
    morphology: complex,
    *,
    initial_phase: float = 0.0,
) -> HarmonicSchedule:
    """Build a continuous-phase piecewise harmonic input on a uniform grid."""

    values = tuple(float(value) for value in frequencies)
    if not values or any(value <= 0.0 for value in values):
        raise ValueError("frequencies must be nonempty and positive")
    if dwell <= 0.0 or step <= 0.0:
        raise ValueError("dwell and step must be positive")
    steps_per_dwell = round(dwell / step)
    if not math.isclose(steps_per_dwell * step, dwell, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("dwell must be an integer multiple of step")

    total_steps = steps_per_dwell * len(values)
    time = np.arange(total_steps + 1, dtype=float) * step
    segment = np.minimum((np.arange(total_steps + 1) // steps_per_dwell), len(values) - 1)
    phase_offsets = np.concatenate(
        ([float(initial_phase)], initial_phase + np.cumsum(np.asarray(values[:-1]) * dwell))
    )
    segment_start = segment.astype(float) * dwell
    phase = phase_offsets[segment] + np.asarray(values)[segment] * (time - segment_start)
    forcing = complex(morphology) * np.exp(1j * phase)
    return HarmonicSchedule(
        time=time,
        phase=phase,
        forcing=forcing.astype(complex),
        segment=segment.astype(np.int64),
        frequencies=values,
        dwell=float(dwell),
    )


def mittag_leffler_relaxation(order: float, damping: float, time: float) -> float:
    """Evaluate ``E_alpha(-damping*time**alpha)`` by a positive spectral integral."""

    if not 0.0 < order <= 1.0:
        raise ValueError("order must satisfy 0 < order <= 1")
    if damping <= 0.0 or time < 0.0:
        raise ValueError("damping must be positive and time nonnegative")
    if time == 0.0:
        return 1.0
    if order == 1.0:
        return math.exp(-damping * time)

    sine = math.sin(math.pi * order)
    cosine = math.cos(math.pi * order)

    def density(rate: float) -> float:
        rate_order = rate**order
        denominator = rate_order**2 + 2.0 * damping * rate_order * cosine + damping**2
        return (
            math.exp(-rate * time) * damping * sine / math.pi * rate ** (order - 1.0) / denominator
        )

    value, _error = quad(
        density,
        0.0,
        math.inf,
        epsabs=1.0e-11,
        epsrel=1.0e-10,
        limit=800,
    )
    return float(value)


def finite_dwell_bound(
    order: float,
    damping: float,
    dwell_time: float,
    initial_magnitude: float,
    previous_input_bound: float,
    new_input_bound: float,
) -> float:
    """Return the theorem's pointwise finite-dwell/history-pollution bound."""

    if min(initial_magnitude, previous_input_bound, new_input_bound) < 0.0:
        raise ValueError("magnitudes and input bounds must be nonnegative")
    coefficient = initial_magnitude + (previous_input_bound + new_input_bound) / damping
    return coefficient * mittag_leffler_relaxation(order, damping, dwell_time)


def piecewise_harmonic_caputo_oracle(
    order: float,
    damping: float,
    frequencies: Sequence[float],
    dwell: float,
    morphology: complex,
    initial: complex,
    time: float,
    *,
    initial_phase: float = 0.0,
) -> complex:
    """Evaluate the full-history solution by positive-spectrum quadrature.

    This is an independent checkpoint oracle for the linear piecewise harmonic
    problem.  It performs no time stepping and includes every completed history
    segment in the convolution.
    """

    values = tuple(float(value) for value in frequencies)
    if not values or any(value <= 0.0 for value in values):
        raise ValueError("frequencies must be nonempty and positive")
    if not 0.0 < order <= 1.0 or damping <= 0.0 or dwell <= 0.0 or time < 0.0:
        raise ValueError("invalid order, damping, dwell, or time")
    if time > len(values) * dwell + 1.0e-12:
        raise ValueError("time lies beyond the supplied frequency schedule")
    if time == 0.0:
        return complex(initial)

    segments: list[tuple[float, float, float, float, float]] = []
    phase_start = float(initial_phase)
    for index, frequency in enumerate(values):
        start = index * dwell
        if start >= time:
            break
        end = min((index + 1) * dwell, time)
        phase_end = phase_start + frequency * (end - start)
        segments.append((start, end, phase_start, phase_end, frequency))
        phase_start += frequency * dwell

    if order == 1.0:
        result = complex(initial) * math.exp(-damping * time)
        for start, end, phase_a, phase_b, frequency in segments:
            result += (
                complex(morphology)
                * (
                    math.exp(-damping * (time - end)) * cmath_exp(phase_b)
                    - math.exp(-damping * (time - start)) * cmath_exp(phase_a)
                )
                / (damping + 1j * frequency)
            )
        return result

    sine = math.sin(math.pi * order)
    cosine = math.cos(math.pi * order)

    def integrand(rate: float) -> complex:
        if rate == 0.0:
            return 0.0j
        rate_order = rate**order
        density = (
            sine
            / math.pi
            * rate_order
            / (rate_order**2 + 2.0 * damping * rate_order * cosine + damping**2)
        )
        convolution = 0.0j
        for start, end, phase_a, phase_b, frequency in segments:
            convolution += (
                complex(morphology)
                * (
                    math.exp(-rate * (time - end)) * cmath_exp(phase_b)
                    - math.exp(-rate * (time - start)) * cmath_exp(phase_a)
                )
                / (rate + 1j * frequency)
            )
        return density * convolution

    convolution, _error = quad_vec(
        integrand,
        0.0,
        math.inf,
        epsabs=2.0e-11,
        epsrel=2.0e-10,
        limit=1000,
    )
    return complex(initial) * mittag_leffler_relaxation(order, damping, time) + complex(convolution)


def cmath_exp(phase: float) -> complex:
    """Return exp(i*phase) without importing cmath in callers."""

    return complex(math.cos(phase), math.sin(phase))


def steady_harmonic_state(
    order: float,
    damping: float,
    morphology: complex,
    frequency: float,
    phase: Sequence[float] | FloatArray,
) -> ComplexArray:
    """Return the infinite-history periodic steady state on supplied phases."""

    from .multirate import harmonic_response

    coefficient = harmonic_response(order, damping, morphology, frequency)
    return coefficient * np.exp(1j * np.asarray(phase, dtype=float))


def demodulate_harmonic(
    time: Sequence[float] | FloatArray,
    state: Sequence[complex] | ComplexArray,
    phase: Sequence[float] | FloatArray,
) -> complex:
    """Estimate a complex harmonic coefficient by time-domain demodulation."""

    grid = np.asarray(time, dtype=float)
    values = np.asarray(state, dtype=complex)
    phases = np.asarray(phase, dtype=float)
    if grid.shape != values.shape or grid.shape != phases.shape or len(grid) < 2:
        raise ValueError("time, state, and phase must share a nontrivial shape")
    duration = float(grid[-1] - grid[0])
    if duration <= 0.0:
        raise ValueError("demodulation duration must be positive")
    baseband = values * np.exp(-1j * phases)
    return complex(np.trapezoid(baseband, grid) / duration)
