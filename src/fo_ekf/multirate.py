"""Data-free formulas for the multi-heart-rate fractional-order certificate.

This module implements the algebra frozen in
``research/control_theory_innovation_spec_2026-07-11.md``.  It is deliberately
limited to the ideal steady harmonic model; it is not an ECG estimator and it
does not hide finite-dwell or finite-memory errors.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class MultiRateRecovery:
    """Parameters recovered from three ideal steady harmonic responses."""

    order: float
    damping: float
    morphology: complex
    invariant: complex


def _three_frequencies(frequencies: Sequence[float]) -> tuple[float, float, float]:
    values = tuple(float(value) for value in frequencies)
    if len(values) != 3:
        raise ValueError("exactly three frequencies are required")
    if not (0.0 < values[0] < values[1] < values[2]):
        raise ValueError("frequencies must satisfy 0 < nu_1 < nu_2 < nu_3")
    return values


def harmonic_response(
    order: float,
    damping: float,
    morphology: complex,
    frequency: float,
) -> complex:
    """Return q / (lambda + (i*nu)**alpha) on the principal branch."""

    if order <= 0.0:
        raise ValueError("order must be positive")
    if damping <= 0.0:
        raise ValueError("damping must be positive")
    if frequency <= 0.0:
        raise ValueError("frequency must be positive")
    if morphology == 0.0:
        raise ValueError("morphology coefficient must be nonzero")
    fractional_frequency = frequency**order * cmath.exp(0.5j * math.pi * order)
    return morphology / (damping + fractional_frequency)


def order_invariant(order: float, frequencies: Sequence[float]) -> float:
    """Return F(alpha) for an ordered frequency triple."""

    if order <= 0.0:
        raise ValueError("order must be positive")
    nu_1, nu_2, nu_3 = _three_frequencies(frequencies)
    log_middle = math.log(nu_2 / nu_1)
    log_upper = math.log(nu_3 / nu_1)
    return math.expm1(log_upper * order) / math.expm1(log_middle * order)


def order_invariant_slope(order: float, frequencies: Sequence[float]) -> float:
    """Return the strictly positive derivative F'(alpha)."""

    if order <= 0.0:
        raise ValueError("order must be positive")
    nu_1, nu_2, nu_3 = _three_frequencies(frequencies)
    log_middle = math.log(nu_2 / nu_1)
    log_upper = math.log(nu_3 / nu_1)
    value = order_invariant(order, frequencies)
    log_slope = log_upper / (-math.expm1(-log_upper * order)) - log_middle / (
        -math.expm1(-log_middle * order)
    )
    return value * log_slope


def response_invariant(responses: Sequence[complex]) -> complex:
    """Return (W_3-W_1)/(W_2-W_1), where W_r=1/Z_r."""

    values = tuple(complex(value) for value in responses)
    if len(values) != 3:
        raise ValueError("exactly three responses are required")
    if any(value == 0.0 for value in values):
        raise ValueError("responses must be nonzero")
    inverse = tuple(1.0 / value for value in values)
    denominator = inverse[1] - inverse[0]
    if denominator == 0.0:
        raise ValueError("inverse responses do not provide separated excitation")
    return (inverse[2] - inverse[0]) / denominator


def invert_order_invariant(
    invariant: float,
    frequencies: Sequence[float],
    *,
    lower: float = 1.0e-6,
    upper: float = 1.0,
    tolerance: float = 1.0e-13,
    max_iterations: int = 200,
) -> float:
    """Invert the strictly increasing F(alpha) by a bracketed bisection."""

    _three_frequencies(frequencies)
    if not 0.0 < lower < upper:
        raise ValueError("order bounds must satisfy 0 < lower < upper")
    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive")

    target = float(invariant)
    value_lower = order_invariant(lower, frequencies)
    value_upper = order_invariant(upper, frequencies)
    slack = 10.0 * tolerance * max(1.0, abs(target))
    if target < value_lower - slack or target > value_upper + slack:
        raise ValueError("invariant lies outside the selected order interval")

    if abs(target - value_lower) <= slack:
        return lower
    if abs(target - value_upper) <= slack:
        return upper

    left, right = lower, upper
    for _ in range(max_iterations):
        middle = 0.5 * (left + right)
        value_middle = order_invariant(middle, frequencies)
        if abs(value_middle - target) <= tolerance or right - left <= tolerance:
            return middle
        if value_middle < target:
            left = middle
        else:
            right = middle
    raise RuntimeError("order-invariant bisection did not converge")


def recover_parameters(
    responses: Sequence[complex],
    frequencies: Sequence[float],
    *,
    order_bounds: tuple[float, float] = (1.0e-6, 1.0),
    imaginary_tolerance: float = 1.0e-9,
) -> MultiRateRecovery:
    """Recover alpha, lambda, and the effective complex morphology coefficient.

    A common lead gain or common phase rotation is absorbed into ``morphology``;
    the recovered order and damping remain unchanged.
    """

    nu_1, nu_2, _ = _three_frequencies(frequencies)
    values = tuple(complex(value) for value in responses)
    invariant = response_invariant(values)
    if abs(invariant.imag) > imaginary_tolerance * max(1.0, abs(invariant.real)):
        raise ValueError("response invariant has a non-negligible imaginary component")

    order = invert_order_invariant(
        invariant.real,
        frequencies,
        lower=order_bounds[0],
        upper=order_bounds[1],
    )
    inverse = tuple(1.0 / value for value in values)
    scale = (inverse[1] - inverse[0]) / (nu_2**order - nu_1**order)
    morphology = cmath.exp(0.5j * math.pi * order) / scale
    intercept = inverse[0] - scale * nu_1**order
    damping_complex = intercept * morphology
    if abs(damping_complex.imag) > imaginary_tolerance * max(1.0, abs(damping_complex.real)):
        raise ValueError("recovered damping is not real within tolerance")
    if damping_complex.real <= 0.0:
        raise ValueError("recovered damping is not positive")
    return MultiRateRecovery(
        order=order,
        damping=damping_complex.real,
        morphology=morphology,
        invariant=invariant,
    )


def inverse_response_error_bound(response_magnitude: float, response_error: float) -> float:
    """Bound |1/Z_hat - 1/Z| from |Z_hat-Z| and |Z|."""

    magnitude = float(response_magnitude)
    error = float(response_error)
    if magnitude <= 0.0 or error < 0.0:
        raise ValueError("response magnitude must be positive and error nonnegative")
    if error >= magnitude:
        raise ValueError("response error must be smaller than the true magnitude")
    return error / (magnitude * (magnitude - error))


def invariant_error_bound(
    inverse_responses: Sequence[complex],
    inverse_response_error: float,
) -> float:
    """Return the rigorous quotient perturbation bound from Theorem 3."""

    values = tuple(complex(value) for value in inverse_responses)
    if len(values) != 3:
        raise ValueError("exactly three inverse responses are required")
    epsilon = float(inverse_response_error)
    if epsilon < 0.0:
        raise ValueError("inverse-response error must be nonnegative")
    denominator = values[1] - values[0]
    numerator = values[2] - values[0]
    if 2.0 * epsilon >= abs(denominator):
        raise ValueError("error is too large for a separated-excitation certificate")
    return (
        2.0
        * epsilon
        * (abs(denominator) + abs(numerator))
        / (abs(denominator) * (abs(denominator) - 2.0 * epsilon))
    )


def gl_tail_mass(order: float, memory_length: int) -> float:
    """Return the exact GL coefficient mass after ``memory_length`` terms."""

    if not 0.0 < order < 1.0:
        raise ValueError("GL tail formula requires 0 < order < 1")
    if memory_length < 0:
        raise ValueError("memory length must be nonnegative")
    log_mass = (
        math.lgamma(memory_length + 1.0 - order)
        - math.lgamma(1.0 - order)
        - math.lgamma(memory_length + 1.0)
    )
    return math.exp(log_mass)
