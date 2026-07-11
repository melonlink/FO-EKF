import cmath
import math
from itertools import pairwise

import pytest

from fo_ekf.multirate import (
    gl_tail_mass,
    harmonic_response,
    invariant_error_bound,
    inverse_response_error_bound,
    order_invariant,
    order_invariant_slope,
    recover_parameters,
    response_invariant,
)


def test_order_invariant_is_strictly_increasing() -> None:
    frequencies = (0.8, 1.3, 2.1)
    orders = [0.05 + 0.01 * index for index in range(96)]
    values = [order_invariant(order, frequencies) for order in orders]

    assert all(right > left for left, right in pairwise(values))
    assert all(order_invariant_slope(order, frequencies) > 0.0 for order in orders)


def test_three_rates_recover_order_damping_and_effective_morphology() -> None:
    order = 0.73
    damping = 0.42
    morphology = 1.7 - 0.6j
    frequencies = (0.7, 1.2, 2.0)
    responses = [
        harmonic_response(order, damping, morphology, frequency) for frequency in frequencies
    ]

    recovered = recover_parameters(responses, frequencies)

    assert recovered.order == pytest.approx(order, abs=2.0e-12)
    assert recovered.damping == pytest.approx(damping, abs=2.0e-12)
    assert recovered.morphology == pytest.approx(morphology, abs=2.0e-12)
    assert recovered.invariant.imag == pytest.approx(0.0, abs=1.0e-12)


def test_common_gain_and_phase_gauge_do_not_change_order_or_damping() -> None:
    order = 0.97
    damping = 0.8
    morphology = 0.9 + 0.4j
    frequencies = (0.9, 1.5, 2.4)
    common_measurement_factor = 3.2 * cmath.exp(0.63j)
    responses = [
        common_measurement_factor * harmonic_response(order, damping, morphology, frequency)
        for frequency in frequencies
    ]

    recovered = recover_parameters(responses, frequencies)

    assert recovered.order == pytest.approx(order, abs=2.0e-11)
    assert recovered.damping == pytest.approx(damping, abs=2.0e-11)
    assert recovered.morphology == pytest.approx(
        common_measurement_factor * morphology,
        abs=2.0e-11,
    )


def test_independent_phase_gauges_are_rejected() -> None:
    order = 0.8
    frequencies = (0.8, 1.4, 2.3)
    responses = [
        harmonic_response(order, 0.6, 1.0 + 0.2j, frequency) * cmath.exp(1j * phase)
        for frequency, phase in zip(frequencies, (0.0, 0.25, -0.31), strict=True)
    ]

    assert abs(response_invariant(responses).imag) > 0.01
    with pytest.raises(ValueError, match="imaginary"):
        recover_parameters(responses, frequencies)


def test_quotient_error_bound_contains_actual_perturbation() -> None:
    order = 0.68
    frequencies = (0.75, 1.25, 2.2)
    responses = tuple(
        harmonic_response(order, 0.55, 1.2 - 0.7j, frequency) for frequency in frequencies
    )
    perturbations = (1.0e-4 + 2.0e-4j, -2.0e-4 + 1.0e-4j, 1.5e-4 - 1.0e-4j)
    perturbed = tuple(
        response + perturbation
        for response, perturbation in zip(responses, perturbations, strict=True)
    )

    response_error = max(abs(value) for value in perturbations)
    epsilon_w = max(
        inverse_response_error_bound(abs(response), response_error) for response in responses
    )
    inverse = tuple(1.0 / response for response in responses)
    bound = invariant_error_bound(inverse, epsilon_w)
    actual = abs(response_invariant(perturbed) - response_invariant(responses))

    assert actual <= bound


def test_gl_tail_mass_matches_coefficient_recurrence() -> None:
    order = 0.37
    memory_length = 250
    coefficient = order
    partial_sum = coefficient
    for index in range(2, memory_length + 1):
        coefficient *= (index - 1.0 - order) / index
        partial_sum += coefficient

    assert gl_tail_mass(order, memory_length) == pytest.approx(
        1.0 - partial_sum,
        rel=2.0e-12,
        abs=2.0e-14,
    )
    assert gl_tail_mass(order, 2 * memory_length) < gl_tail_mass(order, memory_length)
    asymptotic_scaled = gl_tail_mass(order, 10_000) * 10_000**order * math.gamma(1.0 - order)
    assert asymptotic_scaled == pytest.approx(1.0, rel=5.0e-5)
