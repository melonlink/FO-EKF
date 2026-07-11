import cmath

import pytest

from fo_ekf.certificate import (
    disk_order_intervals,
    evaluate_chain_certificate,
    slope_lower_bound,
)
from fo_ekf.multirate import harmonic_response, order_invariant_slope


def test_chain_and_disk_certificates_cover_true_order() -> None:
    order = 0.78
    frequencies = (0.8, 1.4, 2.4)
    truth = tuple(
        harmonic_response(order, 0.65, 1.2 - 0.4j, frequency) for frequency in frequencies
    )
    radius = 2.0e-4
    perturbations = (
        radius * cmath.exp(0.2j),
        radius * cmath.exp(2.1j),
        radius * cmath.exp(-1.4j),
    )
    measured = tuple(value + error for value, error in zip(truth, perturbations, strict=True))

    chain = evaluate_chain_certificate(measured, frequencies, radius)
    disk = disk_order_intervals(measured, frequencies, radius)

    assert chain.status == "PASS"
    assert chain.order_interval is not None
    assert chain.order_interval.contains(order)
    assert any(interval.contains(order) for interval in disk)
    assert sum(interval.width for interval in disk) <= chain.order_interval.width


def test_protocol_failure_is_not_certifiable_even_if_numeric_model_fits() -> None:
    frequencies = (0.8, 1.4, 2.4)
    measured = tuple(
        harmonic_response(0.82, 0.6, 1.0 + 0.1j, frequency) for frequency in frequencies
    )

    result = evaluate_chain_certificate(
        measured,
        frequencies,
        1.0e-8,
        protocol_certified=False,
    )

    assert result.numeric_pass
    assert result.status == "NOT_CERTIFIABLE"
    assert "protocol_G0" in result.reasons


def test_independent_phase_error_is_numerically_rejected() -> None:
    frequencies = (0.8, 1.4, 2.4)
    measured = tuple(
        harmonic_response(0.82, 0.6, 1.0 + 0.1j, frequency) * cmath.exp(1j * phase)
        for frequency, phase in zip(frequencies, (0.0, 0.2, -0.3), strict=True)
    )

    result = evaluate_chain_certificate(measured, frequencies, 1.0e-5)

    assert result.status == "REJECT_NUMERIC"
    assert "imaginary_inconsistency" in result.reasons


def test_close_frequencies_fail_the_precision_gate() -> None:
    frequencies = (1.0, 1.001, 1.002)
    measured = tuple(
        harmonic_response(0.9, 0.7, 1.0 - 0.2j, frequency) for frequency in frequencies
    )

    result = evaluate_chain_certificate(measured, frequencies, 1.0e-4)

    assert result.status == "REJECT_NUMERIC"
    assert (
        "unseparated_inverse_response" in result.reasons or "order_bound_too_wide" in result.reasons
    )


def test_strict_convexity_places_slope_minimum_at_lower_order() -> None:
    frequencies = (0.65, 1.3, 2.7)
    lower = 0.2
    minimum = slope_lower_bound((lower, 1.0), frequencies)

    assert minimum == pytest.approx(order_invariant_slope(lower, frequencies))
    assert all(
        order_invariant_slope(lower + index * 0.01, frequencies) >= minimum for index in range(81)
    )
