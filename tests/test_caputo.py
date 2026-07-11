import math

import numpy as np
import pytest
from scipy.special import erfcx

from fo_ekf.caputo import (
    make_harmonic_schedule,
    mittag_leffler_relaxation,
    piecewise_harmonic_caputo_oracle,
    solve_linear_caputo_pece,
    solve_linear_caputo_product_integration,
)


def test_mittag_leffler_relaxation_special_cases() -> None:
    for time in (0.1, 1.0, 5.0):
        assert mittag_leffler_relaxation(1.0, 0.7, time) == pytest.approx(
            math.exp(-0.7 * time), rel=1.0e-13
        )
        assert mittag_leffler_relaxation(0.5, 0.7, time) == pytest.approx(
            erfcx(0.7 * math.sqrt(time)), rel=3.0e-10, abs=3.0e-11
        )


def test_full_history_solvers_converge_to_homogeneous_reference() -> None:
    order = 0.5
    damping = 0.7
    final_time = 4.0
    errors = []
    for step in (0.04, 0.02, 0.01):
        time = np.arange(round(final_time / step) + 1) * step
        forcing = np.zeros_like(time, dtype=complex)
        state = solve_linear_caputo_product_integration(order, damping, time, forcing, 1.0)
        errors.append(abs(state[-1].real - mittag_leffler_relaxation(order, damping, final_time)))
    assert errors[2] < errors[1] < errors[0]
    assert errors[-1] < 2.0e-5


def test_time_stepper_and_positive_spectrum_oracle_agree_after_switch() -> None:
    order = 0.7
    damping = 0.8
    frequencies = (0.9, 1.7)
    dwell = 4.0
    step = 0.01
    morphology = 1.1 - 0.3j
    initial = 0.2 + 0.1j
    schedule = make_harmonic_schedule(frequencies, dwell, step, morphology)
    implicit = solve_linear_caputo_product_integration(
        order, damping, schedule.time, schedule.forcing, initial
    )
    pece = solve_linear_caputo_pece(order, damping, schedule.time, schedule.forcing, initial)

    for checkpoint in (2.0, 4.0, 6.0, 8.0):
        index = round(checkpoint / step)
        oracle = piecewise_harmonic_caputo_oracle(
            order,
            damping,
            frequencies,
            dwell,
            morphology,
            initial,
            checkpoint,
        )
        assert implicit[index] == pytest.approx(oracle, abs=3.0e-4)
        assert pece[index] == pytest.approx(oracle, abs=4.0e-4)


def test_schedule_keeps_phase_continuous_at_switch() -> None:
    schedule = make_harmonic_schedule((0.7, 1.9, 1.1), 2.0, 0.01, 1.0)
    for switch in (2.0, 4.0):
        index = round(switch / 0.01)
        left_extrapolation = schedule.phase[index - 1] + (
            schedule.phase[index - 1] - schedule.phase[index - 2]
        )
        assert schedule.phase[index] == pytest.approx(left_extrapolation, abs=1.0e-12)
