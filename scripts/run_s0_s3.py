"""Execute the data-free S0-S3 theory-validation sequence.

The script writes only compact, reproducible tables under ``research/results``.
It does not read patient data and does not download any dependency or artifact.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from fo_ekf.caputo import (
    demodulate_harmonic,
    finite_dwell_bound,
    make_harmonic_schedule,
    mittag_leffler_relaxation,
    piecewise_harmonic_caputo_oracle,
    solve_linear_caputo_pece,
    solve_linear_caputo_product_integration,
    steady_harmonic_state,
)
from fo_ekf.certificate import (
    disk_order_intervals,
    evaluate_chain_certificate,
    inverse_response_disks,
)
from fo_ekf.multirate import (
    harmonic_response,
    order_invariant,
    recover_parameters,
    response_invariant,
)

DEFAULT_OUTPUT = Path("research/results/s0_s3_2026-07-11")
SEED = 20260711


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def quantile(values: list[float], probability: float) -> float:
    if not values:
        return math.nan
    return float(np.quantile(np.asarray(values, dtype=float), probability))


def run_s0(output: Path) -> dict[str, Any]:
    orders = (0.5, 0.7, 0.9, 0.97, 1.0)
    dampings = (0.2, 0.8, 2.0)
    morphologies = (1.0 + 0.2j, 0.6 - 1.1j)
    rates = (0.7, 1.0, 1.4, 1.9, 2.6)
    harmonics = (1, 2)
    measurement_factors = (1.0 + 0.0j, 3.2 * np.exp(0.63j))
    triples = tuple(itertools.combinations(range(len(rates)), 3))
    rows: list[dict[str, Any]] = []

    for order, damping, morphology, harmonic, triple, factor in itertools.product(
        orders,
        dampings,
        morphologies,
        harmonics,
        triples,
        measurement_factors,
    ):
        frequencies = tuple(harmonic * rates[index] for index in triple)
        effective_morphology = factor * morphology * (1.0 + 0.15j * (harmonic - 1))
        responses = tuple(
            harmonic_response(order, damping, effective_morphology, frequency)
            for frequency in frequencies
        )
        recovered = recover_parameters(responses, frequencies)
        order_error = abs(recovered.order - order)
        damping_error = abs(recovered.damping - damping)
        morphology_error = abs(recovered.morphology - effective_morphology)
        rows.append(
            {
                "order": order,
                "damping": damping,
                "harmonic": harmonic,
                "triple": "-".join(str(index) for index in triple),
                "nu1": frequencies[0],
                "nu2": frequencies[1],
                "nu3": frequencies[2],
                "common_factor_abs": abs(factor),
                "common_factor_phase": float(np.angle(factor)),
                "order_error": order_error,
                "damping_error": damping_error,
                "morphology_error": morphology_error,
                "invariant_imag_abs": abs(recovered.invariant.imag),
                "pass": int(
                    order_error <= 5.0e-10
                    and damping_error <= 5.0e-9
                    and morphology_error <= 5.0e-9
                ),
            }
        )

    write_csv(output / "s0_analytic_recovery.csv", rows)
    summary = {
        "cases": len(rows),
        "passed": sum(row["pass"] for row in rows),
        "max_order_error": max(row["order_error"] for row in rows),
        "max_damping_error": max(row["damping_error"] for row in rows),
        "max_morphology_error": max(row["morphology_error"] for row in rows),
        "max_invariant_imag_abs": max(row["invariant_imag_abs"] for row in rows),
    }
    summary["status"] = "PASS" if summary["passed"] == summary["cases"] else "FAIL"
    return summary


def run_s1(output: Path) -> dict[str, Any]:
    orders = (0.5, 0.7, 0.9, 0.97, 1.0)
    damping = 0.7
    morphology = 1.2 - 0.4j
    initial = 0.3 + 0.2j
    frequencies = (0.8, 1.4, 2.3)
    dwell = 20.0
    step = 0.02
    checkpoints = (1.0, 2.0, 5.0, 10.0, 20.0)
    dwell_rows: list[dict[str, Any]] = []
    convergence_rows: list[dict[str, Any]] = []

    for order in orders:
        schedule = make_harmonic_schedule(frequencies, dwell, step, morphology)
        implicit = solve_linear_caputo_product_integration(
            order, damping, schedule.time, schedule.forcing, initial
        )
        pece = solve_linear_caputo_pece(order, damping, schedule.time, schedule.forcing, initial)
        for segment_index, frequency in enumerate(frequencies):
            segment_start = segment_index * dwell
            for time_since_switch in checkpoints:
                absolute_time = segment_start + time_since_switch
                sample_index = round(absolute_time / step)
                oracle = piecewise_harmonic_caputo_oracle(
                    order,
                    damping,
                    frequencies,
                    dwell,
                    morphology,
                    initial,
                    absolute_time,
                )
                steady = steady_harmonic_state(
                    order,
                    damping,
                    morphology,
                    frequency,
                    [schedule.phase[sample_index]],
                )[0]
                previous_bound = 0.0 if segment_index == 0 else abs(morphology)
                bound = finite_dwell_bound(
                    order,
                    damping,
                    time_since_switch,
                    abs(initial),
                    previous_bound,
                    abs(morphology),
                )
                actual_error = abs(oracle - steady)
                dwell_rows.append(
                    {
                        "order": order,
                        "segment": segment_index,
                        "frequency": frequency,
                        "time_since_switch": time_since_switch,
                        "absolute_time": absolute_time,
                        "actual_error": actual_error,
                        "theory_bound": bound,
                        "bound_utilization": actual_error / bound if bound else 0.0,
                        "implicit_oracle_error": abs(implicit[sample_index] - oracle),
                        "pece_oracle_error": abs(pece[sample_index] - oracle),
                        "pass": int(actual_error <= bound + 5.0e-10),
                    }
                )

            window_start = segment_start + dwell / 2.0
            left = round(window_start / step)
            right = round((segment_start + dwell) / step)
            coefficient = demodulate_harmonic(
                schedule.time[left : right + 1],
                implicit[left : right + 1],
                schedule.phase[left : right + 1],
            )
            ideal_coefficient = harmonic_response(order, damping, morphology, frequency)
            coefficient_error = abs(coefficient - ideal_coefficient)
            dwell_budget = finite_dwell_bound(
                order,
                damping,
                dwell / 2.0,
                abs(initial),
                0.0 if segment_index == 0 else abs(morphology),
                abs(morphology),
            )
            dwell_rows.append(
                {
                    "order": order,
                    "segment": segment_index,
                    "frequency": frequency,
                    "time_since_switch": "demodulated_last_half",
                    "absolute_time": segment_start + dwell,
                    "actual_error": coefficient_error,
                    "theory_bound": dwell_budget,
                    "bound_utilization": coefficient_error / dwell_budget,
                    "implicit_oracle_error": math.nan,
                    "pece_oracle_error": math.nan,
                    "pass": int(coefficient_error <= dwell_budget + 5.0e-10),
                }
            )

        convergence_schedule: dict[float, tuple[Any, np.ndarray]] = {}
        convergence_frequencies = (0.9, 1.7)
        convergence_dwell = 4.0
        final_time = len(convergence_frequencies) * convergence_dwell
        oracle_final = piecewise_harmonic_caputo_oracle(
            order,
            damping,
            convergence_frequencies,
            convergence_dwell,
            morphology,
            initial,
            final_time,
        )
        for convergence_step in (0.04, 0.02, 0.01):
            schedule_c = make_harmonic_schedule(
                convergence_frequencies,
                convergence_dwell,
                convergence_step,
                morphology,
            )
            state_c = solve_linear_caputo_product_integration(
                order,
                damping,
                schedule_c.time,
                schedule_c.forcing,
                initial,
            )
            convergence_schedule[convergence_step] = (schedule_c, state_c)
            convergence_rows.append(
                {
                    "order": order,
                    "step": convergence_step,
                    "final_oracle_error": abs(state_c[-1] - oracle_final),
                }
            )

    constant_rows: list[dict[str, Any]] = []
    for order, lambda_value, switch_scale, dwell_scale in itertools.product(
        orders,
        (0.25, 1.0, 4.0),
        (0.5, 2.0, 8.0),
        (0.0625, 0.25, 1.0, 4.0, 16.0),
    ):
        characteristic_time = lambda_value ** (-1.0 / order)
        switch_time = switch_scale * characteristic_time
        dwell_time = dwell_scale * characteristic_time
        before = mittag_leffler_relaxation(order, lambda_value, switch_time + dwell_time)
        after = mittag_leffler_relaxation(order, lambda_value, dwell_time)
        input_bound = 1.0
        exact_error = abs(
            -(input_bound / lambda_value) * before + (2.0 * input_bound / lambda_value) * after
        )
        bound = (2.0 * input_bound / lambda_value) * after
        constant_rows.append(
            {
                "order": order,
                "damping": lambda_value,
                "switch_scale": switch_scale,
                "dwell_scale": dwell_scale,
                "exact_error": exact_error,
                "theory_bound": bound,
                "bound_utilization": exact_error / bound if bound else 0.0,
                "pass": int(exact_error <= bound + 1.0e-12),
            }
        )

    write_csv(output / "s1_full_memory_dwell.csv", dwell_rows)
    write_csv(output / "s1_step_convergence.csv", convergence_rows)
    write_csv(output / "s1_constant_switch_stress.csv", constant_rows)
    convergence_by_order: dict[float, list[float]] = defaultdict(list)
    for row in convergence_rows:
        convergence_by_order[float(row["order"])].append(float(row["final_oracle_error"]))
    monotone_convergence = all(
        errors[2] < errors[1] < errors[0] for errors in convergence_by_order.values()
    )
    point_rows = [row for row in dwell_rows if row["time_since_switch"] != "demodulated_last_half"]
    summary = {
        "pointwise_cases": len(point_rows),
        "pointwise_passed": sum(row["pass"] for row in point_rows),
        "max_bound_utilization": max(row["bound_utilization"] for row in point_rows),
        "max_implicit_oracle_error": max(row["implicit_oracle_error"] for row in point_rows),
        "max_pece_oracle_error": max(row["pece_oracle_error"] for row in point_rows),
        "demodulation_cases": len(dwell_rows) - len(point_rows),
        "demodulation_passed": sum(row["pass"] for row in dwell_rows)
        - sum(row["pass"] for row in point_rows),
        "constant_switch_cases": len(constant_rows),
        "constant_switch_passed": sum(row["pass"] for row in constant_rows),
        "max_constant_bound_utilization": max(row["bound_utilization"] for row in constant_rows),
        "monotone_step_convergence": monotone_convergence,
    }
    summary["status"] = (
        "PASS"
        if summary["pointwise_passed"] == summary["pointwise_cases"]
        and summary["constant_switch_passed"] == summary["constant_switch_cases"]
        and summary["monotone_step_convergence"]
        else "FAIL"
    )
    return summary


def _record_s2_triples(
    rows: list[dict[str, Any]],
    scenario: str,
    severity: float,
    frequencies: tuple[float, ...],
    responses: tuple[complex, ...],
    response_error: float,
    protocol_certified: bool,
    expected: str,
) -> None:
    for triple in itertools.combinations(range(len(frequencies)), 3):
        triple_frequencies = tuple(frequencies[index] for index in triple)
        triple_responses = tuple(responses[index] for index in triple)
        result = evaluate_chain_certificate(
            triple_responses,
            triple_frequencies,
            response_error,
            order_bounds=(0.4, 1.0),
            max_order_error=0.05,
            protocol_certified=protocol_certified,
        )
        rows.append(
            {
                "scenario": scenario,
                "severity": severity,
                "triple": "-".join(str(index) for index in triple),
                "expected": expected,
                "protocol_certified": int(protocol_certified),
                "status": result.status,
                "numeric_pass": int(result.numeric_pass),
                "reasons": "|".join(result.reasons),
                "point_order": result.point_order,
                "point_order_error": result.point_order_error,
                "interval_lower": (result.order_interval.lower if result.order_interval else None),
                "interval_upper": (result.order_interval.upper if result.order_interval else None),
                "imaginary_invariant": (
                    result.invariant.imag if result.invariant is not None else None
                ),
                "invariant_error": result.invariant_error,
                "reciprocal_margin": result.reciprocal_margin,
                "denominator_margin": result.denominator_margin,
            }
        )


def run_s2(output: Path) -> dict[str, Any]:
    order = 0.8
    damping = 0.7
    morphology = 1.0 - 0.3j
    heart_rates = np.asarray((50.0, 65.0, 80.0, 100.0, 125.0))
    frequencies = tuple(2.0 * math.pi * heart_rates / 60.0)
    base = tuple(
        harmonic_response(order, damping, morphology, frequency) for frequency in frequencies
    )
    tiny_error = 1.0e-6 * min(abs(value) for value in base)
    centered = np.asarray((-1.0, -0.5, 0.0, 0.5, 1.0))
    alternating = np.asarray((-1.0, 1.0, -1.0, 1.0, -1.0))
    rows: list[dict[str, Any]] = []

    _record_s2_triples(rows, "valid", 0.0, frequencies, base, tiny_error, True, "PASS")
    common_factor = 2.4 * np.exp(0.7j)
    common_gauge = tuple(common_factor * value for value in base)
    _record_s2_triples(
        rows,
        "common_gain_phase",
        0.0,
        frequencies,
        common_gauge,
        tiny_error,
        True,
        "PASS",
    )

    for severity in (0.01, 0.05, 0.2):
        q_values = tuple(
            morphology * np.exp(severity * amplitude + 1j * severity * phase)
            for amplitude, phase in zip(centered, alternating, strict=True)
        )
        responses = tuple(
            harmonic_response(order, damping, q_value, frequency)
            for q_value, frequency in zip(q_values, frequencies, strict=True)
        )
        _record_s2_triples(
            rows,
            "independent_q",
            severity,
            frequencies,
            responses,
            tiny_error,
            False,
            "NOT_CERTIFIABLE",
        )

        lambda_values = tuple(damping * np.exp(severity * value) for value in centered)
        responses = tuple(
            harmonic_response(order, lambda_value, morphology, frequency)
            for lambda_value, frequency in zip(lambda_values, frequencies, strict=True)
        )
        _record_s2_triples(
            rows,
            "independent_lambda",
            severity,
            frequencies,
            responses,
            tiny_error,
            False,
            "NOT_CERTIFIABLE",
        )

        normalized = tuple(
            value * np.exp(severity * pattern)
            for value, pattern in zip(base, centered, strict=True)
        )
        _record_s2_triples(
            rows,
            "per_window_normalization",
            severity,
            frequencies,
            normalized,
            tiny_error,
            False,
            "NOT_CERTIFIABLE",
        )

    for degrees in (1.0, 5.0, 20.0):
        radians = math.radians(degrees)
        shifted = tuple(
            value * np.exp(1j * radians * pattern)
            for value, pattern in zip(base, centered, strict=True)
        )
        _record_s2_triples(
            rows,
            "independent_phase",
            degrees,
            frequencies,
            shifted,
            tiny_error,
            False,
            "NOT_CERTIFIABLE",
        )

    target_order = 0.65
    target_damping = 0.9
    target_morphology = 0.8 + 0.4j
    adversarial = tuple(
        harmonic_response(target_order, target_damping, target_morphology, frequency)
        for frequency in frequencies
    )
    _record_s2_triples(
        rows,
        "adversarial_observation_equivalent_q",
        abs(target_order - order),
        frequencies,
        adversarial,
        tiny_error,
        False,
        "NOT_CERTIFIABLE_NUMERIC_PASS",
    )

    for spacing_bpm in (0.25, 1.0, 3.0, 8.0):
        close_hr = 75.0 + centered * spacing_bpm
        close_frequencies = tuple(2.0 * math.pi * close_hr / 60.0)
        close_responses = tuple(
            harmonic_response(order, damping, morphology, frequency)
            for frequency in close_frequencies
        )
        close_error = 1.0e-3 * min(abs(value) for value in close_responses)
        _record_s2_triples(
            rows,
            "close_heart_rates",
            spacing_bpm,
            close_frequencies,
            close_responses,
            close_error,
            True,
            "CONDITION_DEPENDENT",
        )

    for delay_ms in (1.0, 5.0, 20.0, 80.0):
        delay = delay_ms / 1000.0
        delayed = tuple(
            value * np.exp(-1j * frequency * delay)
            for value, frequency in zip(base, frequencies, strict=True)
        )
        _record_s2_triples(
            rows,
            "uncalibrated_delay",
            delay_ms,
            frequencies,
            delayed,
            tiny_error,
            False,
            "NOT_CERTIFIABLE",
        )
        corrected = tuple(
            value * np.exp(1j * frequency * delay)
            for value, frequency in zip(delayed, frequencies, strict=True)
        )
        _record_s2_triples(
            rows,
            "delay_corrected",
            delay_ms,
            frequencies,
            corrected,
            tiny_error,
            True,
            "PASS",
        )

    short_dwell = 1.0
    short_step = 0.01
    short_schedule = make_harmonic_schedule(frequencies, short_dwell, short_step, morphology)
    short_state = solve_linear_caputo_product_integration(
        order,
        damping,
        short_schedule.time,
        short_schedule.forcing,
        0.0,
    )
    short_coefficients: list[complex] = []
    for index in range(len(frequencies)):
        left = round((index * short_dwell + short_dwell / 2.0) / short_step)
        right = round(((index + 1) * short_dwell) / short_step)
        short_coefficients.append(
            demodulate_harmonic(
                short_schedule.time[left : right + 1],
                short_state[left : right + 1],
                short_schedule.phase[left : right + 1],
            )
        )
    short_error = 1.0e-6 * min(abs(value) for value in short_coefficients)
    _record_s2_triples(
        rows,
        "short_dwell",
        short_dwell,
        frequencies,
        tuple(short_coefficients),
        short_error,
        False,
        "NOT_CERTIFIABLE",
    )

    write_csv(output / "s2_gate_scenarios.csv", rows)
    grouped: dict[tuple[str, float], Counter[str]] = defaultdict(Counter)
    numeric_grouped: dict[tuple[str, float], Counter[int]] = defaultdict(Counter)
    for row in rows:
        key = (str(row["scenario"]), float(row["severity"]))
        grouped[key][str(row["status"])] += 1
        numeric_grouped[key][int(row["numeric_pass"])] += 1
    group_rows = []
    for key in sorted(grouped):
        counts = grouped[key]
        numeric = numeric_grouped[key]
        total = sum(counts.values())
        group_rows.append(
            {
                "scenario": key[0],
                "severity": key[1],
                "cases": total,
                "pass": counts["PASS"],
                "reject_numeric": counts["REJECT_NUMERIC"],
                "not_certifiable": counts["NOT_CERTIFIABLE"],
                "numeric_pass": numeric[1],
                "numeric_pass_rate": numeric[1] / total,
            }
        )
    write_csv(output / "s2_gate_summary.csv", group_rows)

    positive_rows = [
        row for row in rows if row["scenario"] in {"valid", "common_gain_phase", "delay_corrected"}
    ]
    invalid_protocol_rows = [row for row in rows if not row["protocol_certified"]]
    adversarial_rows = [
        row for row in rows if row["scenario"] == "adversarial_observation_equivalent_q"
    ]
    summary = {
        "cases": len(rows),
        "positive_control_cases": len(positive_rows),
        "positive_control_passed": sum(row["status"] == "PASS" for row in positive_rows),
        "invalid_protocol_cases": len(invalid_protocol_rows),
        "invalid_protocol_not_certifiable": sum(
            row["status"] == "NOT_CERTIFIABLE" for row in invalid_protocol_rows
        ),
        "adversarial_cases": len(adversarial_rows),
        "adversarial_numeric_pass": sum(row["numeric_pass"] for row in adversarial_rows),
        "adversarial_not_certifiable": sum(
            row["status"] == "NOT_CERTIFIABLE" for row in adversarial_rows
        ),
    }
    summary["status"] = (
        "PASS"
        if summary["positive_control_passed"] == summary["positive_control_cases"]
        and summary["invalid_protocol_not_certifiable"] == summary["invalid_protocol_cases"]
        and summary["adversarial_numeric_pass"] == summary["adversarial_cases"]
        and summary["adversarial_not_certifiable"] == summary["adversarial_cases"]
        else "FAIL"
    )
    return summary


def run_s3(output: Path) -> dict[str, Any]:
    rng = np.random.default_rng(SEED)
    orders = (0.5, 0.7, 0.9, 0.97)
    dampings = (0.3, 0.8, 1.5)
    morphologies = (1.0 + 0.2j, 0.6 - 1.1j)
    triples = (
        (0.7, 1.0, 1.4),
        (0.7, 1.4, 2.6),
        (0.7, 1.9, 2.6),
        (1.0, 1.4, 1.9),
    )
    relative_errors = (1.0e-4, 1.0e-3, 5.0e-3, 1.0e-2)
    draws = 128
    raw_metrics: dict[tuple[Any, ...], dict[str, list[float] | int]] = {}
    worst_rows: list[dict[str, Any]] = []
    total_violations = 0
    total_trials = 0

    for order, damping, morphology, frequencies, relative_error in itertools.product(
        orders, dampings, morphologies, triples, relative_errors
    ):
        truth = tuple(
            harmonic_response(order, damping, morphology, frequency) for frequency in frequencies
        )
        radius = relative_error * min(abs(value) for value in truth)
        key = (
            order,
            damping,
            round(abs(morphology), 8),
            "-".join(str(value) for value in frequencies),
            relative_error,
        )
        metrics: dict[str, list[float] | int] = {
            "u_w": [],
            "u_r": [],
            "u_alpha": [],
            "chain_width": [],
            "disk_width": [],
            "chain_coverage": 0,
            "disk_true_feasible": 0,
            "numeric_pass": 0,
            "primary_pass": 0,
            "empty_chain": 0,
        }
        max_alpha_row: dict[str, Any] | None = None
        max_alpha_utilization = -1.0

        for draw in range(draws):
            perturbations = tuple(
                radius * math.sqrt(float(rng.random())) * np.exp(2j * math.pi * float(rng.random()))
                for _ in range(3)
            )
            measured = tuple(
                value + perturbation
                for value, perturbation in zip(truth, perturbations, strict=True)
            )
            certificate = evaluate_chain_certificate(
                measured,
                frequencies,
                radius,
                order_bounds=(0.4, 1.0),
                max_order_error=0.05,
            )
            primary_certificate = evaluate_chain_certificate(
                measured,
                frequencies,
                radius,
                order_bounds=(0.4, 1.0),
                max_order_error=0.02,
            )
            true_inverse = tuple(1.0 / value for value in truth)
            measured_inverse = tuple(1.0 / value for value in measured)
            posterior_inverse_bounds = tuple(
                radius / (abs(value) * (abs(value) - radius)) for value in measured
            )
            u_w = max(
                abs(measured_value - true_value) / bound
                for measured_value, true_value, bound in zip(
                    measured_inverse,
                    true_inverse,
                    posterior_inverse_bounds,
                    strict=True,
                )
            )
            true_invariant = order_invariant(order, frequencies)
            measured_invariant = response_invariant(measured)
            if certificate.invariant_error is None or certificate.point_order is None:
                metrics["empty_chain"] = int(metrics["empty_chain"]) + 1
                continue
            u_r = abs(measured_invariant - true_invariant) / certificate.invariant_error
            alpha_error = abs(certificate.point_order - order)
            alpha_bound = certificate.point_order_error or math.inf
            u_alpha = alpha_error / alpha_bound if alpha_bound > 0.0 else 0.0
            chain_contains = bool(
                certificate.order_interval
                and certificate.order_interval.contains(order, tolerance=2.0e-9)
            )
            if chain_contains:
                metrics["chain_coverage"] = int(metrics["chain_coverage"]) + 1
            if certificate.numeric_pass:
                metrics["numeric_pass"] = int(metrics["numeric_pass"]) + 1
            if primary_certificate.numeric_pass:
                metrics["primary_pass"] = int(metrics["primary_pass"]) + 1
            if certificate.order_interval:
                cast_widths = metrics["chain_width"]
                assert isinstance(cast_widths, list)
                cast_widths.append(certificate.order_interval.width)

            centers, radii = inverse_response_disks(measured, radius)
            ratio = order_invariant(order, frequencies)
            disk_residual = abs(centers[2] - centers[0] - ratio * (centers[1] - centers[0])) - (
                radii[2] + (ratio - 1.0) * radii[0] + ratio * radii[1]
            )
            if disk_residual <= 2.0e-10:
                metrics["disk_true_feasible"] = int(metrics["disk_true_feasible"]) + 1

            if draw < 4:
                disk_intervals = disk_order_intervals(
                    measured,
                    frequencies,
                    radius,
                    order_bounds=(0.4, 1.0),
                    grid_size=1025,
                )
                cast_disk_widths = metrics["disk_width"]
                assert isinstance(cast_disk_widths, list)
                cast_disk_widths.append(sum(interval.width for interval in disk_intervals))

            for name, value in (("u_w", u_w), ("u_r", u_r), ("u_alpha", u_alpha)):
                cast_values = metrics[name]
                assert isinstance(cast_values, list)
                cast_values.append(float(value))
            violated = u_w > 1.0 + 2.0e-9 or u_r > 1.0 + 2.0e-9 or u_alpha > 1.0 + 2.0e-8
            total_violations += int(violated)
            total_trials += 1
            if u_alpha > max_alpha_utilization:
                max_alpha_utilization = u_alpha
                max_alpha_row = {
                    "order": order,
                    "damping": damping,
                    "morphology_abs": abs(morphology),
                    "frequencies": key[3],
                    "relative_error": relative_error,
                    "draw": draw,
                    "u_w": u_w,
                    "u_r": u_r,
                    "u_alpha": u_alpha,
                    "alpha_error": alpha_error,
                    "alpha_bound": alpha_bound,
                    "chain_contains": int(chain_contains),
                    "disk_residual": disk_residual,
                }
        raw_metrics[key] = metrics
        if max_alpha_row:
            worst_rows.append(max_alpha_row)

    phase_set = (0.0, 0.5 * math.pi, math.pi, 1.5 * math.pi)
    adversarial_rows: list[dict[str, Any]] = []
    representative_truth = tuple(
        harmonic_response(0.97, 0.8, 1.0 - 0.3j, frequency) for frequency in (0.7, 1.4, 2.6)
    )
    for relative_error in relative_errors:
        radius = relative_error * min(abs(value) for value in representative_truth)
        for phases in itertools.product(phase_set, repeat=3):
            measured = tuple(
                value + radius * np.exp(1j * phase)
                for value, phase in zip(representative_truth, phases, strict=True)
            )
            certificate = evaluate_chain_certificate(
                measured,
                (0.7, 1.4, 2.6),
                radius,
                order_bounds=(0.4, 1.0),
                max_order_error=1.0,
            )
            measured_invariant = response_invariant(measured)
            true_invariant = order_invariant(0.97, (0.7, 1.4, 2.6))
            u_r = (
                abs(measured_invariant - true_invariant) / certificate.invariant_error
                if certificate.invariant_error
                else math.nan
            )
            u_alpha = (
                abs((certificate.point_order or 0.97) - 0.97) / certificate.point_order_error
                if certificate.point_order_error
                else math.nan
            )
            adversarial_rows.append(
                {
                    "relative_error": relative_error,
                    "phase1": phases[0],
                    "phase2": phases[1],
                    "phase3": phases[2],
                    "u_r": u_r,
                    "u_alpha": u_alpha,
                    "chain_contains": int(
                        bool(
                            certificate.order_interval and certificate.order_interval.contains(0.97)
                        )
                    ),
                }
            )

    summary_rows: list[dict[str, Any]] = []
    for key, metrics in raw_metrics.items():
        u_w = metrics["u_w"]
        u_r = metrics["u_r"]
        u_alpha = metrics["u_alpha"]
        chain_width = metrics["chain_width"]
        disk_width = metrics["disk_width"]
        assert isinstance(u_w, list)
        assert isinstance(u_r, list)
        assert isinstance(u_alpha, list)
        assert isinstance(chain_width, list)
        assert isinstance(disk_width, list)
        certifiable = len(u_alpha)
        summary_rows.append(
            {
                "order": key[0],
                "damping": key[1],
                "morphology_abs": key[2],
                "frequencies": key[3],
                "relative_error": key[4],
                "draws": draws,
                "certifiable": certifiable,
                "primary_pass_0p02": metrics["primary_pass"],
                "exploratory_pass_0p05": metrics["numeric_pass"],
                "chain_coverage": metrics["chain_coverage"] / certifiable,
                "disk_true_feasible": metrics["disk_true_feasible"] / certifiable,
                "max_u_w": max(u_w),
                "max_u_r": max(u_r),
                "max_u_alpha": max(u_alpha),
                "p95_u_alpha": quantile(u_alpha, 0.95),
                "median_chain_width": quantile(chain_width, 0.5),
                "median_disk_width_subset": quantile(disk_width, 0.5),
            }
        )

    write_csv(output / "s3_bound_summary.csv", summary_rows)
    write_csv(output / "s3_worst_witnesses.csv", worst_rows)
    write_csv(output / "s3_adversarial_circle.csv", adversarial_rows)
    max_u_w = max(row["max_u_w"] for row in summary_rows)
    max_u_r = max(row["max_u_r"] for row in summary_rows)
    max_u_alpha = max(row["max_u_alpha"] for row in summary_rows)
    adversarial_max_u_r = max(row["u_r"] for row in adversarial_rows)
    adversarial_max_u_alpha = max(row["u_alpha"] for row in adversarial_rows)
    summary = {
        "configurations": len(summary_rows),
        "bounded_random_trials": total_trials,
        "bound_violations": total_violations,
        "min_chain_coverage": min(row["chain_coverage"] for row in summary_rows),
        "min_disk_true_feasible": min(row["disk_true_feasible"] for row in summary_rows),
        "max_u_w": max_u_w,
        "max_u_r": max_u_r,
        "max_u_alpha": max_u_alpha,
        "adversarial_trials": len(adversarial_rows),
        "adversarial_max_u_r": adversarial_max_u_r,
        "adversarial_max_u_alpha": adversarial_max_u_alpha,
        "adversarial_chain_coverage": sum(row["chain_contains"] for row in adversarial_rows)
        / len(adversarial_rows),
    }
    summary["status"] = (
        "PASS"
        if summary["bound_violations"] == 0
        and summary["min_chain_coverage"] == 1.0
        and summary["min_disk_true_feasible"] == 1.0
        and summary["adversarial_chain_coverage"] == 1.0
        else "FAIL"
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    summaries: dict[str, Any] = {"seed": SEED}
    for name, runner in (("S0", run_s0), ("S1", run_s1), ("S2", run_s2), ("S3", run_s3)):
        stage_started = time.perf_counter()
        summaries[name] = runner(args.output)
        summaries[name]["elapsed_seconds"] = time.perf_counter() - stage_started
        print(json.dumps({name: summaries[name]}, ensure_ascii=False, sort_keys=True))
        if summaries[name]["status"] != "PASS":
            summaries["overall_status"] = "FAIL"
            break
    else:
        summaries["overall_status"] = "PASS"
    summaries["elapsed_seconds"] = time.perf_counter() - started
    (args.output / "summary.json").write_text(
        json.dumps(summaries, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"summary": summaries}, ensure_ascii=False, sort_keys=True))
    return 0 if summaries["overall_status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
