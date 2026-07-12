"""Deterministic multi-rate order certificates and rejection gates."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .multirate import (
    invert_order_invariant,
    order_invariant,
    order_invariant_slope,
    response_invariant,
)


@dataclass(frozen=True)
class OrderInterval:
    """Closed interval of fractional orders."""

    lower: float
    upper: float

    @property
    def width(self) -> float:
        return self.upper - self.lower

    def contains(self, value: float, *, tolerance: float = 1.0e-10) -> bool:
        return self.lower - tolerance <= value <= self.upper + tolerance


@dataclass(frozen=True)
class ChainCertificate:
    """Closed-form posterior certificate for one frequency triple."""

    status: str
    reasons: tuple[str, ...]
    numeric_pass: bool
    invariant: complex | None
    invariant_error: float | None
    point_order: float | None
    point_order_error: float | None
    order_interval: OrderInterval | None
    boundary_hit: bool
    reciprocal_margin: float
    denominator_margin: float
    slope_lower_bound: float


def _error_tuple(errors: float | Sequence[float]) -> tuple[float, float, float]:
    if isinstance(errors, (int, float)):
        values = (float(errors),) * 3
    else:
        values = tuple(float(value) for value in errors)
    if len(values) != 3 or any(value < 0.0 for value in values):
        raise ValueError("three nonnegative response-error radii are required")
    return values


def slope_lower_bound(order_bounds: tuple[float, float], frequencies: Sequence[float]) -> float:
    """Return inf F' on the interval; strict convexity places it at the left end."""

    lower, upper = (float(value) for value in order_bounds)
    if not 0.0 < lower < upper:
        raise ValueError("order bounds must satisfy 0 < lower < upper")
    return order_invariant_slope(lower, frequencies)


def _projected_order(
    invariant_real: float,
    frequencies: Sequence[float],
    order_bounds: tuple[float, float],
) -> tuple[float, bool]:
    lower, upper = order_bounds
    invariant_lower = order_invariant(lower, frequencies)
    invariant_upper = order_invariant(upper, frequencies)
    projected = min(max(invariant_real, invariant_lower), invariant_upper)
    boundary_hit = projected != invariant_real or projected in (invariant_lower, invariant_upper)
    return (
        invert_order_invariant(projected, frequencies, lower=lower, upper=upper),
        boundary_hit,
    )


def evaluate_chain_certificate(
    measured_responses: Sequence[complex],
    frequencies: Sequence[float],
    response_errors: float | Sequence[float],
    *,
    order_bounds: tuple[float, float] = (0.4, 1.0),
    max_order_error: float = 0.05,
    protocol_certified: bool = True,
    numerical_tolerance: float = 1.0e-12,
) -> ChainCertificate:
    """Evaluate the posterior chain bound and PASS/REJECT/NOT_CERTIFIABLE state."""

    responses = tuple(complex(value) for value in measured_responses)
    if len(responses) != 3:
        raise ValueError("exactly three measured responses are required")
    errors = _error_tuple(response_errors)
    slope = slope_lower_bound(order_bounds, frequencies)
    reasons: list[str] = []
    reciprocal_margin = min(
        (abs(response) - error) / abs(response) if response != 0.0 else -math.inf
        for response, error in zip(responses, errors, strict=True)
    )
    if any(abs(response) <= error for response, error in zip(responses, errors, strict=True)):
        reasons.append("reciprocal_disk_contains_zero")
        return ChainCertificate(
            status="NOT_CERTIFIABLE",
            reasons=tuple((["protocol_G0"] if not protocol_certified else []) + reasons),
            numeric_pass=False,
            invariant=None,
            invariant_error=None,
            point_order=None,
            point_order_error=None,
            order_interval=None,
            boundary_hit=False,
            reciprocal_margin=reciprocal_margin,
            denominator_margin=-math.inf,
            slope_lower_bound=slope,
        )

    inverse = tuple(1.0 / response for response in responses)
    inverse_errors = tuple(
        error / (abs(response) * (abs(response) - error))
        for response, error in zip(responses, errors, strict=True)
    )
    difference = inverse[1] - inverse[0]
    difference_error = inverse_errors[0] + inverse_errors[1]
    denominator_margin = (abs(difference) - difference_error) / abs(difference)
    if abs(difference) <= difference_error:
        reasons.append("unseparated_inverse_response")
        return ChainCertificate(
            status="NOT_CERTIFIABLE",
            reasons=tuple((["protocol_G0"] if not protocol_certified else []) + reasons),
            numeric_pass=False,
            invariant=None,
            invariant_error=None,
            point_order=None,
            point_order_error=None,
            order_interval=None,
            boundary_hit=False,
            reciprocal_margin=reciprocal_margin,
            denominator_margin=denominator_margin,
            slope_lower_bound=slope,
        )

    invariant = response_invariant(responses)
    invariant_error = (
        inverse_errors[2]
        + abs(invariant) * inverse_errors[1]
        + abs(1.0 - invariant) * inverse_errors[0]
    ) / (abs(difference) - difference_error)
    point_order, boundary_hit = _projected_order(invariant.real, frequencies, order_bounds)
    point_order_error = invariant_error / slope

    invariant_lower = order_invariant(order_bounds[0], frequencies)
    invariant_upper = order_invariant(order_bounds[1], frequencies)
    if abs(invariant.imag) > invariant_error + numerical_tolerance:
        reasons.append("imaginary_inconsistency")
    if invariant.real < invariant_lower - invariant_error - numerical_tolerance or (
        invariant.real > invariant_upper + invariant_error + numerical_tolerance
    ):
        reasons.append("invariant_out_of_range")
    if point_order_error > max_order_error:
        reasons.append("order_bound_too_wide")

    interval: OrderInterval | None = None
    if abs(invariant.imag) <= invariant_error + numerical_tolerance:
        real_half_width = math.sqrt(max(0.0, invariant_error**2 - invariant.imag**2))
        real_lower = max(invariant_lower, invariant.real - real_half_width)
        real_upper = min(invariant_upper, invariant.real + real_half_width)
        if real_lower <= real_upper + numerical_tolerance:
            real_lower = min(max(real_lower, invariant_lower), invariant_upper)
            real_upper = min(max(real_upper, invariant_lower), invariant_upper)
            interval = OrderInterval(
                invert_order_invariant(
                    real_lower,
                    frequencies,
                    lower=order_bounds[0],
                    upper=order_bounds[1],
                ),
                invert_order_invariant(
                    real_upper,
                    frequencies,
                    lower=order_bounds[0],
                    upper=order_bounds[1],
                ),
            )
        else:
            reasons.append("empty_chain_order_set")

    numeric_pass = not reasons and interval is not None
    incompatibility_reasons = {
        "imaginary_inconsistency",
        "invariant_out_of_range",
        "empty_chain_order_set",
    }
    if not protocol_certified:
        status = "NOT_CERTIFIABLE"
        reasons.insert(0, "protocol_G0")
    elif numeric_pass:
        status = "PASS"
    elif not incompatibility_reasons.intersection(reasons):
        status = "NOT_CERTIFIABLE"
    else:
        status = "REJECT_NUMERIC"
    return ChainCertificate(
        status=status,
        reasons=tuple(reasons),
        numeric_pass=numeric_pass,
        invariant=invariant,
        invariant_error=invariant_error,
        point_order=point_order,
        point_order_error=point_order_error,
        order_interval=interval,
        boundary_hit=boundary_hit,
        reciprocal_margin=reciprocal_margin,
        denominator_margin=denominator_margin,
        slope_lower_bound=slope,
    )


def inverse_response_disks(
    measured_responses: Sequence[complex],
    response_errors: float | Sequence[float],
) -> tuple[tuple[complex, ...], tuple[float, ...]]:
    """Map response uncertainty disks exactly through W=1/Z."""

    responses = tuple(complex(value) for value in measured_responses)
    errors = _error_tuple(response_errors)
    if len(responses) != 3:
        raise ValueError("exactly three measured responses are required")
    if any(abs(response) <= error for response, error in zip(responses, errors, strict=True)):
        raise ValueError("a response disk contains zero")
    denominators = tuple(
        abs(response) ** 2 - error**2 for response, error in zip(responses, errors, strict=True)
    )
    centers = tuple(
        response.conjugate() / denominator
        for response, denominator in zip(responses, denominators, strict=True)
    )
    radii = tuple(
        error / denominator for error, denominator in zip(errors, denominators, strict=True)
    )
    return centers, radii


def disk_order_intervals(
    measured_responses: Sequence[complex],
    frequencies: Sequence[float],
    response_errors: float | Sequence[float],
    *,
    order_bounds: tuple[float, float] = (0.4, 1.0),
    grid_size: int = 4097,
) -> tuple[OrderInterval, ...]:
    """Approximate the relaxed three-point disk-consistent order set.

    This enforces only the algebraic relation
    ``W3 - W1 = F(alpha) * (W2 - W1)``. It does not enforce a positive
    real damping value or joint physical consistency across more than one
    triplet/harmonic, so it is not the exact feasible set of the full model.
    """

    if grid_size < 33:
        raise ValueError("grid_size must be at least 33")
    centers, radii = inverse_response_disks(measured_responses, response_errors)

    def residual(order: float) -> float:
        ratio = order_invariant(order, frequencies)
        center_residual = centers[2] - centers[0] - ratio * (centers[1] - centers[0])
        radius = radii[2] + (ratio - 1.0) * radii[0] + ratio * radii[1]
        return abs(center_residual) - radius

    grid = np.linspace(order_bounds[0], order_bounds[1], grid_size)
    values = np.asarray([residual(float(order)) for order in grid])
    feasible = values <= 0.0
    if not np.any(feasible):
        return ()

    intervals: list[OrderInterval] = []
    start: float | None = None
    for index, is_feasible in enumerate(feasible):
        if is_feasible and start is None:
            if index == 0:
                start = float(grid[index])
            else:
                start = float(brentq(residual, float(grid[index - 1]), float(grid[index])))
        is_last = index == len(grid) - 1
        if start is not None and (not is_feasible or is_last):
            if is_feasible and is_last:
                end = float(grid[index])
            else:
                end = float(brentq(residual, float(grid[index - 1]), float(grid[index])))
            intervals.append(OrderInterval(start, end))
            start = None
    return tuple(intervals)
