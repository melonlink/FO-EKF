import math
from pathlib import Path

import numpy as np
import pytest

from fo_ekf import fantasia_pilot as fantasia_module
from fo_ekf.fantasia_pilot import (
    FANTASIA_V1_SHA256,
    FantasiaPilotConfig,
    build_rr_windows,
    canonical_json_hash,
    cardiac_phase,
    coefficient_metrics,
    estimate_window_coefficients,
    fantasia_implementation_manifest,
    fantasia_runtime_environment,
    fit_fractional_model,
    fit_order_one_model,
    pilot_decision_payload,
    predict_model,
)


def test_locked_config_has_five_minute_guard_and_stable_hash() -> None:
    config = FantasiaPilotConfig()
    config.validate()

    assert config.guard_end_seconds - config.identification_end_seconds == 300.0
    assert config.lock_hash() == config.lock_hash()
    assert config.front_end == "identity"
    assert config.normalization == "none"
    assert len(FANTASIA_V1_SHA256) == 12
    assert set(FANTASIA_V1_SHA256) == {
        f"{record}.{extension}" for record in config.records for extension in ("hea", "dat", "ecg")
    }
    assert (
        canonical_json_hash(FANTASIA_V1_SHA256)
        == "a8d8c4a42479d2721f4797d73dcb4e0d6aa7c8170fde100e9183914895256ffb"
    )


def test_production_manifest_binds_sources_and_runtime() -> None:
    manifest = fantasia_implementation_manifest()
    environment = fantasia_runtime_environment()

    assert len(manifest["implementation_sha256"]) == 64
    assert set(manifest["source_sha256"]) == {
        "scripts/run_fantasia_pilot.py",
        "src/fo_ekf/fantasia_pilot.py",
        "src/fo_ekf/wls_execution.py",
    }
    assert all(len(value) == 64 for value in manifest["source_sha256"].values())
    assert {"python", "numpy", "scipy", "wfdb", "python_flint", "flint"} <= set(environment)


def test_sealed_plain_text_writers_emit_lf_bytes(tmp_path: Path) -> None:
    json_path = tmp_path / "payload.json"
    csv_path = tmp_path / "payload.csv"

    fantasia_module._write_json(json_path, {"message": "line one\nline two"})
    fantasia_module._write_csv(csv_path, [{"a": 1, "b": 2}], ("a", "b"))

    for path in (json_path, csv_path):
        payload = path.read_bytes()
        assert b"\r" not in payload
        assert payload.endswith(b"\n")

    with pytest.raises(ValueError, match="LF line endings"):
        fantasia_module._write_text_lf(tmp_path / "invalid.txt", "bad\r\n")


def test_empirical_decision_is_derived_from_counts() -> None:
    negative = pilot_decision_payload(
        record_count=4,
        locked_record_count=2,
        alpha_boundary_hits=4,
        materially_better_than_both=0,
        locked_better_than_constant=1,
        locked_better_than_shuffle=1,
    )
    changed = pilot_decision_payload(
        record_count=4,
        locked_record_count=2,
        alpha_boundary_hits=3,
        materially_better_than_both=1,
        locked_better_than_constant=2,
        locked_better_than_shuffle=2,
    )

    assert negative["distinct_fractional_order_evidence"] == "NOT_SUPPORTED"
    assert negative["confirmatory_rate_response_evidence"] == "NOT_SUPPORTED"
    assert changed["distinct_fractional_order_evidence"] == "INCONCLUSIVE"
    assert changed["confirmatory_rate_response_evidence"] == "INCONCLUSIVE"


def test_empirical_decision_rejects_inconsistent_counts() -> None:
    with pytest.raises(ValueError, match="internally inconsistent"):
        pilot_decision_payload(
            record_count=4,
            locked_record_count=2,
            alpha_boundary_hits=5,
            materially_better_than_both=0,
            locked_better_than_constant=1,
            locked_better_than_shuffle=1,
        )


def test_rr_windows_are_complete_nonoverlapping_and_unknown_breaks_run() -> None:
    samples = list(range(0, 42))
    symbols = ["N"] * len(samples)
    symbols[18] = "?"
    windows = build_rr_windows(
        samples,
        symbols,
        sampling_frequency=1.0,
        segment_start_seconds=0.0,
        segment_end_seconds=41.0,
        split="identification",
    )

    assert len(windows) == 2
    assert windows[0].beat_samples == tuple(range(17))
    assert windows[1].beat_samples == tuple(range(19, 36))
    assert windows[0].end_sample <= windows[1].start_sample
    assert all(len(window.beat_samples) == 17 for window in windows)


def test_identity_wls_recovers_fixed_phase_coefficients_without_normalization() -> None:
    beats = tuple(range(10, 171, 10))
    phase = cardiac_phase(beats)
    expected = (1.2 - 0.7j, -0.3 + 0.2j, 0.08 - 0.04j)
    response = np.full(200, 4.25)
    within = np.full_like(phase, 4.25)
    for harmonic, coefficient in enumerate(expected, start=1):
        within += 2.0 * coefficient.real * np.cos(harmonic * phase)
        within += -2.0 * coefficient.imag * np.sin(harmonic * phase)
    response[beats[0] : beats[-1]] = within

    estimate = estimate_window_coefficients(response, beats, (1, 2, 3))

    assert estimate.dc == pytest.approx(4.25, abs=1.0e-12)
    assert estimate.coefficients == pytest.approx(expected, abs=1.0e-12)
    assert estimate.residual_rmse <= 1.0e-12
    assert estimate.gram_rank == 7
    assert estimate.gram_min_eigenvalue > 0.0


def test_common_fractional_fit_recovers_synthetic_parameters() -> None:
    rates = np.linspace(0.7, 1.8, 24)
    harmonics = (1, 2, 3)
    alpha = 0.68
    damping = 0.42
    tau_star_seconds = 2.0
    morphology = np.asarray((1.2 - 0.5j, -0.4 + 0.2j, 0.1 + 0.03j))
    dimensionless_frequency = (
        2.0 * np.pi * tau_star_seconds * rates[:, None] * np.asarray(harmonics)[None, :]
    )
    observations = morphology[None, :] / (damping + np.power(1j * dimensionless_frequency, alpha))

    fitted = fit_fractional_model(
        rates,
        observations,
        harmonics,
        tau_star_seconds=tau_star_seconds,
    )
    predicted = predict_model(fitted, rates, harmonics, model_name="fractional_common")

    assert fitted["alpha"] == pytest.approx(alpha, abs=2.0e-6)
    assert fitted["lambda"] == pytest.approx(damping, abs=2.0e-5)
    assert fitted["tau_star_seconds"] == tau_star_seconds
    assert fitted["q"] == pytest.approx(morphology, abs=2.0e-5)
    assert coefficient_metrics(observations, predicted)["rmse_mV"] < 1.0e-7


def test_alpha_one_is_exactly_nested_and_metric_uses_complex_energy() -> None:
    rates = np.linspace(0.8, 1.4, 12)
    harmonics = (1, 2)
    damping = 0.8
    morphology = np.asarray((0.7 - 0.3j, 0.2 + 0.1j))
    omega = 2.0 * np.pi * rates[:, None] * np.asarray(harmonics)[None, :]
    observations = morphology[None, :] / (damping + 1j * omega)

    fitted = fit_order_one_model(rates, observations, harmonics)
    predicted = predict_model(fitted, rates, harmonics, model_name="alpha_one_nested")
    metrics = coefficient_metrics(observations, predicted)

    assert fitted["alpha"] == 1.0
    assert fitted["lambda"] == pytest.approx(damping, abs=2.0e-6)
    assert metrics["rmse_mV"] < 1.0e-8
    assert math.isfinite(metrics["nrmse_energy"])
