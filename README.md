# FO-EKF

Fractional-order state estimation for a lightweight, ECG-driven cardiac model.

This repository starts from the research plan in
[`PaperPlan/innovation_points_fractional_order_cardiac_digital_twin.md`](PaperPlan/innovation_points_fractional_order_cardiac_digital_twin.md).
The first milestone is a falsifiable feasibility study: determine whether a fractional-order model and observer improve out-of-subject ECG reconstruction or forecasting over parameter-matched integer-order baselines.

## Branches

- `main`: reviewed baseline and final releases.
- `dev-codex`: research checkpoints and staged work.

## Data

Raw/shared data remain outside this repository under `../PUBLIC` and are treated as read-only. The defaults are:

```text
ECG_DATA_DIR=../PUBLIC/data_PTB-XL
ECG_PROCESSED_DIR=./data/processed
ECG_OUTPUT_DIR=./artifacts
```

Copy `.env.example` to `.env` only when local overrides are needed. Never commit raw ECG files or patient-derived caches.
The current machine stores PTB-XL under `../PUBLIC/data_PTB-XL/raw`; the environment check verifies that location instead of trusting a stale path description.

## Python environment

Python 3.10+ is supported. On Windows PowerShell:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"
.\.venv\Scripts\python.exe scripts\check_environment.py
.\.venv\Scripts\python.exe -m pytest
```

The bootstrap package deliberately contains only path/configuration code. Model equations, fractional discretization, and FO-EKF code will be added only after the feasibility and identifiability review fixes their definitions.

## Feasibility and theory checkpoints

The initial review found that the original broad "first FO-EKF cardiac digital twin" claim is not defensible, but a narrower correlation-consistent, identifiability-gated finite-memory observer is conditionally feasible.

The subsequent control-theory review narrows the primary paper further: a periodicity-consistent mixed integer-fractional ECG surrogate, a gauge-invariant multi-heart-rate certificate for fractional-order identifiability, and an approximation-aware Mittag-Leffler observer-error tube. Generic FO-EKF, finite-memory compensation, and LMI stability remain prior-art tools rather than headline contributions.

- [Feasibility and novelty review](research/feasibility_review_2026-07-11.md)
- [Control-theory innovation specification](research/control_theory_innovation_spec_2026-07-11.md)
- [S0-S3 theory-validation execution report](research/s0_s3_execution_report_2026-07-11.md)
- [Environment snapshot](research/environment_snapshot_2026-07-11.md)
- [Curated literature manifest](literature/manifest.csv)
- [Verified local literature lock](literature/download-lock.json)

The analytical recovery, finite-dwell bound, counterexamples, and deterministic order-error chain have passed their data-free S0-S3 tests. Patient-data calibration remains gated on the preregistered Fantasia protocol: no per-window amplitude normalization, explicit phasor-error radii, common-parameter checks, and leave-one-rate prediction.

## Research artifacts

- `PaperPlan/`: original research concept.
- `literature/`: reproducible literature manifest and local download script.
- `research/`: feasibility, novelty, and experimental-design decisions.
- `src/fo_ekf/`: implementation package (after method freeze).
- `tests/`: tests that do not require the raw datasets.
