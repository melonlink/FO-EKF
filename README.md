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
.\.venv\Scripts\python.exe -m pip install -e ".[dev,cert]"
.\.venv\Scripts\python.exe -m flint.test --quiet
.\.venv\Scripts\python.exe scripts\check_environment.py
.\.venv\Scripts\python.exe -m pytest
```

The `cert` extra pins `python-flint==0.9.0` for rigorous Arb/ACB ball
arithmetic. Ordinary floating-point geometry may propose candidates, but only
the ball verifier may issue an R1 certificate. Python-FLINT's bindings are
MIT-licensed, while its binary wheel bundles LGPL components; redistributed
binary packages must preserve the upstream license notices.

## Feasibility and theory checkpoints

The initial review found that the original broad "first FO-EKF cardiac digital twin" claim is not defensible, but a narrower correlation-consistent, identifiability-gated finite-memory observer is conditionally feasible.

The subsequent control-theory review narrows the primary paper further: a periodicity-consistent mixed integer-fractional ECG surrogate, a gauge-invariant multi-heart-rate certificate for fractional-order identifiability, and a finite-record set-membership uncertainty chain. Approximation-aware observer-error tubes remain an earlier exploration; observer stability under fractional-operator mismatch is deferred and is not a proved contribution of the current theory core. Generic FO-EKF, finite-memory compensation, and LMI stability remain prior-art tools rather than headline contributions.

- [Feasibility and novelty review](research/feasibility_review_2026-07-11.md)
- [Control-theory innovation specification](research/control_theory_innovation_spec_2026-07-11.md)
- [S0-S3 theory-validation execution report](research/s0_s3_execution_report_2026-07-11.md)
- [R0-R4 identifiability and novelty gate](research/r0_r4_execution_report_2026-07-11.md)
- [Finite-record theory core](paper/fo_ecg_theory_core.tex)
- [R2 preregistration protocol](research/r2_protocol_spec_2026-07-11.md)
- [R2 machine-readable template](config/r2_protocol_template.toml)
- [R1 certification and R2 protocol checkpoint](research/r1_cert_r2_protocol_execution_2026-07-11.md)
- [Environment snapshot](research/environment_snapshot_2026-07-11.md)
- [Curated literature manifest](literature/manifest.csv)
- [Verified local literature lock](literature/download-lock.json)

The two-rate quotient theorem and the worst-case ambiguity theorem have passed
their data-free gates. The R1 ball-arithmetic core certifies every issued inner
or outer box, but a projected order set is exact only when no `UNKNOWN` leaf
remains. The R2 schema and propagation validator now fail closed on missing,
leaking, or downward-rounded evidence. Real-ECG order certification nevertheless
remains blocked until every data-to-disk term has an independent calibrated
bound. Missing bounds return `NOT_CERTIFIABLE`; they are never filled from the
target model residual.

## Research artifacts

- `PaperPlan/`: original research concept.
- `literature/`: reproducible literature manifest and local download script.
- `research/`: feasibility, novelty, and experimental-design decisions.
- `src/fo_ekf/`: implementation package (after method freeze).
- `tests/`: tests that do not require the raw datasets.
