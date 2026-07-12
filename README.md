# FO-EKF

Control-theoretic identifiability and fail-closed finite-record certification
for a fractional surface-ECG harmonic surrogate.

This repository starts from the research plan in
[`PaperPlan/innovation_points_fractional_order_cardiac_digital_twin.md`](PaperPlan/innovation_points_fractional_order_cardiac_digital_twin.md).
The repository name is historical. The defensible paper contribution is not a
generic FO-EKF and not a cardiac digital twin: it is a common-gauge multi-rate
identifiability theorem, its no-go boundaries, and a replayable
samples-to-response-disk certificate chain. A design-locked real-data pilot is
retained even though it does not support a distinct fractional order.

## Branches

- `main`: reviewed baseline and final releases.
- `dev-codex`: research checkpoints and staged work.

## Data

Raw/shared data remain outside this repository under `../PUBLIC` and are
treated as read-only. Typical PowerShell overrides are:

```powershell
$env:ECG_DATA_DIR = (Resolve-Path ..\PUBLIC\fantasia)
$env:ECG_PROCESSED_DIR = "$PWD\data\processed"
$env:ECG_OUTPUT_DIR = "$PWD\artifacts"
```

The code reads process environment variables directly; `.env.example` is only
a field template and `.env` is not auto-loaded. Never commit raw WFDB files,
patient-derived caches, or local paths. See [DATA_ATTRIBUTION.md](DATA_ATTRIBUTION.md)
for the Fantasia v1.0.0 DOI, ODC Attribution 1.0 boundary, and required
citations. No data or software download above 100,000,000 bytes may be started
without explicit user approval.

## Python environment

Python 3.11+ is required because the protocol loader uses the standard-library
`tomllib`. The checked environment is recorded in
`requirements-repro-lock.txt`. On Windows PowerShell:

```powershell
py -3.11 -m venv .venv
# Run installs only after checking the download size against the approval rule.
.\.venv\Scripts\python.exe -m pip install -r requirements-repro-lock.txt
.\.venv\Scripts\python.exe -m pip install --no-build-isolation -e . --no-deps
.\.venv\Scripts\python.exe -m flint.test --quiet
.\.venv\Scripts\python.exe scripts\check_environment.py
.\.venv\Scripts\python.exe -m pytest -q
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
- [T1 analytic-selector tube checkpoint](research/t1_analytic_selector_tube_execution_2026-07-12.md)
- [T2 theory-closure checkpoint](research/t2_theory_closure_execution_2026-07-12.md)
- [T3 end-to-end evidence-bundle checkpoint](research/t3_end_to_end_bundle_execution_2026-07-12.md)
- [T4 design-locked Fantasia pilot](research/t4_fantasia_pilot_execution_2026-07-12.md)
- [T5 exact-sample synthetic closure](research/t5_exact_sample_synthetic_chain_execution_2026-07-12.md)
- [T6 two-rate/multi-rate inverse audit](research/t6_two_rate_recovery_execution_2026-07-12.md)
- [Submission manuscript](paper/fo_ecg_certified_identifiability.tex)
- [Anonymous technical supplement](paper/fo_ecg_theory_core.tex)
- [Reproduction instructions](REPRODUCIBILITY.md)
- [Environment snapshot](research/environment_snapshot_2026-07-11.md)
- [Curated literature manifest](literature/manifest.csv)
- [Verified local literature lock](literature/download-lock.json)

The two-rate quotient theorem and the worst-case ambiguity theorem have passed
their data-free gates. The R1 ball-arithmetic core certifies every issued inner
or outer box, but a projected order set is exact only when no `UNKNOWN` leaf
remains. The T1 analytic Gram selector now adds replayable
`ROBUST_INNER / CLOSED_INNER / UNKNOWN` certificates for
`forall alpha, exists nuisance`, including a gap-free alpha-only subdivision
audit. Its exact rational inverse makes the state-direction Krawczyk
contraction term identically zero; all physical claims are still rechecked
against the original Arb disk, annulus, and damping inequalities. The R2
schema and propagation validator now fail closed on missing,
leaking, or downward-rounded evidence. Real-ECG order certification nevertheless
remains blocked until every data-to-disk term has an independent calibrated
bound. Missing bounds return `NOT_CERTIFIABLE`; they are never filled from the
target model residual.

The current T3/T5 chain uses R2 v4 exact integer sample geometry. It binds the
digital payload, deterministic sample indices and integer window endpoints,
joint complex WLS execution, an independent full exact-functional Arb/ACB
interval, the component bounds, the protocol result, and the downstream T1
selector witness. The full numerical radius is added once; any missing link,
operator mismatch, altered payload, or broken seal returns `NOT_CERTIFIABLE`.
Continuous-lock-in weights remain a separate exact-timestamp contract, and a
verified Hermitian LDL* path handles positive-definite Gramians that
Gershgorin cannot certify. This closes the conditional machine-composition
path; it does not manufacture the independent acquisition and model-error
bounds still missing on real ECG.

## Research artifacts

- `PaperPlan/`: original research concept.
- `literature/`: reproducible literature manifest and local download script.
- `research/`: feasibility, novelty, and experimental-design decisions.
- `src/fo_ekf/`: implementation package (after method freeze).
- `tests/`: tests that do not require the raw datasets.
