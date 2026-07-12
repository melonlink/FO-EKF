# Reproducibility guide

This guide distinguishes three claims:

1. theorem and certificate-code checks, which need no ECG data;
2. the exact-sample synthetic certificate chain, which uses eight generated
   integer samples and no external data;
3. the Fantasia engineering pilot, which needs the four read-only source
   records but remains `NOT_CERTIFIABLE` for real-ECG order inference.

No command in this guide downloads data. Do not install software or fetch data
above 100,000,000 bytes without explicit user approval.

## Environment

The evidence release was checked on CPython 3.11.4. Exact package versions are
in `requirements-repro-lock.txt`. If the environment is not already present,
inspect the proposed download size before running an install.

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-repro-lock.txt
.\.venv\Scripts\python.exe -m pip install --no-build-isolation -e . --no-deps
.\.venv\Scripts\python.exe -m flint.test --quiet
```

## Data-free quality gate

Clear path overrides so tests cannot accidentally consume local ECG data:

```powershell
Remove-Item Env:ECG_DATA_DIR -ErrorAction SilentlyContinue
Remove-Item Env:ECG_OUTPUT_DIR -ErrorAction SilentlyContinue
Remove-Item Env:ECG_PROCESSED_DIR -ErrorAction SilentlyContinue
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\ruff.exe check src tests scripts
.\.venv\Scripts\ruff.exe format --check src tests scripts
git diff --check
```

## Exact-sample synthetic chain

The canonical directory contains `artifacts.json`, `summary.json`, and
`SHA256SUMS`. The runner verifies the existing seals, binds the same production
source manifest and code commit, recomputes the complete chain, and refuses
dirty production sources.

```powershell
.\.venv\Scripts\python.exe scripts\run_synthetic_certificate_chain.py
.\.venv\Scripts\python.exe scripts\run_synthetic_certificate_chain.py
```

The two runs must retain internal `summary_sha256`
`ad69630e12a7fde03e51874cc2f8a34961f8d699f86d572cf78231a40a92ca7d`
and a byte-identical three-file checksum set. A positive relation covers only the frozen synthetic
sample functional and declared synthetic bounds; it is not real-ECG evidence.

## Two-rate and multi-rate inverse sweep

```powershell
.\.venv\Scripts\python.exe scripts\run_two_rate_recovery_sweep.py
.\.venv\Scripts\python.exe scripts\run_two_rate_recovery_sweep.py
```

The sealed fixed-grid result contains 1,260 exact two-rate inversions, common
complex-gain checks, positive three-rate consistency checks, and independently
valid third-rate tamper controls. Both executions must preserve internal seal
`7e58edc544419a4196571c4b538ca429d528bc506cfe34a6502551404be4ba75`.

## Fantasia pilot

The source is Fantasia Database v1.0.0 (DOI `10.13026/C2RG61`). The generator
checks all 12 selected source files against the official v1.0.0 SHA-256 values.
Set the read-only source directory explicitly:

```powershell
$env:ECG_DATA_DIR = (Resolve-Path ..\PUBLIC\fantasia)
.\.venv\Scripts\python.exe scripts\replay_fantasia_wls.py `
  --bundle-dir research\results\fantasia_pilot_2026-07-12
```

The replay command is read-only and must report 1,859/1,859 matches. To test a
fresh regeneration without overwriting the tracked canonical bundle, use the
generator's default ignored output directory:

```powershell
.\.venv\Scripts\python.exe scripts\run_fantasia_pilot.py
.\.venv\Scripts\python.exe scripts\replay_fantasia_wls.py `
  --bundle-dir artifacts\fantasia_pilot_candidate
```

Compare the candidate `artifacts/fantasia_pilot_candidate/SHA256SUMS` with the
tracked bundle. Environment metadata and production-source hashes are embedded
in `design_lock.json` and `summary.json`. Numerical equality of the WLS path
does not supply the missing acquisition, annotation, anti-alias, memory,
front-end, harmonic-tail, or model-discrepancy bounds.

For the sealed bundle, SHA-256 of the `SHA256SUMS` file is
`c6d4018f84e90a2abe7130b384e9b0ca8c5ac7859bc04bf6b19d98964d80518a`.

## Figures and manuscripts

```powershell
.\.venv\Scripts\python.exe scripts\make_paper_figures.py
Push-Location paper
latexmk -pdf -interaction=nonstopmode -halt-on-error fo_ecg_certified_identifiability.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error fo_ecg_theory_core.tex
Pop-Location
```

The PDF release is accepted only after all pages are rasterized and visually
checked. The paper must continue to report `NOT_SUPPORTED` for a distinct
fractional-order advantage in the four-record pilot and `NOT_CERTIFIABLE` for
real-ECG set inference.
