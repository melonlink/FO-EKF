# FO-EKF environment snapshot

## Repository

- Root: `D:\PROJECT\AI_PROJECT\ECG_Proj\FO-EKF`
- Remote: `https://github.com/melonlink/FO-EKF.git`
- Default branch: `main`
- Staged-research branch: `dev-codex`
- Bootstrap commit: `900fa66 chore: bootstrap FO-EKF research workspace`

Branch policy:

- reviewed baselines and final versions go to `main`;
- coherent research checkpoints go to `dev-codex`;
- small exploratory changes are committed locally before checkpoint publication.

## Python research environment

- Project environment: `.venv`
- Interpreter: Python 3.11.4
- Core verified packages:
  - NumPy 2.4.6
  - SciPy 1.17.1
  - pandas 2.3.3
  - Matplotlib 3.11.0
  - WFDB 4.3.1
  - pytest 8.4.2
  - Ruff 0.15.21
- Certification extra:
  - python-flint 0.9.0
  - Arb/ACB rigorous ball arithmetic
  - official self-test: 227 tests and 3,092 doctests passed
- Installation target: editable project with `.[dev,cert]`
- Current repository validation: 107 tests passed, including 21 R1-CERT
  tests and 53 R2 protocol/adversarial tests; full Ruff lint and the scoped
  format check for all files changed in this stage passed. The repository-wide
  format check still flags the pre-existing `src/fo_ekf/certificate.py`, which
  this checkpoint deliberately leaves untouched.
- Theory manuscript: `pdflatex` built 12 pages with no unresolved
  references, warnings, or overfull/underfull boxes.

The project environment is intentionally lightweight. It does not duplicate the existing PyTorch installation.
The certification wheel was 9.15 MiB; the earlier 0.8.0 wheel used during
version discovery was 10.4 MB. Both were below the 100 MB approval threshold.

## Existing PyTorch environment

- Interpreter: `D:\TOOLS\AI\anaconda3\envs\myPyTorch\python.exe`
- Python: 3.10.19
- PyTorch: 2.9.0+cu130
- CUDA available: yes
- CUDA runtime reported by PyTorch: 13.0
- Device: NVIDIA RTX A2000 8GB Laptop GPU

This environment is reserved for later deep baselines or differentiable optimization. The initial FO-EKF implementation should use NumPy/SciPy unless a PyTorch-specific experiment is justified.

## Data boundary

- Actual PTB-XL root: `D:\PROJECT\AI_PROJECT\ECG_Proj\PUBLIC\data_PTB-XL`
- Verified metadata: `raw\ptbxl_database.csv`
- Shared source data remain read-only and outside Git.
- Derived files default to `data/processed`; outputs default to `artifacts`.
- `ECG_DATA_DIR`, `ECG_PROCESSED_DIR`, and `ECG_OUTPUT_DIR` can override these paths.

## Download boundary

- Any individual dataset, software package, model, or file expected to exceed 100 MB requires explicit user approval before download.
- The literature downloader enforces a hard 100 MiB per-file limit and refuses incomplete PDFs.
- Current literature set: 15 PDFs, 73,553,418 bytes total; largest individual file 24,449,147 bytes.
