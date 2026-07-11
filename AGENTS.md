# FO-EKF repository guardrails

## Repository policy

- `main` contains reviewed, mature baselines and final versions.
- `dev-codex` contains coherent research checkpoints and staged integration work.
- Small exploratory changes should be committed locally before they are grouped into a remote checkpoint.
- Do not promote work from `dev-codex` to `main` until the research claim, implementation, tests, and evidence have been reviewed.

## Data boundary

- Shared ECG data are under `../PUBLIC` and are read-only for this project.
- Never copy raw WFDB signals or other large source datasets into this repository.
- Resolve input data through `ECG_DATA_DIR`; write derived data to `ECG_PROCESSED_DIR` and reports/results to `ECG_OUTPUT_DIR`.
- Patient-level splits and the official PTB-XL `strat_fold` must be preserved; external validation must exclude overlapping PTB/PTB-XL records.

## Evidence policy

- Do not describe a method as a cardiac digital twin until subject-specific synchronization and prospective or held-out validation are demonstrated.
- Do not claim observability, identifiability, stability, real-time performance, or physiological meaning without an explicit theorem/test and reproducible evidence.
- Literature PDFs may be downloaded locally from the manifest, but only metadata and legally redistributable small assets belong in Git.

