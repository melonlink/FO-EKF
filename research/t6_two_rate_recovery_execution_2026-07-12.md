# T6 deterministic two-rate and multi-rate inverse audit

- Date: 2026-07-12
- Branch target: `dev-codex`
- Production runner: `scripts/run_two_rate_recovery_sweep.py`
- Formal artifact directory: `research/results/two_rate_recovery_2026-07-12/`
- Bound production commit: `82f88c384c383e0cb3c3945e887fc5d9e0c18bcb`
- Formal artifact status: **sealed and replayed**
- Data/software downloads: none

## Verdict

The production two-rate global inverse completed all 1260 exact, noiseless cases without a
failure.  The same cases also passed the common-complex-gain invariance check and the exact
three-rate reference-pair consistency check.  A deterministic perturbation of only the third
response was detected in all 1260 cases, both by the geometric cross-ratio deviation and by the
reference-pair model check.

This is a floating-point implementation audit of the analytical exact-model results.  It does
not replace the global-inverse or consistency proofs, does not establish robustness to noise,
and is not evidence from real ECG or a cardiac digital twin.

## Frozen design

The runner uses a fixed Cartesian grid and no pseudo-random generator:

- 21 orders from 0.05 through 1.0, including both endpoints;
- five damping values from `1e-2` through `1e2` on a base-10 logarithmic grid;
- rate ratios `1.25`, `1.75`, `3`, and `6`;
- three fixed `(base-frequency, nonzero complex morphology)` condition pairs
  (paired by index, not crossed with each other);
- four nonzero common complex gains, cycled deterministically across the cases;
- geometric triples `(nu, r*nu, r^2*nu)` for the three-rate checks;
- a legal alternate pair with damping increased by 5%, morphology chosen to keep the first
  response exact, and only the third response replaced.

Every two-rate inversion uses the full declared order bracket `[1e-6, 1]`; no local initial
guess or true-parameter neighbourhood is supplied to the bisection.

## Final sealed result

| Audit | Result | Worst recorded error |
|---|---:|---:|
| Two-rate global inverse | 1260 / 1260 | order `4.07e-13`; damping `8.86e-10` |
| Recovered morphology | 1260 / 1260 | absolute error `2.19e-11` |
| Common-gain invariance | 1260 / 1260 | order difference `4.55e-13`; damping difference `1.09e-9` |
| Gained morphology covariance | 1260 / 1260 | absolute error `1.13e-10` |
| Three-rate positive consistency | 1260 / 1260 | order spread `3.98e-13`; damping spread `7.98e-10` |
| Geometric cross-ratio identity | 1260 / 1260 | absolute error `3.87e-12` |
| Third-rate tamper detection | 1260 / 1260 | minimum cross-ratio deviation `4.50e-4` |

Every tampered first/third pair remained a valid exact response pair on its own, but recovered
the deliberately altered damping.  Therefore all 1260 cases returned an explicit inconsistent
multi-rate audit; none relied on rejection outside the physical quotient image, and there were
no missed tamper cases.

## Determinism and frozen provenance

Each execution independently renders the 1260-case table twice and records that the two byte
strings match.  The formal runner additionally enforces a two-commit evidence freeze:

1. `pyproject.toml`, `requirements-repro-lock.txt`, the runner, and `multirate.py` form the
   production-source manifest.  Default generation refuses to run if any of them is dirty.
2. The first clean execution records the full production-code commit and seals `cases.csv`, the
   canonical summary, and `SHA256SUMS`.
3. After those results are committed separately, a default rerun validates the existing
   summary seal, artifact hash, complete `SHA256SUMS`, and source manifest before reusing the
   recorded production-code commit.  The result commit therefore does not rewrite its own
   provenance.
4. If any production-source hash changes, the runner refuses to overwrite the old evidence and
   requires it to be archived before a new clean generation.

The earlier files bound to the pre-production `bd989f1...` base commit were deliberately
removed and are not evidence.  The final artifact binds
`82f88c384c383e0cb3c3945e887fc5d9e0c18bcb`. Two consecutive post-freeze runs retained
byte-identical `cases.csv`, `summary.json`, and `SHA256SUMS`.

- internal summary seal: `7e58edc544419a4196571c4b538ca429d528bc506cfe34a6502551404be4ba75`;
- `cases.csv`: `e07ab743991daedaa0ae883d278efbcb43b507a5878e1585f80b72af1a7db834`;
- `summary.json`: `a208a0466dde1b138681b31eed8a022c12898f262609c2520210336e4e60791f`;
- `SHA256SUMS`: `6fc7cd9e3c875d1b660d706f3c51c0b1752a92cc092f47ae1288f37dd6abe6db`.

## Reproduction

From the repository root:

```powershell
git status --short -- pyproject.toml requirements-repro-lock.txt scripts/run_two_rate_recovery_sweep.py src/fo_ekf/multirate.py
.\.venv\Scripts\python.exe scripts\run_two_rate_recovery_sweep.py
.\.venv\Scripts\python.exe scripts\run_two_rate_recovery_sweep.py
.\.venv\Scripts\python.exe -m pytest -q tests\test_two_rate_recovery_sweep.py
.\.venv\Scripts\ruff.exe check scripts\run_two_rate_recovery_sweep.py tests\test_two_rate_recovery_sweep.py
.\.venv\Scripts\ruff.exe format --check scripts\run_two_rate_recovery_sweep.py tests\test_two_rate_recovery_sweep.py
```

The first status command must print nothing. A default rerun validates the existing evidence,
retains the recorded production commit in `code_commit`, and must produce byte-identical files.

The pytest suite checks the case count, grid endpoints, inverse-error ceilings, all positive
consistency outcomes, all tamper detections, and both artifact hashes.  It also simulates a
later result commit and verifies reuse of commit A, rejects dirty production sources, rejects a
changed source manifest, and fails closed on summary, CSV, and `SHA256SUMS` tampering.
