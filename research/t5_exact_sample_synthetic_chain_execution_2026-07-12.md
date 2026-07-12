# T5 exact-sample synthetic certificate chain execution

- Date: 2026-07-12
- Branch target: `dev-codex`
- Runner: `scripts/run_synthetic_certificate_chain.py`
- Machine-readable result: `research/results/synthetic_certificate_chain_2026-07-12/summary.json`
- Replay artifacts: `artifacts.json` and `SHA256SUMS` in the same directory
- Bound production commit: `82f88c384c383e0cb3c3945e887fc5d9e0c18bcb`
- Data/software downloads: none

## Verdict

The deterministic synthetic chain closes end to end:

`8 signed-int64 samples at indices 80..87 -> replayed WLS execution -> CERTIFIED_INTERVAL -> R2 v4 exact-index component bound -> PASS_DETERMINISTIC protocol -> R2/T1 bundle -> replayed PRESERVED_ROBUST`.

This is the first repository artifact in which the response center is not merely declared. The selected digital payload is hash-bound and replayed, and the exact-functional interval certificate independently encloses conversion, phase/exponential evaluation, direct Gram/right-hand-side accumulation, verified joint solution, and front-end division. The bundle reports `mathematical_exact_sample_wls_coverage=true` and `uncovered_numeric_term=null`.

## Replayed stages

| Stage | Result |
|---|---|
| WLS execution | `MATCH` |
| Exact integer-payload interval | `MATCH / CERTIFIED_INTERVAL` |
| R2 v4 exact-index component certificate | valid `CERTIFIED_BOUND` |
| Frozen protocol | `PASS_DETERMINISTIC` |
| R2--T1 bundle | `PRESERVED_ROBUST` |
| Bundle replay | valid / `replay_verified` |

The window is 20--22 s with R-peak anchors at samples 80, 84, and 88 (`fs=4 Hz`). Both WLS and the R2 Gram proof use the same integer phase operator `2*pi*N*(j-a)/(b-a)`; binary64 seconds are retained only for time-domain error terms.

The bundle radius uses

`component sampling floor + full exact-functional numerical floor`.

The older solve-only WLS floor is not added again, because the full interval radius already contains Gram/right-hand-side accumulation, solve, and front-end division. The published numerical upper is a complex-modulus disk radius and is included once; no extra factor of `sqrt(2)` is applied.

The generator directly emits the retained-harmonic steady-state response for `alpha=0.7`, `lambda=0.5`, `q=1-0.2i`, and `omega=2*pi`. Consequently history/dwell transients, timestamp/interpolation/anti-alias errors, measurement noise, and model residual are zero by construction—not values estimated as zero from ECG. Only ADC quantization and the exact-functional numerical enclosure are nonzero.

The known response is approximately `0.0999000589 - 0.2436352405i`. The replayed center differs by about `6.51e-11`; the published disk radius is about `1.00e-9`, contains the known response, and has radius-to-response-magnitude ratio about `3.80e-9`. This makes the synthetic case an informative positive closure rather than only a loose composition test.

## Fail-closed controls

Two adversarial controls were executed through the same public bundle API:

1. Changing one external digital sample produced `NOT_CERTIFIABLE / wls_execution_record_not_replayable`.
2. Changing the published interval radius without a valid seal produced `NOT_CERTIFIABLE / wls_interval_certificate_not_replayable`.

The protocol validator now requires every `PASS_DETERMINISTIC` row to bind the execution, window, sample-index, digital, time, normalized-weight, and interval-certificate hashes, declare exact sample-functional coverage, and publish the target exact-functional numerical radius.

## Interpretation boundary

This execution proves implementation-level closure for one frozen synthetic integer-sample functional and its directly constructed R2/T1 bounds. It does not close acquisition, annotation, anti-alias, interpolation, front-end calibration, model-discrepancy, or physiological-validity errors on real ECG. Therefore it is not real-data R2 evidence and is not a cardiac digital-twin claim.

## Reproduction

From the repository root:

```powershell
.\.venv\Scripts\python.exe scripts\run_synthetic_certificate_chain.py
```

The runner refuses dirty production sources. Evidence records both the code commit and a stable production-source manifest hash. After the evidence files are committed, rerunning at the result commit reuses the recorded code commit only when the source manifest is unchanged, avoiding a result-commit self-reference. For byte-identical reproduction, check out the recorded code commit or retain identical production sources, then run the command twice; the two summary hashes must match.

Final seals:

- internal summary seal: `ad69630e12a7fde03e51874cc2f8a34961f8d699f86d572cf78231a40a92ca7d`;
- `summary.json`: `8806830dda55b13c5b92d9f88ef69feb314c301a570eea778879c35dc9e4f7f2`;
- `artifacts.json`: `8764174da8a2364d2219e03ae9eba6a41bdea6ed677663fc56f13a3501c102da`;
- `SHA256SUMS`: `5fb63cd1077b92c965101b22db135a13bb2b366c95861132a77c7c82824e48d6`.

Three consecutive post-freeze executions retained byte-identical
`artifacts.json`, `summary.json`, and `SHA256SUMS`.
