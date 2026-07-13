# Paper artifacts

The submission-oriented manuscript is
[fo_ecg_certified_identifiability.tex](fo_ecg_certified_identifiability.tex).
The longer [fo_ecg_theory_core.tex](fo_ecg_theory_core.tex) is the technical
supplement containing the complete theorem and certificate derivations.

Build from this directory with:

    latexmk -pdf -interaction=nonstopmode -halt-on-error fo_ecg_certified_identifiability.tex
    latexmk -pdf -interaction=nonstopmode -halt-on-error fo_ecg_theory_core.tex

Regenerate the method and data-derived figures from the repository root with:

    .\.venv\Scripts\python.exe scripts\make_paper_figures.py

The generator creates the structural-identifiability geometry, the
point-estimate/certified-set workflow comparison, the core-method evidence
figure, and the Fantasia pilot summary.  Data-derived panels read only the
sealed exact-data bundles and versioned legacy results under
`research/results`; raw ECG remains outside the repository.

The manuscript intentionally preserves two negative conclusions:

- the four-record pilot does not support a distinct fractional order; all four
  fitted orders terminate at the nested alpha = 1 boundary;
- the real records remain NOT_CERTIFIABLE because independent acquisition,
  annotation, memory, anti-alias, and model-discrepancy bounds are incomplete.

These statements must not be replaced by a cardiac-digital-twin or
physiological-order claim.
