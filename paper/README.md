# Paper artifacts

The submission-oriented manuscript is
[fo_ecg_certified_identifiability.tex](fo_ecg_certified_identifiability.tex).
The longer [fo_ecg_theory_core.tex](fo_ecg_theory_core.tex) is the technical
supplement containing the complete theorem and certificate derivations.

Build from this directory with:

    latexmk -pdf -interaction=nonstopmode -halt-on-error fo_ecg_certified_identifiability.tex
    latexmk -pdf -interaction=nonstopmode -halt-on-error fo_ecg_theory_core.tex

Regenerate the data-derived figure from the repository root with:

    .\.venv\Scripts\python.exe scripts\make_paper_figures.py

The figure reads only the locked derived bundle under
research/results/fantasia_pilot_2026-07-12. Raw ECG remains outside the
repository.

The manuscript intentionally preserves two negative conclusions:

- the four-record pilot does not support a distinct fractional order; all four
  fitted orders terminate at the nested alpha = 1 boundary;
- the real records remain NOT_CERTIFIABLE because independent acquisition,
  annotation, memory, anti-alias, and model-discrepancy bounds are incomplete.

These statements must not be replaced by a cardiac-digital-twin or
physiological-order claim.
