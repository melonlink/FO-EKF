# Data attribution and redistribution boundary

The source ECG data are not part of this repository. The four-record
engineering pilot uses the **Fantasia Database, version 1.0.0**, obtained
through PhysioNet:

- dataset DOI: https://doi.org/10.13026/C2RG61
- authoritative landing page: https://physionet.org/content/fantasia/1.0.0/
- file license: Open Data Commons Attribution License v1.0

Required scholarly attribution includes both:

1. Iyengar N, Peng C-K, Morin R, Goldberger AL, Lipsitz LA. “Age-related
   alterations in the fractal scaling of cardiac interbeat interval
   dynamics.” *American Journal of Physiology*, 271 (1996), R1078–R1084.
2. Goldberger AL et al. “PhysioBank, PhysioToolkit, and PhysioNet: Components
   of a new research resource for complex physiologic signals.”
   *Circulation*, 101 (2000), e215–e220.

The tracked pilot bundle contains no waveform samples. It does contain source
hashes, fitted coefficients, and selected beat-sample indices derived from
the Fantasia annotations inside serialized WLS execution records. Those
derived materials do not remove the source-dataset attribution and license
boundary. The repository's MIT license applies to original code and documents;
it does not replace the terms applicable to the source database or its
derived content.

Do not commit `.dat`, `.hea`, `.ecg`, or other raw WFDB source files. Keep the
external dataset read-only and resolve it through `ECG_DATA_DIR`.
