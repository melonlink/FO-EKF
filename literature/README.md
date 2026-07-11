# Literature workspace

`manifest.csv` is the source of truth for references used in the feasibility review. PDFs are downloaded to `papers/` only when an authoritative open-access or author-repository source is available. The PDFs remain local and are ignored by Git; this repository tracks citation and access metadata, while `download-lock.json` records the locally verified size and SHA-256.

Download the curated set from the repository root:

```powershell
.\.venv\Scripts\python.exe scripts\download_literature.py --dry-run
.\.venv\Scripts\python.exe scripts\download_literature.py
```

The helper refuses any individual response larger than 100 MiB. Larger data, models, or software require explicit user approval before download. A manifest entry with `download=0` records metadata only; it must not be fetched from an unofficial mirror.
