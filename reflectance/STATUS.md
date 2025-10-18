# Status

## Chromatogram & DAD Snapshot
- Phase-2 concentration tables from `Initial_Calibration/stats_for_paper/robust_concentration_vs_dose_thesis_phase_2_results` were duplicated into this repository on 2025-10-16T10:25:49Z.
- `canonical_dataset/manifest.json` records the source directories and SHA-256 checksums for every frozen CSV/JSON in `canonical_dataset/`.
- Chromatogram/DAD processing code is now retired; reflectance analyses should use the canonical tables exclusively.
- Run `python tests/verify_canonical_dataset.py` after any changes to confirm the frozen data remain byte-identical to the archived snapshot.
