# Phase 8 Report Workspace

Pulls together the Phase 7 summaries and supporting tables into a comprehensive write-up.

## Contents
- `inputs/phase8_input_manifest.json` &ndash; references Phase 7 manifest and staged summary tables (modality/kind summaries, baseline metrics).
- `tables/` &ndash; houses aggregated tables created during the write-up (e.g., `comparison_overview.csv`).
- `plots/` &ndash; exported figures used in the report.
- `report/phase8_overview.ipynb` &ndash; notebook skeleton that loads metrics and houses the narrative.

## Workflow
1. Load Phase 7 artefacts using the manifest.
2. Build any additional aggregated tables (e.g., merged dependence vs. baseline metrics). `analysis/phase8/build_report_tables.py` currently generates `tables/comparison_overview.csv` and updates the manifest.
3. Generate report assets with `analysis/phase8/generate_report_assets.py` (creates plots and `report/summary.md`, updates manifest entries).
4. Draft the narrative/plots in the notebook; export to PDF/HTML if needed (place outputs in `report/`).
4. If Phase 7 data changes, rerun `analysis/phase7/build_summary.py` before refreshing the report.
