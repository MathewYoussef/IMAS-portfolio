# Restructure Proposal — Agent 1

## Current Snapshot
- Repository is reflectance-heavy: `reflectance/` weighs ~108 MB and carries all NumPy staging data, analyses, and canonical tables; other thesis blocks are absent.
- Lightweight shared tooling remains at root (`imas_pipeline`, `scripts`, `environment.yml`), and a standalone GA smoother lives under `projects/`.
- Duplicate staging arrays (9 032 hash groups ≈ 26 MB) exist between `reflectance/denoised_full_run_staging` and `reflectance/from_dose_concentration_chromo_dad_work_upstream/denoised_full_run_staging`.
- `.DS_Store` slipped into version control; no `.gitignore` guards it.

## Target Top-Level Layout
```text
hub/
  README.md (portfolio overview)
  config/
    environment.yml
  packages/
    imas_pipeline/
  scripts/
    run_catalogue.sh
    run_demo.sh
  data/ (shared tiny fixtures only)
  docs/
    release_notes.md (future)
  ops/ (Agent deliverables stay under ops/)
Reflectance/
  analysis/phase{3..9}/
  data/
    raw/
      manifests/
    processed/
      denoised_full_run_staging/
    final/
      canonical_dataset/
  docs/
    AGENTS.md
    STATUS.md
  tests/
archive/
  Reflectance/
    legacy_equations/
    upstream_drops/
Initial_Calibration/           (placeholder for incoming block)
Act_of_God_Mamba_Results/      (placeholder; expect notebooks + plots)
Supplements/                   (placeholder; move supplement docs here)
subprojects/
  ga_smoother/
ops/ (operations workspace — unchanged)
```

## Key Moves & Justification
- **Promote `reflectance/` → `Reflectance/`** to align with thesis block naming and prepare for sibling blocks. Internal paths follow normalised tiers (`analysis/`, `data/raw|processed|final/`, `docs/`, `tests/`).
- **Deduplicate staging arrays** by keeping the canonical copy under `Reflectance/data/processed/denoised_full_run_staging/` and archiving the redundant upstream export. Introduce manifests linking back to original acquisition metadata before pruning.
- **Archive older artefacts**: move `reflectance/archive/` and the upstream DAD drop under `archive/Reflectance/` to reduce cognitive load while maintaining provenance.
- **Establish a `hub/`** to host shared tooling (`imas_pipeline`), configuration, scripts, and portfolio-level docs. This reduces noise at root and sets expectations for future shared assets.
- **Extract `projects/ga_smoother`** to `subprojects/ga_smoother/` (or external repo) because it has its own install flow, demo data, and is reusable beyond the reflectance block.
- **Reserve block stubs** for `Initial_Calibration`, `Act_of_God_Mamba_Results`, and `Supplements` so later agents can drop their assets without further root churn.

## Risks & Mitigations
- **Path-breaking changes**: Python modules and manifests reference the old `reflectance/` paths. Mitigate by staging a refactor PR with search-and-replace plus regression tests (`Reflectance/tests/verify_canonical_dataset.py`).
- **Deduplication safety**: Removing upstream copies before manifests are rewritten could lose acquisition traceability. Mitigate by snapshotting metadata to `/archive/Reflectance/upstream_drops/manifests/` and validating row counts vs. canonical tables.
- **Subproject extraction**: `ga_smoother` dependencies diverge from the main env. Coordinate with its maintainer before moving to `subprojects/` or a standalone repo; consider packaging it (pip installable) for reuse.
- **Missing blocks**: The target tree reserves folders that are currently empty. Follow-up tasks must bring in their assets (from other repos or file shares) and stub README/manifest files so the structure is meaningful.

## Next Steps for Follow-on Agents
1. Update `.gitignore` to cover `.DS_Store`, `ops/output/`, and virtualenv artefacts.
2. Script the reflectance path migration and re-run the verification suite to prove no regressions.
3. Design manifests for raw vs processed data to backfill the `Reflectance/data/raw/` tier (even if raw currently lives elsewhere).
4. Coordinate with stakeholders to ingest Initial Calibration, Mamba results, and Supplements content into their reserved blocks following the same tiering pattern.
5. Decide whether `ga_smoother` graduates to an external repository or remains as a `subprojects/` module with its own CI.
