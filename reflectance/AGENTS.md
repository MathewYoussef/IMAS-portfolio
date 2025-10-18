# Repository Guidelines

## Project Structure & Module Organization
Root-level Python scripts build and audit the reflectance dataset. Active tooling such as `build_canonical_dataset.py`, `compute_precision_weighted_concentrations.py`, and `ensure_physical_dose_columns.py` operate on archived CSVs. Canonical tables live in `canonical_dataset/`, with checksums recorded in `canonical_dataset/manifest.json`. Legacy pipelines remain under `archive/`, and exploratory notebooks and exports are grouped by phase in `analysis/`. Tests reside in `tests/`.

## Build, Test, and Development Commands
Create a Python 3.10+ virtual environment (`python3 -m venv .venv && source .venv/bin/activate`) and install the required libraries manually (`pip install pandas numpy`). Rebuild the canonical bundle with `python build_canonical_dataset.py --output canonical_dataset`; the command validates source hashes before overwriting outputs. Refresh precision-weighted summaries via `python compute_precision_weighted_concentrations.py`. Use `python ensure_physical_dose_columns.py` if you extend dose metadata columns.

## Coding Style & Naming Conventions
Follow PEP 8 with 4-space indentation, snake_case modules and functions, and short docstrings that describe the data contract. Maintain type hints on public functions (see `build_canonical_dataset.py`) and prefer pure, idempotent helpers. CLI entry points should fail fast with explicit stderr guidance when a workflow is deprecated or misused.

## Testing Guidelines
Run `python tests/verify_canonical_dataset.py` after any data or script change; it compares each canonical artefact against `archive/aggregated_reflectance/` under a tight float tolerance. Extend the script when introducing new tables by defining deterministic sort keys and allowed schema deltas. Name additional tests `test_*` and keep fixtures self-contained to avoid polluting the canonical snapshots.

## Commit & Pull Request Guidelines
Write commit subjects in the imperative mood with ≤72 characters, optionally followed by a short body that notes regenerated artefacts (e.g., “Rebuild canonical dataset after dose fix”). Reference issue IDs or lab notebook pages when available. Pull requests should summarise the motivation, list the exact commands executed, attach verification output, and call out any canonical dataset diffs or manifest updates.

## Data Integrity & Security
Treat `archive/aggregated_reflectance/` as the immutable audit source. Regenerate rather than edit files in place, and capture provenance for any archival refresh. When canonical outputs change, update `canonical_dataset/manifest.json` and double-check downstream analyses for copies of stale data. Keep sensitive lab identifiers out of commits unless stored in encrypted archives.
