# Phase 9 Equation Workspace

This directory holds the curated ridge-based concentration equations and supporting assets. The canonical path is `analysis/phase9/equations/`; a legacy copy (`phase9/equations/`) has been archived under `archive/legacy_phase9_equations/` and should be removed in an upcoming writable session to avoid drift.

## Contents
- `model_inventory.json` – curated coefficients/intercepts plus metrics and reporting-role flags.
- `equation_reference.md` – human-readable tables for each equation.
- `evaluate_equations.py` – sanity-check tool that recomputes predictions using the inventory and compares against canonical data.
- `feature_subset_sweep.py` – LOOCV ablation runner for predefined single- and two-feature subsets.
- `validation_harness.py` – replay utility that regenerates per-dose predictions/residuals for curated equations.
- `PROGRESS.md` – ongoing build log documenting decisions and work performed in this directory.
- `PLAN.md` – high-level plan captured during the initial scoping conversation.
- `outputs/` – stores generated sweep outputs (summary + per-fold CSVs).

## Next Steps
1. Remove the duplicate `phase9/` directory once auditing is complete (its contents already exist under `archive/legacy_phase9_equations/`).
2. Run `python analysis/phase9/equations/evaluate_equations.py` locally to sanity-check the equations after any data refresh.
3. Use `PROGRESS.md` to log all future changes/jobs in this workspace for auditability.
