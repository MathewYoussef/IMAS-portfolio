# Scytonemin Thesis Audit Hub

This site is the reviewer-facing hub for the Scytonemin thesis audit, scoped to this repository. It explains how evidence is organized, where generated artifacts will land, and how to verify claims.

## How to Audit in 10 Minutes

1. Review Agent 1 reports in `ops/output/inventory/` and `ops/output/proposals/` to understand current layout and proposed moves.
2. Jump to a block page (Reflectance, Initial Calibration, Mamba‑SSM, Supplements) to see expected figures/tables and their `scaffold/<block>/` destinations.
3. Check the Claims section to trace thesis statements to notebooks, datasets, and tests (to be wired by Agents 3/4).

## Responsibilities & Handoffs

- Agent 2 (Scaffold) — This docs site, `mkdocs.yml`, `Makefile`, and `scaffold/**`.
- Agent 3 (Data & Analysis) — Catalogs data, publishes processed tables and generated figures into scaffold.
- Agent 4 (Media) — Produces web-weight media and page embeddings.
- Manuscript source — `thesis.docx.md` (one level up) provides authoritative prose excerpts while drafting.
