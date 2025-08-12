#!/usr/bin/env bash
set -euo pipefail
python -m imas_pipeline.catalogue --root data/sample --out catalog.csv
