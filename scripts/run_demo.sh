#!/usr/bin/env bash
set -euo pipefail
python -m imas_pipeline.demo_netcdf_ingest --input data/sample/sample.csv --out demo_output.png
echo "Demo complete → see demo_output.png"
