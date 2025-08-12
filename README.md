# IMAS Portfolio

Compact, reproducible Python examples for IMAS/IMOS-style data work:
- Data discovery, cataloguing, preprocessing
- NetCDF/xarray time series handling
- Simple ML demos (CNN/GA) with toy data

## Quick start (conda)
```bash
conda env create -f environment.yml
conda activate imas
python -m imas_pipeline.demo_netcdf_ingest --help
IMAS-portfolio/
  imas_pipeline/        # importable package for scripts
    __init__.py
    demo_netcdf_ingest.py
    utils_io.py
  data/
    sample/             # tiny, anonymized examples (no large files)
  scripts/
    run_demo.sh
Check it exists:
```bash
ls -l README.md
cat > .gitignore << 'EOF'
__pycache__/
*.pyc
# editors
.vscode/
.idea/
# data (keep tiny samples only)
data/raw/
data/interim/
*.nc
*.h5
*.tif
# env / secrets / OS
.env
.DS_Store
