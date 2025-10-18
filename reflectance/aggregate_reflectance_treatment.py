#!/usr/bin/env python3
"""
Retired dose-level reflectance aggregation script.

Dose-level outputs are now built as part of the canonical dataset
(`python build_canonical_dataset.py`).  The legacy implementation is preserved at
`archive/pipeline/aggregate_reflectance_treatment.py` for audit purposes only.
"""

import sys


def main() -> None:
    sys.stderr.write(
        "aggregate_reflectance_treatment.py is retired.\n"
        "Use the canonical dataset in canonical_dataset/ for current analyses, "
        "or run archive/pipeline/aggregate_reflectance_treatment.py if you must "
        "recreate the legacy artefacts.\n"
    )
    sys.exit(1)


if __name__ == "__main__":
    main()
