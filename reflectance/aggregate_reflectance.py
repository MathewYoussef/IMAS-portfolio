#!/usr/bin/env python3
"""
Retired entry point for the legacy reflectance aggregation pipeline.

The pipeline has been superseded by the canonical dataset build that works
directly from `archive/aggregated_reflectance/`.  This stub remains so older
invocations fail fast with a clear message instead of silently regenerating
outdated tables.
"""

import sys


def main() -> None:
    sys.stderr.write(
        "aggregate_reflectance.py has been retired.\n"
        "The raw outputs it produced are frozen under archive/aggregated_reflectance/, "
        "and all current analyses should consume canonical_dataset/.\n"
        "If you truly need the legacy behaviour, run the archived script at "
        "archive/pipeline/aggregate_reflectance.py.\n"
    )
    sys.exit(1)


if __name__ == "__main__":
    main()
