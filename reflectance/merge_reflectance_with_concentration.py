#!/usr/bin/env python3
"""
Legacy reflectance↔concentration crosswalk builder.

Sample-level joins are no longer supported because destructive sampling prevents
perfect alignment.  Treatment-level alignment is documented in
canonical_dataset/dose_level_canonical_summary.csv.  The original script is
preserved in archive/tools/merge_reflectance_with_concentration.py if you need
to audit historical mappings.
"""

import sys


def main() -> None:
    sys.stderr.write(
        "merge_reflectance_with_concentration.py is retired.\n"
        "Refer to canonical_dataset/ for harmonised dose-level tables or use the archived "
        "script in archive/tools/ if you must rebuild the legacy crosswalk.\n"
    )
    sys.exit(1)


if __name__ == "__main__":
    main()
