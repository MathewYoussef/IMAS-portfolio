#!/usr/bin/env python3
"""
Retired helper for back-filling UVA/UVB columns.

All active tables already store physical dose metadata; the legacy fixer lives
under archive/pipeline/ensure_physical_dose_columns.py for completeness.
"""

import sys


def main() -> None:
    sys.stderr.write(
        "ensure_physical_dose_columns.py is retired.\n"
        "No action required—canonical_dataset/ already contains UVA/UVB columns.\n"
        "If you truly need the legacy routine, run the archived copy in "
        "archive/pipeline/ensure_physical_dose_columns.py.\n"
    )
    sys.exit(1)


if __name__ == "__main__":
    main()
