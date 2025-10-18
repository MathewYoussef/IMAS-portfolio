#!/usr/bin/env python3
"""
Retired DAD dose-level aggregation script.

The canonical dataset already includes the trimmed DAD summaries; this stub is
kept only to maintain backwards compatibility for tooling that imports the old
module name.
"""

import sys


def main() -> None:
    sys.stderr.write(
        "aggregate_dad_treatment.py is retired.\n"
        "Canonical outputs live in canonical_dataset/; the legacy script is archived at "
        "archive/pipeline/aggregate_dad_treatment.py for historical reproduction.\n"
    )
    sys.exit(1)


if __name__ == "__main__":
    main()
