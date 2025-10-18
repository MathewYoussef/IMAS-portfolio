#!/usr/bin/env python3
"""
Retired precision-weighted concentration fusion script.

Latent concentrations are already bundled in canonical_dataset/; this shim
exists so accidental invocations produce a clear error instead of recomputing
from scratch.  The historical implementation is preserved under
archive/pipeline/compute_precision_weighted_concentrations.py.
"""

import sys


def main() -> None:
    sys.stderr.write(
        "compute_precision_weighted_concentrations.py is retired.\n"
        "Use canonical_dataset/precision_weighted_concentrations*.csv for current work, "
        "or consult archive/pipeline/compute_precision_weighted_concentrations.py if you "
        "must rerun the legacy workflow.\n"
    )
    sys.exit(1)


if __name__ == "__main__":
    main()
