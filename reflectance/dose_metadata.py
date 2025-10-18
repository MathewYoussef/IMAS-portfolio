#!/usr/bin/env python3
"""
Centralised dose metadata used across reflectance and concentration pipelines.

Each dose label (``dose_1``–``dose_6``) is mapped to its measured UVA and UVB
irradiance.  Downstream tables should always carry these physical columns so
joins are based on explicit intensities instead of ambiguous ordinals.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple


@dataclass(frozen=True)
class DoseRecord:
    label: str
    uva_mw_cm2: float
    uvb_mw_cm2: float

    @property
    def as_tuple(self) -> Tuple[float, float]:
        return (self.uva_mw_cm2, self.uvb_mw_cm2)


_DOSE_SEQUENCE: List[DoseRecord] = [
    DoseRecord("dose_1", 0.0, 0.0),
    DoseRecord("dose_2", 0.647, 0.246),
    DoseRecord("dose_3", 1.095, 0.338),
    DoseRecord("dose_4", 1.692, 0.584),
    DoseRecord("dose_5", 2.488, 0.768),
    DoseRecord("dose_6", 3.185, 0.707),
]


def iter_dose_records() -> Iterable[DoseRecord]:
    """Yield the canonical dose records in ascending UVA order."""
    return iter(_DOSE_SEQUENCE)


def dose_mapping() -> Dict[str, DoseRecord]:
    """Return a mapping from dose label to the corresponding metadata record."""
    return {record.label: record for record in _DOSE_SEQUENCE}


def attach_dose_metadata(label: str) -> DoseRecord:
    """Return the dose metadata for ``label`` or raise if the key is unknown."""
    mapping = dose_mapping()
    if label not in mapping:
        raise KeyError(f"Unknown dose label: {label!r}")
    return mapping[label]
