"""Monitors for the constraints a geodesic must satisfy."""

from __future__ import annotations

import numpy as np


def max_abs(values) -> float:
    return float(np.max(np.abs(values)))


def relative_drift(values) -> float:
    """max |q(l) - q(0)| / max(|q(0)|, 1)."""
    v = np.asarray(values, float)
    return float(np.max(np.abs(v - v[0])) / max(abs(v[0]), 1.0))
