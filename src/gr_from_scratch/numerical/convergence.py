"""Convergence-order estimation."""

from __future__ import annotations

import numpy as np


def observed_order(step_sizes, errors) -> float:
    """Least-squares slope of log(error) against log(h)."""
    h, e = np.asarray(step_sizes, float), np.asarray(errors, float)
    mask = e > 0
    return float(np.polyfit(np.log(h[mask]), np.log(e[mask]), 1)[0])


def richardson_order(e_coarse: float, e_fine: float, ratio: float = 2.0) -> float:
    """Pairwise order estimate log(e_coarse / e_fine) / log(ratio)."""
    return float(np.log(e_coarse / e_fine) / np.log(ratio))
