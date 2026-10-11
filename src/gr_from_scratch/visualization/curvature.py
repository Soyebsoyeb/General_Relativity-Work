"""Curvature-invariant plots."""

from __future__ import annotations

import numpy as np

from .style import apply_style


def plot_kretschmann_radial(ax, r, k_values, analytic=None, **kw):
    apply_style()
    ax.loglog(r, k_values, **kw)
    if analytic is not None:
        ax.loglog(r, analytic, "--", color="0.4", lw=0.8)
    ax.set_xlabel(r"$r/M$")
    ax.set_ylabel(r"$K\,M^4$")
    return ax


def loglog_slope(x, y) -> float:
    return float(np.polyfit(np.log(x), np.log(y), 1)[0])
