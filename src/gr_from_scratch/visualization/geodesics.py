"""Orbit and ray plots in the equatorial plane."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .style import apply_style


def plot_equatorial_paths(ax, trajectories, labels=None, horizon: float | None = 2.0, **kw):
    """Plot paths given as arrays with columns (t, r, theta, phi) in the x-y plane."""
    apply_style()
    for i, tr in enumerate(trajectories):
        r, ph = tr[:, 1], tr[:, 3]
        ax.plot(r * np.cos(ph), r * np.sin(ph), label=None if labels is None else labels[i], **kw)
    if horizon:
        ax.add_patch(plt.Circle((0, 0), horizon, color="k", zorder=5))
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x/M$")
    ax.set_ylabel(r"$y/M$")
    return ax
