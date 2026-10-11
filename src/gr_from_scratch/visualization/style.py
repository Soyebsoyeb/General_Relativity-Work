"""Matplotlib defaults: plain, print-oriented figures."""

from __future__ import annotations

import matplotlib as mpl
from cycler import cycler

PALETTE = ["#1b1b1b", "#c0392b", "#2c6fbb", "#e69f00", "#2e8b57", "#8e5ea2", "#7f7f7f"]


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "stix",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "legend.fontsize": 8,
            "legend.frameon": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "axes.spines.top": True,
            "axes.spines.right": True,
            "axes.prop_cycle": cycler(color=PALETTE),
            "lines.linewidth": 1.2,
            "figure.dpi": 120,
            "savefig.dpi": 220,
            "savefig.bbox": "tight",
            "axes.grid": False,
            "figure.constrained_layout.use": True,
        }
    )
