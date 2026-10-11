"""Geometric sections: horizons and ergosurfaces."""

from __future__ import annotations

import numpy as np

from .style import apply_style


def kerr_sections(ax, spin: float, mass: float = 1.0, color="k", label=None):
    """Meridional (x-z) section of the Kerr horizons (solid) and ergosurface (dashed).

    Boyer-Lindquist r is mapped to the Kerr-Schild-like embedding
    x = sqrt(r^2 + a^2) sin(theta), z = r cos(theta) for display only.
    """
    apply_style()
    th = np.linspace(0, np.pi, 400)
    rp = mass + np.sqrt(mass**2 - spin**2)
    rm = mass - np.sqrt(mass**2 - spin**2)
    rE = mass + np.sqrt(mass**2 - spin**2 * np.cos(th) ** 2)

    def xz(r, th):
        return np.sqrt(r**2 + spin**2) * np.sin(th), r * np.cos(th)

    for rad in (rp, rm):
        x, z = xz(rad * np.ones_like(th), th)
        ax.plot(x, z, color=color, lw=1.0, label=label if rad == rp else None)
        ax.plot(-x, z, color=color, lw=1.0)
    x, z = xz(rE, th)
    ax.plot(x, z, "--", color=color, lw=0.9)
    ax.plot(-x, z, "--", color=color, lw=0.9)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\sqrt{r^2+a^2}\sin\theta\ /M$")
    ax.set_ylabel(r"$r\cos\theta\ /M$")
    return ax
