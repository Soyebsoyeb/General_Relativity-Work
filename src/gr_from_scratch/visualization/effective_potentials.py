"""Effective-potential plots for Schwarzschild geodesics."""

from __future__ import annotations

from .style import apply_style


def v_timelike(r, L, M=1.0):
    return (1 - 2 * M / r) * (1 + L**2 / r**2)


def v_null(r, M=1.0):
    """V/L^2 for photons."""
    return (1 - 2 * M / r) / r**2


def plot_timelike(ax, r, L_values, M=1.0):
    apply_style()
    for L in L_values:
        ax.plot(r, v_timelike(r, L, M), label=rf"$L={L:g}M$")
    ax.set_xlabel(r"$r/M$")
    ax.set_ylabel(r"$V_{\rm eff}=E^2-\dot r^2$")
    ax.axvline(2 * M, color="0.5", lw=0.6)
    ax.legend()
    return ax
