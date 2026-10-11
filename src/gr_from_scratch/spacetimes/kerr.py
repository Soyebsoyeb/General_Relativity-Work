"""Kerr spacetime in Boyer-Lindquist coordinates.

The metric and its closed-form inverse are *entered*. Whether they solve the
vacuum equations is checked numerically at sample points in
``tests/test_kerr.py``; no symbolic proof is attempted here.
"""

from __future__ import annotations

import sympy as sp

from ..core.metric import Metric
from ..core.simplification import identity
from .schwarzschild import M, phi, r, t, theta

a = sp.Symbol("a", real=True)


def sigma(mass=M, spin=a):
    return r**2 + spin**2 * sp.cos(theta) ** 2


def delta(mass=M, spin=a):
    return r**2 - 2 * mass * r + spin**2


def kerr(mass=M, spin=a, symbolic_simplification: bool = False) -> Metric:
    S, D = sigma(mass, spin), delta(mass, spin)
    s2 = sp.sin(theta) ** 2
    A = (r**2 + spin**2) ** 2 - spin**2 * D * s2
    g = sp.Matrix(
        [
            [-(1 - 2 * mass * r / S), 0, 0, -2 * mass * r * spin * s2 / S],
            [0, S / D, 0, 0],
            [0, 0, S, 0],
            [-2 * mass * r * spin * s2 / S, 0, 0, A * s2 / S],
        ]
    )
    ginv = sp.Matrix(
        [
            [-A / (S * D), 0, 0, -2 * mass * r * spin / (S * D)],
            [0, D / S, 0, 0],
            [0, 0, 1 / S, 0],
            [-2 * mass * r * spin / (S * D), 0, 0, (D - spin**2 * s2) / (S * D * s2)],
        ]
    )
    m = Metric([t, r, theta, phi], g, inverse=ginv, name="Kerr")
    if not symbolic_simplification:
        # Skip per-component simplification: expressions are evaluated at points.
        m.simplifier = identity
    return m


def horizons(mass=M, spin=a):
    d = sp.sqrt(mass**2 - spin**2)
    return mass - d, mass + d


def ergosurface(mass=M, spin=a, polar=theta):
    """Outer boundary of the ergoregion: g_tt = 0."""
    return mass + sp.sqrt(mass**2 - spin**2 * sp.cos(polar) ** 2)


def kretschmann_closed_form(mass=M, spin=a):
    """Known closed form for K = R_{abcd} R^{abcd} (used only as an external check)."""
    c2 = sp.cos(theta) ** 2
    S = r**2 + spin**2 * c2
    return 48 * mass**2 * (r**2 - spin**2 * c2) * (S**2 - 16 * r**2 * spin**2 * c2) / S**6


# --- closed-form equatorial circular-orbit data (Bardeen, Press & Teukolsky 1972) ---
# These are literature formulas. tests/test_kerr.py checks them against the
# geodesic equation built from this package's own Christoffel symbols.


def circular_orbit(radius: float, spin: float, mass: float = 1.0, prograde: bool = True):
    """Return (E, L, Omega) for the equatorial circular geodesic at ``radius``."""
    import numpy as np

    s = 1.0 if prograde else -1.0
    sm = np.sqrt(mass)
    den = radius**0.75 * np.sqrt(radius**1.5 - 3 * mass * radius**0.5 + 2 * s * spin * sm)
    energy = (radius**1.5 - 2 * mass * radius**0.5 + s * spin * sm) / den
    ang = s * sm * (radius**2 - 2 * s * spin * np.sqrt(mass * radius) + spin**2) / den
    omega = s * sm / (radius**1.5 + s * spin * sm)
    return float(energy), float(ang), float(omega)


def isco_radius(spin: float, mass: float = 1.0, prograde: bool = True) -> float:
    """ISCO radius in Boyer-Lindquist coordinates (|spin| <= mass)."""
    import numpy as np

    a = abs(spin) / mass
    z1 = 1 + (1 - a**2) ** (1 / 3) * ((1 + a) ** (1 / 3) + (1 - a) ** (1 / 3))
    z2 = np.sqrt(3 * a**2 + z1**2)
    sign = -1.0 if prograde else 1.0
    return float(mass * (3 + z2 + sign * np.sqrt((3 - z1) * (3 + z1 + 2 * z2))))
