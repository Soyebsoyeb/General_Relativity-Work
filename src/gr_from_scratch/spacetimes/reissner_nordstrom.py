"""Reissner-Nordstrom spacetime and its electromagnetic potential.

Units: G = c = 1 and Gaussian electromagnetic units, so that
G_{mu nu} = 8 pi T_{mu nu} with T_{mu nu} = (1/4 pi)(F_{mu a} F_nu^a - (1/4) g F^2).
The charge Q is then measured in units of length (sqrt(G) M scale).
"""

from __future__ import annotations

import sympy as sp

from ..core.metric import Metric
from .schwarzschild import M, phi, r, t, theta

Q = sp.Symbol("Q", positive=True)


def lapse(mass=M, charge=Q):
    return 1 - 2 * mass / r + charge**2 / r**2


def reissner_nordstrom(mass=M, charge=Q) -> Metric:
    ff = lapse(mass, charge)
    return Metric(
        [t, r, theta, phi],
        sp.diag(-ff, 1 / ff, r**2, r**2 * sp.sin(theta) ** 2),
        name="Reissner-Nordstrom",
    )


def vector_potential(charge=Q):
    """A_mu = (-Q/r, 0, 0, 0)."""
    return [-charge / r, 0, 0, 0]


def horizons(mass=M, charge=Q):
    d = sp.sqrt(mass**2 - charge**2)
    return mass - d, mass + d
