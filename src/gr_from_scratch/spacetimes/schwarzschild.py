"""Schwarzschild spacetime in Schwarzschild and ingoing Eddington-Finkelstein coordinates.

``schwarzschild()`` merely *enters* the metric. That it solves the vacuum
equations is established by ``field_equations.vacuum_equations``, which both
verifies the entered form and re-derives it from the static spherical ansatz.
"""

from __future__ import annotations

import sympy as sp

from ..core.metric import Metric

t = sp.Symbol("t", real=True)
r = sp.Symbol("r", positive=True)
theta = sp.Symbol("theta", positive=True)
phi = sp.Symbol("phi", real=True)
M = sp.Symbol("M", positive=True)
v = sp.Symbol("v", real=True)


def f(mass=M, radius=r):
    return 1 - 2 * mass / radius


def schwarzschild(mass=M) -> Metric:
    ff = f(mass)
    return Metric(
        [t, r, theta, phi],
        sp.diag(-ff, 1 / ff, r**2, r**2 * sp.sin(theta) ** 2),
        name="Schwarzschild",
    )


def tortoise(radius=r, mass=M):
    """r* = r + 2M ln(r/2M - 1), valid for r > 2M."""
    return radius + 2 * mass * sp.log(radius / (2 * mass) - 1)


def eddington_finkelstein(mass=M) -> Metric:
    """Ingoing coordinates v = t + r*: ds^2 = -f dv^2 + 2 dv dr + r^2 dOmega^2."""
    s = schwarzschild(mass)
    return s.transform([v, r, theta, phi], [v - tortoise(r, mass), r, theta, phi])


def horizon_radius(mass=M):
    return 2 * mass
