"""Flat spacetimes in several coordinate systems (progressively harder tests)."""

from __future__ import annotations

import sympy as sp

from ..core.metric import Metric

t, x, y, z = sp.symbols("t x y z", real=True)
r = sp.Symbol("r", positive=True)
theta = sp.Symbol("theta", positive=True)
phi = sp.Symbol("phi", real=True)


def minkowski(dim: int = 4) -> Metric:
    coords = [t, x, y, z][:dim]
    return Metric(coords, sp.diag(-1, *([1] * (dim - 1))), name=f"Minkowski{dim}D")


def euclidean_polar() -> Metric:
    """ds^2 = dr^2 + r^2 dphi^2: the flat plane in polar coordinates."""
    return Metric([r, phi], sp.diag(1, r**2), name="polar")


def euclidean_spherical() -> Metric:
    return Metric([r, theta, phi], sp.diag(1, r**2, r**2 * sp.sin(theta) ** 2), name="spherical3d")


def minkowski_spherical() -> Metric:
    return Metric(
        [t, r, theta, phi],
        sp.diag(-1, 1, r**2, r**2 * sp.sin(theta) ** 2),
        name="Minkowski-spherical",
    )


def unit_sphere() -> Metric:
    """S^2 of unit radius: constant positive curvature, R = 2."""
    return Metric([theta, phi], sp.diag(1, sp.sin(theta) ** 2), name="S2")


def polar_to_cartesian():
    """(x, y) as functions of (r, phi)."""
    return [r * sp.cos(phi), r * sp.sin(phi)]


def spherical_to_cartesian():
    """(t, x, y, z) as functions of (t, r, theta, phi)."""
    return [t, r * sp.sin(theta) * sp.cos(phi), r * sp.sin(theta) * sp.sin(phi), r * sp.cos(theta)]
