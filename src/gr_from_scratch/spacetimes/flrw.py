"""Friedmann-Lemaitre-Robertson-Walker spacetime, ds^2 = -dt^2 + a(t)^2 dSigma_k^2."""

from __future__ import annotations

import sympy as sp

from ..core.metric import Metric

t = sp.Symbol("t", real=True)
r = sp.Symbol("r", positive=True)
theta = sp.Symbol("theta", positive=True)
phi = sp.Symbol("phi", real=True)
k = sp.Symbol("k", real=True)
a = sp.Function("a", positive=True)


def flrw(curvature=k, scale=a) -> Metric:
    A = scale(t)
    return Metric(
        [t, r, theta, phi],
        sp.diag(-1, A**2 / (1 - curvature * r**2), A**2 * r**2, A**2 * r**2 * sp.sin(theta) ** 2),
        name="FLRW",
    )
