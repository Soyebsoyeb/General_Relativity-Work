"""Symbolic test of conservation laws along geodesics.

For a vector field xi and the quantity Q = g_{mu nu} xi^mu u^nu, with du/dl
given by the geodesic equation,

    dQ/dl = (1/2) (L_xi g)_{ab} u^a u^b,

which vanishes for every geodesic exactly when xi is a Killing vector. The
function below evaluates dQ/dl directly from the definition, using the
Euler-Lagrange accelerations, so it does not assume the Killing property.
"""

from __future__ import annotations

import sympy as sp

from ..core.metric import Metric
from .variational_principles import euler_lagrange_acceleration, velocity_symbols


def conservation_rate(metric: Metric, xi_up) -> sp.Expr:
    n, x = metric.dim, metric.coords
    u = velocity_symbols(metric)
    acc = euler_lagrange_acceleration(metric, u)
    xi_dn = [sum(metric.g[a, b] * xi_up[b] for b in range(n)) for a in range(n)]
    dq = 0
    for m in range(n):
        dq += sum(sp.diff(xi_dn[m], x[a]) * u[a] * u[m] for a in range(n)) + xi_dn[m] * acc[m]
    return metric.simplifier(sp.expand(dq))


def is_conserved_along_geodesics(metric: Metric, xi_up) -> bool:
    return conservation_rate(metric, xi_up) == 0
