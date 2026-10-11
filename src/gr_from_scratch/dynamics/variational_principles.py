"""Geodesics from a variational principle, derived without Christoffel symbols.

The Lagrangian L = (1/2) g_{mu nu} xdot^mu xdot^nu has Euler-Lagrange equations

    d/dl (g_{mu nu} u^nu) - (1/2) d_mu g_{ab} u^a u^b = 0.

Expanding the total derivative gives g_{mu nu} a^nu = (1/2) d_mu g_{ab} u^a u^b
- d_a g_{mu nu} u^a u^nu, which is solved for the acceleration with g^{-1}. The
test suite compares the result with -Gamma u u from the connection.
"""

from __future__ import annotations

import sympy as sp

from ..core.metric import Metric


def velocity_symbols(metric: Metric) -> list[sp.Symbol]:
    return [sp.Symbol(f"u{i}") for i in range(metric.dim)]


def euler_lagrange_acceleration(metric: Metric, u=None) -> list[sp.Expr]:
    n, x, g = metric.dim, metric.coords, metric.g
    u = u or velocity_symbols(metric)
    rhs = []
    for m in range(n):
        half = sum(sp.diff(g[a, b], x[m]) * u[a] * u[b] for a in range(n) for b in range(n)) / 2
        conn = sum(sp.diff(g[m, v], x[a]) * u[a] * u[v] for a in range(n) for v in range(n))
        rhs.append(half - conn)
    acc = []
    for s in range(n):
        acc.append(metric.simplifier(sum(metric.ginv[s, m] * rhs[m] for m in range(n))))
    return acc
