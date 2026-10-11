"""Killing vectors and the conserved quantities they generate.

A vector field xi is Killing iff nabla_mu xi_nu + nabla_nu xi_mu = 0, which is
equivalent to L_xi g = 0. Along a geodesic with tangent u, g(xi, u) is
constant. This module *verifies* supplied candidates; it does not solve the
Killing equation for an unknown field.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from ..core import index_operations as ops
from ..core.tensor import Tensor
from .connection import Connection


def killing_residual(conn: Connection, xi_up) -> np.ndarray:
    m = conn.metric
    n = m.dim
    xi_dn = [sum(m.g[a, b] * xi_up[b] for b in range(n)) for a in range(n)]
    out = ops.zeros(n, 2)
    for a in range(n):
        for b in range(n):
            nab = sp.diff(xi_dn[b], m.coords[a]) - sum(conn.gamma[s, a, b] * xi_dn[s] for s in range(n))
            nba = sp.diff(xi_dn[a], m.coords[b]) - sum(conn.gamma[s, b, a] * xi_dn[s] for s in range(n))
            out[a, b] = m.simplifier(nab + nba)
    return out


def is_killing(conn: Connection, xi_up) -> bool:
    return all(e == 0 for e in killing_residual(conn, xi_up).ravel())


def conserved_quantity(metric, xi_up, u):
    """g_{mu nu} xi^mu u^nu (constant along geodesics if xi is Killing)."""
    return sum(metric.g[a, b] * xi_up[a] * u[b] for a in range(metric.dim) for b in range(metric.dim))


def as_vector(metric, xi_up) -> Tensor:
    return Tensor(metric, xi_up, "u")
