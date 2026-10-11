"""Riemann tensor.

R^r_{s m n} = d_m G^r_{n s} - d_n G^r_{m s} + G^r_{m l} G^l_{n s} - G^r_{n l} G^l_{m s}

Derived from [nabla_m, nabla_n] V^r = R^r_{s m n} V^s for zero torsion.
"""

from __future__ import annotations

import sympy as sp

from ..core import index_operations as ops
from ..core.tensor import Tensor
from .connection import Connection


def riemann(conn: Connection, simplify: bool = True) -> Tensor:
    """R^rho_{sigma mu nu} with variance ``"uddd"``."""
    m = conn.metric
    n, x, G = m.dim, m.coords, conn.gamma
    out = ops.zeros(n, 4)
    for r in range(n):
        for s in range(n):
            for mu in range(n):
                for nu in range(mu + 1, n):
                    val = sp.diff(G[r, nu, s], x[mu]) - sp.diff(G[r, mu, s], x[nu])
                    for lam in range(n):
                        val += G[r, mu, lam] * G[lam, nu, s] - G[r, nu, lam] * G[lam, mu, s]
                    val = m.simplifier(val) if simplify else val
                    out[r, s, mu, nu] = val
                    out[r, s, nu, mu] = -val
    return Tensor(m, out, "uddd")


def riemann_lowered(riem: Tensor) -> Tensor:
    """R_{rho sigma mu nu}."""
    return riem.lower_index(0)
