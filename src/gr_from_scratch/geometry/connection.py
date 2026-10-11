"""Levi-Civita connection.

Derivation (sketch; full version in ``docs/derivations/01_connection.md``)
    Metric compatibility  nabla_l g_{mn} = d_l g_{mn} - G^s_{lm} g_{sn} - G^s_{ln} g_{ms} = 0
    and zero torsion  G^r_{mn} = G^r_{nm}  are written three times with the
    indices cyclically permuted. Adding two copies and subtracting the third
    isolates G_{s m n} = g_{s r} G^r_{mn}, and raising the index gives

        G^r_{mn} = 1/2 g^{rs} (d_m g_{ns} + d_n g_{ms} - d_s g_{mn}).

The connection is stored as a plain array ``gamma[rho, mu, nu]``. It is
deliberately **not** a ``Tensor``: its components do not transform
homogeneously (see ``core.coordinates.transform_connection``).
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from ..core import index_operations as ops
from ..core.metric import Metric


class Connection:
    def __init__(self, metric: Metric, gamma: np.ndarray) -> None:
        self.metric = metric
        self.gamma = gamma

    def __getitem__(self, idx):
        return self.gamma[idx]

    def nonzero(self) -> dict:
        s = self.metric.simplifier
        return {i: s(self.gamma[i]) for i in np.ndindex(self.gamma.shape) if s(self.gamma[i]) != 0}


def christoffel(metric: Metric) -> Connection:
    n, x, g, gi, simp = metric.dim, metric.coords, metric.g, metric.ginv, metric.simplifier
    dg = [[[sp.diff(g[a, b], x[c]) for c in range(n)] for b in range(n)] for a in range(n)]
    gam = ops.zeros(n, 3)
    for r in range(n):
        for m in range(n):
            for v in range(m, n):
                tot = sp.S.Zero
                for s in range(n):
                    if gi[r, s] != 0:
                        # d_m g_{v s} + d_v g_{m s} - d_s g_{m v}
                        tot += gi[r, s] * (dg[v][s][m] + dg[m][s][v] - dg[m][v][s])
                val = simp(tot / 2)
                gam[r, m, v] = val
                gam[r, v, m] = val
    return Connection(metric, gam)


def torsion(conn: Connection) -> np.ndarray:
    """T^rho_{mu nu} = Gamma^rho_{mu nu} - Gamma^rho_{nu mu}."""
    return conn.gamma - np.transpose(conn.gamma, (0, 2, 1))


def metric_compatibility_residual(conn: Connection) -> np.ndarray:
    """nabla_l g_{mn} from the definition, as an array [l, m, n]."""
    m = conn.metric
    n = m.dim
    out = ops.zeros(n, 3)
    for lam in range(n):
        for a in range(n):
            for b in range(n):
                val = sp.diff(m.g[a, b], m.coords[lam])
                for s in range(n):
                    val -= conn.gamma[s, lam, a] * m.g[s, b] + conn.gamma[s, lam, b] * m.g[a, s]
                out[lam, a, b] = m.simplifier(val)
    return out
