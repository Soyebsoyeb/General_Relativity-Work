"""Geodesic deviation (tidal acceleration).

D^2 xi^mu / dtau^2 = - R^mu_{nu rho sigma} u^nu xi^rho u^sigma.

In an orthonormal frame carried by the observer, with E_ij = R^i_{0 j 0}
(frame components), the spatial relative acceleration is
d^2 xi^i / dtau^2 = - E_ij xi^j.
"""

from __future__ import annotations

import numpy as np

from ..core.tensor import Tensor
from ..geometry.tetrads import Coframe


def deviation_acceleration(riem: Tensor, u, xi) -> list:
    n = riem.metric.dim
    out = []
    for m in range(n):
        tot = 0
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    tot += riem.c[m, a, b, c] * u[a] * xi[b] * u[c]
        out.append(riem.metric.simplifier(-tot))
    return out


def tidal_tensor(coframe: Coframe) -> np.ndarray:
    """E_ij = R^i_{0 j 0} in the orthonormal frame (rows/columns 1..n-1)."""
    R = coframe.riemann_frame()
    n = coframe.n
    out = np.empty((n - 1, n - 1), dtype=object)
    for i in range(1, n):
        for j in range(1, n):
            out[i - 1, j - 1] = coframe.metric.simplifier(R[i, 0, j, 0])
    return out
