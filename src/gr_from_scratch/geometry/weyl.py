"""Weyl (conformal) curvature tensor in n >= 3 dimensions.

C_{abcd} = R_{abcd}
           - 2/(n-2) (g_{a[c} R_{d]b} - g_{b[c} R_{d]a})
           + 2/((n-1)(n-2)) R g_{a[c} g_{d]b}

Dimension dependence: the Weyl tensor vanishes identically for n = 3, where
the Riemann tensor is fixed by the Ricci tensor. The code does not special-case
this; the n = 3 test checks it as an outcome.
"""

from __future__ import annotations

import sympy as sp

from ..core import index_operations as ops
from ..core.tensor import Tensor
from .ricci import ricci_scalar, ricci_tensor
from .riemann import riemann_lowered


class DimensionError(ValueError):
    pass


def weyl_tensor(riem: Tensor) -> Tensor:
    m = riem.metric
    n = m.dim
    if n < 3:
        raise DimensionError("the Weyl tensor is not defined for dimension < 3")
    Rdown = riemann_lowered(riem)
    ric = ricci_tensor(riem)
    R = ricci_scalar(ric)
    g = m.g
    out = ops.zeros(n, 4)
    half = sp.Rational(1, 2)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    t1 = half * (g[a, c] * ric.c[d, b] - g[a, d] * ric.c[c, b])
                    t2 = half * (g[b, c] * ric.c[d, a] - g[b, d] * ric.c[c, a])
                    t3 = half * (g[a, c] * g[d, b] - g[a, d] * g[c, b])
                    out[a, b, c, d] = (
                        Rdown.c[a, b, c, d]
                        - sp.Rational(2, n - 2) * (t1 - t2)
                        + sp.Rational(2, (n - 1) * (n - 2)) * R * t3
                    )
    return Tensor(m, out, "dddd").simplify()
