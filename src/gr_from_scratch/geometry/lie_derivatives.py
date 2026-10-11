"""Lie derivative of a tensor along a vector field (coordinate formula)."""

from __future__ import annotations

import numpy as np
import sympy as sp

from ..core import index_operations as ops
from ..core.tensor import Tensor


def lie_derivative(t: Tensor, xi) -> Tensor:
    """(L_xi T)^{a...}_{b...} = xi^c d_c T - sum_up d_c xi^a T^{..c..} + sum_down d_b xi^c T_{..c..}."""
    m = t.metric
    n, x = m.dim, m.coords
    out = ops.zeros(n, t.rank)
    for idx in np.ndindex(*(n,) * t.rank):
        val = sum(xi[c] * sp.diff(t.c[idx], x[c]) for c in range(n))
        for slot, v in enumerate(t.variance):
            for c in range(n):
                j = idx[:slot] + (c,) + idx[slot + 1 :]
                if v == "u":
                    val -= sp.diff(xi[idx[slot]], x[c]) * t.c[j]
                else:
                    val += sp.diff(xi[c], x[idx[slot]]) * t.c[j]
        out[idx] = val
    return Tensor(m, out, t.variance).simplify()
