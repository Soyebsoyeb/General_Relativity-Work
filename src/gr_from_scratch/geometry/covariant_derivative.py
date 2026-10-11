"""Covariant derivative of tensors of arbitrary rank and mixed variance.

For each slot the connection adds a term that compensates the change of the
local basis:

    nabla_c T^a  = d_c T^a  + Gamma^a_{c l} T^l       (upper slot)
    nabla_c T_a  = d_c T_a  - Gamma^l_{c a} T_l       (lower slot)

The new derivative index is appended as the last slot.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from ..core import index_operations as ops
from ..core.tensor import Tensor
from .connection import Connection


def covariant_derivative(t: Tensor, conn: Connection) -> Tensor:
    m = t.metric
    n, rank = m.dim, t.rank
    out = ops.zeros(n, rank + 1)
    for idx in np.ndindex(*(n,) * rank):
        for c in range(n):
            val = sp.diff(t.c[idx], m.coords[c])
            for slot, v in enumerate(t.variance):
                for lam in range(n):
                    j = idx[:slot] + (lam,) + idx[slot + 1 :]
                    if v == "u":
                        val += conn.gamma[idx[slot], c, lam] * t.c[j]
                    else:
                        val -= conn.gamma[lam, c, idx[slot]] * t.c[j]
            out[idx + (c,)] = val
    return Tensor(m, out, t.variance + "d").simplify()


def divergence(t: Tensor, conn: Connection, slot: int) -> Tensor:
    """nabla_c T^{... c ...} contracted on the given upper slot."""
    if t.variance[slot] != "u":
        raise ValueError("divergence needs an upper index to contract")
    return covariant_derivative(t, conn).contract(slot, t.rank)
