"""Coordinate transformations of tensors and of the connection.

The connection is *not* a tensor. ``transform_connection`` implements its
inhomogeneous transformation law, which tests use to show that Christoffel
symbols can be made to vanish at a point by a change of coordinates while the
curvature cannot.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from . import index_operations as ops
from .metric import Metric, as_symbols
from .tensor import Tensor


def jacobian(old_in_new, new_coords) -> sp.Matrix:
    """J[mu, a] = d x^mu / d x'^a."""
    new = as_symbols(new_coords)
    return sp.Matrix(len(old_in_new), len(new), lambda m, a: sp.diff(old_in_new[m], new[a]))


def transform_tensor(t: Tensor, new_metric: Metric, old_in_new) -> Tensor:
    """Components of ``t`` in the coordinates of ``new_metric``.

    ``old_in_new`` expresses the old coordinates through the new ones. Upper
    indices transform with J^-1 = dx'/dx, lower indices with J = dx/dx'.
    """
    old = t.metric.coords
    sub = dict(zip(old, old_in_new, strict=True))
    jac = jacobian(old_in_new, new_metric.coords)
    jinv = jac.inv().applyfunc(new_metric.simplifier)
    arr = np.empty(t.c.shape, dtype=object)
    for idx in np.ndindex(t.c.shape):
        arr[idx] = sp.sympify(t.c[idx]).subs(sub)
    for axis, v in enumerate(t.variance):
        arr = ops.apply_matrix(arr, jinv if v == "u" else jac.T, axis)
    return Tensor(new_metric, arr, t.variance).simplify()


def transform_connection(gamma, old_coords, new_coords, old_in_new, simplifier):
    """Christoffel symbols in new coordinates via the inhomogeneous law.

    Gamma'^a_{bc} = (dx'^a/dx^r)(dx^m/dx'^b)(dx^n/dx'^c) Gamma^r_{mn}
                    + (dx'^a/dx^r) d^2 x^r / dx'^b dx'^c
    """
    new = as_symbols(new_coords)
    n = len(new)
    sub = dict(zip(old_coords, old_in_new, strict=True))
    jac = jacobian(old_in_new, new)
    jinv = jac.inv().applyfunc(simplifier)
    out = ops.zeros(n, 3)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                tot = sp.S.Zero
                for r in range(n):
                    second = sp.diff(old_in_new[r], new[b], new[c])
                    tot += jinv[a, r] * second
                    for m in range(n):
                        for k in range(n):
                            g = sp.sympify(gamma[r, m, k]).subs(sub)
                            tot += jinv[a, r] * jac[m, b] * jac[k, c] * g
                out[a, b, c] = simplifier(tot)
    return out
