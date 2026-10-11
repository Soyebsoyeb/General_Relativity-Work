"""Scalar curvature invariants.

No single scalar detects every singularity: a vanishing invariant does not
imply regularity (null singularities and pp-wave spacetimes have all
polynomial scalar invariants equal to zero), and a divergent invariant
certifies a curvature singularity only in the sense that no coordinate
change can remove it. The Kretschmann scalar is the standard first check.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from ..core import index_operations as ops
from ..core.tensor import Tensor
from .ricci import ricci_scalar, ricci_tensor
from .riemann import riemann_lowered
from .weyl import weyl_tensor


def kretschmann(riem: Tensor):
    """K = R_{abcd} R^{abcd}."""
    Rd = riemann_lowered(riem)
    return Rd.full_contraction(Rd)


def ricci_squared(riem: Tensor):
    """R_{ab} R^{ab}."""
    ric = ricci_tensor(riem)
    return ric.full_contraction(ric)


def weyl_squared(riem: Tensor):
    """C_{abcd} C^{abcd}."""
    C = weyl_tensor(riem)
    return C.full_contraction(C)


def euler_density(riem: Tensor):
    """Gauss-Bonnet density R^2 - 4 R_{ab}R^{ab} + R_{abcd}R^{abcd}."""
    R = ricci_scalar(ricci_tensor(riem))
    return riem.metric.simplifier(R**2 - 4 * ricci_squared(riem) + kretschmann(riem))


def pontryagin_density(riem: Tensor):
    """(1/2) eps^{abcd} R^m_{n a b} R^n_{m c d}  in four dimensions.

    The orientation sign follows eps^{0123} = -1/sqrt(-g). Only the square or
    the vanishing of this scalar is convention independent.
    """
    m = riem.metric
    if m.dim != 4:
        raise ValueError("Pontryagin density is implemented for dimension 4 only")
    eps = ops.levi_civita_symbol(4)
    sg = m.sqrt_minus_g
    tot = sp.S.Zero
    R = riem.c
    for a, b, c, d in np.ndindex(4, 4, 4, 4):
        e = eps[a, b, c, d]
        if e == 0:
            continue
        inner = sum(R[mm, nn, a, b] * R[nn, mm, c, d] for mm in range(4) for nn in range(4))
        tot += -e / sg * inner
    return m.simplifier(tot / 2)


def kretschmann_at(riem: Tensor, point: dict) -> float:
    """Numerical K = R_{abcd} R^{abcd} at a coordinate point.

    For metrics whose symbolic components were not simplified (Kerr). The
    components are evaluated at 30 digits, then contracted in double precision.
    """
    m = riem.metric
    n = m.dim
    R = np.empty((n,) * 4)
    for idx in np.ndindex(R.shape):
        R[idx] = float(sp.N(sp.sympify(riem.c[idx]).subs(point), 30))
    g = np.array(m.g.subs(point).evalf(30).tolist(), dtype=float)
    gi = np.array(m.ginv.subs(point).evalf(30).tolist(), dtype=float)
    low = np.einsum("ae,ebcd->abcd", g, R)  # R_{abcd}
    up = np.einsum("ap,bq,cr,ds,pqrs->abcd", gi, gi, gi, gi, low)
    return float(np.einsum("abcd,abcd->", low, up))
