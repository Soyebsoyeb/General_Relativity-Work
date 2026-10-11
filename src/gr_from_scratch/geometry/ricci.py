"""Ricci tensor and scalar curvature."""

from __future__ import annotations

from ..core.tensor import Tensor


def ricci_tensor(riem: Tensor) -> Tensor:
    """R_{sigma nu} = R^rho_{sigma rho nu}."""
    return riem.contract(0, 2).simplify()


def ricci_scalar(ricci: Tensor):
    """R = g^{sigma nu} R_{sigma nu}."""
    m = ricci.metric
    tot = 0
    for a in range(m.dim):
        for b in range(m.dim):
            tot += m.ginv[a, b] * ricci.c[a, b]
    return m.simplifier(tot)
