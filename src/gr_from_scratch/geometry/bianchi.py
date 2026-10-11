"""Residuals of the algebraic and differential Bianchi identities.

Each function returns an array that should vanish. None of them uses the
expected answer; the differential identities rely on the generic covariant
derivative, so they test the connection terms as well.
"""

from __future__ import annotations

import numpy as np

from ..core import index_operations as ops
from ..core.tensor import Tensor
from .connection import Connection
from .covariant_derivative import covariant_derivative
from .riemann import riemann_lowered


def pair_and_antisymmetry_residuals(riem: Tensor) -> dict[str, np.ndarray]:
    R = riemann_lowered(riem).c
    return {
        "antisym_first_pair": R + np.transpose(R, (1, 0, 2, 3)),
        "antisym_second_pair": R + np.transpose(R, (0, 1, 3, 2)),
        "pair_exchange": R - np.transpose(R, (2, 3, 0, 1)),
    }


def algebraic_bianchi_residual(riem: Tensor) -> np.ndarray:
    """R_{a[bcd]}, i.e. the cyclic sum over the last three indices."""
    R = riemann_lowered(riem).c
    return R + np.transpose(R, (0, 2, 3, 1)) + np.transpose(R, (0, 3, 1, 2))


def differential_bianchi_residual(riem: Tensor, conn: Connection) -> np.ndarray:
    """nabla_e R_{abcd} + nabla_a R_{becd} + nabla_b R_{eacd}  (cyclic in e,a,b).

    ``covariant_derivative`` puts the derivative index last: DR[a,b,c,d,e].
    """
    DR = covariant_derivative(riemann_lowered(riem), conn).c
    # term1[a,b,c,d,e] = DR[a,b,c,d,e]
    # term2[a,b,c,d,e] = nabla_a R_{b e c d} = DR[b,e,c,d,a]
    # term3[a,b,c,d,e] = nabla_b R_{e a c d} = DR[e,a,c,d,b]
    t2 = np.transpose(DR, (4, 0, 2, 3, 1))
    t3 = np.transpose(DR, (1, 4, 2, 3, 0))
    return DR + t2 + t3


def contracted_bianchi_residual(einstein_up: Tensor, conn: Connection) -> np.ndarray:
    """nabla_mu G^{mu nu}; ``einstein_up`` must have variance ``"uu"``."""
    D = covariant_derivative(einstein_up, conn)
    return D.contract(0, 2).c


def residual_is_zero(arr: np.ndarray, simplifier) -> bool:
    return all(simplifier(e) == 0 for e in np.asarray(arr, dtype=object).ravel())


__all__ = [
    "algebraic_bianchi_residual",
    "contracted_bianchi_residual",
    "differential_bianchi_residual",
    "ops",
    "pair_and_antisymmetry_residuals",
    "residual_is_zero",
]
