"""Einstein tensor G_{mu nu} = R_{mu nu} - (1/2) R g_{mu nu}."""

from __future__ import annotations

import numpy as np
import sympy as sp

from ..core.tensor import Tensor


def einstein_tensor(ricci: Tensor, scalar) -> Tensor:
    m = ricci.metric
    g = np.array(m.g.tolist(), dtype=object)
    return Tensor(m, ricci.c - sp.Rational(1, 2) * scalar * g, "dd").simplify()


def field_equation_lhs(ricci: Tensor, scalar, cosmological_constant=0) -> Tensor:
    """G_{mu nu} + Lambda g_{mu nu} (vanishes in vacuum with Lambda = 0)."""
    m = ricci.metric
    g = np.array(m.g.tolist(), dtype=object)
    G = einstein_tensor(ricci, scalar)
    return Tensor(m, G.c + cosmological_constant * g, "dd").simplify()
