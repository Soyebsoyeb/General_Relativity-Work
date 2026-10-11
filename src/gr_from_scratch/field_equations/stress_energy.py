"""Stress-energy tensors (G = c = 1, Gaussian units for electromagnetism)."""

from __future__ import annotations

import numpy as np
import sympy as sp

from ..core import index_operations as ops
from ..core.metric import Metric
from ..core.tensor import Tensor


def perfect_fluid(metric: Metric, rho, pressure, u_up) -> Tensor:
    """T_{mu nu} = (rho + p) u_mu u_nu + p g_{mu nu}; ``u_up`` must satisfy u.u = -1."""
    n = metric.dim
    u_dn = [sum(metric.g[a, b] * u_up[b] for b in range(n)) for a in range(n)]
    T = ops.zeros(n, 2)
    for a in range(n):
        for b in range(n):
            T[a, b] = (rho + pressure) * u_dn[a] * u_dn[b] + pressure * metric.g[a, b]
    return Tensor(metric, T, "dd").simplify()


def field_strength(metric: Metric, A_dn) -> Tensor:
    """F_{mu nu} = d_mu A_nu - d_nu A_mu."""
    n = metric.dim
    F = ops.zeros(n, 2)
    for a in range(n):
        for b in range(n):
            F[a, b] = sp.diff(A_dn[b], metric.coords[a]) - sp.diff(A_dn[a], metric.coords[b])
    return Tensor(metric, F, "dd").simplify()


def electromagnetic(metric: Metric, F: Tensor) -> Tensor:
    """T_{mu nu} = (1/4 pi) (F_{mu a} F_nu^a - (1/4) g_{mu nu} F_{ab} F^{ab})."""
    n = metric.dim
    Fup2 = F.raise_index(1)  # F_{mu}^{a}
    F2 = F.full_contraction(F)
    T = ops.zeros(n, 2)
    for a in range(n):
        for b in range(n):
            s = sum(F.c[a, c] * Fup2.c[b, c] for c in range(n))
            T[a, b] = (s - sp.Rational(1, 4) * metric.g[a, b] * F2) / (4 * sp.pi)
    return Tensor(metric, T, "dd").simplify()


def maxwell_divergence(metric: Metric, F: Tensor) -> list:
    """(1/sqrt|g|) d_mu (sqrt|g| F^{mu nu}); vanishes for source-free Maxwell fields."""
    Fuu = F.raise_index(0).raise_index(1)
    sg = metric.sqrt_minus_g
    out = []
    for nu in range(metric.dim):
        tot = sum(sp.diff(sg * Fuu.c[mu, nu], metric.coords[mu]) for mu in range(metric.dim))
        out.append(metric.simplifier(tot / sg))
    return out


def scalar_field(metric: Metric, phi_expr, potential=0) -> Tensor:
    """T_{mu nu} = d_mu phi d_nu phi - g_{mu nu} ((1/2)(d phi)^2 + V)."""
    n = metric.dim
    d = [sp.diff(phi_expr, x) for x in metric.coords]
    kin = sum(metric.ginv[a, b] * d[a] * d[b] for a in range(n) for b in range(n))
    T = ops.zeros(n, 2)
    for a in range(n):
        for b in range(n):
            T[a, b] = d[a] * d[b] - metric.g[a, b] * (kin / 2 + potential)
    return Tensor(metric, T, "dd").simplify()


__all__ = ["electromagnetic", "field_strength", "maxwell_divergence", "perfect_fluid", "scalar_field", "np"]
