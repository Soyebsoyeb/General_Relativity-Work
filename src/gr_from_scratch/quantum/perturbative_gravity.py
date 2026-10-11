"""Quadratic action of general relativity about flat space.

Status: *partially implemented, symbolic only.* It expands sqrt(-g) R to
second order in h and compares its Euler-Lagrange equations with
(i) the linearised Einstein tensor and (ii) the Fierz-Pauli Lagrangian.
There is no graviton quantisation, no propagator and no amplitude code.

With g = eta + eps h, S = int sqrt(-g) R gives
    delta S / delta g_{mn} = - sqrt(-g) G^{mn},
so at order eps the Euler-Lagrange derivative of the eps^2 Lagrangian with
respect to h_{mn} is  - G^{(1) mn}[h].
"""

from __future__ import annotations

import sympy as sp
from sympy.calculus.euler import euler_equations

from ..field_equations.linearized_gravity import ETA, linearized_einstein, trace
from ..field_equations.perturbation_series import order_coefficient, perturbative_geometry


def symmetric_perturbation(x, name="h", depends_on=None):
    """Independent component functions h_{mn}, m <= n, of the coordinates ``depends_on``.

    ``depends_on`` defaults to all of ``x``. Restricting it (for example to
    (t, x)) keeps the symbolic expansion cheap while all ten components and
    the full four-dimensional tensor structure are retained.
    """
    n = len(x)
    args = tuple(depends_on) if depends_on is not None else tuple(x)
    comps = {}
    h = sp.zeros(n)
    for m in range(n):
        for v in range(m, n):
            f = sp.Function(f"{name}{m}{v}")(*args)
            comps[(m, v)] = f
            h[m, v] = h[v, m] = f
    return h, comps


def quadratic_lagrangian(h, x):
    """Coefficient of eps^2 in sqrt(-g) R (no integration by parts applied)."""
    geo = perturbative_geometry(h, x, order=2)
    return sp.expand(order_coefficient(sp.expand(geo["sqrt_g"] * geo["R"]), 2))


def fierz_pauli_lagrangian(h, x):
    """L_FP = -(1/2) d_l h_mn d^l h^mn + d_m h^mn d^l h_ln - d_m h^mn d_n h + (1/2) d_l h d^l h."""
    n = len(x)
    ei = ETA.inv()

    def up2(a, b):
        return ei[a, a] * ei[b, b] * h[a, b]

    t1 = -sp.Rational(1, 2) * sum(
        ei[l, l] * sp.diff(h[m, v], x[l]) * sp.diff(up2(m, v), x[l])
        for l in range(n)
        for m in range(n)
        for v in range(n)
    )
    div = [sum(ei[m, m] * sp.diff(h[m, v], x[m]) for m in range(n)) for v in range(n)]  # d^m h_{mv}
    divu = [ei[v, v] * div[v] for v in range(n)]
    t2 = sum(div[v] * divu[v] for v in range(n))
    tr = trace(h)
    t3 = -sum(divu[v] * sp.diff(tr, x[v]) for v in range(n))
    t4 = sp.Rational(1, 2) * sum(ei[l, l] * sp.diff(tr, x[l]) ** 2 for l in range(n))
    return sp.expand(t1 + t2 + t3 + t4)


def euler_lagrange_components(L, comps, variables):
    funcs = list(comps.values())
    eqs = euler_equations(L, funcs, list(variables))
    return {key: eq.lhs for key, eq in zip(comps.keys(), eqs, strict=True)}


def compare_with_linearized_einstein(x, scale=1, depends_on=None):
    """Residuals E-L(L2)_{mn} - w_{mn} * scale * (-G^{(1) mn}), w = 1 diagonal, 2 off-diagonal.

    The factor 2 for m != n counts the two equal components h_{mn} = h_{nm}.
    """
    h, comps = symmetric_perturbation(x, depends_on=depends_on)
    L2 = quadratic_lagrangian(h, x)
    el = euler_lagrange_components(L2, comps, depends_on if depends_on is not None else x)
    G1 = linearized_einstein(h.tolist(), x)
    ei = ETA.inv()
    out = {}
    for (m, v), val in el.items():
        w = 1 if m == v else 2
        Gup = ei[m, m] * ei[v, v] * G1[m, v]
        out[(m, v)] = sp.simplify(sp.expand(val - w * scale * (-Gup)))
    return out, L2, h, comps
