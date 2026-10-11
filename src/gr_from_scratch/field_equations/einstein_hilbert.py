"""Selected variational identities for the Einstein-Hilbert action.

What is implemented (each is a test, not an assertion):

* delta sqrt|g| = (1/2) sqrt|g| g^{mn} h_{mn}                       for g -> g + eps h
* delta g^{mn} = - h^{mn}
* delta R = - R^{mn} h_{mn} + nabla_m nabla_n h^{mn} - box h
* total-derivative structure: sqrt|g| nabla_m V^m = d_m (sqrt|g| V^m),
  V^m = nabla_n h^{mn} - nabla^m h
* delta [ sqrt|g| (R - 2 Lambda) ] = - sqrt|g| (G^{mn} + Lambda g^{mn}) h_{mn} + d_m (sqrt|g| V^m)
* a reduced (minisuperspace) action for static spherical metrics whose
  Euler-Lagrange equations reproduce the Einstein tensor components.

What is not implemented: a general-purpose symbolic variational engine. The
identities are verified on user-supplied backgrounds and perturbations h.
The reduced action relies on the principle of symmetric criticality, which
holds for the compact symmetry group SO(3) used here (Palais 1979).
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from sympy.calculus.euler import euler_equations

from ..core.metric import Metric
from ..core.simplification import identity
from ..core.tensor import Tensor
from ..geometry.covariant_derivative import covariant_derivative
from .vacuum_equations import compute_geometry

eps = sp.Symbol("epsilon", real=True)


def _vol(metric: Metric, g=None):
    """sqrt|det g| with the sign of det fixed from the background metric.

    Using Abs would introduce non-differentiable artefacts; the background
    determinant has a definite sign on the patch, which stays true for small eps.
    """
    pt = {c: sp.Rational(37, 100) + sp.Rational(11, 100) * i for i, c in enumerate(metric.coords)}
    sign = -1 if sp.N(metric.det.subs(pt)) < 0 else 1
    gg = metric.g if g is None else g
    return sp.sqrt(sign * gg.det())


def _perturbed(metric: Metric, h) -> Metric:
    g = metric.g + eps * sp.Matrix(h)
    return Metric(metric.coords, g, simplifier=identity, name="perturbed")


def _d_eps0(expr):
    return sp.simplify(sp.diff(expr, eps).subs(eps, 0))


def det_variation_residual(metric: Metric, h):
    """d/d(eps) sqrt|g + eps h| at 0 minus (1/2) sqrt|g| g^{mn} h_{mn}."""
    pm = _perturbed(metric, h)
    lhs = _d_eps0(_vol(metric, pm.g))
    sg = _vol(metric)
    rhs = sg * sum(metric.ginv[a, b] * h[a][b] for a in range(metric.dim) for b in range(metric.dim)) / 2
    return sp.simplify(lhs - rhs)


def inverse_variation_residual(metric: Metric, h):
    """d/d(eps) (g + eps h)^{-1} at 0 plus g^{-1} h g^{-1}."""
    pm = _perturbed(metric, h)
    lhs = pm.ginv.applyfunc(_d_eps0)
    rhs = -(metric.ginv * sp.Matrix(h) * metric.ginv)
    return (lhs - rhs).applyfunc(sp.simplify)


def _scalar_curvature(metric: Metric):
    return compute_geometry(metric).scalar


def ricci_scalar_variation_terms(metric: Metric, h):
    """Return (-R^{mn} h_{mn}, nabla_m nabla_n h^{mn}, box h) for the background."""
    geom = compute_geometry(metric)
    n = metric.dim
    ht = Tensor(metric, np.array(sp.Matrix(h).tolist(), dtype=object), "dd")
    Dh = covariant_derivative(ht, geom.conn)
    DDh = covariant_derivative(Dh, geom.conn)  # [mu, nu, alpha, beta] = nabla_beta nabla_alpha h_{mu nu}
    gi = metric.ginv
    div_div = sum(
        gi[b, m] * gi[a, v] * DDh.c[m, v, a, b] for m in range(n) for v in range(n) for a in range(n) for b in range(n)
    )
    box_tr = sum(
        gi[a, b] * gi[m, v] * DDh.c[m, v, a, b] for m in range(n) for v in range(n) for a in range(n) for b in range(n)
    )
    ricci_up = geom.ricci.raise_index(0).raise_index(1)
    rh = sum(ricci_up.c[a, b] * h[a][b] for a in range(n) for b in range(n))
    return metric.simplifier(-rh), metric.simplifier(div_div), metric.simplifier(box_tr)


def ricci_scalar_variation_residual(metric: Metric, h):
    """d R/d eps - ( -R^{mn} h_{mn} + nabla nabla h - box h ); should vanish."""
    pm = _perturbed(metric, h)
    lhs = _d_eps0(_scalar_curvature(pm))
    a, b, c = ricci_scalar_variation_terms(metric, h)
    return sp.simplify(lhs - (a + b - c))


def boundary_term_residual(metric: Metric, h):
    """nabla_m V^m - (1/sqrt|g|) d_m (sqrt|g| V^m) for V^m = nabla_n h^{mn} - nabla^m h."""
    geom = compute_geometry(metric)
    n = metric.dim
    gi = metric.ginv
    ht = Tensor(metric, np.array(sp.Matrix(h).tolist(), dtype=object), "dd")
    Dh = covariant_derivative(ht, geom.conn)  # [m, n, a]
    trace = sum(gi[a, b] * h[a][b] for a in range(n) for b in range(n))
    # V_m = g^{ab} nabla_a h_{m b} - d_m trace
    Vd = [
        sum(gi[a, b] * Dh.c[m, b, a] for a in range(n) for b in range(n)) - sp.diff(trace, metric.coords[m])
        for m in range(n)
    ]
    Vu = [sum(gi[m, k] * Vd[k] for k in range(n)) for m in range(n)]
    Vt = Tensor(metric, np.array(Vu, dtype=object), "u")
    cov_div = covariant_derivative(Vt, geom.conn).contract(0, 1).c.item()
    sg = _vol(metric)
    dens_div = sum(sp.diff(sg * Vu[m], metric.coords[m]) for m in range(n)) / sg
    return sp.simplify(cov_div - dens_div)


def full_action_variation_residual(metric: Metric, h, cosmological_constant=0):
    """d/d eps [sqrt|g|(R - 2 Lambda)] - ( -sqrt|g|(G^{mn} + Lambda g^{mn}) h_{mn} + d_m(sqrt|g| V^m) )."""
    n = metric.dim
    geom = compute_geometry(metric)
    pm = _perturbed(metric, h)
    L_eps = _vol(metric, pm.g) * (_scalar_curvature(pm) - 2 * cosmological_constant)
    lhs = _d_eps0(L_eps)
    sg = _vol(metric)
    Gup = geom.einstein.raise_index(0).raise_index(1)
    bulk = -sg * sum(
        (Gup.c[a, b] + cosmological_constant * metric.ginv[a, b]) * h[a][b] for a in range(n) for b in range(n)
    )
    gi = metric.ginv
    ht = Tensor(metric, np.array(sp.Matrix(h).tolist(), dtype=object), "dd")
    Dh = covariant_derivative(ht, geom.conn)
    trace = sum(gi[a, b] * h[a][b] for a in range(n) for b in range(n))
    Vd = [
        sum(gi[a, b] * Dh.c[m, b, a] for a in range(n) for b in range(n)) - sp.diff(trace, metric.coords[m])
        for m in range(n)
    ]
    Vu = [sum(gi[m, k] * Vd[k] for k in range(n)) for m in range(n)]
    boundary = sum(sp.diff(sg * Vu[m], metric.coords[m]) for m in range(n))
    return sp.simplify(lhs - (bulk + boundary))


def reduced_static_spherical_equations():
    """Euler-Lagrange equations of the reduced Lagrangian sqrt|g| R for
    ds^2 = -A(r) dt^2 + B(r) dr^2 + r^2 dOmega^2, versus the Einstein tensor.

    With S = int sqrt|g| R and dS = - int sqrt|g| G^{mn} dg_{mn} (up to boundary
    terms), varying A (g_tt = -A) and B (g_rr = B) gives
        E_A =  sqrt|g| G^{tt},        E_B = - sqrt|g| G^{rr}.
    The residuals returned here should be zero.
    """
    r, th, ph, t = sp.symbols("r theta phi t", positive=True)
    A, B = sp.Function("A", positive=True)(r), sp.Function("B", positive=True)(r)
    m = Metric([t, r, th, ph], sp.diag(-A, B, r**2, r**2 * sp.sin(th) ** 2))
    geom = compute_geometry(m)
    sg = r**2 * sp.sin(th) * sp.sqrt(A * B)
    L = sp.simplify(sg * geom.scalar)
    eqs = euler_equations(L, [A, B], r)
    Gup = geom.einstein.raise_index(0).raise_index(1).c
    resA = sp.simplify(eqs[0].lhs - sg * Gup[0, 0])
    resB = sp.simplify(eqs[1].lhs + sg * Gup[1, 1])
    return resA, resB
