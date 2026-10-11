"""Linearised gravity about Minkowski space.

Metric g = eta + h with |h| << 1; indices are moved with eta = diag(-1,1,1,1).
Curvature is expanded to first order in h, in the package convention:

    R_{abcd} = 1/2 ( d_b d_c h_{ad} + d_a d_d h_{bc} - d_a d_c h_{bd} - d_b d_d h_{ac} )
    R_{bd}   = 1/2 ( d_a d_b h^a_d + d_a d_d h^a_b - box h_{bd} - d_b d_d h )
    G_{ab}   = R_{ab} - (1/2) eta_{ab} R

Gauge: h -> h + d_a xi_b + d_b xi_a leaves R_{abcd} unchanged (checked
symbolically with arbitrary functions). The statement that this removes
unphysical degrees of freedom is a counting argument given in
``docs/derivations/09_linearized_gravity.md``.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

ETA = sp.diag(-1, 1, 1, 1)


def _nest(h):
    """Accept a nested list or a Matrix and return a nested list."""
    return sp.Matrix(h).tolist()


def linearized_riemann(h, x) -> np.ndarray:
    """R^{(1)}_{abcd} as a rank-4 object array."""
    h = _nest(h)
    n = len(x)
    out = np.empty((n,) * 4, dtype=object)
    for a, b, c, d in np.ndindex(n, n, n, n):
        out[a, b, c, d] = sp.Rational(1, 2) * (
            sp.diff(h[a][d], x[b], x[c])
            + sp.diff(h[b][c], x[a], x[d])
            - sp.diff(h[b][d], x[a], x[c])
            - sp.diff(h[a][c], x[b], x[d])
        )
    return out


def linearized_ricci(h, x) -> sp.Matrix:
    n = len(x)
    Rm = linearized_riemann(h, x)
    eta_inv = ETA.inv()
    return sp.Matrix(n, n, lambda b, d: sum(eta_inv[a, c] * Rm[a, b, c, d] for a in range(n) for c in range(n)))


def linearized_scalar(h, x):
    eta_inv = ETA.inv()
    ric = linearized_ricci(h, x)
    return sum(eta_inv[a, b] * ric[a, b] for a in range(len(x)) for b in range(len(x)))


def linearized_einstein(h, x) -> sp.Matrix:
    h = _nest(h)
    ric, R = linearized_ricci(h, x), linearized_scalar(h, x)
    return ric - sp.Rational(1, 2) * ETA * R


def gauge_transform(h, xi_up, x):
    """h_{ab} -> h_{ab} + d_a xi_b + d_b xi_a, with xi_b = eta_{bc} xi^c."""
    h = _nest(h)
    n = len(x)
    xi_dn = [sum(ETA[b, c] * xi_up[c] for c in range(n)) for b in range(n)]
    return sp.Matrix(n, n, lambda a, b: h[a][b] + sp.diff(xi_dn[b], x[a]) + sp.diff(xi_dn[a], x[b]))


def trace(h):
    h = _nest(h)
    eta_inv = ETA.inv()
    return sum(eta_inv[a, b] * h[a][b] for a in range(4) for b in range(4))


def trace_reverse(h):
    """hbar_{ab} = h_{ab} - (1/2) eta_{ab} h."""
    h = _nest(h)
    tr = trace(h)
    return sp.Matrix(4, 4, lambda a, b: h[a][b] - sp.Rational(1, 2) * ETA[a, b] * tr)


def box(expr, x):
    return -sp.diff(expr, x[0], 2) + sum(sp.diff(expr, x[i], 2) for i in range(1, len(x)))


def harmonic_gauge_residual(h, x):
    """d^a hbar_{ab}; harmonic (Lorenz) gauge sets this to zero."""
    hb = trace_reverse(h)
    eta_inv = ETA.inv()
    return [sum(eta_inv[a, c] * sp.diff(hb[c, b], x[a]) for a in range(4) for c in range(4)) for b in range(4)]


def plane_wave_tt(x, amp_plus, amp_cross, omega, k=None):
    """TT plane wave along +z: h_xx = -h_yy = A+ cos(w t - k z), h_xy = Ax cos(w t - k z)."""
    t, _, _, z = x
    k = omega if k is None else k
    ph = sp.cos(omega * t - k * z)
    h = sp.zeros(4, 4)
    h[1, 1], h[2, 2] = amp_plus * ph, -amp_plus * ph
    h[1, 2] = h[2, 1] = amp_cross * ph
    return h


def tt_conditions(h, x):
    """Return (trace, h_{0a} components, d^a h_{ab} divergence) for a TT check."""
    eta_inv = ETA.inv()
    div = [sum(eta_inv[a, c] * sp.diff(h[c, b], x[a]) for a in range(4) for c in range(4)) for b in range(4)]
    return sp.simplify(trace(h)), [sp.simplify(h[0, a]) for a in range(4)], [sp.simplify(d) for d in div]
