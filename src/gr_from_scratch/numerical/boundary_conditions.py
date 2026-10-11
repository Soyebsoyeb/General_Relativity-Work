"""Boundary data for the radial wave equation psi'' = (V - omega^2) psi in tortoise coordinate.

Horizon (r* -> -infinity): purely ingoing, psi = exp(-i w r*).
Infinity (r* -> +infinity): psi = A_in u^-(w r*) + A_out u^+(w r*), where u^+/-
are the flat-space radial solutions of the centrifugal problem,
u^+(x) = i^(l+1) x h_l^(1)(x) -> e^{ix}, u^-(x) = conj(u^+(x)) -> e^{-ix}.
With l = None plane waves e^{+/-ix} are used. The Wronskian of u^+ and u^- is
2i, so |A_in|^2 - |A_out|^2 is conserved by the true evolution and its
deviation from 1 measures the matching error (~ M ln r* / r*^3 for the
Bessel matching, ~ l(l+1)/(w^2 r*^2) for plane waves).
"""

from __future__ import annotations

import numpy as np


def ingoing_horizon_data(omega: float, rstar: float):
    psi = np.exp(-1j * omega * rstar)
    return psi, -1j * omega * psi


def _riccati_hankel(l: int, x: float):
    """u^+(x) and du^+/dx with u^+ = i^(l+1) x h_l^(1)(x)."""
    from scipy.special import spherical_jn, spherical_yn

    h = spherical_jn(l, x) + 1j * spherical_yn(l, x)
    dh = spherical_jn(l, x, derivative=True) + 1j * spherical_yn(l, x, derivative=True)
    c = 1j ** (l + 1)
    return c * x * h, c * (h + x * dh)


def asymptotic_amplitudes(omega: float, rstar: float, psi: complex, dpsi: complex, l: int | None = None):
    if l is None:
        a_in = 0.5 * np.exp(1j * omega * rstar) * (psi + 1j * dpsi / omega)
        a_out = 0.5 * np.exp(-1j * omega * rstar) * (psi - 1j * dpsi / omega)
        return a_in, a_out
    x = omega * rstar
    up, dup = _riccati_hankel(l, x)
    um, dum = np.conj(up), np.conj(dup)
    # psi = a_in um + a_out up ; dpsi/omega = a_in dum + a_out dup
    mat = np.array([[um, up], [dum, dup]])
    a_in, a_out = np.linalg.solve(mat, np.array([psi, dpsi / omega]))
    return a_in, a_out
