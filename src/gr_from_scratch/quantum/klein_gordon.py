"""Scalar field on a fixed curved background: (box - m^2 - xi R) phi = 0.

Classification: a *classical* field on curved spacetime. The scattering
numbers below are classical transmission coefficients (greybody factors).
They are not statements about quantised gravity; by themselves they say
nothing about Hawking radiation spectra unless combined with the usual
thermal-emission argument, which is outside this module.

Derived result (checked symbolically in the tests). With
phi = exp(-i w t) P(theta) psi(r) / r on Schwarzschild, where
P'' + cot(theta) P' = -l(l+1) P, the equation reduces to

    d^2 psi / dr*^2 + (w^2 - V) psi = 0,   V = f ( l(l+1)/r^2 + 2M/r^3 + m^2 + xi R ),

with R = 0 in vacuum and dr*/dr = 1/f.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.special import lambertw

from ..core.metric import Metric
from ..numerical.boundary_conditions import asymptotic_amplitudes, ingoing_horizon_data


def box_scalar(metric: Metric, phi, sqrt_g=None):
    """(1/sqrt|g|) d_mu ( sqrt|g| g^{mu nu} d_nu phi )."""
    sg = sqrt_g if sqrt_g is not None else metric.sqrt_minus_g
    n = metric.dim
    tot = sum(
        sp.diff(sg * metric.ginv[a, b] * sp.diff(phi, metric.coords[b]), metric.coords[a])
        for a in range(n)
        for b in range(n)
    )
    return tot / sg


def schwarzschild_radial_reduction():
    """Return the residual of the claimed reduction to the radial equation (should be 0)."""
    from ..spacetimes.schwarzschild import M, f, r, schwarzschild, t, theta

    w, lam, m2 = sp.symbols("omega lambda m2", positive=True)
    psi = sp.Function("psi")(r)
    P = sp.Function("P")(theta)
    Phi = sp.exp(-sp.I * w * t) * P * psi / r
    met = schwarzschild()
    op = box_scalar(met, Phi, sqrt_g=r**2 * sp.sin(theta)) - m2 * Phi
    # angular equation: P'' = -cot(theta) P' - lambda P, with lambda = l(l+1)
    op = op.subs(sp.Derivative(P, (theta, 2)), -sp.cot(theta) * sp.diff(P, theta) - lam * P).doit()
    ff = f(M, r)
    V = ff * (lam / r**2 + 2 * M / r**3 + m2)
    claimed = ff * sp.diff(ff * sp.diff(psi, r), r) + (w**2 - V) * psi
    lhs = sp.simplify(op * r * ff / (sp.exp(-sp.I * w * t) * P))
    return sp.simplify(lhs - claimed)


def radius_from_tortoise(rstar, mass: float = 1.0):
    """Invert r* = r + 2M ln(r/2M - 1) for r > 2M."""
    z = rstar / (2 * mass) - 1.0
    if z < 20:
        x = float(np.real(lambertw(np.exp(z))))
    else:  # x = z - ln x, fixed-point iteration
        x = z
        for _ in range(12):
            x = z - np.log(x)
    return 2 * mass * (1 + x)


def potential(r: float, l: int, mass: float = 1.0, m2: float = 0.0) -> float:
    f = 1 - 2 * mass / r
    return f * (l * (l + 1) / r**2 + 2 * mass / r**3 + m2)


@dataclass
class Scattering:
    l: int
    omega: float
    transmission: float
    reflection: float
    flux_residual: float  # 1 - T - R = (|A_in|^2 - |A_out|^2 - 1)/|A_in|^2, should be ~0
    a_in: complex
    a_out: complex


def scattering(
    l: int, omega: float, mass: float = 1.0, rstar_min: float = -80.0, rstar_max: float = 3000.0, rtol: float = 1e-10
) -> Scattering:
    """Transmission and reflection for a massless scalar wave, unit ingoing amplitude at the horizon."""

    def rhs(rs, y):
        r = radius_from_tortoise(rs, mass)
        return [y[1], (potential(r, l, mass) - omega**2) * y[0]]

    psi0, dpsi0 = ingoing_horizon_data(omega, rstar_min)
    sol = solve_ivp(rhs, (rstar_min, rstar_max), [psi0, dpsi0], method="DOP853", rtol=rtol, atol=1e-13)
    a_in, a_out = asymptotic_amplitudes(omega, rstar_max, sol.y[0, -1], sol.y[1, -1], l=l)
    T = 1.0 / abs(a_in) ** 2
    Rr = abs(a_out) ** 2 / abs(a_in) ** 2
    return Scattering(l, omega, T, Rr, 1.0 - T - Rr, a_in, a_out)
