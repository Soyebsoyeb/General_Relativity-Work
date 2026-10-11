"""Relativistic observables in Schwarzschild spacetime (G = c = 1, lengths in units of M).

Every observable is computed two independent ways where possible:

* ``*_quad``   a one-dimensional quadrature of the exact orbit equation, and
* ``*_geodesic`` integration of the full geodesic equation from the package's
  Christoffel symbols.
They are compared against the weak-field analytic expressions.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from ..spacetimes.schwarzschild import M as _Msym
from ..spacetimes.schwarzschild import schwarzschild
from .geodesics import GeodesicSystem, event_function

_SYSTEMS: dict = {}


def schwarzschild_system(mass: float = 1.0) -> GeodesicSystem:
    if mass not in _SYSTEMS:
        _SYSTEMS[mass] = GeodesicSystem(schwarzschild(), {_Msym: mass})
    return _SYSTEMS[mass]


# ---- light deflection ---------------------------------------------------------
def deflection_weak_field(b: float, mass: float = 1.0, order: int = 2) -> float:
    """4M/b (+ 15 pi/4 (M/b)^2): series of the exact deflection in M/b."""
    x = mass / b
    return 4 * x + (15 * np.pi / 4 * x**2 if order >= 2 else 0.0)


def deflection_quad(b: float, mass: float = 1.0) -> float:
    """Exact deflection of a null ray with impact parameter b > 3 sqrt(3) M.

    d phi / du = 1 / sqrt(1/b^2 - u^2 + 2 M u^3), u = 1/r. The turning point
    u0 is the smallest positive root of 1/b^2 - u^2 + 2 M u^3, and the
    substitution u = u0 (1 - s^2) removes the endpoint singularity.
    """
    roots = np.roots([2 * mass, -1.0, 0.0, 1.0 / b**2])
    u0 = min(r.real for r in roots if abs(r.imag) < 1e-12 and r.real > 0)

    def integrand(s):
        u = u0 * (1 - s * s)
        return 2 * u0 * s / np.sqrt(1 / b**2 - u * u + 2 * mass * u**3)

    half, _ = quad(integrand, 0.0, 1.0, epsabs=1e-14, epsrel=1e-13, limit=200)
    return 2 * half - np.pi


def deflection_geodesic(b: float, mass: float = 1.0, r0: float | None = None, rtol: float = 1e-12) -> float:
    """Deflection from integrating the null geodesic between radius r0 on both sides.

    The finite start radius truncates the asymptotic region; the resulting error
    grows with b/r0, so r0 defaults to max(1e4, 100 b).

    Subtracts the flat-space angle pi - 2 asin(b / r0) swept between the same radii.
    """
    r0 = max(1e4, 100 * b) if r0 is None else r0
    gs = schwarzschild_system(mass)
    f0 = 1 - 2 * mass / r0
    ut, uphi = 1.0 / f0, b / r0**2
    ur = -np.sqrt(1.0 - f0 * b**2 / r0**2)
    x0 = np.array([0.0, r0, np.pi / 2, 0.0])
    u0 = np.array([ut, ur, 0.0, uphi])
    ev = event_function(1, r0, direction=+1, terminal=True)
    tr = gs.integrate(x0, u0, 10 * r0, kappa=0.0, rtol=rtol, atol=1e-14, events=[ev])
    dphi = tr.x[-1, 3]
    return dphi - (np.pi - 2 * np.arcsin(b / r0))


# ---- periapsis precession ----------------------------------------------------------
def periapsis_advance_leading(p: float, mass: float = 1.0) -> float:
    """6 pi M / p, the leading weak-field result for a bound orbit with semi-latus rectum p."""
    return 6 * np.pi * mass / p


def _orbit_constants(p: float, e: float, mass: float):
    """E, L for a bound orbit with periapsis p/(1+e) and apoapsis p/(1-e).

    u = 1/r obeys (du/dphi)^2 = A + 2 M B u - u^2 + 2 M u^3 with A = (E^2-1)/L^2,
    B = 1/L^2; the two turning points fix A and B (linear equations).
    """
    up, ua = (1 + e) / p, (1 - e) / p
    mat = np.array([[1.0, 2 * mass * up], [1.0, 2 * mass * ua]])
    rhs = np.array([up**2 - 2 * mass * up**3, ua**2 - 2 * mass * ua**3])
    A, B = np.linalg.solve(mat, rhs)
    L = 1 / np.sqrt(B)
    E = np.sqrt(1 + A * L**2)
    return E, L, up, ua


def periapsis_advance_quad(p: float, e: float, mass: float = 1.0) -> float:
    """Exact advance per orbit.

    P(u) = 2M (u - ua)(u - up)(u - u3) with u3 = 1/(2M) - ua - up (sum of roots).
    With u = c + h sin(psi): dphi = d psi / sqrt(2M (u3 - u)).
    """
    E, L, up, ua = _orbit_constants(p, e, mass)
    u3 = 1 / (2 * mass) - ua - up
    c, h = (up + ua) / 2, (up - ua) / 2
    val, _ = quad(
        lambda ps: 1 / np.sqrt(2 * mass * (u3 - (c + h * np.sin(ps)))),
        -np.pi / 2,
        np.pi / 2,
        epsabs=1e-14,
        epsrel=1e-13,
    )
    return 2 * val - 2 * np.pi


def periapsis_advance_geodesic(p: float, e: float, mass: float = 1.0, n_orbits: int = 3, rtol: float = 1e-12):
    """Advance per orbit from the geodesic integrator; returns (mean, per-orbit array)."""
    gs = schwarzschild_system(mass)
    E, L, up, ua = _orbit_constants(p, e, mass)
    rp = 1 / up
    f = 1 - 2 * mass / rp
    x0 = np.array([0.0, rp, np.pi / 2, 0.0])
    u0 = np.array([E / f, 0.0, 0.0, L / rp**2])
    period_guess = 2 * np.pi * np.sqrt(((rp + 1 / ua) / 2) ** 3 / mass)
    ev = event_function(5, 0.0, direction=+1)  # u^r crosses zero upward: periapsis
    tr = gs.integrate(x0, u0, 1.5 * period_guess * (n_orbits + 1), kappa=1.0, rtol=rtol, atol=1e-13, events=[ev])
    phis = np.array([y[3] for y in tr.sol.y_events[0]])
    if len(phis) < 2:
        raise RuntimeError("fewer than two periapsis passages detected")
    adv = np.diff(phis) - 2 * np.pi
    return float(adv.mean()), adv


def newtonian_period(a: float, mass: float = 1.0) -> float:
    return 2 * np.pi * np.sqrt(a**3 / mass)


# ---- redshift and time delay ----------------------------------------------------------
def redshift_static_exact(r_emit: float, r_obs: float, mass: float = 1.0) -> float:
    """nu_obs / nu_emit for static emitter and observer: sqrt(f(r_emit)/f(r_obs))."""
    return float(np.sqrt((1 - 2 * mass / r_emit) / (1 - 2 * mass / r_obs)))


def redshift_weak_field(r_emit: float, r_obs: float, mass: float = 1.0) -> float:
    """1 - Delta Phi with Phi = -M/r."""
    return float(1 - (mass / r_emit - mass / r_obs))


def redshift_geodesic(r_emit: float, r_obs: float, mass: float = 1.0, rtol: float = 1e-12) -> float:
    """nu_obs / nu_emit from an outgoing radial photon, with nu = -g(u, u_static) at both ends.

    Uses the integrated u^t, so it tests conservation of E along the path.
    """
    gs = schwarzschild_system(mass)
    f_e = 1 - 2 * mass / r_emit
    x0 = np.array([0.0, r_emit, np.pi / 2, 0.0])
    u0 = np.array([1.0 / f_e, 1.0, 0.0, 0.0])  # null: f ut^2 = ur^2/f
    ev = event_function(1, r_obs, direction=+1, terminal=True)
    tr = gs.integrate(x0, u0, 10 * r_obs, kappa=0.0, rtol=rtol, atol=1e-14, events=[ev])
    f_o = 1 - 2 * mass / tr.x[-1, 1]
    nu_e = np.sqrt(f_e) * tr.u[0, 0]
    nu_o = np.sqrt(f_o) * tr.u[-1, 0]
    return float(nu_o / nu_e)


def shapiro_round_trip_excess(r1: float, r2: float, r0: float, mass: float = 1.0) -> float:
    """Exact coordinate-time excess of a round trip r1 -> r0 -> r2 -> r0 -> r1 over flat space.

    A null ray with closest approach r0 has b^2 = r0^2 / f(r0) and
    dt/dr = 1 / (f sqrt(1 - f b^2 / r^2)). With r = r0 + s^2 the endpoint
    singularity is removed. The excess is measured against the straight-line
    flat-space path with the same closest approach.
    """
    f0 = 1 - 2 * mass / r0
    b2 = r0**2 / f0

    def leg(rmax):
        smax = np.sqrt(rmax - r0)

        def integrand(s):
            r = r0 + s * s
            f = 1 - 2 * mass / r
            return 2 * s / (f * np.sqrt(max(1 - f * b2 / r**2, 1e-300)))

        v, _ = quad(integrand, 0.0, smax, epsabs=1e-12, epsrel=1e-12, limit=400)
        return v - np.sqrt(rmax**2 - r0**2)

    return 2 * (leg(r1) + leg(r2))


def shapiro_weak_field(r1: float, r2: float, r0: float, mass: float = 1.0) -> float:
    """Leading weak-field round-trip delay 4M [1 + ln(4 r1 r2 / r0^2)] (Shapiro 1964)."""
    return 4 * mass * (1 + np.log(4 * r1 * r2 / r0**2))
