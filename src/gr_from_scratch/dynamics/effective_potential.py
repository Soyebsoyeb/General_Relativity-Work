"""Radial motion and effective potentials for static spherically symmetric metrics.

For equatorial motion with conserved E = -g_tt u^t and L = g_pp u^phi, the
normalisation g(u,u) = -kappa gives

    g_rr rdot^2 = -kappa + E^2 / (-g_tt) - L^2 / g_pp.

If g_tt g_rr = -1 (Schwarzschild, Reissner-Nordstrom) this is
rdot^2 = E^2 - V(r) with V = -g_tt (kappa + L^2 / g_pp).
"""

from __future__ import annotations

import sympy as sp

from ..core.metric import Metric

E, L = sp.symbols("E L", real=True)


def radial_equation(metric: Metric, kappa=1):
    """Return (rdot_squared, V_eff) for equatorial motion, assuming g_tt g_rr = -1."""
    g = metric.g
    gtt, grr = g[0, 0], g[1, 1]
    gpp = g[3, 3].subs(metric.coords[2], sp.pi / 2)
    if sp.simplify(gtt * grr + 1) != 0:
        raise ValueError("effective-potential form needs g_tt g_rr = -1")
    rdot2 = sp.simplify((-kappa + E**2 / (-gtt) - L**2 / gpp) / grr)
    V = sp.simplify(E**2 - rdot2)
    return rdot2, V


def circular_orbit_analysis(metric: Metric, kappa=1):
    """Circular-orbit conditions V' = 0, stability V'' and the ISCO V' = V'' = 0.

    Returns a dict with L^2(r) for circular orbits, E^2(r), and the list of
    radii satisfying the marginal-stability condition.
    """
    r = metric.coords[1]
    _, V = radial_equation(metric, kappa)
    dV = sp.diff(V, r)
    d2V = sp.diff(V, r, 2)
    L2_sol = sp.solve(dV, L**2) if kappa != 0 else []
    out = {"V": V}
    if kappa == 1:
        L2 = sp.simplify(L2_sol[0])
        E2 = sp.simplify(V.subs(L**2, L2))
        stability = sp.simplify(d2V.subs(L**2, L2))
        out.update(L2=L2, E2=E2, stability=stability, isco=sp.solve(sp.numer(sp.together(stability)), r))
    else:
        # null: V = -g_tt L^2/g_pp is proportional to L^2; extremise the ratio
        ratio = sp.simplify(V / L**2)
        out["photon_sphere"] = sp.solve(sp.diff(ratio, r), r)
    return out
