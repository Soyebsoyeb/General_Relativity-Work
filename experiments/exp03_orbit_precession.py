"""Experiment 03: perihelion precession of bound Schwarzschild orbits.

Question
    How does the periapsis advance per orbit from the full geodesic equation
    compare with an independent exact quadrature and with 6 pi M / p?
Equations
    (du/dphi)^2 = A + 2 M B u - u^2 + 2 M u^3, u = 1/r, turning points fix A, B.
    Exact advance = 2 int dphi - 2 pi (sine substitution, see observables.py).
    Weak field: delta phi = 6 pi M / p, p = semi-latus rectum.
Parameters
    M = 1, eccentricity e = 0.3 for the scan, p in [20, 2000].
Numerics
    DOP853, rtol = 1e-12, atol = 1e-13, three orbits, periapsis detected by an
    event on u^r. Quadrature: scipy.quad, rel. tol. 1e-13.
Interpretation
    Geodesic and quadrature should agree to integrator tolerance; both differ
    from 6 pi M / p by a relative amount O(M/p) (absolute O((M/p)^2)), which
    shows up as a slope of about -1 in the log-log relative deviation.
Limitations
    Test-particle motion only; no self-force, no spin, no radiation reaction.
"""

import matplotlib.pyplot as plt
import numpy as np
from _common import save_fig, save_results

from gr_from_scratch.dynamics import observables as ob
from gr_from_scratch.visualization.geodesics import plot_equatorial_paths
from gr_from_scratch.visualization.style import apply_style

apply_style()
e = 0.3
ps = np.logspace(np.log10(20), np.log10(2000), 9)
quad_v = np.array([ob.periapsis_advance_quad(p, e) for p in ps])
geo_v = np.array([ob.periapsis_advance_geodesic(p, e, n_orbits=2)[0] for p in ps])
lead = np.array([ob.periapsis_advance_leading(p) for p in ps])
agree = np.abs(geo_v - quad_v) / quad_v
dev = np.abs(quad_v - lead) / lead
slope = float(np.polyfit(np.log(ps), np.log(dev), 1)[0])

# rosette
p0, e0 = 12.0, 0.4
gs = ob.schwarzschild_system()
E, L, up, ua = ob._orbit_constants(p0, e0, 1.0)
rp = 1 / up
f = 1 - 2 / rp
tr = gs.integrate(
    np.array([0, rp, np.pi / 2, 0.0]),
    np.array([E / f, 0, 0, L / rp**2]),
    5 * 2 * np.pi * np.sqrt(((rp + 1 / ua) / 2) ** 3),
    rtol=1e-11,
    atol=1e-13,
    max_step=2.0,
)

fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
plot_equatorial_paths(ax[0], [tr.x], lw=0.7)
ax[0].set_title(rf"$p={p0:g}M,\ e={e0}$")
ax[1].loglog(ps, quad_v, "-", label="exact quadrature")
ax[1].loglog(ps, geo_v, "o", ms=4, mfc="none", label="geodesic integration")
ax[1].loglog(ps, lead, "--", color="0.5", label=r"$6\pi M/p$")
ax[1].set(xlabel=r"$p/M$", ylabel=r"advance per orbit [rad]")
ax[1].legend()
save_fig(fig, "exp03_orbit_precession")
save_results(
    "exp03_orbit_precession",
    {
        "e": e,
        "p": ps,
        "quad": quad_v,
        "geodesic": geo_v,
        "leading": lead,
        "max_geodesic_vs_quad_rel": float(agree.max()),
        "slope_of_deviation_from_leading": slope,
    },
)
print("agreement", agree.max(), "dev slope", slope)
