"""Experiment 04: deflection of light by a Schwarzschild mass.

Question
    How accurate are the first- and second-order weak-field deflection formulas,
    and does geodesic integration agree with exact quadrature?
Equations
    Exact: alpha = 2 int_0^{u0} du / sqrt(1/b^2 - u^2 + 2 M u^3) - pi.
    Weak field: alpha = 4M/b + (15 pi / 4)(M/b)^2 + O((M/b)^3).
Parameters
    M = 1, b from 6 to 2000 (critical b_c = 3 sqrt(3) M ~ 5.196).
Numerics
    Geodesic: DOP853 rtol 1e-12, started and stopped at r0 = max(1e4, 100 b)
    with the flat-space angle pi - 2 asin(b/r0) subtracted. Quadrature: scipy.quad.
Interpretation
    Residual of the first-order formula scales as (M/b)^2 and of the second
    as (M/b)^3; the geodesic differs from quadrature by a finite-start-radius
    truncation error, which falls rapidly as r0 grows (see results JSON for
    the measured maximum relative difference).
Limitations
    Point-mass lens, test rays, no plasma, no lens-plane approximation.
"""

import matplotlib.pyplot as plt
import numpy as np
from _common import save_fig, save_results

from gr_from_scratch.dynamics import observables as ob
from gr_from_scratch.visualization.geodesics import plot_equatorial_paths
from gr_from_scratch.visualization.style import apply_style

apply_style()
bs = np.logspace(np.log10(6), np.log10(2000), 14)
ex = np.array([ob.deflection_quad(b) for b in bs])
geo = np.array([ob.deflection_geodesic(b) for b in bs[::2]])
w1 = np.array([ob.deflection_weak_field(b, order=1) for b in bs])
w2 = np.array([ob.deflection_weak_field(b, order=2) for b in bs])
geo_err = np.abs(geo - ex[::2]) / ex[::2]
big = bs > 30
s1 = float(np.polyfit(np.log(bs[big]), np.log(np.abs(ex - w1)[big]), 1)[0])
s2 = float(np.polyfit(np.log(bs[big]), np.log(np.abs(ex - w2)[big]), 1)[0])

gs = ob.schwarzschild_system()
paths = []
for b in (5.5, 6.5, 8, 11, 16):
    r0 = 60.0
    f0 = 1 - 2 / r0
    x0 = np.array([0, r0, np.pi / 2, 0.0])
    u0 = np.array([1 / f0, -np.sqrt(1 - f0 * b**2 / r0**2), 0, b / r0**2])
    ev = ob.event_function(1, r0, direction=+1, terminal=True)
    paths.append(gs.integrate(x0, u0, 400, kappa=0, rtol=1e-11, atol=1e-13, events=[ev], max_step=0.5).x)

fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.7), gridspec_kw={"width_ratios": [1.1, 1, 1]})
plot_equatorial_paths(ax[0], paths, lw=0.8)
ax[0].set_title("null rays, $r_0=60M$")
ax[1].loglog(bs, ex, "-", label="exact")
ax[1].loglog(bs, w1, "--", color="0.5", label=r"$4M/b$")
ax[1].loglog(bs, w2, ":", color="C1", label="+ second order")
ax[1].loglog(bs[::2], geo, "o", ms=3, mfc="none", color="C2", label="geodesic")
ax[1].set(xlabel=r"$b/M$", ylabel=r"$\alpha$ [rad]")
ax[1].legend()
ax[2].loglog(bs, np.abs(ex - w1) / ex, label="1st order")
ax[2].loglog(bs, np.abs(ex - w2) / ex, label="2nd order")
ax[2].loglog(bs[::2], geo_err + 1e-16, "o", ms=3, mfc="none", label="geodesic vs exact")
ax[2].set(xlabel=r"$b/M$", ylabel="relative error")
ax[2].legend()
save_fig(fig, "exp04_light_deflection")
save_results(
    "exp04_light_deflection",
    {
        "b": bs,
        "exact": ex,
        "geodesic_b": bs[::2],
        "geodesic": geo,
        "max_geodesic_vs_exact_rel": float(geo_err.max()),
        "slope_err_first_order": s1,
        "slope_err_second_order": s2,
    },
)
print("geo err", geo_err.max(), "slopes", s1, s2)
