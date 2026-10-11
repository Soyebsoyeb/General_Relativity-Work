"""Experiment 06: where does the Newtonian description break down?

Question
    How large is the relativistic correction to (i) the circular-orbit period
    and (ii) periapsis motion, as a function of M/r, and where does it reach 1%?
Equations
    Circular orbits: coordinate period T = 2 pi r^{3/2}/sqrt(M) exactly (Kepler
    form holds in Schwarzschild time); proper period tau = T sqrt(1 - 3M/r).
    Relative proper-period correction 1 - sqrt(1 - 3M/r) ~ 3M/(2r).
    Precession per orbit 6 pi M / p (weak field), computed exactly in exp03.
    Newtonian limit: g_00 = -(1 + 2 Phi), Phi = -M/r, d^2 x/dt^2 = -grad Phi.
Numerics
    Proper period obtained by integrating the circular geodesic through
    delta phi = 2 pi (event detection). DOP853 rtol 1e-12.
Interpretation
    Corrections are first order in M/r; the Newtonian result is accurate to 1%
    only for r above roughly a few hundred M for the period and for precession
    to be negligible over one orbit. Over many orbits the accumulated angular
    error grows linearly with orbit number.
Limitations
    Test-mass limit; no post-Newtonian two-body dynamics.
"""

import matplotlib.pyplot as plt
import numpy as np
from _common import save_fig, save_results

from gr_from_scratch.dynamics import observables as ob
from gr_from_scratch.dynamics.geodesics import event_function
from gr_from_scratch.visualization.style import apply_style

apply_style()
gs = ob.schwarzschild_system()
radii = np.logspace(np.log10(8), np.log10(2000), 10)
num, exact = [], []
for r in radii:
    ut = 1 / np.sqrt(1 - 3 / r)
    x0 = np.array([0, r, np.pi / 2, 0.0])
    u0 = np.array([ut, 0, 0, ut / r**1.5])
    ev = event_function(3, 2 * np.pi, direction=+1, terminal=True)
    tr = gs.integrate(x0, u0, 3 * 2 * np.pi * r**1.5, rtol=1e-12, atol=1e-13, events=[ev])
    tau = tr.sol.t_events[0][0]
    num.append(tau / (2 * np.pi * r**1.5) - 1)
    exact.append(np.sqrt(1 - 3 / r) - 1)
num, exact = np.array(num), np.array(exact)
slope = float(np.polyfit(np.log(radii), np.log(np.abs(num)), 1)[0])
r1pct = float(radii[np.argmin(np.abs(np.abs(num) - 0.01))])
r_one = float(1.5 / 0.01 * 1.0)  # leading-order estimate 3/(2r) = 0.01

p0, e0 = 15.0, 0.3
E, L, up, ua = ob._orbit_constants(p0, e0, 1.0)
rp = 1 / up
f = 1 - 2 / rp
a_sm = (rp + 1 / ua) / 2
T = 2 * np.pi * np.sqrt(a_sm**3)
tr = gs.integrate(
    np.array([0, rp, np.pi / 2, 0.0]), np.array([E / f, 0, 0, L / rp**2]), 3 * T, rtol=1e-11, atol=1e-13, max_step=1.0
)
phi = np.linspace(0, 2 * np.pi, 400)
rN = p0 / (1 + e0 * np.cos(phi))

fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
ax[0].loglog(radii, np.abs(num), "o", ms=4, label="geodesic integration")
ax[0].loglog(radii, np.abs(exact), "-", lw=0.8, label=r"$1-\sqrt{1-3M/r}$")
ax[0].loglog(radii, 1.5 / radii, "--", color="0.5", lw=0.8, label=r"$3M/2r$")
ax[0].axhline(0.01, color="C1", lw=0.6)
ax[0].set(xlabel=r"$r/M$", ylabel=r"$|\tau/\tau_N - 1|$")
ax[0].legend()
ax[1].plot(tr.x[:, 1] * np.cos(tr.x[:, 3]), tr.x[:, 1] * np.sin(tr.x[:, 3]), lw=0.6, label="geodesic, 3 orbits")
ax[1].plot(rN * np.cos(phi), rN * np.sin(phi), "--", color="C1", lw=0.9, label="Newtonian ellipse")
ax[1].plot([0], [0], "ko", ms=3)
ax[1].set_aspect("equal")
ax[1].set(xlabel=r"$x/M$", ylabel=r"$y/M$", title=rf"$p={p0:g}M,\ e={e0}$")
ax[1].legend(fontsize=7)
save_fig(fig, "exp06_newtonian_limit")
save_results(
    "exp06_newtonian_limit",
    {
        "r": radii,
        "period_correction_numeric": num,
        "period_correction_exact": exact,
        "max_abs_numeric_minus_exact": float(np.abs(num - exact).max()),
        "loglog_slope": slope,
        "leading_order_r_for_1pct": r_one,
    },
)
print("slope", slope, "max diff", np.abs(num - exact).max())
