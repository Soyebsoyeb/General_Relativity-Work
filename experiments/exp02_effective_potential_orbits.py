"""Experiment 02: effective potentials, photon sphere and ISCO.

Question
    Do the circular-orbit radii derived symbolically from the effective potential
    (photon sphere, ISCO) agree with the textbook values 3M and 6M, and what
    does the potential landscape look like for different angular momenta?
Equations
    Timelike: rdot^2 = E^2 - V, V = (1 - 2M/r)(1 + L^2/r^2).
    Null:     rdot^2 = E^2 - L^2 (1 - 2M/r)/r^2.
    Circular: V' = 0;  ISCO: V' = V'' = 0.
Assumptions
    Equatorial geodesics, G = c = 1, M = 1.
Numerics
    Symbolic differentiation and solving; no integration.
Limitations
    Exact; the only approximation is plotting resolution.
"""

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from _common import save_fig, save_results

from gr_from_scratch.dynamics.effective_potential import circular_orbit_analysis
from gr_from_scratch.spacetimes.schwarzschild import M, schwarzschild
from gr_from_scratch.visualization.effective_potentials import plot_timelike, v_null
from gr_from_scratch.visualization.style import apply_style

apply_style()
met = schwarzschild()
tl = circular_orbit_analysis(met, 1)
nl = circular_orbit_analysis(met, 0)
isco, photon = tl["isco"][0], nl["photon_sphere"][0]
assert isco == 6 * M and photon == 3 * M

rr = np.linspace(2.02, 24, 600)
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
plot_timelike(ax[0], rr, [3.2, 3.464, 4.0, 5.0])
ax[0].set_ylim(0.85, 1.12)
ax[0].set_title("timelike")
ax[1].plot(rr, v_null(rr))
ax[1].axvline(3, color="0.5", ls="--", lw=0.7)
ax[1].set(xlabel=r"$r/M$", ylabel=r"$V/L^2$", title="null (photon sphere at $3M$)")
rc = np.linspace(6, 40, 200)
E2 = sp.lambdify(tl["E2"].free_symbols - {M} and list(tl["E2"].free_symbols), tl["E2"], "numpy")
sym = [s for s in tl["E2"].free_symbols if s.name == "r"][0]
E2f = sp.lambdify(sym, tl["E2"].subs(M, 1), "numpy")
ax[2].plot(rc, 1 - np.sqrt(E2f(rc)))
ax[2].set(xlabel=r"$r/M$", ylabel=r"$1-E/\mu$", title="binding energy")
ax[2].axhline(1 - np.sqrt(8 / 9), color="0.5", ls=":", lw=0.7)
save_fig(fig, "exp02_effective_potential_orbits")
save_results(
    "exp02_effective_potential_orbits",
    {
        "isco_symbolic": str(isco),
        "photon_sphere_symbolic": str(photon),
        "binding_energy_at_isco": float(1 - sp.sqrt(tl["E2"].subs({M: 1, sym: 6}))),
    },
)
print(isco, photon)
