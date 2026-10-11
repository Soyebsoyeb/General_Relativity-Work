"""Experiment 10: massless scalar field scattering off a Schwarzschild black hole.

Classification
    Classical field on a fixed curved background (not quantised gravity).
Question
    What are the transmission coefficients T_l(w), does flux conservation hold
    numerically, and does the l = 0 result approach the horizon-area limit?
Equations
    d^2 psi/dr*^2 + (w^2 - V) psi = 0,  V = f (l(l+1)/r^2 + 2M/r^3),  f = 1 - 2M/r.
    Boundary conditions: psi = exp(-i w r*) at r*_min = -80 M (ingoing at the
    horizon); at r*_max psi is matched to A_in u^- + A_out u^+, the Riccati-Hankel
    solutions of the centrifugal problem (numerical/boundary_conditions.py).
    T = 1/|A_in|^2, R = |A_out|^2/|A_in|^2, flux check 1 - T - R (relative form;
    the absolute form |A_in|^2 - |A_out|^2 - 1 is dominated by round-off when T << 1).
    Low-frequency limit, l = 0: sigma_abs = pi T/w^2 -> 16 pi M^2 (horizon area),
    i.e. T/(16 M^2 w^2) -> 1.
Parameters
    M = 1; w in [0.05, 1.0]; r*_max = max(3000, 80/w).
Numerics
    DOP853 in r* with complex state, rtol 1e-10; r(r*) by Lambert W.
Tolerances
    Flux residual reported; T/(16 w^2) extrapolation reported at several w.
Limitations
    Matching at finite r*_max is quantified through the flux residual and by
    repeating one point with doubled r*_max. No Hawking-emission claim is made.
"""

import matplotlib.pyplot as plt
import numpy as np
from _common import save_fig, save_results

from gr_from_scratch.quantum.klein_gordon import scattering
from gr_from_scratch.visualization.style import apply_style

apply_style()
ws = np.linspace(0.05, 1.0, 24)
Tl, flux = {}, []
for l in range(4):
    vals = []
    for w in ws:
        s = scattering(l, float(w), rstar_max=max(3000.0, 80 / w))
        vals.append(s.transmission)
        flux.append(abs(s.flux_residual))
    Tl[l] = np.array(vals)
sigma = np.pi / ws**2 * sum((2 * l + 1) * Tl[l] for l in range(4))

low = [0.005, 0.01, 0.02, 0.05]
ratio = []
for w in low:
    s = scattering(0, w, rstar_max=max(3000.0, 120 / w))
    ratio.append(s.transmission / (16 * w * w))
s1 = scattering(2, 0.3, rstar_max=3000.0)
s2 = scattering(2, 0.3, rstar_max=6000.0)

fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.7))
for l in range(4):
    ax[0].plot(ws, Tl[l], label=rf"$\ell={l}$")
ax[0].set(xlabel=r"$M\omega$", ylabel=r"$T_\ell$")
ax[0].legend()
ax[1].plot(ws, sigma / (16 * np.pi), label=r"$\ell\le3$ sum")
ax[1].axhline(1, color="0.5", ls="--", lw=0.7)
ax[1].set(xlabel=r"$M\omega$", ylabel=r"$\sigma_{\rm abs}/(16\pi M^2)$", title="truncated at $\\ell=3$")
ax[2].plot(low, ratio, "o-", ms=4)
ax[2].axhline(1, color="0.5", ls="--", lw=0.7)
ax[2].set(xlabel=r"$M\omega$", ylabel=r"$T_0/(16M^2\omega^2)$", title=r"$\ell=0$ low-frequency limit")
save_fig(fig, "exp10_scalar_scattering")
save_results(
    "exp10_scalar_scattering",
    {
        "omega": ws,
        "T": {str(k): v for k, v in Tl.items()},
        "max_flux_residual": max(flux),
        "low_omega": low,
        "T0_over_16w2": ratio,
        "rstar_max_sensitivity_l2_w0.3": abs(s1.transmission - s2.transmission),
    },
)
print("flux", max(flux), "ratio", ratio, "sens", abs(s1.transmission - s2.transmission))
