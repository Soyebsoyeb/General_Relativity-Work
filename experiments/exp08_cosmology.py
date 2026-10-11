"""Experiment 08: FLRW expansion against analytic solutions.

Question
    Does numerical integration of the second-order Friedmann (acceleration)
    equation reproduce known closed-form a(t), and does the first-order Friedmann
    constraint stay satisfied?
Models (H0 = 1, a(t0) = 1 normalisation; matching initial data from the analytic solution)
    matter-dominated flat:      a = (3t/2)^{2/3}
    radiation-dominated flat:   a = (2t)^{1/2}
    de Sitter (Om_L = 1):       a = exp(t)
    flat LCDM (Om_m = 0.3):     a = (Om_m/Om_L)^{1/3} sinh^{2/3}(3 sqrt(Om_L) t / 2)
    closed dust (Om_m = 2):     cycloid, a = (a_max/2)(1 - cos eta), parametric in eta
Numerics
    DOP853 rtol 1e-11 for (a, adot). Constraint residual (adot/a)^2 - E^2(a).
Interpretation
    Relative error is bounded by the integrator tolerance; constraint
    residual is conserved up to accumulated round-off, as the Bianchi identity
    guarantees consistency of the two Friedmann equations with continuity.
Limitations
    Homogeneous, isotropic background only; no perturbations; the closed model is
    shown only up to eta = 2.6 (before recollapse); the analytic reference inverts
    t(eta) by root finding.
"""

import matplotlib.pyplot as plt
import numpy as np
from _common import save_fig, save_results

from gr_from_scratch.numerical import friedmann as fr
from gr_from_scratch.visualization.style import apply_style

apply_style()


def run(cosmo, a_fn, t0, t1):
    a0 = a_fn(t0)
    h = 1e-6 * max(t0, 1.0)
    ad0 = (a_fn(t0 + h) - a_fn(t0 - h)) / (2 * h)
    t, a, ad, con = fr.integrate(cosmo, a0, ad0, t1 - t0, n=400)
    t = t + t0
    ref = a_fn(t)
    return t, a, ref, con


cases = {}
cases["matter (flat)"] = run(fr.Cosmology(om_m=1.0), fr.a_matter_flat, 0.2, 3.0)
cases["radiation (flat)"] = run(fr.Cosmology(om_r=1.0), fr.a_radiation_flat, 0.05, 1.5)
cases["de Sitter"] = run(fr.Cosmology(om_l=1.0), lambda t: fr.a_de_sitter(t), 0.0, 3.0)
cases[r"flat $\Lambda$CDM"] = run(fr.Cosmology(om_m=0.3, om_l=0.7), lambda t: fr.a_lcdm_flat(t, 0.3), 0.1, 2.5)

om = 2.0
eta = np.linspace(0.3, 2.6, 400)
tc, ac = fr.closed_dust(eta, om)
cc = fr.Cosmology(om_m=om)
ad0 = (ac[1] - ac[0]) / (tc[1] - tc[0])
# exact adot from adot = a'(eta)/t'(eta)
ok = abs(1 - om)
amax = om / ok
da = 0.5 * amax * np.sin(eta)
dt = om / (2 * ok**1.5) * (1 - np.cos(eta))
tn, an, adn, conn = fr.integrate(cc, ac[0], (da / dt)[0], tc[-1] - tc[0], n=400)
tn += tc[0]
from scipy.optimize import brentq

eta_n = np.array(
    [brentq(lambda x, tt=tt: fr.closed_dust(x, om)[0] - tt, 1e-9, 2 * np.pi - 1e-9, xtol=1e-14) for tt in tn]
)
ref_closed = fr.closed_dust(eta_n, om)[1]
cases[r"closed dust ($\Omega_m=2$)"] = (tn, an, ref_closed, conn)

fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.7))
summary = {}
for name, (t, a, ref, con) in cases.items():
    rel = np.abs(a - ref) / ref
    summary[name] = {"max_rel_error": float(rel.max()), "max_constraint": float(np.abs(con).max())}
    ax[0].plot(t, a, label=name, lw=0.9)
    ax[1].semilogy(t, rel + 1e-16, lw=0.9)
    ax[2].semilogy(t, np.abs(con) + 1e-18, lw=0.9)
ax[0].set(xlabel=r"$H_0 t$", ylabel=r"$a(t)$", yscale="log")
ax[0].legend(fontsize=6)
ax[1].set(xlabel=r"$H_0 t$", ylabel="relative error vs analytic")
ax[2].set(xlabel=r"$H_0 t$", ylabel=r"$|(\dot a/a)^2-E^2(a)|$")
save_fig(fig, "exp08_cosmology")
save_results("exp08_cosmology", summary)
for k, v in summary.items():
    print(k, v)
