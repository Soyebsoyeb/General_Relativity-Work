"""Experiment 01: curvature invariants of Schwarzschild and Kerr.

Question
    Does the package's Riemann tensor reproduce K = 48 M^2 / r^6 for
    Schwarzschild, and the independently known closed form for Kerr?
Equations
    K = R_{abcd} R^{abcd};  Schwarzschild K = 48 G^2 M^2 / (c^4 r^6).
    Kerr closed form: K = 48 M^2 (r^2 - a^2 c^2)(S^2 - 16 r^2 a^2 c^2)/S^6,
    S = r^2 + a^2 c^2, c = cos(theta).
Assumptions
    G = c = M = 1. Kerr in Boyer-Lindquist coordinates, a/M = 0.6, theta = pi/3.
Numerics
    Schwarzschild: exact symbolic contraction. Kerr: components evaluated at
    30 digits at each sample radius, contracted in double precision.
Tolerance
    Relative error reported; double precision contraction, so ~1e-12 expected.
Limitations
    The Kerr check is at sample points, not a symbolic proof.
"""

import numpy as np
import sympy as sp
from _common import save_fig, save_results

import matplotlib.pyplot as plt
from gr_from_scratch.geometry.connection import christoffel
from gr_from_scratch.geometry.invariants import kretschmann, kretschmann_at
from gr_from_scratch.geometry.riemann import riemann
from gr_from_scratch.spacetimes import kerr as K
from gr_from_scratch.spacetimes.schwarzschild import M, r, schwarzschild
from gr_from_scratch.visualization.curvature import loglog_slope
from gr_from_scratch.visualization.style import apply_style

apply_style()
Ks = sp.simplify(kretschmann(riemann(christoffel(schwarzschild()))))
assert sp.simplify(Ks - 48 * M**2 / r**6) == 0

rr = np.logspace(np.log10(2.05), 2, 60)
k_s = np.array([float(Ks.subs({M: 1, r: x})) for x in rr])
slope = loglog_slope(rr, k_s)

spin = sp.Rational(3, 5)
km = K.kerr().subs({K.M: 1, K.a: spin})
riem = riemann(christoffel(km), simplify=False)
th0 = sp.pi / 3
radii = [2.2, 3.0, 5.0, 8.0, 15.0]
num, ref = [], []
for x in radii:
    pt = {K.r: sp.Float(x, 30), K.theta: th0}
    num.append(kretschmann_at(riem, pt))
    ref.append(float(K.kretschmann_closed_form(1, spin).subs(pt)))
rel = [abs(a - b) / abs(b) for a, b in zip(num, ref, strict=True)]

fig, ax = plt.subplots(1, 2, figsize=(7, 2.8))
ax[0].loglog(rr, k_s, label="symbolic contraction")
ax[0].loglog(rr, 48 / rr**6, "--", color="0.5", lw=0.8, label=r"$48M^2/r^6$")
ax[0].axvline(2, color="0.6", lw=0.6)
ax[0].set(xlabel=r"$r/M$", ylabel=r"$K M^4$", title="Schwarzschild")
ax[0].legend()
rk = np.logspace(np.log10(1.9), 1.3, 200)
ax[1].loglog(
    rk,
    [abs(float(K.kretschmann_closed_form(1, spin).subs({K.r: x, K.theta: th0}))) for x in rk],
    color="0.5",
    label="closed form",
)
ax[1].loglog(radii, np.abs(num), "o", ms=4, color="C1", label="engine, point values")
ax[1].set(xlabel=r"$r/M$", ylabel=r"$|K| M^4$", title=r"Kerr, $a/M=0.6$, $\theta=\pi/3$")
ax[1].legend()
save_fig(fig, "exp01_curvature_invariants")
save_results(
    "exp01_curvature_invariants",
    {
        "schwarzschild_loglog_slope": slope,
        "kerr_radii": radii,
        "kerr_engine": num,
        "kerr_closed_form": ref,
        "kerr_max_relative_error": max(rel),
    },
)
print("slope", slope, "kerr max rel err", max(rel))
