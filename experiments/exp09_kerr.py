"""Experiment 09: Kerr geometry, frame dragging and the ISCO.

Question
    (a) Where are the horizons and ergosurface? (b) Do the literature
    circular-orbit formulas satisfy the geodesic equation built from this
    package's own metric? (c) Is the literature ISCO radius the point where
    the radial effective potential becomes marginally stable?
Equations
    Boyer-Lindquist metric (see spacetimes/kerr.py), M = 1.
    r_pm = M +- sqrt(M^2 - a^2);  ergosurface r_E = M + sqrt(M^2 - a^2 cos^2 theta).
    Frame dragging: omega(r) = -g_{t phi}/g_{phi phi} (ZAMO angular velocity).
    Radial equation, equatorial: g^{rr} u_r^2 = W(r;E,L) with
        W = -1 - g^{tt} E^2 + 2 g^{t phi} E L - g^{phi phi} L^2.
    Circular: W = W' = 0; stable if W'' < 0; ISCO: W'' = 0 along the family of
    circular orbits (E(r), L(r)) from the closed forms.
Numerics
    Residuals evaluated in double precision from the engine's inverse metric
    and Christoffel symbols; ISCO by brentq on W''.
Limitations
    The metric is entered, not derived; vacuum character is tested in
    tests/test_kerr.py at sample points only.
"""

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from _common import save_fig, save_results
from scipy.optimize import brentq

from gr_from_scratch.dynamics.geodesics import GeodesicSystem
from gr_from_scratch.spacetimes import kerr as K
from gr_from_scratch.visualization.spacetime_diagrams import kerr_sections
from gr_from_scratch.visualization.style import apply_style

apply_style()
spins = [0.0, 0.5, 0.9, 0.998]

# W(r; E, L, a) from the engine's inverse metric at theta = pi/2
km = K.kerr()
Es, Ls = sp.symbols("E L", real=True)
gi = km.ginv.subs({K.theta: sp.pi / 2, K.M: 1})
W = -1 - gi[0, 0] * Es**2 + 2 * gi[0, 3] * Es * Ls - gi[3, 3] * Ls**2
Wf = sp.lambdify((K.r, K.a, Es, Ls), W, "numpy")
W1f = sp.lambdify((K.r, K.a, Es, Ls), sp.diff(W, K.r), "numpy")
W2f = sp.lambdify((K.r, K.a, Es, Ls), sp.diff(W, K.r, 2), "numpy")

rows = []
for a in spins:
    for pro in (True, False):
        r_lit = K.isco_radius(a, 1.0, pro)

        def w2(r, a=a, pro=pro):
            E, L, _ = K.circular_orbit(r, a, 1.0, pro)
            return W2f(r, a, E, L)

        # narrow bracket around the literature value; a sign change confirms a root
        r_num = brentq(w2, 0.97 * r_lit, 1.03 * r_lit, xtol=1e-14)
        r_test = 1.5 * r_lit
        E, L, _ = K.circular_orbit(r_test, a, 1.0, pro)
        rows.append(
            {
                "a": a,
                "prograde": pro,
                "isco_literature": r_lit,
                "isco_from_W2": r_num,
                "W_residual": float(Wf(r_test, a, E, L)),
                "Wprime_residual": float(W1f(r_test, a, E, L)),
            }
        )

# geodesic-equation residual at a sample radius (independent route: Christoffels)
gs = GeodesicSystem(K.kerr().subs({K.M: 1, K.a: sp.Rational(3, 5)}))
r0 = 7.0
E, L, Om = K.circular_orbit(r0, 0.6, 1.0, True)
x0 = np.array([0, r0, np.pi / 2, 0.0])
g = gs.g(x0)
ut = 1 / np.sqrt(-(g[0, 0] + 2 * g[0, 3] * Om + g[3, 3] * Om**2))
u = np.array([ut, 0, 0, Om * ut])
acc_r = abs(gs.acceleration(x0, u)[1])
E_eng = -(g[0, 0] * u[0] + g[0, 3] * u[3])
L_eng = g[3, 0] * u[0] + g[3, 3] * u[3]

fig, ax = plt.subplots(1, 3, figsize=(8.2, 3.1))
for a, c in zip(spins[1:], ["C0", "C1", "C2"], strict=True):
    kerr_sections(ax[0], a, color=c, label=f"a={a}")
ax[0].legend(fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.22), columnspacing=1.0, handlelength=1.2)
ax[0].set_title("horizons (solid), ergosurface (dashed)", fontsize=8)
rr = np.linspace(1.6, 10, 300)
for a in [0.5, 0.9]:
    ge = km.g.subs({K.theta: sp.pi / 2, K.M: 1, K.a: a})
    om = sp.lambdify(K.r, -ge[0, 3] / ge[3, 3], "numpy")
    ax[1].plot(rr, om(rr), label=f"a={a}")
ax[1].set(xlabel=r"$r/M$", ylabel=r"$\omega_{\rm ZAMO}\,M$", title="frame dragging, equator")
ax[1].legend()
aa = np.linspace(0, 0.998, 200)
ax[2].plot(aa, [K.isco_radius(x, 1, True) for x in aa], label="prograde")
ax[2].plot(aa, [K.isco_radius(x, 1, False) for x in aa], label="retrograde")
for row in rows:
    ax[2].plot(row["a"], row["isco_from_W2"], "k.", ms=4)
ax[2].set(xlabel=r"$a/M$", ylabel=r"$r_{\rm ISCO}/M$", title="ISCO (dots: from $W''=0$)")
ax[2].legend()
save_fig(fig, "exp09_kerr")
bad = max(abs(r["isco_literature"] - r["isco_from_W2"]) for r in rows if np.isfinite(r["isco_from_W2"]))
save_results(
    "exp09_kerr",
    {
        "isco_table": rows,
        "max_isco_discrepancy": bad,
        "circular_geodesic_radial_acceleration": acc_r,
        "E_literature_minus_engine": E - E_eng,
        "L_literature_minus_engine": L - L_eng,
    },
)
print("max ISCO discrepancy", bad, "acc_r", acc_r, E - E_eng, L - L_eng)
for r in rows:
    print(r)
