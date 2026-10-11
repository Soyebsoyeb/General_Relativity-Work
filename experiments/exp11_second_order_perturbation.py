"""Experiment 11: first- versus second-order curvature of a plane gravitational wave.

Question
    For g = eta + eps h with h a TT plane wave, the linearised Einstein tensor
    vanishes. What does the second-order truncation give, and how does it
    compare with the exact Einstein tensor of the same metric?
Metric
    ds^2 = -dt^2 + dz^2 + (1 + eps f) dx^2 + (1 - eps f) dy^2,
    f = A sin(w (t - z)). This family is an exact metric (Rosen-type), not a
    vacuum solution beyond first order.
Truncation
    Orders 1 and 2 in eps (explicitly stated; nothing beyond eps^2 is kept).
Numerics
    Symbolic: perturbative_geometry(order=2); exact G from the full engine;
    point evaluation at t = 0.3, z = 0.1, A = w = 1.
Interpretation
    G^(1) = 0 for a vacuum wave, G^(2) is non-zero: the wave sources curvature
    at second order. The truncation error of G_tt scales as eps^4, not eps^3: eps -> -eps
    is the relabelling x <-> y, so G_tt is even in eps and the odd order vanishes.
Limitations
    This is a series in eps, not a nonlinear evolution. The exact metric used
    for comparison is a specific closed family.
"""

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from _common import save_fig, save_results

from gr_from_scratch.core.metric import Metric
from gr_from_scratch.field_equations.perturbation_series import eps, order_coefficient, perturbative_geometry
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry
from gr_from_scratch.visualization.style import apply_style

apply_style()
t, x, y, z = sp.symbols("t x y z", real=True)
f = sp.sin(t - z)
h = sp.zeros(4, 4)
h[1, 1], h[2, 2] = f, -f
geo = perturbative_geometry(h, [t, x, y, z], order=2)
G1 = geo["G"].applyfunc(lambda e: sp.simplify(order_coefficient(e, 1)))
G2 = geo["G"].applyfunc(lambda e: sp.simplify(order_coefficient(e, 2)))

exact = Metric([t, x, y, z], sp.diag(-1, 1 + eps * f, 1 - eps * f, 1), simplifier=sp.simplify)
Gx = compute_geometry(exact).einstein.c
pt = {t: 0.3, z: 0.1}
comp = (0, 0)
G2_val = float(G2[comp].subs(pt))
eps_vals = np.logspace(-3, -0.5, 10)
err = []
for ev in eps_vals:
    ex = float(sp.N(Gx[comp].subs(pt).subs(eps, ev)))
    err.append(abs(ex - ev**2 * G2_val))
slope = float(np.polyfit(np.log(eps_vals), np.log(err), 1)[0])

fig, ax = plt.subplots(figsize=(3.4, 2.6))
ax.loglog(eps_vals, err, "o-", ms=4, label=f"slope {slope:.2f}")
ax.loglog(eps_vals, err[0] * (eps_vals / eps_vals[0]) ** 4, "--", color="0.5", lw=0.8, label=r"$\propto\epsilon^4$")
ax.set(xlabel=r"$\epsilon$", ylabel=r"$|G_{tt}^{\rm exact}-\epsilon^2 G^{(2)}_{tt}|$")
ax.legend()
save_fig(fig, "exp11_second_order_perturbation")
save_results(
    "exp11_second_order_perturbation",
    {
        "G1_is_zero": bool(all(e == 0 for e in G1)),
        "G2_tt": str(G2[0, 0]),
        "G2_nonzero_components": {str(k): str(G2[k]) for k in [(0, 0), (0, 3), (3, 3), (1, 1)] if G2[k] != 0},
        "truncation_error_slope": slope,
    },
)
print("G1 zero:", all(e == 0 for e in G1), "G2_tt:", G2[0, 0], "slope", slope)
