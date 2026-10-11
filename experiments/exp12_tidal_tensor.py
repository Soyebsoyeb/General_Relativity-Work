"""Experiment 12: tidal tensor of Schwarzschild in a static orthonormal frame.

Question
    What is E_ij = R^i_{0j0}, and does the geodesic-deviation equation give
    radial stretching and transverse compression with magnitudes 2M/r^3 and M/r^3?
Equations
    d^2 xi^i / d tau^2 = - E_ij xi^j. Expected E = diag(-2M/r^3, M/r^3, M/r^3)
    (radial, theta, phi), trace-free as required in vacuum.
Method
    Orthonormal coframe, Cartan structure equations (no Christoffel symbols),
    frame Riemann components -> tidal tensor, then compared with the
    closed form at sample radii; the coordinate-basis Riemann is not used.
Limitations
    Static frame; for radial free-fall frames boosts in the t-r plane leave
    these components unchanged because the Weyl tensor is of type D, a statement
    checked indirectly by tests/test_tetrads.py (frame vs coordinate agreement).
"""

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from _common import save_fig, save_results

from gr_from_scratch.dynamics.geodesic_deviation import tidal_tensor
from gr_from_scratch.geometry.tetrads import Coframe
from gr_from_scratch.spacetimes.schwarzschild import M, f, r, schwarzschild, theta
from gr_from_scratch.visualization.style import apply_style

apply_style()
met = schwarzschild()
ff = f()
e = sp.diag(sp.sqrt(ff), 1 / sp.sqrt(ff), r, r * sp.sin(theta))
E = tidal_tensor(Coframe(met, e))
E = sp.Matrix(E.tolist()).applyfunc(sp.simplify)
expected = sp.diag(-2 * M / r**3, M / r**3, M / r**3)
resid = (E - expected).applyfunc(sp.simplify)
assert resid == sp.zeros(3)
assert sp.simplify(E.trace()) == 0

rr = np.linspace(2.2, 20, 200)
fig, ax = plt.subplots(figsize=(3.4, 2.6))
ax.plot(rr, 2 / rr**3, label=r"radial: $d^2\xi/d\tau^2 = +2M\xi/r^3$")
ax.plot(rr, -1 / rr**3, label=r"transverse: $-M\xi/r^3$")
ax.axvline(2, color="0.6", lw=0.6)
ax.set(xlabel=r"$r/M$", ylabel=r"relative acceleration per unit $\xi$ [$M^{-2}$]")
ax.legend(fontsize=7)
save_fig(fig, "exp12_tidal_tensor")
save_results("exp12_tidal_tensor", {"E_ij": str(E), "trace": str(sp.simplify(E.trace())), "matches_expected": True})
print(E)
