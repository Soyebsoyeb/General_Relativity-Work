"""Experiment 07: free test particles in a plane gravitational wave.

Question
    How does a ring of free particles deform under + and x polarisations, and
    how well does the linearised geodesic-deviation equation reproduce the
    exact proper separation in the TT metric?
Equations
    h_+ = A sin(w t), h_x = A sin(w t) (separate runs), TT gauge, z = 0.
    Proper separation L = sqrt((1+h+) dx^2 + (1-h+) dy^2 + 2 hx dx dy).
    Deviation: xi'' = (1/2) hddot_ij xi^j, xi(0) = d0, xi'(0) = (1/2) hdot(0) d0.
Parameters
    w = 1 (time unit 1/w), amplitudes A in {1e-3 ... 0.2}; ring radius 1 for display.
Numerics
    DOP853 rtol 1e-11 for the deviation equation; exact separation is analytic.
Interpretation
    The two descriptions differ at O(A^2): the measured maximum difference over
    two periods scales with slope ~2 in A.
Limitations
    Linearised theory; amplitudes up to 0.2 are drawn for visibility only and
    exceed the regime where |h| << 1 is a good approximation.
"""

import matplotlib.pyplot as plt
import numpy as np
from _common import save_fig, save_results

from gr_from_scratch.dynamics.gravitational_waves import integrate_deviation, proper_separation
from gr_from_scratch.visualization.style import apply_style

apply_style()
w = 1.0
amps = np.array([1e-3, 3e-3, 1e-2, 3e-2, 1e-1])
th = np.linspace(0, 2 * np.pi, 33)[:-1]
diffs = []
for A in amps:
    worst = 0.0
    for d0 in [(1, 0), (0, 1), (np.cos(np.pi / 4), np.sin(np.pi / 4))]:
        ts, xi = integrate_deviation(d0, A, 0.0, w, 4 * np.pi, n=801)
        lin = np.linalg.norm(xi, axis=1)
        ex = proper_separation(ts, d0, A, 0.0, w)
        worst = max(worst, np.abs(lin - ex).max())
    diffs.append(worst)
slope = float(np.polyfit(np.log(amps), np.log(diffs), 1)[0])

A = 0.2
phases = [0, np.pi / 2, np.pi, 3 * np.pi / 2]
fig, ax = plt.subplots(2, 4, figsize=(7.4, 4.0), sharex=True, sharey=True)
for row, (ap, ac, name) in enumerate([(A, 0, r"$+$"), (0, A, r"$\times$")]):
    for col, ph in enumerate(phases):
        xs, ys = [], []
        for t_ in th:
            d0 = (np.cos(t_), np.sin(t_))
            ts, xi = integrate_deviation(d0, ap, ac, w, ph if ph > 0 else 1e-9, n=3)
            xs.append(xi[-1, 0])
            ys.append(xi[-1, 1])
        ax[row, col].plot(np.append(xs, xs[0]), np.append(ys, ys[0]), "-o", ms=2, lw=0.8)
        ax[row, col].plot(np.cos(np.append(th, th[0])), np.sin(np.append(th, th[0])), color="0.8", lw=0.6)
        ax[row, col].set_aspect("equal")
        ax[row, col].set_title(rf"{name}, $\omega t={ph / np.pi:g}\pi$", fontsize=8)
ax[0, 0].set_xlim(-1.4, 1.4)
ax[0, 0].set_ylim(-1.4, 1.4)
save_fig(fig, "exp07_gravitational_waves")
fig2, ax2 = plt.subplots(figsize=(3.4, 2.6))
ax2.loglog(amps, diffs, "o-", ms=4, label=f"slope {slope:.2f}")
ax2.loglog(amps, diffs[0] * (amps / amps[0]) ** 2, "--", color="0.5", lw=0.8, label=r"$\propto A^2$")
ax2.set(xlabel="amplitude $A$", ylabel=r"max $|L_{\rm lin}-L_{\rm exact}|$")
ax2.legend()
save_fig(fig2, "exp07_gw_linear_vs_exact")
save_results("exp07_gravitational_waves", {"amplitudes": amps, "max_difference": diffs, "loglog_slope": slope})
print("slope", slope)
