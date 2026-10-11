"""Experiment 05: convergence of the geodesic integration and constraint monitoring.

Question
    Does a fixed-step RK4 applied to the geodesic equation converge at fourth
    order, and how do adaptive tolerances map to achieved error, normalisation
    drift and drift of E and L?
Setup
    Schwarzschild, M = 1, bound orbit p = 20, e = 0.5, one radial period of
    proper time (Newtonian estimate used only to size the interval).
    Reference: DOP853 with rtol = 1e-13, atol = 1e-14.
Metrics
    Error = Euclidean norm of (x, u) at the final time versus the reference.
    Constraint = |g(u,u) + 1|. Drift = max |Q(l) - Q(0)| for E = -g(d_t, u), L = g(d_phi, u).
Tolerances
    Adaptive rtol from 1e-4 to 1e-12, atol = 1e-3 * rtol.
Limitations
    The reference is itself numerical; errors below ~1e-12 are not meaningful.
    Fixed-step RK4 is for order measurement only, not recommended for production.
"""

import matplotlib.pyplot as plt
import numpy as np
from _common import save_fig, save_results

from gr_from_scratch.dynamics import observables as ob
from gr_from_scratch.numerical.convergence import observed_order
from gr_from_scratch.numerical.integration import rk4_fixed
from gr_from_scratch.visualization.style import apply_style

apply_style()
gs = ob.schwarzschild_system()
p, e = 20.0, 0.5
E, L, up, ua = ob._orbit_constants(p, e, 1.0)
rp = 1 / up
f = 1 - 2 / rp
x0 = np.array([0, rp, np.pi / 2, 0.0])
u0 = np.array([E / f, 0, 0, L / rp**2])
a = (rp + 1 / ua) / 2
T = 2 * np.pi * np.sqrt(a**3)
ref = gs.integrate(x0, u0, T, rtol=1e-13, atol=1e-14).sol.y[:, -1]

steps = np.array([200, 400, 800, 1600, 3200])
errs = []
for n in steps:
    _, y = rk4_fixed(gs.rhs, np.concatenate([x0, u0]), T, int(n))
    errs.append(np.linalg.norm(y[-1] - ref))
h = T / steps
order = observed_order(h, errs)

xi_t, xi_p = [1, 0, 0, 0], [0, 0, 0, 1]
Et, Lp = gs.conserved_function(xi_t), gs.conserved_function(xi_p)
rtols = 10.0 ** -np.arange(4, 13)
ad_err, ad_con, ad_E, ad_L = [], [], [], []
for rt in rtols:
    tr = gs.integrate(x0, u0, T, rtol=rt, atol=1e-3 * rt)
    ad_err.append(np.linalg.norm(tr.sol.y[:, -1] - ref))
    ad_con.append(np.abs(tr.constraint).max())
    ad_E.append(max(abs(Et(x, u) - Et(tr.x[0], tr.u[0])) for x, u in zip(tr.x, tr.u, strict=True)))
    ad_L.append(max(abs(Lp(x, u) - Lp(tr.x[0], tr.u[0])) for x, u in zip(tr.x, tr.u, strict=True)))

fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8))
ax[0].loglog(h, errs, "o-", ms=4, label=f"RK4, observed order {order:.2f}")
ax[0].loglog(h, errs[0] * (h / h[0]) ** 4, "--", color="0.5", lw=0.8, label=r"$\propto h^4$")
ax[0].set(xlabel=r"step $h$ (proper time)", ylabel="final-state error")
ax[0].legend()
ax[1].loglog(rtols, ad_err, "o-", ms=4, label="final-state error")
ax[1].loglog(rtols, np.array(ad_con) + 1e-17, "s-", ms=3, label=r"max $|g(u,u)+1|$")
ax[1].loglog(rtols, np.array(ad_E) + 1e-17, "^-", ms=3, label=r"drift in $E$")
ax[1].loglog(rtols, np.array(ad_L) + 1e-17, "v-", ms=3, label=r"drift in $L$")
ax[1].set(xlabel="DOP853 rtol", ylabel="error / drift")
ax[1].legend(fontsize=7)
save_fig(fig, "exp05_integrator_convergence")
save_results(
    "exp05_integrator_convergence",
    {
        "rk4_steps": steps,
        "rk4_errors": errs,
        "rk4_observed_order": order,
        "rtol": rtols,
        "adaptive_error": ad_err,
        "adaptive_constraint": ad_con,
        "adaptive_E_drift": ad_E,
        "adaptive_L_drift": ad_L,
    },
)
print("RK4 order", order)
