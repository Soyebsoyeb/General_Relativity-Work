"""Free test particles in a plane gravitational wave (linearised theory, TT gauge).

Metric at z = 0:  g_xx = 1 + h+, g_yy = 1 - h+, g_xy = hx,  with h = A sin(w t).
Two descriptions are compared:

1. Exact proper separation of particles fixed at TT coordinates (they are
   geodesics of the TT metric because Gamma^i_{00} = 0).
2. The linearised geodesic-deviation equation  xi'' = (1/2) hddot_ij xi^j.

They agree to O(A^2); the difference is the neglected second-order term.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp


def h_components(t, a_plus, a_cross, omega):
    s = np.sin(omega * t)
    return a_plus * s, a_cross * s


def proper_separation(t, delta0, a_plus, a_cross, omega):
    """sqrt(g_ij D^i D^j) for TT-fixed coordinate separation ``delta0`` = (Dx, Dy)."""
    hp, hx = h_components(t, a_plus, a_cross, omega)
    dx, dy = delta0
    return np.sqrt((1 + hp) * dx**2 + (1 - hp) * dy**2 + 2 * hx * dx * dy)


def integrate_deviation(delta0, a_plus, a_cross, omega, t_end, n=2001, rtol=1e-11):
    """Solve xi'' = (1/2) hddot_ij xi^j with xi(0) = delta0 and xi'(0) = (1/2) hdot(0) xi(0)."""
    d0 = np.asarray(delta0, float)

    def hddot(t):
        c = -(omega**2) * np.sin(omega * t)
        return np.array([[a_plus * c, a_cross * c], [a_cross * c, -a_plus * c]])

    hdot0 = omega * np.array([[a_plus, a_cross], [a_cross, -a_plus]])
    v0 = 0.5 * hdot0 @ d0

    def rhs(t, y):
        return np.concatenate([y[2:], 0.5 * hddot(t) @ y[:2]])

    ts = np.linspace(0, t_end, n)
    sol = solve_ivp(rhs, (0, t_end), np.concatenate([d0, v0]), t_eval=ts, rtol=rtol, atol=1e-14, method="DOP853")
    return ts, sol.y[:2].T
