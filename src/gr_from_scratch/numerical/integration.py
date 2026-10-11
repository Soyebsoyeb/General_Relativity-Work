"""Fixed-step Runge-Kutta used to measure convergence order independently of SciPy's adaptive solvers."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np


def rk4_fixed(rhs: Callable, y0, t_end: float, n_steps: int):
    """Classical RK4 on [0, t_end] with ``n_steps`` equal steps. Returns (t, y)."""
    y0 = np.asarray(y0, dtype=float)
    h = t_end / n_steps
    t = np.linspace(0.0, t_end, n_steps + 1)
    y = np.empty((n_steps + 1, y0.size))
    y[0] = y0
    for i in range(n_steps):
        k1 = rhs(t[i], y[i])
        k2 = rhs(t[i] + h / 2, y[i] + h / 2 * k1)
        k3 = rhs(t[i] + h / 2, y[i] + h / 2 * k2)
        k4 = rhs(t[i] + h, y[i] + h * k3)
        y[i + 1] = y[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return t, y
