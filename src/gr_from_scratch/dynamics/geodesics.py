"""Numerical geodesics.

The geodesic equation in first-order form, with y = (x^mu, u^mu):

    dx^mu/dl = u^mu,        du^mu/dl = - Gamma^mu_{ab}(x) u^a u^b.

``l`` is proper time for timelike curves (g(u,u) = -1) and an affine parameter
for null curves (g(u,u) = 0). The Christoffel symbols come from this
package's own ``geometry.connection`` and are lambdified once.

Validation hooks
    ``constraint``: g(u,u) + kappa along the path (should stay at round-off).
    ``conserved``: values of g(xi, u) for supplied Killing vectors.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import OptimizeResult

from ..core.metric import Metric
from ..geometry.connection import christoffel


@dataclass
class Trajectory:
    lam: np.ndarray
    x: np.ndarray  # shape (N, n)
    u: np.ndarray  # shape (N, n)
    constraint: np.ndarray  # g(u,u) + kappa
    sol: OptimizeResult  # scipy's OdeResult is a subclass

    @property
    def success(self) -> bool:
        return bool(self.sol.success)

    @property
    def events(self):
        return self.sol.t_events, self.sol.y_events


class GeodesicSystem:
    def __init__(self, metric: Metric, params: dict | None = None) -> None:
        self.metric = metric
        self.n = metric.dim
        self.params = dict(params or {})
        conn = christoffel(metric)
        coords = metric.coords
        nz = []
        for idx in np.ndindex(*conn.gamma.shape):
            e = sp.sympify(conn.gamma[idx]).subs(self.params)
            if e != 0:
                nz.append((idx, e))
        self._idx = np.array([i for i, _ in nz], dtype=int)
        self._gam = sp.lambdify(coords, [e for _, e in nz], "numpy", cse=True)
        self._g = sp.lambdify(coords, metric.g.subs(self.params), "numpy")

    # -- geometry at a point ----------------------------------------------
    def g(self, x) -> np.ndarray:
        return np.array(self._g(*x), dtype=float)

    def christoffel(self, x) -> np.ndarray:
        out = np.zeros((self.n,) * 3)
        vals = self._gam(*x)
        for (r, m, v), val in zip(self._idx, vals, strict=True):
            out[r, m, v] = val
        return out

    def norm(self, x, u) -> float:
        return float(u @ self.g(x) @ u)

    def acceleration(self, x, u) -> np.ndarray:
        a = np.zeros(self.n)
        vals = self._gam(*x)
        for (r, m, v), val in zip(self._idx, vals, strict=True):
            a[r] -= val * u[m] * u[v]
        return a

    def rhs(self, _lam, y) -> np.ndarray:
        x, u = y[: self.n], y[self.n :]
        return np.concatenate([u, self.acceleration(x, u)])

    # -- initial data --------------------------------------------------------
    def future_directed_u(self, x, u_spatial, kappa: float) -> np.ndarray:
        """Complete (u^0, u^i) so that g(u,u) = -kappa with u^0 > 0.

        Assumes coordinate 0 is the time coordinate. Solves the quadratic
        g_00 u0^2 + 2 g_0i u^i u0 + g_ij u^i u^j = -kappa.
        """
        g = self.g(x)
        us = np.asarray(u_spatial, dtype=float)
        a, b, c = g[0, 0], g[0, 1:] @ us, us @ g[1:, 1:] @ us
        if abs(a) < 1e-14:
            roots = [-(c + kappa) / (2 * b)]
        else:
            disc = b * b - a * (c + kappa)
            if disc < 0:
                raise ValueError("no real u^0: the requested spatial velocity is not allowed here")
            roots = [(-b + s * np.sqrt(disc)) / a for s in (1, -1)]
        pos = [r for r in roots if r > 0]
        if not pos:
            raise ValueError("no future-directed solution")
        return np.concatenate([[max(pos)], us])

    # -- integration -------------------------------------------------------------
    def integrate(
        self,
        x0,
        u0,
        lam_max: float,
        *,
        kappa: float = 1.0,
        rtol: float = 1e-10,
        atol: float = 1e-12,
        method: str = "DOP853",
        events: Sequence[Callable] | None = None,
        t_eval=None,
        max_step: float = np.inf,
    ) -> Trajectory:
        y0 = np.concatenate([x0, u0]).astype(float)
        sol = solve_ivp(
            self.rhs,
            (0.0, lam_max),
            y0,
            method=method,
            rtol=rtol,
            atol=atol,
            events=events,
            t_eval=t_eval,
            max_step=max_step,
        )
        x, u = sol.y[: self.n].T, sol.y[self.n :].T
        cons = np.array([self.norm(xi, ui) + kappa for xi, ui in zip(x, u, strict=True)])
        return Trajectory(sol.t, x, u, cons, sol)

    def conserved_function(self, xi_up) -> Callable:
        """Return f(x, u) = g_{mu nu}(x) xi^mu(x) u^nu for a Killing vector xi."""
        xi = sp.lambdify(self.metric.coords, [sp.sympify(c).subs(self.params) for c in xi_up], "numpy")

        def f(x, u):
            return float(np.asarray(xi(*x), dtype=float) @ self.g(x) @ u)

        return f


def event_function(index: int, value: float = 0.0, direction: float = 0.0, terminal: bool = False, n: int = 4):
    """Event: component ``index`` of y = (x, u) crosses ``value``.

    Use ``index = n + i`` for the velocity component u^i.
    """

    def ev(_lam, y):
        return y[index] - value

    ev.terminal = terminal  # type: ignore[attr-defined]
    ev.direction = direction  # type: ignore[attr-defined]
    return ev
