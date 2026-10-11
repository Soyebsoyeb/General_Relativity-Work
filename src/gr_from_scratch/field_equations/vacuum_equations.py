"""Einstein-equation residuals and a derivation of the static spherical vacuum solution."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import sympy as sp

from ..core.metric import Metric
from ..core.tensor import Tensor
from ..geometry.connection import Connection, christoffel
from ..geometry.einstein import einstein_tensor
from ..geometry.ricci import ricci_scalar, ricci_tensor
from ..geometry.riemann import riemann


@dataclass
class Geometry:
    metric: Metric
    conn: Connection
    riemann: Tensor
    ricci: Tensor
    scalar: sp.Expr
    einstein: Tensor


def compute_geometry(metric: Metric) -> Geometry:
    conn = christoffel(metric)
    riem = riemann(conn)
    ric = ricci_tensor(riem)
    R = ricci_scalar(ric)
    return Geometry(metric, conn, riem, ric, R, einstein_tensor(ric, R))


def field_equation_residual(
    geom: Geometry, T: Tensor | None = None, cosmological_constant=0, coupling=8 * sp.pi
) -> Tensor:
    """E_{mu nu} = G_{mu nu} + Lambda g_{mu nu} - coupling * T_{mu nu}."""
    m = geom.metric
    g = np.array(m.g.tolist(), dtype=object)
    out = geom.einstein.c + cosmological_constant * g
    if T is not None:
        out = out - coupling * T.c
    return Tensor(m, out, "dd").simplify()


def derive_static_spherical_vacuum():
    """Solve the vacuum equations for ds^2 = -A(r) dt^2 + B(r) dr^2 + r^2 dOmega^2.

    Returns a dict containing the mixed Einstein components and the solution.
    Steps:
      1. G^t_t = -(1/r^2) d/dr [ r (1 - 1/B) ]   (checked symbolically)
      2. G^t_t = 0  =>  1/B = 1 - 2M/r  (M an integration constant)
      3. G^r_r - G^t_t = (1/(r B)) (A'/A + B'/B)  =>  A B = const; rescale t so const = 1.
    Assumptions: staticity and spherical symmetry (no appeal to Birkhoff's theorem).
    """
    r, th, ph, t = sp.symbols("r theta phi t", positive=True)
    M = sp.Symbol("M", real=True)
    A, B = sp.Function("A")(r), sp.Function("B")(r)
    m = Metric([t, r, th, ph], sp.diag(-A, B, r**2, r**2 * sp.sin(th) ** 2), name="static-spherical")
    geom = compute_geometry(m)
    Gmixed = (geom.einstein.raise_index(0)).c  # G^mu_nu
    Gtt, Grr = sp.simplify(Gmixed[0, 0]), sp.simplify(Gmixed[1, 1])
    claim1 = sp.simplify(Gtt + sp.diff(r * (1 - 1 / B), r) / r**2)
    u = sp.Function("u")(r)
    sol_u = sp.dsolve(sp.Eq(sp.simplify(Gtt.subs(B, 1 / u).doit()), 0), u)
    claim3 = sp.simplify(Grr - Gtt - (sp.diff(A, r) / A + sp.diff(B, r) / B) / (r * B))
    return {
        "metric": m,
        "G^t_t": Gtt,
        "G^r_r": Grr,
        "check_Gtt_form": claim1,
        "u=1/B solution": sol_u,
        "check_Grr_minus_Gtt_form": claim3,
        "M": M,
    }
