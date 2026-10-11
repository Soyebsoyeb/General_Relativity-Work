"""Higher-curvature corrections and the range of validity of the effective theory.

Status: one worked symbolic example and one dimensional estimate.

Example. For the reduced static spherical Lagrangian
    L = sqrt|g| ( R + a R^2 + b R_{mn}R^{mn} + c R_{abcd}R^{abcd} + d I3 ),
I3 = R_{ab}^{cd} R_{cd}^{ef} R_{ef}^{ab}, the Euler-Lagrange equations are
evaluated on Schwarzschild. The three quadratic terms contribute nothing on
a Ricci-flat background in four dimensions (R^2 and Ric^2 are quadratic in
the Ricci tensor, and Riem^2 differs from the Gauss-Bonnet density, which is
topological, by such terms). The cubic term does not vanish, which is where
vacuum corrections first appear. No claim is made about the actual size of
any coefficient.
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy import constants as const
from sympy.calculus.euler import euler_equations

from ..core.metric import Metric
from ..field_equations.vacuum_equations import compute_geometry
from ..geometry.invariants import euler_density, kretschmann, ricci_squared


def reduced_corrections_on_schwarzschild():
    """Return dict of E-L residual of each term of the reduced Lagrangian on Schwarzschild."""
    r, th, ph, t = sp.symbols("r theta phi t", positive=True)
    M = sp.Symbol("M", positive=True)
    A, B = sp.Function("A", positive=True)(r), sp.Function("B", positive=True)(r)
    m = Metric([t, r, th, ph], sp.diag(-A, B, r**2, r**2 * sp.sin(th) ** 2))
    geom = compute_geometry(m)
    n = 4
    # R_{ab}^{cd} as a 16x16 matrix for the cubic invariant
    Rdn = geom.riemann.lower_index(0)
    Rmix = Rdn.raise_index(2).raise_index(3)
    mat = sp.Matrix(n * n, n * n, lambda i, j: Rmix.c[i // n, i % n, j // n, j % n])
    I3 = sp.simplify((mat * mat * mat).trace())
    invs = {
        "R^2": sp.simplify(geom.scalar**2),
        "Ric^2": sp.simplify(ricci_squared(geom.riemann)),
        "Riem^2": sp.simplify(kretschmann(geom.riemann)),
        "Gauss-Bonnet": sp.simplify(euler_density(geom.riemann)),
        "Riem^3": I3,
    }
    sg = r**2 * sp.sqrt(A * B)  # the sin(theta) factor is common to all terms and dropped
    f = 1 - 2 * M / r
    out = {}
    for name, inv in invs.items():
        L = sg * inv
        eqs = euler_equations(L, [A, B], r)
        out[name] = [sp.simplify(e.lhs.subs({A: f, B: 1 / f}).doit()) for e in eqs]
    return out


def curvature_expansion_parameter(mass_kg: float, radius_in_horizon_units: float = 1.0) -> float:
    """eps = l_P^2 sqrt(K) at r = radius_in_horizon_units * 2GM/c^2 (SI inputs).

    sqrt(K) = 4 sqrt(3) G M / (c^2 r^3) for Schwarzschild. eps << 1 is the
    regime in which a truncated derivative expansion is self-consistent.
    """
    G, c, hbar = const.G, const.c, const.hbar
    lp2 = hbar * G / c**3
    rg = G * mass_kg / c**2
    rr = radius_in_horizon_units * 2 * rg
    return float(lp2 * 4 * np.sqrt(3) * G * mass_kg / (c**2 * rr**3))
