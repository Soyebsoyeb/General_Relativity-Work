"""Orthonormal frames and Cartan's structure equations.

Conventions
    Coframe e^a = e^a_mu dx^mu with  g_{mu nu} = eta_{ab} e^a_mu e^b_nu.
    Frame vectors E_a = E_a^mu d_mu are the inverse: E_a^mu e^b_mu = delta_a^b.

First structure equation (zero torsion):   d e^a + omega^a_b ^ e^b = 0
Second structure equation:                 Omega^a_b = d omega^a_b + omega^a_c ^ omega^c_b

Writing d e^a = (1/2) C^a_{bc} e^b ^ e^c and omega^a_b = omega^a_{bc} e^c, the
first equation is C^a_{bc} = omega^a_{bc} - omega^a_{cb}. With omega_{abc}
antisymmetric in (a, b) (metric compatibility of the frame) the unique
solution is

    omega_{abc} = (1/2) (C_{abc} + C_{bca} - C_{cab}).

Nothing in this module uses Christoffel symbols, so agreement with the
coordinate calculation is an independent check.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from ..core import index_operations as ops
from ..core.metric import Metric
from ..core.tensor import Tensor
from .differential_forms import Form


class Coframe:
    def __init__(self, metric: Metric, e, eta=None) -> None:
        n = metric.dim
        self.metric = metric
        self.n = n
        self.e = sp.ImmutableMatrix(e)  # e[a, mu]
        self.eta = sp.ImmutableMatrix(eta if eta is not None else sp.diag(-1, *([1] * (n - 1))))
        self.E = sp.ImmutableMatrix(self.e.inv().applyfunc(metric.simplifier))  # E[mu, a]
        resid = (self.e.T * self.eta * self.e - metric.g).applyfunc(metric.simplifier)
        if resid != sp.zeros(n):
            raise ValueError("coframe is not orthonormal for this metric")

    def forms(self) -> list[Form]:
        return [Form.one_form(self.metric.coords, [self.e[a, m] for m in range(self.n)]) for a in range(self.n)]

    def structure_functions(self) -> np.ndarray:
        """C^a_{bc} defined by d e^a = (1/2) C^a_{bc} e^b ^ e^c."""
        n, s = self.n, self.metric.simplifier
        C = ops.zeros(n, 3)
        for a, form in enumerate(self.forms()):
            D = form.d().component_array()  # (d e^a)_{mu nu}
            for b in range(n):
                for c in range(n):
                    C[a, b, c] = s(sum(self.E[m, b] * self.E[v, c] * D[m, v] for m in range(n) for v in range(n)))
        return C

    def spin_connection_frame(self) -> np.ndarray:
        """omega^a_{b c} (frame components, first index up)."""
        n, s = self.n, self.metric.simplifier
        C = self.structure_functions()
        Cd = ops.zeros(n, 3)  # lower first index with eta
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    Cd[a, b, c] = sum(self.eta[a, d] * C[d, b, c] for d in range(n))
        w = ops.zeros(n, 3)
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    w[a, b, c] = sp.Rational(1, 2) * (Cd[a, b, c] + Cd[b, c, a] - Cd[c, a, b])
        etainv = self.eta.inv()
        out = ops.zeros(n, 3)
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    out[a, b, c] = s(sum(etainv[a, d] * w[d, b, c] for d in range(n)))
        return out

    def spin_connection_coords(self) -> np.ndarray:
        """omega^a_{b mu} = omega^a_{b c} e^c_mu."""
        n, s = self.n, self.metric.simplifier
        wf = self.spin_connection_frame()
        out = ops.zeros(n, 3)
        for a in range(n):
            for b in range(n):
                for mu in range(n):
                    out[a, b, mu] = s(sum(wf[a, b, c] * self.e[c, mu] for c in range(n)))
        return out

    def first_structure_residual(self) -> np.ndarray:
        """(d e^a)_{mu nu} + omega^a_{b mu} e^b_nu - omega^a_{b nu} e^b_mu."""
        n, s = self.n, self.metric.simplifier
        w = self.spin_connection_coords()
        out = ops.zeros(n, 3)
        for a, form in enumerate(self.forms()):
            D = form.d().component_array()
            for mu in range(n):
                for nu in range(n):
                    val = D[mu, nu] + sum(w[a, b, mu] * self.e[b, nu] - w[a, b, nu] * self.e[b, mu] for b in range(n))
                    out[a, mu, nu] = s(val)
        return out

    def curvature_two_form_components(self) -> np.ndarray:
        """Omega^a_{b mu nu} from the second structure equation."""
        n, x, s = self.n, self.metric.coords, self.metric.simplifier
        w = self.spin_connection_coords()
        out = ops.zeros(n, 4)
        for a in range(n):
            for b in range(n):
                for mu in range(n):
                    for nu in range(n):
                        val = sp.diff(w[a, b, nu], x[mu]) - sp.diff(w[a, b, mu], x[nu])
                        val += sum(w[a, c, mu] * w[c, b, nu] - w[a, c, nu] * w[c, b, mu] for c in range(n))
                        out[a, b, mu, nu] = s(val)
        return out

    def riemann_frame(self) -> np.ndarray:
        """R^a_{bcd} = Omega^a_{b mu nu} E_c^mu E_d^nu."""
        n, s = self.n, self.metric.simplifier
        Om = self.curvature_two_form_components()
        out = ops.zeros(n, 4)
        for a, b, c, d in np.ndindex(n, n, n, n):
            out[a, b, c, d] = s(sum(Om[a, b, m, v] * self.E[m, c] * self.E[v, d] for m in range(n) for v in range(n)))
        return out

    def project_coordinate_riemann(self, riem: Tensor) -> np.ndarray:
        """R^a_{bcd} = e^a_rho E_b^sigma E_c^mu E_d^nu R^rho_{sigma mu nu} (coordinate result)."""
        n, s = self.n, self.metric.simplifier
        R = riem.c
        out = ops.zeros(n, 4)
        for a, b, c, d in np.ndindex(n, n, n, n):
            tot = sp.S.Zero
            for r in range(n):
                if self.e[a, r] == 0:
                    continue
                for sg in range(n):
                    if self.E[sg, b] == 0:
                        continue
                    for m in range(n):
                        if self.E[m, c] == 0:
                            continue
                        for v in range(n):
                            tot += self.e[a, r] * self.E[sg, b] * self.E[m, c] * self.E[v, d] * R[r, sg, m, v]
            out[a, b, c, d] = s(tot)
        return out
