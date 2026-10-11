"""Metric tensor on a coordinate patch."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import sympy as sp

from .simplification import simplify


class MetricError(ValueError):
    pass


def as_symbols(coords) -> tuple[sp.Symbol, ...]:
    if isinstance(coords, str):
        return tuple(sp.symbols(coords, real=True))
    return tuple(coords)


class Metric:
    """A (pseudo-)Riemannian metric g_{mu nu} given by symbolic components.

    ``inverse`` may be supplied (it is checked against ``g``); this is the
    practical route for metrics such as Kerr whose inverse SymPy finds slowly.
    ``simplifier`` is applied to every derived component.
    """

    def __init__(
        self,
        coords,
        g,
        *,
        name: str = "",
        inverse=None,
        simplifier: Callable = simplify,
    ) -> None:
        self.coords = as_symbols(coords)
        self.dim = len(self.coords)
        self.g = sp.ImmutableMatrix(g)
        if self.g.shape != (self.dim, self.dim):
            raise MetricError(f"expected a {self.dim}x{self.dim} matrix, got {self.g.shape}")
        if self.g != self.g.T:
            raise MetricError("metric components must be symmetric")
        self.name = name
        self.simplifier = simplifier
        if inverse is not None:
            self.ginv = sp.ImmutableMatrix(inverse)
            resid = (self.g * self.ginv - sp.eye(self.dim)).applyfunc(simplify)
            if resid != sp.zeros(self.dim):
                raise MetricError("supplied inverse does not satisfy g g^-1 = 1")
        else:
            self.ginv = sp.ImmutableMatrix(self.g.inv()).applyfunc(simplifier)
        self._det = None

    @property
    def det(self):
        if self._det is None:
            self._det = self.simplifier(self.g.det())
        return self._det

    @property
    def sqrt_minus_g(self):
        """sqrt(|det g|), which equals sqrt(-det g) for Lorentzian signature."""
        return sp.sqrt(sp.simplify(sp.Abs(self.det)))

    def signature_at(self, point: dict) -> tuple[int, int]:
        """(number of negative, number of positive) eigenvalues at ``point``."""
        vals = np.linalg.eigvalsh(np.array(self.g.subs(point).evalf(), dtype=float))
        return int((vals < 0).sum()), int((vals > 0).sum())

    def is_lorentzian_at(self, point: dict) -> bool:
        return self.signature_at(point) == (1, self.dim - 1)

    def subs(self, mapping: dict) -> Metric:
        """Substitute parameter values, keeping the coordinates."""
        keep = {k: v for k, v in mapping.items() if k not in self.coords}
        return Metric(
            self.coords,
            self.g.subs(keep),
            name=self.name,
            inverse=self.ginv.subs(keep),
            simplifier=self.simplifier,
        )

    def transform(self, new_coords, old_in_new: Sequence) -> Metric:
        """Metric in coordinates x' given x^mu = x^mu(x').

        g'_{ab} = (dx^mu/dx'^a)(dx^nu/dx'^b) g_{mu nu}(x(x')).
        """
        new = as_symbols(new_coords)
        if len(new) != self.dim:
            raise MetricError("coordinate transformation must preserve dimension")
        sub = dict(zip(self.coords, old_in_new, strict=True))
        jac = sp.Matrix(self.dim, self.dim, lambda m, a: sp.diff(old_in_new[m], new[a]))
        gnew = (jac.T * self.g.subs(sub) * jac).applyfunc(self.simplifier)
        return Metric(new, gnew, name=self.name, simplifier=self.simplifier)

    def __repr__(self) -> str:
        return f"Metric({self.name or 'unnamed'}, coords={self.coords})"
