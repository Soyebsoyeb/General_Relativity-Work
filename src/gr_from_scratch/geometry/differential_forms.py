"""Differential forms: wedge product and exterior derivative.

A p-form is stored as ``{(i_1 < ... < i_p): coefficient}`` meaning
omega = sum coeff dx^{i_1} ^ ... ^ dx^{i_p}. Zero coefficients are omitted.
"""

from __future__ import annotations

import itertools

import numpy as np
import sympy as sp

from ..core import index_operations as ops


def _sort_with_sign(idx: tuple[int, ...]):
    """Sort an index tuple, returning (sorted tuple, sign) or (None, 0) if repeated."""
    if len(set(idx)) < len(idx):
        return None, 0
    sign, lst = 1, list(idx)
    for i in range(len(lst)):
        for j in range(len(lst) - 1 - i):
            if lst[j] > lst[j + 1]:
                lst[j], lst[j + 1] = lst[j + 1], lst[j]
                sign = -sign
    return tuple(lst), sign


class Form:
    def __init__(self, coords, degree: int, comps: dict | None = None) -> None:
        self.coords = tuple(coords)
        self.degree = degree
        self.comps: dict[tuple[int, ...], sp.Expr] = {}
        for k, v in (comps or {}).items():
            key, sign = _sort_with_sign(tuple(k))
            if key is not None:
                v = sign * sp.sympify(v)
                self.comps[key] = self.comps.get(key, 0) + v
        self.comps = {k: v for k, v in self.comps.items() if v != 0}

    @classmethod
    def one_form(cls, coords, components) -> Form:
        return cls(coords, 1, {(i,): c for i, c in enumerate(components)})

    def __add__(self, other: Form) -> Form:
        if self.degree != other.degree:
            raise ValueError("cannot add forms of different degree")
        c = dict(self.comps)
        for k, v in other.comps.items():
            c[k] = c.get(k, 0) + v
        return Form(self.coords, self.degree, c)

    def __mul__(self, s) -> Form:
        return Form(self.coords, self.degree, {k: sp.sympify(s) * v for k, v in self.comps.items()})

    __rmul__ = __mul__

    def wedge(self, other: Form) -> Form:
        out: dict = {}
        for ka, va in self.comps.items():
            for kb, vb in other.comps.items():
                key, sign = _sort_with_sign(ka + kb)
                if key is not None:
                    out[key] = out.get(key, 0) + sign * va * vb
        return Form(self.coords, self.degree + other.degree, out)

    def d(self) -> Form:
        out: dict = {}
        for k, v in self.comps.items():
            for j, xj in enumerate(self.coords):
                dv = sp.diff(v, xj)
                if dv != 0:
                    key, sign = _sort_with_sign((j,) + k)
                    if key is not None:
                        out[key] = out.get(key, 0) + sign * dv
        return Form(self.coords, self.degree + 1, out)

    def simplified(self, simplifier=sp.simplify) -> Form:
        return Form(self.coords, self.degree, {k: simplifier(v) for k, v in self.comps.items()})

    def is_zero(self, simplifier=sp.simplify) -> bool:
        return all(simplifier(v) == 0 for v in self.comps.values())

    def component_array(self) -> np.ndarray:
        """Fully antisymmetric array omega_{i_1...i_p} (omega = (1/p!) omega_I dx^I)."""
        n = len(self.coords)
        arr = ops.zeros(n, self.degree)
        for k, v in self.comps.items():
            for perm in itertools.permutations(range(self.degree)):
                _, sign = _sort_with_sign(tuple(perm))
                arr[tuple(k[p] for p in perm)] = sign * v
        return arr
