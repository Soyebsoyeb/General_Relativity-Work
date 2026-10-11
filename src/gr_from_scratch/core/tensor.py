"""Tensors with explicit index variance on top of a ``Metric``.

A ``Tensor`` stores its components (an object array of SymPy expressions) and
a variance string. The metric is only used to move indices. Coordinate
transformations of tensors live in ``coordinates.py``.

Limitations
    Components are always stored explicitly; there is no abstract index
    algebra and no canonicalisation under index symmetries.
"""

from __future__ import annotations

import itertools

import numpy as np
import sympy as sp

from . import index_operations as ops
from .metric import Metric


class TensorError(ValueError):
    pass


class Tensor:
    def __init__(self, metric: Metric, components, variance: str) -> None:
        arr = np.array(components, dtype=object)
        if set(variance) - {"u", "d"}:
            raise TensorError("variance must contain only 'u' and 'd'")
        if arr.ndim != len(variance) or any(s != metric.dim for s in arr.shape):
            raise TensorError(
                f"components of shape {arr.shape} do not match variance {variance!r} in dimension {metric.dim}"
            )
        self.metric = metric
        self.c = arr
        self.variance = variance

    @property
    def rank(self) -> int:
        return len(self.variance)

    def __getitem__(self, idx):
        return self.c[idx]

    # -- constructors ------------------------------------------------------
    @classmethod
    def from_function(cls, metric: Metric, variance: str, fn) -> Tensor:
        arr = ops.zeros(metric.dim, len(variance))
        for idx in itertools.product(range(metric.dim), repeat=len(variance)):
            arr[idx] = fn(*idx)
        return cls(metric, arr, variance)

    @classmethod
    def metric_tensor(cls, metric: Metric) -> Tensor:
        return cls(metric, np.array(metric.g.tolist(), dtype=object), "dd")

    # -- algebra -----------------------------------------------------------
    def _like(self, comps, variance=None) -> Tensor:
        return Tensor(self.metric, comps, variance or self.variance)

    def _check_compatible(self, other: Tensor) -> None:
        if self.variance != other.variance:
            raise TensorError(f"variance mismatch: {self.variance!r} vs {other.variance!r}")

    def __add__(self, other: Tensor) -> Tensor:
        self._check_compatible(other)
        return self._like(self.c + other.c)

    def __sub__(self, other: Tensor) -> Tensor:
        self._check_compatible(other)
        return self._like(self.c - other.c)

    def __neg__(self) -> Tensor:
        return self._like(-self.c)

    def __mul__(self, scalar) -> Tensor:
        return self._like(self.c * sp.sympify(scalar))

    __rmul__ = __mul__

    def outer(self, other: Tensor) -> Tensor:
        comps = np.multiply.outer(self.c, other.c)
        return Tensor(self.metric, comps, self.variance + other.variance)

    def contract(self, i: int, j: int) -> Tensor:
        if self.variance[i] == self.variance[j]:
            raise TensorError("contraction needs one upper and one lower index")
        keep = "".join(v for k, v in enumerate(self.variance) if k not in (i, j))
        return Tensor(self.metric, ops.contract_axes(self.c, i, j), keep)

    def transpose(self, perm: tuple[int, ...]) -> Tensor:
        return Tensor(self.metric, ops.permute(self.c, perm), "".join(self.variance[p] for p in perm))

    def symmetrize(self, axes: tuple[int, ...]) -> Tensor:
        return self._like(ops.symmetrize(self.c, axes))

    def antisymmetrize(self, axes: tuple[int, ...]) -> Tensor:
        return self._like(ops.antisymmetrize(self.c, axes))

    # -- index gymnastics --------------------------------------------------
    def raise_index(self, axis: int) -> Tensor:
        if self.variance[axis] == "u":
            raise TensorError(f"index {axis} is already contravariant")
        comps = ops.apply_matrix(self.c, self.metric.ginv, axis)
        v = self.variance[:axis] + "u" + self.variance[axis + 1 :]
        return Tensor(self.metric, comps, v).simplify()

    def lower_index(self, axis: int) -> Tensor:
        if self.variance[axis] == "d":
            raise TensorError(f"index {axis} is already covariant")
        comps = ops.apply_matrix(self.c, self.metric.g, axis)
        v = self.variance[:axis] + "d" + self.variance[axis + 1 :]
        return Tensor(self.metric, comps, v).simplify()

    def to_variance(self, target: str) -> Tensor:
        if len(target) != self.rank:
            raise TensorError("target variance has the wrong rank")
        out = self
        for k, (have, want) in enumerate(zip(self.variance, target, strict=True)):
            if have != want:
                out = out.raise_index(k) if want == "u" else out.lower_index(k)
        return out

    # -- utilities ---------------------------------------------------------
    def simplify(self) -> Tensor:
        f = self.metric.simplifier
        out = np.empty(self.c.shape, dtype=object)
        for idx in np.ndindex(self.c.shape):
            out[idx] = f(self.c[idx])
        return self._like(out)

    def diff(self, symbol) -> np.ndarray:
        out = np.empty(self.c.shape, dtype=object)
        for idx in np.ndindex(self.c.shape):
            out[idx] = sp.diff(self.c[idx], symbol)
        return out

    def subs(self, mapping: dict) -> Tensor:
        out = np.empty(self.c.shape, dtype=object)
        for idx in np.ndindex(self.c.shape):
            out[idx] = sp.sympify(self.c[idx]).subs(mapping)
        return self._like(out)

    def is_zero(self) -> bool:
        return all(self.metric.simplifier(e) == 0 for e in self.c.ravel())

    def nonzero_components(self) -> dict:
        return {
            idx: self.metric.simplifier(self.c[idx])
            for idx in np.ndindex(self.c.shape)
            if self.metric.simplifier(self.c[idx]) != 0
        }

    def full_contraction(self, other: Tensor) -> sp.Expr:
        """Contract all indices of ``self`` with ``other`` (variances must be opposite)."""
        if self.rank != other.rank:
            raise TensorError("rank mismatch")
        a, b = self, other
        for k, (u, v) in enumerate(zip(a.variance, b.variance, strict=True)):
            if u == v:  # move b's index so that it is opposite to a's
                b = b.raise_index(k) if u == "d" else b.lower_index(k)
        total = sp.S.Zero
        for idx in np.ndindex(a.c.shape):
            total += a.c[idx] * b.c[idx]
        return self.metric.simplifier(total)

    def __repr__(self) -> str:
        return f"Tensor(variance={self.variance!r}, dim={self.metric.dim})"
