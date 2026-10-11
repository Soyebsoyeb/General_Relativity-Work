"""Array-level index manipulation on object arrays of SymPy expressions."""

from __future__ import annotations

import itertools
import math

import numpy as np
import sympy as sp


def zeros(dim: int, rank: int) -> np.ndarray:
    out = np.empty((dim,) * rank, dtype=object)
    out.fill(sp.S.Zero)
    return out


def apply_matrix(arr: np.ndarray, mat, axis: int) -> np.ndarray:
    """out[..., a, ...] = sum_m mat[a, m] arr[..., m, ...] along ``axis``."""
    m = np.array(sp.Matrix(mat).tolist(), dtype=object)
    return np.moveaxis(np.tensordot(m, arr, axes=([1], [axis])), 0, axis)


def contract_axes(arr: np.ndarray, i: int, j: int) -> np.ndarray:
    """Sum over the diagonal of axes ``i`` and ``j`` (no metric involved)."""
    if i == j:
        raise ValueError("cannot contract an axis with itself")
    dim = arr.shape[i]
    moved = np.moveaxis(arr, (i, j), (0, 1))
    out = zeros(dim, arr.ndim - 2)
    for k in range(dim):
        out = out + moved[k, k]
    return out


def permute(arr: np.ndarray, perm: tuple[int, ...]) -> np.ndarray:
    """out[i_0, ..., i_n] = arr[i_perm(0), ...] i.e. ``np.transpose`` semantics."""
    return np.transpose(arr, perm)


def _perm_sign(p: tuple[int, ...]) -> int:
    sign, seen = 1, [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, length = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                length += 1
            if length % 2 == 0:
                sign = -sign
    return sign


def symmetrize(arr: np.ndarray, axes: tuple[int, ...]) -> np.ndarray:
    """Average over all permutations of ``axes``."""
    out = zeros(arr.shape[0], arr.ndim)
    perms = list(itertools.permutations(range(len(axes))))
    for p in perms:
        order = list(range(arr.ndim))
        for k, axis in enumerate(axes):
            order[axis] = axes[p[k]]
        out = out + np.transpose(arr, order)
    return out / math.factorial(len(axes))


def antisymmetrize(arr: np.ndarray, axes: tuple[int, ...]) -> np.ndarray:
    """Signed average over all permutations of ``axes``: T_{[a b ...]}."""
    out = zeros(arr.shape[0], arr.ndim)
    for p in itertools.permutations(range(len(axes))):
        order = list(range(arr.ndim))
        for k, axis in enumerate(axes):
            order[axis] = axes[p[k]]
        out = out + _perm_sign(p) * np.transpose(arr, order)
    return out / math.factorial(len(axes))


def levi_civita_symbol(dim: int) -> np.ndarray:
    """The permutation symbol (a tensor density, not a tensor)."""
    out = zeros(dim, dim)
    for p in itertools.permutations(range(dim)):
        out[p] = _perm_sign(p)
    return out
