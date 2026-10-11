"""Simplification helpers and honest zero-testing.

Symbolic simplification is not a decision procedure: failing to reduce an
expression to zero does not prove it is non-zero. ``check_zero`` therefore
reports how a conclusion was reached:

``exact``    the expression simplified to the literal 0.
``numeric``  simplification failed but random sampling gave residuals below
             tolerance. This is evidence, not proof.
``nonzero``  a sample exceeded tolerance (the expression is certainly not 0).
"""

from __future__ import annotations

import random
from dataclasses import dataclass

import numpy as np
import sympy as sp

_TRANSCENDENTAL = (sp.sin, sp.cos, sp.tan, sp.exp, sp.log, sp.sinh, sp.cosh, sp.tanh)


def simplify(expr):
    """Cheap rational simplification first, full ``simplify`` only if needed."""
    expr = sp.sympify(expr)
    if expr.is_number:
        return sp.simplify(expr)
    out = sp.cancel(sp.together(expr))
    if out.has(*_TRANSCENDENTAL) or out.has(sp.Pow):
        out = sp.simplify(out)
    if out != 0 and out.has(sp.sin, sp.cos, sp.tan):
        # multiple-angle identities (sin 2x tan x + cos 2x - 1 = 0) need expansion first
        alt = sp.simplify(sp.expand_trig(out))
        if alt == 0 or sp.count_ops(alt) < sp.count_ops(out):
            out = alt
    return out


def identity(expr):
    """No-op simplifier for expressions that will only be evaluated numerically."""
    return expr


@dataclass(frozen=True)
class ZeroCheck:
    status: str  # "exact", "numeric" or "nonzero"
    residual: float

    @property
    def ok(self) -> bool:
        return self.status in ("exact", "numeric")

    def __bool__(self) -> bool:
        return self.ok


def check_zero(expr, domain=None, tol=1e-9, samples=8, seed=0, simplifier=simplify) -> ZeroCheck:
    """Decide whether ``expr`` vanishes identically.

    ``domain`` maps free symbols to ``(lo, hi)`` sampling intervals; symbols
    without an entry are sampled from (0.5, 3.0).
    """
    expr = sp.sympify(expr)
    if expr == 0:
        return ZeroCheck("exact", 0.0)
    s = simplifier(expr)
    if s == 0:
        return ZeroCheck("exact", 0.0)
    rng = random.Random(seed)
    domain = domain or {}
    worst = 0.0
    syms = sorted(s.free_symbols, key=lambda x: x.name)
    for _ in range(samples):
        pt = {}
        for sym in syms:
            lo, hi = domain.get(sym, (0.5, 3.0))
            pt[sym] = sp.Float(rng.uniform(lo, hi), 30)
        worst = max(worst, abs(complex(sp.N(s.subs(pt), 30))))
    return ZeroCheck("numeric" if worst < tol else "nonzero", worst)


def check_zero_array(arr, **kw) -> ZeroCheck:
    """Apply ``check_zero`` to every entry; the weakest status wins."""
    rank = {"exact": 0, "numeric": 1, "nonzero": 2}
    status, worst = "exact", 0.0
    for e in np.asarray(arr, dtype=object).ravel():
        r = check_zero(e, **kw)
        worst = max(worst, r.residual)
        if rank[r.status] > rank[status]:
            status = r.status
    return ZeroCheck(status, worst)
