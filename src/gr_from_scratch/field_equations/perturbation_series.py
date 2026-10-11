"""Perturbation of geometric quantities in powers of epsilon, g = eta + eps h.

Everything is expanded with explicit truncation, so the order is a parameter
that appears in every result. ``order=2`` keeps terms through eps^2. Results
must be read as truncated series, not as solutions of the nonlinear equations.
"""

from __future__ import annotations

import sympy as sp

eps = sp.Symbol("epsilon")


def truncate(expr, order: int):
    expr = sp.expand(expr)
    return sum(eps**k * expr.coeff(eps, k) for k in range(order + 1))


def perturbative_geometry(h, x, order: int = 2, eta=None):
    """Truncated series of Gamma, Riemann, Ricci, R, sqrt(-g) and G for g = eta + eps h."""
    n = len(x)
    eta = eta if eta is not None else sp.diag(-1, *([1] * (n - 1)))
    eta_inv = eta.inv()
    H = sp.Matrix(h)
    # g^{-1} = sum_k (-eps eta^{-1} h)^k eta^{-1}, truncated
    ginv, term = sp.zeros(n), eta_inv
    for _ in range(order + 1):
        ginv += term
        term = (-eps * eta_inv * H) * term
    ginv = ginv.applyfunc(lambda e: truncate(e, order))
    g = eta + eps * H
    gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for m in range(n):
            for v in range(n):
                s = (
                    sum(
                        ginv[r, s_] * (sp.diff(g[v, s_], x[m]) + sp.diff(g[m, s_], x[v]) - sp.diff(g[m, v], x[s_]))
                        for s_ in range(n)
                    )
                    / 2
                )
                gam[r][m][v] = truncate(s, order)
    ric = sp.zeros(n)
    for s_ in range(n):
        for v in range(n):
            tot = 0
            for r in range(n):
                tot += sp.diff(gam[r][v][s_], x[r]) - sp.diff(gam[r][r][s_], x[v])
                for lam in range(n):
                    tot += gam[r][r][lam] * gam[lam][v][s_] - gam[r][v][lam] * gam[lam][r][s_]
            ric[s_, v] = truncate(tot, order)
    R = truncate(sum(ginv[a, b] * ric[a, b] for a in range(n) for b in range(n)), order)
    # sqrt(-det g) = exp( (1/2) tr log(1 + eps eta^-1 h) ) expanded to the requested order
    X = eta_inv * H
    tr1 = X.trace()
    logdet = eps * tr1
    if order >= 2:
        logdet += -(eps**2) * (X * X).trace() / 2
    if order >= 3:
        logdet += eps**3 * (X**3).trace() / 3
    sqrtg = truncate(sum((logdet / 2) ** k / sp.factorial(k) for k in range(order + 1)), order)
    G = (ric - sp.Rational(1, 2) * R * g).applyfunc(lambda e: truncate(e, order))
    return {"gamma": gam, "ricci": ric, "R": R, "sqrt_g": sqrtg, "G": G, "ginv": ginv}


def order_coefficient(expr, k: int):
    return sp.expand(expr).coeff(eps, k)
