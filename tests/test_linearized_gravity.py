import numpy as np
import pytest
import sympy as sp

from gr_from_scratch.core.metric import Metric
from gr_from_scratch.field_equations import linearized_gravity as lg
from gr_from_scratch.field_equations.perturbation_series import eps, order_coefficient, perturbative_geometry
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry

t, x, y, z = sp.symbols("t x y z", real=True)
X = [t, x, y, z]


def simp(e):
    return sp.simplify(e)


def test_linearized_riemann_is_gauge_invariant_with_arbitrary_functions():
    h = sp.Matrix(4, 4, lambda a, b: sp.Function(f"h{min(a, b)}{max(a, b)}")(*X))
    xi = [sp.Function(f"xi{i}")(*X) for i in range(4)]
    h2 = lg.gauge_transform(h, xi, X)
    R1, R2 = lg.linearized_riemann(h, X), lg.linearized_riemann(h2, X)
    assert all(simp(a - b) == 0 for a, b in zip(R1.ravel(), R2.ravel(), strict=True))
    # but the perturbation itself is not invariant
    assert any(simp(a - b) != 0 for a, b in zip(h, h2, strict=True))


def test_linearized_einstein_matches_first_order_of_full_engine():
    """Independent check: expand the full Einstein tensor of eta + eps h to first order."""
    h = sp.zeros(4, 4)
    h[0, 0] = sp.sin(x) * t
    h[1, 1] = z**2
    h[0, 2] = h[2, 0] = sp.cos(t) * x
    h[3, 3] = y * z
    full = Metric(X, sp.diag(-1, 1, 1, 1) + eps * h, simplifier=lambda e: e)
    G = compute_geometry(full).einstein.c
    G1_full = sp.Matrix(4, 4, lambda a, b: sp.diff(G[a, b], eps).subs(eps, 0))
    G1 = lg.linearized_einstein(h, X)
    assert all(simp(a - b) == 0 for a, b in zip(G1_full, G1, strict=True))


def test_linearized_riemann_matches_full_engine_first_order():
    h = sp.zeros(4, 4)
    h[1, 1] = sp.sin(t) * z
    h[2, 3] = h[3, 2] = x * y
    full = Metric(X, sp.diag(-1, 1, 1, 1) + eps * h, simplifier=lambda e: e)
    Rfull = compute_geometry(full).riemann.lower_index(0).c
    R1 = lg.linearized_riemann(h, X)
    for idx in np.ndindex(4, 4, 4, 4):
        d = sp.diff(sp.sympify(Rfull[idx]), eps).subs(eps, 0)
        assert simp(d - R1[idx]) == 0


def test_tt_plane_wave_is_vacuum_solution_only_on_light_cone():
    w = sp.Symbol("omega", positive=True)
    A, B = sp.symbols("A_p A_x", real=True)
    h_on = lg.plane_wave_tt(X, A, B, w)
    tr, h0, div = lg.tt_conditions(h_on, X)
    assert tr == 0 and all(c == 0 for c in h0) and all(c == 0 for c in div)
    assert all(simp(e) == 0 for e in lg.linearized_einstein(h_on, X))
    k = sp.Symbol("k", positive=True)
    h_off = lg.plane_wave_tt(X, A, B, w, k=k)
    G = lg.linearized_einstein(h_off, X)
    assert any(simp(e) != 0 for e in G)
    assert all(simp(e.subs(k, w)) == 0 for e in G)


def test_harmonic_gauge_reduces_einstein_to_wave_operator():
    """In Lorenz gauge, G^(1)_{ab} = -(1/2) box hbar_{ab}."""
    w = sp.Symbol("omega", positive=True)
    h = sp.zeros(4, 4)
    # harmonic-gauge perturbation: only h_xy(t, z) -> traceless, divergence-free
    h[1, 2] = h[2, 1] = sp.sin(w * t - 2 * w * z)
    assert all(simp(e) == 0 for e in lg.harmonic_gauge_residual(h, X))
    hb = lg.trace_reverse(h)
    G = lg.linearized_einstein(h, X)
    for a in range(4):
        for b in range(4):
            assert simp(G[a, b] + lg.box(hb[a, b], X) / 2) == 0


def test_second_order_series_agrees_with_taylor_of_exact_einstein_tensor():
    f = sp.sin(t - z)
    h = sp.zeros(4, 4)
    h[1, 1], h[2, 2] = f, -f
    geo = perturbative_geometry(h, X, order=2)
    exact = Metric(X, sp.diag(-1, 1 + eps * f, 1 - eps * f, 1), simplifier=sp.simplify)
    Gx = compute_geometry(exact).einstein.c
    for idx in [(0, 0), (0, 3), (3, 3), (1, 1), (2, 2)]:
        for k in (0, 1, 2):
            taylor = sp.simplify(sp.diff(Gx[idx], eps, k).subs(eps, 0) / sp.factorial(k))
            assert simp(order_coefficient(geo["G"][idx], k) - taylor) == 0
    # first order vanishes (vacuum wave), second order does not: gravitational waves source curvature
    assert all(simp(order_coefficient(geo["G"][i, j], 1)) == 0 for i in range(4) for j in range(4))
    assert simp(order_coefficient(geo["G"][0, 0], 2)) != 0


@pytest.mark.slow
def test_quadratic_action_euler_lagrange_equals_linearized_einstein():
    from gr_from_scratch.quantum.perturbative_gravity import (
        compare_with_linearized_einstein,
        euler_lagrange_components,
        fierz_pauli_lagrangian,
    )

    res, L2, h, comps = compare_with_linearized_einstein(X, depends_on=[t, x])
    assert all(v == 0 for v in res.values())
    fp = fierz_pauli_lagrangian(h, X)
    el_fp = euler_lagrange_components(fp, comps, [t, x])
    el_2 = euler_lagrange_components(L2, comps, [t, x])
    assert all(sp.simplify(sp.expand(el_2[k] - el_fp[k] / 2)) == 0 for k in el_2)
