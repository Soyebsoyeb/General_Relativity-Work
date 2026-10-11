import sympy as sp

from gr_from_scratch.core.simplification import simplify
from gr_from_scratch.core.tensor import Tensor
from gr_from_scratch.dynamics.conserved_quantities import conservation_rate, is_conserved_along_geodesics
from gr_from_scratch.dynamics.effective_potential import circular_orbit_analysis, radial_equation
from gr_from_scratch.dynamics.variational_principles import euler_lagrange_acceleration, velocity_symbols
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry, derive_static_spherical_vacuum
from gr_from_scratch.geometry.invariants import kretschmann
from gr_from_scratch.geometry.killing import is_killing, killing_residual
from gr_from_scratch.geometry.lie_derivatives import lie_derivative
from gr_from_scratch.spacetimes.schwarzschild import M, eddington_finkelstein, phi, r, schwarzschild, theta


def test_vacuum_equations_hold_and_riemann_is_nonzero(schw_geom):
    assert schw_geom.einstein.is_zero()
    assert schw_geom.ricci.is_zero()
    assert not schw_geom.riemann.is_zero()


def test_static_spherical_derivation_recovers_schwarzschild():
    d = derive_static_spherical_vacuum()
    assert d["check_Gtt_form"] == 0 and d["check_Grr_minus_Gtt_form"] == 0
    sol = d["u=1/B solution"]  # u = 1/B = 1 + C1/r from G^t_t = 0
    C1 = sp.Symbol("C1")
    f = 1 - 2 * M / r
    assert simplify(sol.rhs.subs(C1, -2 * M) - f) == 0
    # the integration constant is fixed by Newtonian matching, C1 = -2M; the full solution then
    # solves every Einstein equation, not only G^t_t = 0
    assert compute_geometry(schwarzschild()).einstein.is_zero()


def test_eddington_finkelstein_is_regular_vacuum_and_same_invariant():
    ef = eddington_finkelstein()
    f = 1 - 2 * M / r
    assert simplify(ef.g[0, 0] + f) == 0 and ef.g[0, 1] == 1 and ef.g[1, 1] == 0
    assert simplify(ef.det + r**4 * sp.sin(theta) ** 2) == 0  # finite and non-degenerate at r = 2M
    g = compute_geometry(ef)
    assert g.ricci.is_zero()
    assert simplify(kretschmann(g.riemann) - 48 * M**2 / r**6) == 0


def test_killing_vectors_and_conservation():
    s = schwarzschild()
    from gr_from_scratch.geometry.connection import christoffel

    conn = christoffel(s)
    assert is_killing(conn, [1, 0, 0, 0])
    assert is_killing(conn, [0, 0, 0, 1])
    # rotation generator mixing theta and phi
    ph = phi
    rot = [0, 0, sp.sin(ph), sp.cos(ph) / sp.tan(theta)]
    assert is_killing(conn, rot)
    # negative controls: radial translation and a scaled Killing field are not Killing
    assert not is_killing(conn, [0, 1, 0, 0])
    assert not is_killing(conn, [r, 0, 0, 0])


def test_conservation_law_from_definition():
    s = schwarzschild()
    assert is_conserved_along_geodesics(s, [1, 0, 0, 0])
    assert is_conserved_along_geodesics(s, [0, 0, 0, 1])
    assert not is_conserved_along_geodesics(s, [0, 1, 0, 0])
    assert conservation_rate(s, [0, 1, 0, 0]) != 0


def test_lie_derivative_of_metric_equals_killing_residual():
    s = schwarzschild()
    from gr_from_scratch.geometry.connection import christoffel

    conn = christoffel(s)
    xi = [sp.Symbol("c1") * 0 + r * sp.sin(theta), 0, 1 / r, sp.cos(theta)]  # arbitrary field
    L = lie_derivative(Tensor.metric_tensor(s), xi)
    res = killing_residual(conn, xi)
    assert all(simplify(a - b) == 0 for a, b in zip(L.c.ravel(), res.ravel(), strict=True))


def test_euler_lagrange_matches_christoffel_geodesic(schw_geom):
    m = schw_geom.metric
    u = velocity_symbols(m)
    el = euler_lagrange_acceleration(m, u)
    G = schw_geom.conn.gamma
    for s_ in range(4):
        geo = -sum(G[s_, a, b] * u[a] * u[b] for a in range(4) for b in range(4))
        assert simplify(el[s_] - geo) == 0


def test_effective_potential_photon_sphere_and_isco():
    s = schwarzschild()
    _, V = radial_equation(s, 1)
    f = 1 - 2 * M / r
    L = sp.Symbol("L", real=True)
    assert simplify(V - f * (1 + L**2 / r**2)) == 0
    tl = circular_orbit_analysis(s, 1)
    nl = circular_orbit_analysis(s, 0)
    assert tl["isco"] == [6 * M]
    assert nl["photon_sphere"] == [3 * M]
    # ISCO energy E = sqrt(8/9)
    E2_isco = simplify(tl["E2"].subs(r, 6 * M))
    assert E2_isco == sp.Rational(8, 9)
