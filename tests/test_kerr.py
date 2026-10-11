import numpy as np
import pytest
import sympy as sp

from gr_from_scratch.dynamics.geodesics import GeodesicSystem
from gr_from_scratch.geometry.connection import christoffel
from gr_from_scratch.geometry.invariants import kretschmann_at
from gr_from_scratch.geometry.killing import killing_residual
from gr_from_scratch.geometry.riemann import riemann
from gr_from_scratch.spacetimes import kerr as K

SPIN = sp.Rational(3, 5)


@pytest.fixture(scope="module")
def kerr_num():
    m = K.kerr().subs({K.M: 1, K.a: SPIN})
    conn = christoffel(m)
    return m, conn, riemann(conn, simplify=False)


def test_inverse_metric_consistency():
    m = K.kerr()  # constructor verifies g g^-1 = 1 and raises otherwise
    assert (m.g * m.ginv).applyfunc(sp.simplify) == sp.eye(4)


@pytest.mark.parametrize("pt", [(3.7, 0.9), (2.4, 1.3)])
def test_kerr_is_ricci_flat_at_sample_points(kerr_num, pt):
    m, conn, R = kerr_num
    sub = {K.r: sp.Float(pt[0], 30), K.theta: sp.Float(pt[1], 30)}
    vals = np.array(
        [
            [[[sp.N(R.c[i, j, k, l].subs(sub), 30) for l in range(4)] for k in range(4)] for j in range(4)]
            for i in range(4)
        ],
        dtype=object,
    )
    ricci = sum(vals[i, :, i, :] for i in range(4))
    assert max(abs(complex(x)) for x in ricci.ravel()) < 1e-25


@pytest.mark.parametrize("pt", [(3.7, 0.9), (6.0, 1.4)])
def test_kretschmann_matches_independent_closed_form(kerr_num, pt):
    m, conn, R = kerr_num
    sub = {K.r: sp.Float(pt[0], 30), K.theta: sp.Float(pt[1], 30)}
    num = kretschmann_at(R, sub)
    ref = float(K.kretschmann_closed_form(1, SPIN).subs(sub))
    assert abs(num - ref) / abs(ref) < 1e-12


def test_riemann_is_not_trivially_zero_but_ricci_is(kerr_num):
    m, conn, R = kerr_num
    sub = {K.r: sp.Float(3.7, 30), K.theta: sp.Float(0.9, 30)}
    assert abs(kretschmann_at(R, sub)) > 1e-6


def test_horizons_and_ergosurface():
    a = sp.Rational(3, 5)
    rm, rp = K.horizons(1, a)
    assert (rm, rp) == (sp.Rational(1, 5), sp.Rational(9, 5))
    assert K.delta(1, a).subs(K.r, rp) == 0
    # g_tt vanishes on the ergosurface for every theta; it touches the horizon at the poles
    m = K.kerr()
    for th in (sp.Rational(1, 3), sp.Rational(1, 2), 1):
        re = K.ergosurface(1, a, th)
        gtt = m.g[0, 0].subs({K.M: 1, K.a: a, K.theta: th, K.r: re})
        assert abs(float(sp.N(gtt))) < 1e-12
    assert K.ergosurface(1, a, 0) == rp
    assert K.ergosurface(1, a, sp.pi / 2) == 2  # equator: r = 2M


def test_killing_vectors_numerically(kerr_num):
    m, conn, R = kerr_num
    sub = {K.r: sp.Float(3.7, 30), K.theta: sp.Float(0.9, 30)}
    for xi in ([1, 0, 0, 0], [0, 0, 0, 1]):
        res = killing_residual(conn, xi)
        assert max(abs(float(sp.N(sp.sympify(e).subs(sub)))) for e in res.ravel()) < 1e-25
    res = killing_residual(conn, [0, 1, 0, 0])
    assert max(abs(float(sp.N(sp.sympify(e).subs(sub)))) for e in res.ravel()) > 1e-3


@pytest.mark.parametrize("prograde", [True, False])
def test_circular_orbit_formulas_satisfy_geodesic_equation(kerr_num, prograde):
    m = kerr_num[0]
    gs = GeodesicSystem(m)
    r0, a = 7.0, 0.6
    E, L, Om = K.circular_orbit(r0, a, 1.0, prograde)
    x0 = np.array([0, r0, np.pi / 2, 0.0])
    g = gs.g(x0)
    ut = 1 / np.sqrt(-(g[0, 0] + 2 * g[0, 3] * Om + g[3, 3] * Om**2))
    u = np.array([ut, 0, 0, Om * ut])
    assert abs(gs.acceleration(x0, u)[1]) < 1e-13
    assert abs(-(g[0, 0] * u[0] + g[0, 3] * u[3]) - E) < 1e-12
    assert abs((g[3, 0] * u[0] + g[3, 3] * u[3]) - L) < 1e-12
    assert abs(gs.norm(x0, u) + 1) < 1e-12


def test_isco_known_values():
    assert abs(K.isco_radius(0.0) - 6.0) < 1e-12
    assert abs(K.isco_radius(1.0, prograde=True) - 1.0) < 1e-6
    assert abs(K.isco_radius(1.0, prograde=False) - 9.0) < 1e-12
    assert K.isco_radius(0.5, prograde=True) < 6 < K.isco_radius(0.5, prograde=False)


def test_circular_orbit_perturbation_stability_straddles_isco(kerr_num):
    """Perturbed orbits stay bounded outside the ISCO and plunge inside it (a = 0.6, prograde)."""
    gs = GeodesicSystem(kerr_num[0])
    r_isco = K.isco_radius(0.6, 1.0, True)

    def run(r0, T):
        E, L, Om = K.circular_orbit(r0, 0.6, 1.0, True)
        x0 = np.array([0, r0, np.pi / 2, 0.0])
        g = gs.g(x0)
        ut = 1 / np.sqrt(-(g[0, 0] + 2 * g[0, 3] * Om + g[3, 3] * Om**2))
        u = np.array([ut, 0, 0, Om * ut])
        # same E, L; add a small radial kick
        u[1] = 1e-3
        ev = lambda lam, y: y[1] - 1.0 * 1.6  # reaches near the horizon region
        ev.terminal = True
        tr = gs.integrate(x0, u, T, kappa=1.0, rtol=1e-10, atol=1e-12, events=[ev])
        return tr.x[:, 1]

    outside = run(1.25 * r_isco, 3000.0)
    inside = run(0.9 * r_isco, 3000.0)
    assert outside.max() < 2.0 * r_isco and outside.min() > 0.5 * r_isco
    assert inside.min() < 0.7 * r_isco  # runs away inward


def test_frame_dragging_far_field_is_lense_thirring():
    m = K.kerr()
    gm = m.g.subs({K.theta: sp.pi / 2, K.M: 1, K.a: sp.Rational(3, 5)})
    r = 1000
    om = float((-gm[0, 3] / gm[3, 3]).subs(K.r, r))
    assert abs(om - 2 * 0.6 / r**3) / om < 1e-4 and om > 0
