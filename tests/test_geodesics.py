import numpy as np
import pytest

from gr_from_scratch.dynamics import observables as ob
from gr_from_scratch.dynamics.geodesics import GeodesicSystem, event_function
from gr_from_scratch.numerical.convergence import observed_order, richardson_order
from gr_from_scratch.numerical.integration import rk4_fixed
from gr_from_scratch.spacetimes import minkowski as mk


@pytest.fixture(scope="module")
def gs():
    return ob.schwarzschild_system()


def test_minkowski_geodesic_is_straight_line():
    g = GeodesicSystem(mk.minkowski())
    x0, u0 = np.zeros(4), np.array([1.2, 0.3, -0.2, 0.1])
    tr = g.integrate(x0, u0, 5.0, kappa=-(u0 @ g.g(x0) @ u0))
    np.testing.assert_allclose(tr.x[-1], x0 + 5.0 * u0, atol=1e-10)


def test_flat_polar_coordinates_give_straight_line():
    g = GeodesicSystem(mk.euclidean_polar())
    x0 = np.array([2.0, 0.5])
    u0 = np.array([0.3, 0.2])
    tr = g.integrate(x0, u0, 3.0, kappa=-1.0, rtol=1e-12, atol=1e-13)
    cart0 = np.array([2 * np.cos(0.5), 2 * np.sin(0.5)])
    vel = np.array([0.3 * np.cos(0.5) - 2 * 0.2 * np.sin(0.5), 0.3 * np.sin(0.5) + 2 * 0.2 * np.cos(0.5)])
    r, ph = tr.x[-1]
    np.testing.assert_allclose([r * np.cos(ph), r * np.sin(ph)], cart0 + 3.0 * vel, atol=1e-9)


def test_future_directed_normalisation_and_errors(gs):
    x = np.array([0, 10.0, np.pi / 2, 0])
    for kappa in (1.0, 0.0):
        u = gs.future_directed_u(x, [0.05, 0, 0.01], kappa)
        assert u[0] > 0
        assert abs(gs.norm(x, u) + kappa) < 1e-13
    # inside the horizon no timelike curve is static: u^0 cannot be real for zero spatial velocity
    with pytest.raises(ValueError):
        gs.future_directed_u(np.array([0, 1.5, np.pi / 2, 0]), [0, 0, 0], 1.0)


def test_timelike_orbit_normalisation_and_conserved_quantities(gs):
    E, L, up, ua = ob._orbit_constants(20.0, 0.5, 1.0)
    rp = 1 / up
    f = 1 - 2 / rp
    x0 = np.array([0, rp, np.pi / 2, 0.0])
    u0 = np.array([E / f, 0, 0, L / rp**2])
    tr = gs.integrate(x0, u0, 600.0, rtol=1e-12, atol=1e-13)
    assert tr.success
    assert np.abs(tr.constraint).max() < 1e-9
    Et, Lp = gs.conserved_function([1, 0, 0, 0]), gs.conserved_function([0, 0, 0, 1])
    e = np.array([Et(x, u) for x, u in zip(tr.x, tr.u, strict=True)])
    l = np.array([Lp(x, u) for x, u in zip(tr.x, tr.u, strict=True)])
    assert np.abs(e - e[0]).max() < 1e-9 and np.abs(l - l[0]).max() < 1e-9
    assert abs(-e[0] - E) < 1e-12 and abs(l[0] - L) < 1e-12  # E = -g(d_t,u)


def test_null_geodesic_stays_null(gs):
    x0 = np.array([0, 30.0, np.pi / 2, 0.0])
    u = gs.future_directed_u(x0, [-1.0, 0, 0.03], 0.0)
    tr = gs.integrate(x0, u, 40.0, kappa=0.0, rtol=1e-12, atol=1e-13)
    assert np.abs(tr.constraint).max() < 1e-9


def test_rk4_has_fourth_order_on_schwarzschild(gs):
    E, L, up, ua = ob._orbit_constants(20.0, 0.5, 1.0)
    rp = 1 / up
    f = 1 - 2 / rp
    y0 = np.array([0, rp, np.pi / 2, 0, E / f, 0, 0, L / rp**2])
    T = 100.0
    ref = gs.integrate(y0[:4], y0[4:], T, rtol=1e-13, atol=1e-14).sol.y[:, -1]
    steps = np.array([100, 200, 400, 800])
    errs = [np.linalg.norm(rk4_fixed(gs.rhs, y0, T, int(n))[1][-1] - ref) for n in steps]
    assert 3.6 < observed_order(T / steps, errs) < 4.4
    assert 3.5 < richardson_order(errs[0], errs[1]) < 4.5


def test_adaptive_error_decreases_with_tolerance(gs):
    E, L, up, ua = ob._orbit_constants(20.0, 0.5, 1.0)
    rp = 1 / up
    f = 1 - 2 / rp
    x0, u0 = np.array([0, rp, np.pi / 2, 0.0]), np.array([E / f, 0, 0, L / rp**2])
    ref = gs.integrate(x0, u0, 200.0, rtol=1e-13, atol=1e-14).sol.y[:, -1]
    errs = [np.linalg.norm(gs.integrate(x0, u0, 200.0, rtol=t, atol=t * 1e-3).sol.y[:, -1] - ref) for t in (1e-6, 1e-9)]
    assert errs[1] < errs[0] / 10


def test_event_detection_finds_periapsis(gs):
    E, L, up, ua = ob._orbit_constants(30.0, 0.4, 1.0)
    rp = 1 / up
    f = 1 - 2 / rp
    ev = event_function(5, 0.0, direction=+1)
    tr = gs.integrate(
        np.array([0, rp, np.pi / 2, 0.0]),
        np.array([E / f, 0, 0, L / rp**2]),
        1500.0,
        events=[ev],
        rtol=1e-12,
        atol=1e-13,
    )
    ys = tr.sol.y_events[0]
    assert len(ys) >= 2 and np.allclose(ys[:, 1], rp, rtol=1e-8)


# ---- observables ----------------------------------------------------------------
@pytest.mark.parametrize("b", [20.0, 100.0])
def test_deflection_geodesic_matches_quadrature(b):
    assert abs(ob.deflection_geodesic(b) - ob.deflection_quad(b)) / ob.deflection_quad(b) < 1e-8


def test_deflection_series_coefficients():
    b = 2000.0
    exact = ob.deflection_quad(b)
    assert abs(exact - 4 / b) / exact < 2e-3
    second = (exact - 4 / b) / (1 / b) ** 2
    assert abs(second - 15 * np.pi / 4) / (15 * np.pi / 4) < 5e-3


@pytest.mark.filterwarnings("ignore:The occurrence of roundoff error")
def test_strong_field_deflection_is_logarithmic_with_known_constant():
    """alpha -> -ln(b/b_c - 1) + ln[216 (7 - 4 sqrt 3)] - pi as b -> b_c = 3 sqrt 3 M
    (strong-deflection limit, Bozza 2002)."""
    bc = 3 * np.sqrt(3)
    d1, d2 = 1e-4, 1e-6
    a1, a2 = ob.deflection_quad(bc * (1 + d1)), ob.deflection_quad(bc * (1 + d2))
    assert abs((a2 - a1) - np.log(d1 / d2)) < 2e-3
    const = np.log(216 * (7 - 4 * np.sqrt(3))) - np.pi
    assert abs(a2 + np.log(d2) - const) < 1e-3


@pytest.mark.parametrize("p,e", [(100.0, 0.3), (30.0, 0.5)])
def test_precession_geodesic_matches_quadrature(p, e):
    mean, _ = ob.periapsis_advance_geodesic(p, e, n_orbits=2)
    assert abs(mean - ob.periapsis_advance_quad(p, e)) / mean < 1e-8


def test_precession_approaches_leading_order():
    devs = [abs(ob.periapsis_advance_quad(p, 0.3) / ob.periapsis_advance_leading(p) - 1) for p in (200.0, 2000.0)]
    assert devs[1] < devs[0] / 5 and devs[1] < 0.01


def test_gravitational_redshift_geodesic_vs_exact_and_weak_field():
    exact = ob.redshift_static_exact(10.0, 1e4)
    assert abs(ob.redshift_geodesic(10.0, 1e4) - exact) < 1e-9
    weak = ob.redshift_weak_field(1e3, 1e5)
    assert abs(ob.redshift_static_exact(1e3, 1e5) - weak) < 1e-5


def test_shapiro_excess_matches_weak_field_with_second_order_remainder():
    r1 = 1e6
    d = {r0: ob.shapiro_round_trip_excess(r1, r1, r0) - ob.shapiro_weak_field(r1, r1, r0) for r0 in (100.0, 300.0)}
    assert all(v > 0 for v in d.values())
    # remainder ~ c M^2 / r0 with c of order 4 pi^2: scaled values agree within 10 %
    assert abs(d[100.0] * 100 - d[300.0] * 300) / (d[100.0] * 100) < 0.10


def test_circular_orbit_proper_period_correction():
    gs_ = ob.schwarzschild_system()
    r = 50.0
    ut = 1 / np.sqrt(1 - 3 / r)
    ev = event_function(3, 2 * np.pi, direction=+1, terminal=True)
    tr = gs_.integrate(
        np.array([0, r, np.pi / 2, 0.0]),
        np.array([ut, 0, 0, ut / r**1.5]),
        3 * 2 * np.pi * r**1.5,
        rtol=1e-12,
        atol=1e-13,
        events=[ev],
    )
    tau = tr.sol.t_events[0][0]
    assert abs(tau - 2 * np.pi * r**1.5 * np.sqrt(1 - 3 / r)) < 1e-8
