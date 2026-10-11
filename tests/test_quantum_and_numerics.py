import numpy as np
import pytest
import sympy as sp

from gr_from_scratch.dynamics.gravitational_waves import integrate_deviation, proper_separation
from gr_from_scratch.numerical.boundary_conditions import _riccati_hankel, asymptotic_amplitudes
from gr_from_scratch.numerical.convergence import observed_order
from gr_from_scratch.numerical.integration import rk4_fixed
from gr_from_scratch.quantum import eft
from gr_from_scratch.quantum import klein_gordon as kg


def test_klein_gordon_radial_reduction_symbolic():
    assert kg.schwarzschild_radial_reduction() == 0


def test_tortoise_inversion_roundtrip():
    for r in (2.001, 2.5, 10.0, 100.0, 5e3):
        rs = r + 2 * np.log(r / 2 - 1)
        assert abs(kg.radius_from_tortoise(rs) - r) / r < 1e-10


def test_scattering_flux_conservation_and_ordering():
    s0, s1, s2 = (kg.scattering(l, 0.3, rstar_max=3000.0) for l in (0, 1, 2))
    for s in (s0, s1, s2):
        assert abs(s.flux_residual) < 1e-7
        assert 0 <= s.transmission <= 1 and abs(s.transmission + s.reflection - 1) < 1e-7
    assert s0.transmission > s1.transmission > s2.transmission  # centrifugal barrier


def test_scattering_insensitive_to_outer_boundary():
    a = kg.scattering(2, 0.3, rstar_max=3000.0).transmission
    b = kg.scattering(2, 0.3, rstar_max=6000.0).transmission
    assert abs(a - b) < 1e-6


def test_low_frequency_absorption_approaches_horizon_area():
    """T_0 / (16 M^2 w^2) -> 1 so that sigma_abs -> 16 pi M^2 = horizon area."""
    ratios = [kg.scattering(0, w, rstar_max=120 / w).transmission / (16 * w * w) for w in (0.02, 0.005)]
    assert abs(ratios[1] - 1) < abs(ratios[0] - 1) < 0.2
    assert abs(ratios[1] - 1) < 0.05


def test_low_frequency_power_law_for_l1():
    """T_l ~ w^(2l+2): slope 4 for l = 1."""
    w = np.array([0.02, 0.04])
    T = [kg.scattering(1, float(x), rstar_max=max(3000.0, 120 / x)).transmission for x in w]
    slope = np.log(T[1] / T[0]) / np.log(2)
    assert 3.6 < slope < 4.6


def test_high_frequency_transmission_tends_to_one():
    assert kg.scattering(0, 1.0, rstar_max=3000.0).transmission > 0.9999


def test_asymptotic_amplitude_extraction_on_exact_flat_solutions():
    w, rs = 0.7, 40.0
    psi = 1.0 * np.exp(-1j * w * rs) + 0.5 * np.exp(1j * w * rs)
    dpsi = -1j * w * 1.0 * np.exp(-1j * w * rs) + 1j * w * 0.5 * np.exp(1j * w * rs)
    a_in, a_out = asymptotic_amplitudes(w, rs, psi, dpsi)
    assert abs(a_in - 1) < 1e-12 and abs(a_out - 0.5) < 1e-12
    for l in (0, 1, 3):
        up, dup = _riccati_hankel(l, w * rs)
        um, dum = np.conj(up), np.conj(dup)
        psi = 2 * um + 0.3j * up
        dpsi = w * (2 * dum + 0.3j * dup)
        a_in, a_out = asymptotic_amplitudes(w, rs, psi, dpsi, l=l)
        assert abs(a_in - 2) < 1e-10 and abs(a_out - 0.3j) < 1e-10


def test_riccati_hankel_asymptotics_and_normalisation():
    x = 500.0
    for l in (0, 2, 4):
        up, _ = _riccati_hankel(l, x)
        assert abs(up - np.exp(1j * x)) < (l * (l + 1)) / x + 1e-12


def test_rk4_order_on_harmonic_oscillator():
    rhs = lambda t, y: np.array([y[1], -y[0]])
    ns = np.array([20, 40, 80, 160])
    errs = [abs(rk4_fixed(rhs, [1.0, 0.0], 2.0, int(n))[1][-1, 0] - np.cos(2.0)) for n in ns]
    assert 3.8 < observed_order(2.0 / ns, errs) < 4.2


def test_gw_deviation_equals_exact_proper_separation_to_second_order():
    diffs = []
    for A in (1e-3, 1e-2):
        worst = 0.0
        for d0 in [(1, 0), (0, 1), (np.cos(0.7), np.sin(0.7))]:
            ts, xi = integrate_deviation(d0, A, 0.4 * A, 1.0, 8.0, n=401)
            lin = np.linalg.norm(xi, axis=1)
            ex = proper_separation(ts, d0, A, 0.4 * A, 1.0)
            worst = max(worst, np.abs(lin - ex).max())
        diffs.append(worst)
    assert 1.8 < np.log(diffs[1] / diffs[0]) / np.log(10) < 2.2


def test_gw_ring_area_preserved_to_second_order():
    """TT waves are trace-free, so the ring area changes only at O(h^2).

    The ring is a 400-gon; its area is compared with the same polygon at A = 0
    so that the discretisation error cancels.
    """
    th = np.linspace(0, 2 * np.pi, 400, endpoint=False)

    def polygon_area(A):
        pts = np.array([integrate_deviation((np.cos(t_), np.sin(t_)), A, 0.0, 1.0, 2.0, n=3)[1][-1] for t_ in th])
        return 0.5 * np.abs(np.sum(pts[:, 0] * np.roll(pts[:, 1], -1) - np.roll(pts[:, 0], -1) * pts[:, 1]))

    A = 1e-3
    assert abs(polygon_area(A) / polygon_area(0.0) - 1) < 5 * A**2
    # negative control: a non-trace-free (isotropic) stretch would change the area at O(h)
    assert abs(((1 + A) ** 2) - 1) > 5 * A**2


@pytest.mark.slow
def test_quadratic_curvature_terms_do_not_modify_schwarzschild_but_cubic_does():
    res = eft.reduced_corrections_on_schwarzschild()
    for name in ("R^2", "Ric^2", "Riem^2", "Gauss-Bonnet"):
        assert all(e == 0 for e in res[name]), name
    assert any(e != 0 for e in res["Riem^3"])


def test_eft_expansion_parameter_scales_as_inverse_mass_squared():
    solar = eft.curvature_expansion_parameter(1.989e30)
    big = eft.curvature_expansion_parameter(1.989e36)
    assert abs(np.log10(solar / big) - 12) < 1e-9
    assert solar < 1e-60
    planck_mass = 2.176e-8
    assert 0.1 < eft.curvature_expansion_parameter(planck_mass) < 100  # O(1) at the Planck scale


def test_newtonian_limit_poisson_equation_from_linearized_einstein():
    """h_00 = -2 Phi, h_ij = -2 Phi delta_ij  =>  G^(1)_00 = 2 laplacian(Phi) = 8 pi rho."""
    from gr_from_scratch.field_equations.linearized_gravity import linearized_einstein

    t, x, y, z = sp.symbols("t x y z", real=True)
    Phi = sp.Function("Phi")(x, y, z)
    h = sp.diag(-2 * Phi, -2 * Phi, -2 * Phi, -2 * Phi)
    G = linearized_einstein(h, [t, x, y, z])
    lap = sum(sp.diff(Phi, v, 2) for v in (x, y, z))
    assert sp.simplify(G[0, 0] - 2 * lap) == 0
    # spatial components vanish for static isotropic weak field in vacuum only if lap Phi = 0 consistency
    assert sp.simplify(G[0, 1]) == 0


def test_newtonian_limit_of_coordinate_acceleration():
    """d^2 r/dt^2 of a particle momentarily at rest is -Gamma^r_tt = -M f / r^2 (exact),
    i.e. Newton's -M/r^2 up to the relative correction 2M/r."""
    from gr_from_scratch.dynamics.observables import schwarzschild_system

    gs = schwarzschild_system()
    prev = None
    for r in (1e2, 1e3, 1e4):
        x = np.array([0, r, np.pi / 2, 0.0])
        a_coord = -gs.christoffel(x)[1, 0, 0]
        ratio = a_coord / (-1.0 / r**2)
        assert abs(ratio - (1 - 2 / r)) < 1e-12
        err = abs(ratio - 1)
        assert prev is None or err < prev
        prev = err
    assert prev < 3e-4
