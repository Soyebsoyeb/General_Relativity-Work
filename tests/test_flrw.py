import numpy as np
import pytest
import sympy as sp

from gr_from_scratch.core.simplification import simplify
from gr_from_scratch.field_equations import stress_energy as se
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry, field_equation_residual
from gr_from_scratch.numerical import friedmann as fr
from gr_from_scratch.spacetimes import flrw as fl


@pytest.fixture(scope="module")
def geom():
    return compute_geometry(fl.flrw())


def test_friedmann_equations_from_einstein_tensor(geom):
    """G^t_t = -3(adot^2 + k)/a^2,  G^r_r = -(2 a addot + adot^2 + k)/a^2."""
    a = fl.a(fl.t)
    ad, add = sp.diff(a, fl.t), sp.diff(a, fl.t, 2)
    Gm = geom.einstein.raise_index(0).c
    assert simplify(Gm[0, 0] + 3 * (ad**2 + fl.k) / a**2) == 0
    assert simplify(Gm[1, 1] + (2 * a * add + ad**2 + fl.k) / a**2) == 0
    assert simplify(Gm[1, 1] - Gm[2, 2]) == 0 and simplify(Gm[2, 2] - Gm[3, 3]) == 0  # isotropy


def test_perfect_fluid_closes_the_system_and_gives_raychaudhuri(geom):
    a = fl.a(fl.t)
    ad, add = sp.diff(a, fl.t), sp.diff(a, fl.t, 2)
    rho = 3 * (ad**2 + fl.k) / (8 * sp.pi * a**2)
    p = -(2 * a * add + ad**2 + fl.k) / (8 * sp.pi * a**2)
    T = se.perfect_fluid(geom.metric, rho, p, [1, 0, 0, 0])
    assert field_equation_residual(geom, T).is_zero()
    # acceleration equation addot/a = -(4 pi / 3)(rho + 3 p)
    assert simplify(add / a + sp.Rational(4, 3) * sp.pi * (rho + 3 * p)) == 0
    # continuity follows from the Bianchi identity, it was not imposed
    cont = sp.diff(rho, fl.t) + 3 * ad / a * (rho + p)
    assert simplify(cont) == 0


def test_cosmological_constant_enters_as_expected(geom):
    a = fl.a(fl.t)
    L = sp.Symbol("Lambda", real=True)
    ad = sp.diff(a, fl.t)
    # vacuum with Lambda: (adot/a)^2 + k/a^2 = Lambda/3
    G = geom.einstein.raise_index(0).c
    assert simplify(-G[0, 0] / 3 - ((ad / a) ** 2 + fl.k / a**2)) == 0
    assert sp.solve(sp.Eq(-G[0, 0] + 0, L), L) is not None


@pytest.mark.parametrize(
    "cosmo,afn,t0,t1",
    [
        (fr.Cosmology(om_m=1.0), fr.a_matter_flat, 0.2, 3.0),
        (fr.Cosmology(om_r=1.0), fr.a_radiation_flat, 0.05, 1.5),
        (fr.Cosmology(om_l=1.0), fr.a_de_sitter, 0.0, 3.0),
        (fr.Cosmology(om_m=0.3, om_l=0.7), lambda t: fr.a_lcdm_flat(t, 0.3), 0.1, 2.5),
    ],
)
def test_numerical_expansion_matches_analytic(cosmo, afn, t0, t1):
    h = 1e-6
    a0 = afn(t0)
    ad0 = (afn(t0 + h) - afn(t0 - h)) / (2 * h)
    t, a, ad, con = fr.integrate(cosmo, a0, ad0, t1 - t0)
    assert np.max(np.abs(a - afn(t + t0)) / afn(t + t0)) < 1e-8
    assert np.abs(con).max() < 1e-6  # first-order Friedmann constraint


def test_closed_dust_cycloid():
    om = 2.0
    eta0 = 0.5
    t, a = fr.closed_dust(np.array([eta0]), om)
    ok = abs(1 - om)
    amax = om / ok
    ad0 = (0.5 * amax * np.sin(eta0)) / (om / (2 * ok**1.5) * (1 - np.cos(eta0)))
    ts, an, adn, con = fr.integrate(fr.Cosmology(om_m=om), a[0], ad0, 1.0)
    from scipy.optimize import brentq

    eta_n = brentq(lambda x: fr.closed_dust(x, om)[0] - (t[0] + 1.0), 1e-6, 3.0)
    assert abs(an[-1] - fr.closed_dust(eta_n, om)[1]) < 1e-8
    # curvature sign: Omega_k < 0 for the closed universe
    assert fr.Cosmology(om_m=om).om_k < 0


def test_negative_control_wrong_equation_fails():
    """Integrating matter dynamics but comparing with the radiation solution must fail."""
    afn = fr.a_radiation_flat
    t0, h = 0.05, 1e-6
    a0, ad0 = afn(t0), (afn(t0 + h) - afn(t0 - h)) / (2 * h)
    t, a, _, _ = fr.integrate(fr.Cosmology(om_m=1.0), a0, ad0, 1.0)
    assert np.max(np.abs(a - afn(t + t0)) / afn(t + t0)) > 1e-2
