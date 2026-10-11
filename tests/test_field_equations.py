import pytest
import sympy as sp

from gr_from_scratch.core.metric import Metric
from gr_from_scratch.core.simplification import simplify
from gr_from_scratch.field_equations import einstein_hilbert as eh
from gr_from_scratch.field_equations import stress_energy as se
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry, field_equation_residual
from gr_from_scratch.geometry.covariant_derivative import divergence
from gr_from_scratch.spacetimes import flrw as fl
from gr_from_scratch.spacetimes import reissner_nordstrom as rn
from gr_from_scratch.spacetimes.schwarzschild import M, r, theta


@pytest.fixture(scope="module")
def rn_geom():
    return compute_geometry(rn.reissner_nordstrom())


def test_reissner_nordstrom_solves_einstein_maxwell(rn_geom):
    m = rn_geom.metric
    F = se.field_strength(m, rn.vector_potential())
    T = se.electromagnetic(m, F)
    res = field_equation_residual(rn_geom, T)
    assert res.is_zero()
    assert not T.is_zero()  # genuinely non-vacuum
    assert all(simplify(x) == 0 for x in se.maxwell_divergence(m, F))


def test_einstein_maxwell_fails_for_wrong_charge_coefficient():
    """Negative controls: change Q^2/r^2 -> 2 Q^2/r^2 in f, or the potential normalisation."""
    t, th, ph = sp.symbols("t theta phi", real=True)
    Q = rn.Q
    f_bad = 1 - 2 * M / r + 2 * Q**2 / r**2
    bad = Metric([t, r, theta, ph], sp.diag(-f_bad, 1 / f_bad, r**2, r**2 * sp.sin(theta) ** 2))
    g = compute_geometry(bad)
    F = se.field_strength(bad, rn.vector_potential())
    assert not field_equation_residual(g, se.electromagnetic(bad, F)).is_zero()
    good = compute_geometry(rn.reissner_nordstrom())
    F2 = se.field_strength(good.metric, [2 * c for c in rn.vector_potential()])
    assert not field_equation_residual(good, se.electromagnetic(good.metric, F2)).is_zero()


def test_reissner_nordstrom_horizons():
    rm, rp = rn.horizons()
    f = rn.lapse()
    assert simplify(f.subs(r, rp)) == 0 and simplify(f.subs(r, rm)) == 0


def test_perfect_fluid_conservation_is_the_continuity_equation():
    """nabla_mu T^{mu nu} = 0 on FLRW reduces to rho' + 3 H (rho + p) = 0."""
    g = compute_geometry(fl.flrw())
    m = g.metric
    rho, p = sp.Function("rho")(fl.t), sp.Function("p")(fl.t)
    T = se.perfect_fluid(m, rho, p, [1, 0, 0, 0]).raise_index(0).raise_index(1)
    div = divergence(T, g.conn, 0).c
    a = fl.a(fl.t)
    expected = sp.diff(rho, fl.t) + 3 * sp.diff(a, fl.t) / a * (rho + p)
    assert simplify(div[0] - expected) == 0
    assert all(simplify(div[i]) == 0 for i in (1, 2, 3))


def test_scalar_field_stress_energy_conserved_on_shell_in_flat_space():
    t, x, y, z = sp.symbols("t x y z", real=True)
    m = Metric([t, x, y, z], sp.diag(-1, 1, 1, 1))
    phi = sp.cos(t - z)  # solves box phi = 0
    T = se.scalar_field(m, phi).raise_index(0).raise_index(1)
    g = compute_geometry(m)
    assert all(simplify(c) == 0 for c in divergence(T, g.conn, 0).c)


# ---- variational identities (3D background with non-vanishing Einstein tensor) ------
@pytest.fixture(scope="module")
def background():
    t, x, y = sp.symbols("t x y", real=True)
    m = Metric([t, x, y], sp.diag(-(1 + x**2), 1, 2 + y**2))
    h = [[y, x * y, 0], [x * y, x**2, 1], [0, 1, t * y]]
    assert not compute_geometry(m).einstein.is_zero()
    return m, h


def test_variation_of_determinant_and_inverse(background):
    m, h = background
    assert eh.det_variation_residual(m, h) == 0
    assert eh.inverse_variation_residual(m, h) == sp.zeros(3)


def test_boundary_term_is_a_total_divergence(background):
    assert eh.boundary_term_residual(*background) == 0


def test_ricci_scalar_variation(background):
    assert eh.ricci_scalar_variation_residual(*background) == 0


def test_ricci_scalar_variation_negative_control(background):
    m, h = background
    a, b, c = eh.ricci_scalar_variation_terms(m, h)
    lhs = eh._d_eps0(eh._scalar_curvature(eh._perturbed(m, h)))
    assert sp.simplify(lhs - (a + b + c)) != 0  # wrong sign of the box term


@pytest.mark.slow
def test_full_einstein_hilbert_variation_with_cosmological_constant(background):
    m, h = background
    assert eh.full_action_variation_residual(m, h, sp.Symbol("Lambda")) == 0


def test_reduced_action_reproduces_einstein_tensor_components():
    assert eh.reduced_static_spherical_equations() == (0, 0)


def test_cosmological_constant_de_sitter_static_patch():
    """G + Lambda g = 0 for the static de Sitter patch f = 1 - Lambda r^2/3."""
    t, ph = sp.symbols("t phi", real=True)
    L = sp.Symbol("Lambda", positive=True)
    f = 1 - L * r**2 / 3
    ds = Metric([t, r, theta, ph], sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(theta) ** 2))
    g = compute_geometry(ds)
    assert field_equation_residual(g, None, cosmological_constant=L).is_zero()
    assert not g.einstein.is_zero()
