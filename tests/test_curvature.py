import numpy as np
import pytest
import sympy as sp

from gr_from_scratch.core.metric import Metric
from gr_from_scratch.core.simplification import simplify
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry
from gr_from_scratch.geometry import bianchi, invariants
from gr_from_scratch.geometry.weyl import DimensionError, weyl_tensor
from gr_from_scratch.spacetimes import flrw
from gr_from_scratch.spacetimes import minkowski as mk
from gr_from_scratch.spacetimes.schwarzschild import M, r


def zero_arr(arr):
    return all(simplify(e) == 0 for e in np.asarray(arr, dtype=object).ravel())


@pytest.mark.parametrize("make", [mk.minkowski, mk.euclidean_polar, mk.euclidean_spherical, mk.minkowski_spherical])
def test_flat_spaces_have_zero_riemann(make):
    g = compute_geometry(make())
    assert g.riemann.is_zero()
    assert g.scalar == 0


def test_unit_sphere_curvature(sphere_geom):
    th = mk.theta
    R = sphere_geom.riemann
    # R^theta_{phi theta phi} = sin^2(theta) for the unit sphere with the package convention
    assert simplify(R.c[0, 1, 0, 1] - sp.sin(th) ** 2) == 0
    assert sphere_geom.scalar == 2
    assert simplify(sphere_geom.ricci.c[0, 0] - 1) == 0
    assert sphere_geom.einstein.is_zero()  # G vanishes identically in 2D


def test_riemann_symmetries_schwarzschild(schw_geom):
    for name, res in bianchi.pair_and_antisymmetry_residuals(schw_geom.riemann).items():
        assert zero_arr(res), name


def test_ricci_symmetric(schw_geom):
    c = schw_geom.ricci.c
    assert zero_arr(c - c.T)


def test_algebraic_bianchi_on_generic_3d_metric():
    x, y, z = sp.symbols("x y z", real=True)
    m = Metric(
        [x, y, z],
        sp.Matrix([[1 + z**2, x * y / 5, 0], [x * y / 5, 2 + x**2, z / 7], [0, z / 7, 1 + y**2]]),
        simplifier=lambda e: e,
    )
    g = compute_geometry(m)
    res = bianchi.algebraic_bianchi_residual(g.riemann)
    val = {x: 0.3, y: -0.4, z: 0.7}
    assert max(abs(float(sp.N(e.subs(val)))) for e in res.ravel()) < 1e-12
    # negative control: the residual is not trivially zero for a deliberately wrong tensor
    wrong = np.array(g.riemann.lower_index(0).c, dtype=object).copy()
    wrong[0, 1, 0, 1] += 1
    cyc = wrong + np.transpose(wrong, (0, 2, 3, 1)) + np.transpose(wrong, (0, 3, 1, 2))
    assert any(simplify(e) != 0 for e in cyc.ravel())


def test_differential_bianchi_sphere_and_schwarzschild(sphere_geom, schw_geom):
    for g in (sphere_geom, schw_geom):
        res = bianchi.differential_bianchi_residual(g.riemann, g.conn)
        assert zero_arr(res)


def test_contracted_bianchi_identity_on_flrw():
    """nabla_mu G^{mu nu} = 0 identically, shown on FLRW where G != 0."""
    m = flrw.flrw()
    g = compute_geometry(m)
    assert not g.einstein.is_zero()
    res = bianchi.contracted_bianchi_residual(g.einstein.raise_index(0).raise_index(1), g.conn)
    assert zero_arr(res)


def test_ricci_flat_schwarzschild_riemann_nonzero(schw_geom):
    assert schw_geom.ricci.is_zero()
    assert schw_geom.scalar == 0
    assert schw_geom.einstein.is_zero()
    assert not schw_geom.riemann.is_zero()


def test_kretschmann_schwarzschild(schw_geom):
    k = invariants.kretschmann(schw_geom.riemann)
    assert simplify(k - 48 * M**2 / r**6) == 0


def test_weyl_equals_riemann_in_vacuum_and_decomposition_identity(schw_geom):
    C = weyl_tensor(schw_geom.riemann)
    assert zero_arr(C.c - schw_geom.riemann.lower_index(0).c)
    # tracelessness
    assert zero_arr(np.einsum("abcd,ac->bd", C.c, np.array(schw_geom.metric.ginv.tolist(), dtype=object)))


def test_weyl_decomposition_of_kretschmann_in_4d():
    """K = C^2 + 2 Ric^2 - R^2/3 on a non-vacuum 4D metric (FLRW with curvature)."""
    g = compute_geometry(flrw.flrw())
    K = invariants.kretschmann(g.riemann)
    C2 = invariants.weyl_squared(g.riemann)
    ric2 = invariants.ricci_squared(g.riemann)
    assert simplify(K - (C2 + 2 * ric2 - g.scalar**2 / 3)) == 0


def test_weyl_vanishes_for_flrw_conformally_flat():
    g = compute_geometry(flrw.flrw())
    assert weyl_tensor(g.riemann).is_zero()


def test_weyl_vanishes_identically_in_three_dimensions():
    x, y, z = sp.symbols("x y z", real=True)
    m = Metric([x, y, z], sp.diag(1 + x**2, 2 + y**2 * z**2, 3 + x * z), simplifier=lambda e: e)
    g = compute_geometry(m)
    C = weyl_tensor(g.riemann)
    val = {x: 0.3, y: 0.9, z: 0.5}
    assert max(abs(float(sp.N(e.subs(val)))) for e in C.c.ravel()) < 1e-12
    # non-trivial: the curvature itself is not zero at the sample point
    assert max(abs(float(sp.N(e.subs(val)))) for e in g.riemann.c.ravel()) > 1e-3


def test_weyl_undefined_in_two_dimensions(sphere_geom):
    with pytest.raises(DimensionError):
        weyl_tensor(sphere_geom.riemann)


def test_gauss_bonnet_density_equals_kretschmann_in_vacuum(schw_geom):
    assert simplify(invariants.euler_density(schw_geom.riemann) - invariants.kretschmann(schw_geom.riemann)) == 0


def test_pontryagin_vanishes_for_static_spherical(schw_geom):
    assert invariants.pontryagin_density(schw_geom.riemann) == 0


def test_einstein_tensor_is_ricci_minus_half_R_g():
    g = compute_geometry(flrw.flrw())
    m = g.metric
    expected = g.ricci.c - sp.Rational(1, 2) * g.scalar * np.array(m.g.tolist(), dtype=object)
    assert zero_arr(g.einstein.c - expected)
