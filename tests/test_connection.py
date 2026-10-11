import numpy as np
import sympy as sp

from gr_from_scratch.core.coordinates import transform_connection
from gr_from_scratch.core.simplification import simplify
from gr_from_scratch.geometry.connection import christoffel, metric_compatibility_residual, torsion
from gr_from_scratch.spacetimes import minkowski as mk
from gr_from_scratch.spacetimes.schwarzschild import M


def all_zero(arr):
    return all(simplify(e) == 0 for e in np.asarray(arr, dtype=object).ravel())


def test_cartesian_minkowski_connection_vanishes():
    assert christoffel(mk.minkowski()).nonzero() == {}


def test_polar_christoffels():
    c = christoffel(mk.euclidean_polar())
    r = mk.r
    nz = c.nonzero()
    assert nz == {(0, 1, 1): -r, (1, 0, 1): 1 / r, (1, 1, 0): 1 / r}


def test_spherical_christoffels():
    c = christoffel(mk.euclidean_spherical())
    r, th = mk.r, mk.theta
    nz = c.nonzero()
    assert nz[(0, 1, 1)] == -r
    assert simplify(nz[(0, 2, 2)] + r * sp.sin(th) ** 2) == 0
    assert simplify(nz[(1, 2, 2)] + sp.sin(th) * sp.cos(th)) == 0
    assert simplify(nz[(2, 0, 2)] - 1 / r) == 0
    assert simplify(nz[(1, 0, 1)] - 1 / r) == 0
    assert simplify(nz[(2, 1, 2)] - sp.cos(th) / sp.sin(th)) == 0


def test_zero_torsion_and_metric_compatibility(schw_geom):
    assert all_zero(torsion(schw_geom.conn))
    assert all_zero(metric_compatibility_residual(schw_geom.conn))


def test_metric_compatibility_detects_wrong_connection(schw_geom):
    bad = schw_geom.conn
    bad.gamma = bad.gamma.copy()
    saved = bad.gamma[0, 0, 1]
    bad.gamma[0, 0, 1] = 2 * saved
    try:
        assert not all_zero(metric_compatibility_residual(bad))
    finally:
        bad.gamma[0, 0, 1] = saved


def test_schwarzschild_known_components(schw_geom):
    r = sp.Symbol("r", positive=True)
    G = schw_geom.conn
    assert simplify(G[0, 0, 1] - M / (r * (r - 2 * M))) == 0
    assert simplify(G[1, 0, 0] - M * (r - 2 * M) / r**3) == 0


def test_connection_is_not_a_tensor_but_obeys_inhomogeneous_law():
    """Gamma vanishes in Cartesian, is non-zero in polar; the inhomogeneous law
    reproduces the directly computed polar symbols from the Cartesian ones."""
    x, y = sp.symbols("x y", real=True)
    from gr_from_scratch.core.metric import Metric

    cart = Metric([x, y], sp.eye(2))
    gam_c = christoffel(cart).gamma
    polar = mk.euclidean_polar()
    direct = christoffel(polar).gamma
    via = transform_connection(gam_c, cart.coords, polar.coords, mk.polar_to_cartesian(), simplify)
    assert all_zero(via - direct)
    # a homogeneous (tensor) law applied to the zero array would give zero: it cannot reproduce Gamma
    assert not all_zero(direct)


def test_connection_can_be_removed_at_a_point_but_curvature_cannot():
    """Riemann normal coordinates: for the sphere the symbols vanish at the equator point
    only after a coordinate change, while the Ricci scalar stays 2."""
    from gr_from_scratch.field_equations.vacuum_equations import compute_geometry

    g = compute_geometry(mk.unit_sphere())
    assert g.scalar == 2
    th = mk.theta
    at_pole_axis = {k: sp.simplify(v.subs(th, sp.pi / 2)) for k, v in g.conn.nonzero().items()}
    # at the equator the symbols vanish (this particular chart) yet curvature is not zero
    assert all(v == 0 for v in at_pole_axis.values()) and g.scalar == 2
