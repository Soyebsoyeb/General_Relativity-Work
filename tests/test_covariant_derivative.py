import numpy as np
import sympy as sp

from gr_from_scratch.core.simplification import simplify
from gr_from_scratch.core.tensor import Tensor
from gr_from_scratch.geometry.connection import christoffel
from gr_from_scratch.geometry.covariant_derivative import covariant_derivative, divergence
from gr_from_scratch.spacetimes import minkowski as mk


def zeros(arr):
    return all(simplify(e) == 0 for e in np.asarray(arr, dtype=object).ravel())


def test_constant_cartesian_vector_has_zero_derivative_in_polar():
    """V = d_x has polar components (cos phi, -sin phi / r); nabla V must vanish."""
    m = mk.euclidean_polar()
    r, ph = mk.r, mk.phi
    V = Tensor(m, [sp.cos(ph), -sp.sin(ph) / r], "u")
    assert zeros(covariant_derivative(V, christoffel(m)).c)


def test_naive_partial_derivative_does_not_vanish():
    m = mk.euclidean_polar()
    r, ph = mk.r, mk.phi
    V = [sp.cos(ph), -sp.sin(ph) / r]
    partial = [[sp.diff(V[a], x) for x in m.coords] for a in range(2)]
    assert not zeros(np.array(partial, dtype=object))


def test_gradient_hessian_symmetric_and_sign_of_lower_slot():
    m = mk.euclidean_polar()
    r, ph = mk.r, mk.phi
    f = r**2 * sp.cos(ph) * sp.sin(ph)  # = x y
    df = Tensor(m, [sp.diff(f, r), sp.diff(f, ph)], "d")
    H = covariant_derivative(df, christoffel(m))
    assert zeros(H.c - np.transpose(H.c))
    # Cartesian Hessian of x y: [[0,1],[1,0]] -> check trace = Laplacian = 0
    assert simplify(sum(m.ginv[a, b] * H.c[a, b] for a in range(2) for b in range(2))) == 0


def test_metric_compatibility_via_generic_derivative(schw_geom):
    g = Tensor.metric_tensor(schw_geom.metric)
    assert zeros(covariant_derivative(g, schw_geom.conn).c)
    gi = Tensor(schw_geom.metric, np.array(schw_geom.metric.ginv.tolist(), dtype=object), "uu")
    assert zeros(covariant_derivative(gi, schw_geom.conn).c)


def test_leibniz_rule():
    m = mk.euclidean_spherical()
    conn = christoffel(m)
    r, th, ph = mk.r, mk.theta, mk.phi
    A = Tensor(m, [r * th, sp.sin(ph), r], "u")
    B = Tensor(m, [sp.cos(th), r**2, ph], "d")
    lhs = covariant_derivative(A.outer(B), conn)  # (A^a B_b ; c)
    DA, DB = covariant_derivative(A, conn), covariant_derivative(B, conn)
    n = 3
    rhs = np.empty((n, n, n), dtype=object)
    for a, b, c in np.ndindex(n, n, n):
        rhs[a, b, c] = DA.c[a, c] * B.c[b] + A.c[a] * DB.c[b, c]
    assert zeros(lhs.c - rhs)


def test_contraction_commutes_with_derivative():
    m = mk.euclidean_spherical()
    conn = christoffel(m)
    r, th, ph = mk.r, mk.theta, mk.phi
    T = Tensor(m, np.array([[r, th, ph], [th * r, 1, r**2], [sp.sin(ph), r, th]], dtype=object), "ud")
    a = covariant_derivative(T.contract(0, 1), conn).c  # d (trace)
    b = covariant_derivative(T, conn).contract(0, 1).c
    assert zeros(a - b)


def test_divergence_of_vector_matches_density_formula(schw_geom):
    m = schw_geom.metric
    r, th = m.coords[1], m.coords[2]
    V = Tensor(m, [0, r**3 * sp.sin(th), th, 0], "u")
    div = divergence(V, schw_geom.conn, 0).c.item()
    sg = r**2 * sp.sin(th)  # sqrt(-g) for Schwarzschild
    dens = sum(sp.diff(sg * V.c[i], m.coords[i]) for i in range(4)) / sg
    assert simplify(div - dens) == 0
