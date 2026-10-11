import numpy as np
import pytest
import sympy as sp

from gr_from_scratch.core import index_operations as ops
from gr_from_scratch.core.coordinates import transform_tensor
from gr_from_scratch.core.metric import Metric, MetricError
from gr_from_scratch.core.simplification import check_zero
from gr_from_scratch.core.tensor import Tensor, TensorError
from gr_from_scratch.spacetimes import minkowski as mk


def test_metric_rejects_asymmetric_and_bad_shape():
    with pytest.raises(MetricError):
        Metric("x y", sp.Matrix([[1, 1], [0, 1]]))
    with pytest.raises(MetricError):
        Metric("x y", sp.eye(3))


def test_supplied_inverse_is_checked():
    x, y = sp.symbols("x y", real=True)
    with pytest.raises(MetricError):
        Metric([x, y], sp.diag(1, 2), inverse=sp.diag(1, 3))


def test_det_inverse_and_signature():
    m = mk.minkowski()
    assert m.det == -1
    assert m.ginv == m.g
    assert m.signature_at({}) == (1, 3)
    assert m.is_lorentzian_at({})
    assert not mk.euclidean_polar().is_lorentzian_at({mk.r: 1})


def test_raise_lower_roundtrip_and_metric_identity():
    m = mk.euclidean_spherical()
    g = Tensor.metric_tensor(m)
    gup = g.raise_index(0).raise_index(1)
    assert gup.variance == "uu"
    assert sp.Matrix(gup.c.tolist()) == m.ginv
    delta = g.raise_index(0)
    assert sp.Matrix(delta.c.tolist()) == sp.eye(3)
    assert np.all(gup.lower_index(0).lower_index(1).c == g.c)


def test_variance_errors():
    m = mk.minkowski()
    t = Tensor(m, np.zeros((4, 4), dtype=int), "dd")
    with pytest.raises(TensorError):
        t.contract(0, 1)
    with pytest.raises(TensorError):
        t.lower_index(0)
    with pytest.raises(TensorError):
        Tensor(m, np.zeros((4,), dtype=int), "dd")


def test_contraction_trace_and_outer():
    m = mk.euclidean_polar()
    t = Tensor(m, [[1, 2], [3, 4]], "ud")
    assert t.contract(0, 1).c.item() == 5
    v = Tensor(m, [1, 2], "u")
    w = Tensor(m, [3, 4], "d")
    assert v.outer(w).contract(0, 1).c.item() == 11


def test_symmetrization_projectors():
    rng = np.random.default_rng(0)
    arr = np.array(rng.integers(-5, 5, (3, 3, 3)), dtype=object)
    s = ops.symmetrize(arr, (0, 1, 2))
    a = ops.antisymmetrize(arr, (0, 1, 2))
    assert np.all(s == np.transpose(s, (1, 0, 2)))
    assert np.all(a == -np.transpose(a, (1, 0, 2)))
    assert np.all(ops.antisymmetrize(a, (0, 1, 2)) == a)  # projector
    assert np.all(ops.symmetrize(a, (0, 1)) == 0 * a)


def test_tensor_transformation_law_polar_to_cartesian():
    """The metric and a vector field transform covariantly under x = r cos(phi), y = r sin(phi)."""
    x, y = sp.symbols("x y", real=True)
    cart = Metric([x, y], sp.eye(2))
    polar = mk.euclidean_polar()
    old_in_new = mk.polar_to_cartesian()  # (x, y) as functions of (r, phi)
    g_c = Tensor.metric_tensor(cart)
    g_p = transform_tensor(g_c, polar, old_in_new)
    assert sp.Matrix(g_p.c.tolist()).applyfunc(sp.simplify) == polar.g
    # rotation generator -y d_x + x d_y becomes d_phi
    v = Tensor(cart, [-y, x], "u")
    vp = transform_tensor(v, polar, old_in_new)
    assert [sp.simplify(c) for c in vp.c] == [0, 1]


def test_transform_matches_metric_transform():
    cart = Metric(sp.symbols("x y", real=True), sp.eye(2))
    polar_direct = mk.euclidean_polar()
    polar_via = cart.transform([mk.r, mk.phi], mk.polar_to_cartesian())
    assert polar_via.g == polar_direct.g


def test_check_zero_reports_method():
    x = sp.Symbol("x", real=True)
    assert check_zero(sp.sin(x) ** 2 + sp.cos(x) ** 2 - 1).status == "exact"
    assert check_zero(x - 1).status == "nonzero"
    ok = check_zero(sp.Float("1e-30") * x)
    assert ok.ok
