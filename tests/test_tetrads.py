import numpy as np
import sympy as sp

from gr_from_scratch.core.simplification import simplify
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry
from gr_from_scratch.geometry.differential_forms import Form
from gr_from_scratch.geometry.tetrads import Coframe
from gr_from_scratch.spacetimes import flrw as fl
from gr_from_scratch.spacetimes import minkowski as mk
from gr_from_scratch.spacetimes.schwarzschild import f, r, schwarzschild, theta


def zero(arr):
    return all(simplify(e) == 0 for e in np.asarray(arr, dtype=object).ravel())


def test_exterior_derivative_squares_to_zero():
    x, y, z = sp.symbols("x y z", real=True)
    w = Form.one_form([x, y, z], [x * y * z, sp.sin(x) * z**2, y + x**3])
    assert w.d().d().is_zero()
    h = Form(w.coords, 0, {(): x * sp.exp(y) * z})
    assert h.d().d().is_zero()


def test_wedge_antisymmetry_and_leibniz():
    x, y, z = sp.symbols("x y z", real=True)
    a = Form.one_form([x, y, z], [y, z, x])
    b = Form.one_form([x, y, z], [z**2, x, y * z])
    ab, ba = a.wedge(b), b.wedge(a)
    assert (ab + ba).is_zero()
    assert a.wedge(a).is_zero()
    # d(a ^ b) = da ^ b - a ^ db
    lhs = a.wedge(b).d()
    rhs = a.d().wedge(b) + (a.wedge(b.d())) * (-1)
    assert (lhs + rhs * (-1)).is_zero()


def test_schwarzschild_frame_vs_coordinate_curvature():
    s = schwarzschild()
    ff = f()
    e = sp.diag(sp.sqrt(ff), 1 / sp.sqrt(ff), r, r * sp.sin(theta))
    cf = Coframe(s, e)
    assert zero(cf.first_structure_residual())
    geom = compute_geometry(s)
    Rf = cf.riemann_frame()
    Rc = cf.project_coordinate_riemann(geom.riemann)
    assert zero(Rf - Rc)
    # nonzero: the frame Riemann is not trivially zero
    assert not zero(Rf)


def test_spin_connection_antisymmetric_in_frame_indices():
    s = schwarzschild()
    ff = f()
    e = sp.diag(sp.sqrt(ff), 1 / sp.sqrt(ff), r, r * sp.sin(theta))
    cf = Coframe(s, e)
    w = cf.spin_connection_frame()  # omega^a_{bc}
    eta = cf.eta
    low = np.empty(w.shape, dtype=object)
    for a, b, c in np.ndindex(w.shape):
        low[a, b, c] = sum(eta[a, d] * w[d, b, c] for d in range(4))
    assert zero(low + np.transpose(low, (1, 0, 2)))


def test_flrw_frame_vs_coordinate_curvature_with_spatial_curvature():
    k = fl.k
    a = sp.Function("a", positive=True)(fl.t)
    s = fl.flrw()
    e = sp.diag(1, a / sp.sqrt(1 - k * fl.r**2), a * fl.r, a * fl.r * sp.sin(fl.theta))
    cf = Coframe(s, e)
    geom = compute_geometry(s)
    assert zero(cf.riemann_frame() - cf.project_coordinate_riemann(geom.riemann))


def test_sphere_frame_curvature_is_one():
    s = mk.unit_sphere()
    e = sp.diag(1, sp.sin(mk.theta))
    cf = Coframe(s, e, eta=sp.eye(2))
    R = cf.riemann_frame()
    assert simplify(R[0, 1, 0, 1] - 1) == 0  # K = +1
