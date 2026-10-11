import pytest
import sympy as sp

from gr_from_scratch.field_equations.vacuum_equations import compute_geometry
from gr_from_scratch.spacetimes.schwarzschild import schwarzschild


@pytest.fixture(scope="session")
def schw_geom():
    return compute_geometry(schwarzschild())


@pytest.fixture(scope="session")
def sphere_geom():
    from gr_from_scratch.spacetimes.minkowski import unit_sphere

    return compute_geometry(unit_sphere())


def zero(x):
    return sp.simplify(x) == 0
