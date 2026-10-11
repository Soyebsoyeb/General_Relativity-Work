# Module overview

The package lives in `src/gr_from_scratch/`. This page is a map, not a reference manual; each module has docstrings that state its conventions.

| Package | Contents |
|---|---|
| `core` | `Metric`, `Tensor`, index operations, coordinate transformations of tensors and of the connection, and `check_zero` |
| `geometry` | connection, covariant derivative, Riemann, Ricci, Weyl, Einstein, Bianchi identities, invariants, Killing vectors, Lie derivative, differential forms, tetrads |
| `dynamics` | geodesics, variational derivation, conservation laws, effective potentials, geodesic deviation, relativistic observables, gravitational-wave test particles |
| `field_equations` | stress-energy tensors, vacuum equations, action variation, linearised gravity, perturbation series |
| `spacetimes` | Minkowski, Schwarzschild, Reissner-Nordstrom, Kerr, FLRW |
| `numerical` | fixed-step Runge-Kutta, convergence tools, constraint monitors, boundary data, Friedmann dynamics |
| `quantum` | scalar field on Schwarzschild, quadratic gravitational action, effective-field-theory example |
| `visualization` | plotting style and helpers |

## A short example

```python
from gr_from_scratch.spacetimes.schwarzschild import schwarzschild
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry
from gr_from_scratch.geometry.invariants import kretschmann

g = compute_geometry(schwarzschild())
g.einstein.is_zero()       # True: the vacuum equations hold
g.riemann.is_zero()        # False: vacuum, but not flat
kretschmann(g.riemann)     # 48*M**2/r**6
```

## Design choices worth knowing

- The connection is stored as a plain array, deliberately not as a `Tensor`, because its components do not transform like a tensor.
- Every derived tensor passes through a simplifier. For metrics where simplification is too slow (Kerr), expressions are left unsimplified and evaluated at points.
- Tests never use the same formula to produce both the answer and the expectation.
