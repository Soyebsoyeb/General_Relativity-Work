# Module overview

| Package | Contents |
|---|---|
| `core` | `Metric`, `Tensor`, `index_operations`, `coordinates` (tensor and connection transformation laws), `simplification` (`check_zero`), `conventions` |
| `geometry` | `connection`, `covariant_derivative`, `riemann`, `ricci`, `weyl`, `einstein`, `bianchi`, `invariants`, `killing`, `lie_derivatives`, `differential_forms`, `tetrads` |
| `dynamics` | `geodesics`, `variational_principles`, `conserved_quantities`, `effective_potential`, `geodesic_deviation`, `observables`, `gravitational_waves` |
| `field_equations` | `stress_energy`, `vacuum_equations`, `einstein_hilbert`, `linearized_gravity`, `perturbation_series` |
| `spacetimes` | `minkowski`, `schwarzschild`, `reissner_nordstrom`, `kerr`, `flrw` |
| `numerical` | `integration` (fixed-step RK4), `convergence`, `constraints`, `boundary_conditions`, `friedmann` |
| `quantum` | `klein_gordon`, `perturbative_gravity`, `eft` |
| `visualization` | `style`, `geodesics`, `curvature`, `effective_potentials`, `spacetime_diagrams` |

Minimal example:

```python
from gr_from_scratch.spacetimes.schwarzschild import schwarzschild
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry
from gr_from_scratch.geometry.invariants import kretschmann

g = compute_geometry(schwarzschild())
g.einstein.is_zero()  # True: vacuum
kretschmann(g.riemann)  # 48*M**2/r**6
```
