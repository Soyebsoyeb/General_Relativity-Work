# Implementation status

| Area | Status | Notes |
|---|---|---|
| Metric, tensors, index operations, coordinate transformations | complete | explicit components; no abstract-index algebra |
| Christoffel, covariant derivative (any rank), Riemann, Ricci, Weyl, Einstein | complete | Weyl for $n\ge3$ |
| Bianchi identities, curvature symmetries, Kretschmann, Euler density, Pontryagin density | complete | Pontryagin 4D only |
| Killing vectors | partial | **verification** of supplied fields and conserved quantities; no solver for the Killing equation |
| Lie derivative | complete | coordinate formula |
| Geodesic equation: symbolic (Euler–Lagrange) and numerical | complete | adaptive DOP853, events, constraint and conservation monitors, fixed-step RK4 for order tests |
| Geodesic deviation, tidal tensor | complete | static frame in the experiment |
| Einstein–Hilbert variation | partial | identities checked on supplied backgrounds; reduced static-spherical action; no boundary (GHY) term; no general variational engine |
| Schwarzschild (derivation from ansatz, EF coordinates, invariants, orbits, observables) | complete | Birkhoff not proved |
| Reissner–Nordström | complete | exact Einstein–Maxwell verification |
| Kerr | partial | vacuum character verified **numerically** at sample points; equatorial circular orbits and ISCO validated; no symbolic proof, no off-equatorial orbits |
| Newtonian limit | complete | analytic and numeric |
| Light deflection, periapsis advance, redshift, Shapiro delay | complete | dual computation (quadrature + geodesics) except Shapiro (quadrature vs weak-field formula) |
| Linearised gravity, gauge, TT waves, deviation | complete | flat background only |
| Second-order perturbation series | partial | explicit $\epsilon^2$ truncation; specific families |
| FLRW | complete | analytic comparisons for five models |
| Differential forms, tetrads, Cartan equations | complete | frames supplied by user; no spinors, no first-order action |
| Scalar field on Schwarzschild | partial | classical scattering; no quantisation |
| Quadratic gravitational action, Fierz–Pauli | partial | $h$ restricted to two coordinates in tests |
| EFT corrections | partial | one symbolic example + dimensional estimate |
| String theory | explanatory text only | not implemented |
| **Ray tracing / imaging (`numerical/ray_tracing.py`)** | **planned** | null geodesic bundles are in experiment 04; no image-plane code |
| Graviton quantisation, amplitudes, loop calculations | **not implemented** | not claimed |
| Numerical relativity (3+1 evolution) | **not implemented** | not claimed |
| Notebooks | complete | 12 notebooks, executed |

**This repository does not solve quantum gravity.**
