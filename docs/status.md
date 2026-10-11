# Implementation status

This page lists what is finished, what is partly done, and what does not exist. It is meant to be an accurate account, not a sales list.

## Complete

| Area | Notes |
|---|---|
| Metric, tensors, index operations, coordinate transformations | components are stored explicitly |
| Christoffel symbols and covariant derivative of any rank | any mixture of upper and lower indices |
| Riemann, Ricci, Weyl and Einstein tensors | Weyl tensor for dimension 3 or higher |
| Bianchi identities, curvature symmetries | tested on explicit metrics |
| Curvature invariants | Kretschmann, Ricci squared, Weyl squared, Euler density; Pontryagin density in four dimensions |
| Lie derivative | coordinate formula |
| Geodesics | symbolic (from the connection and from the Euler-Lagrange equations) and numerical |
| Geodesic deviation and the tidal tensor | |
| Schwarzschild | derivation from the ansatz, horizon-regular coordinates, orbits, observables |
| Reissner-Nordstrom | exact Einstein-Maxwell check |
| Newtonian limit | analytic and numerical |
| Linearised gravity, gauge freedom, plane waves | flat background only |
| FLRW cosmology | compared with five analytic solutions |
| Differential forms, tetrads, Cartan equations | frames are supplied by the user |

## Partial

| Area | What exists | What is missing |
|---|---|---|
| Killing vectors | verifying a supplied vector, and conserved quantities | a solver for the Killing equation |
| Einstein-Hilbert variation | identities checked on supplied backgrounds; a reduced action for static spherical metrics | a general variational engine; boundary terms |
| Kerr | metric, horizons, ergosurface, circular orbits, ISCO; vacuum property checked at sample points | a symbolic vacuum proof; off-equatorial orbits |
| Second-order perturbations | explicit truncated series for specific families | a general framework |
| Scalar field on Schwarzschild | classical scattering | quantisation |
| Quadratic gravitational action | verified for perturbations depending on two coordinates | all four coordinates (a cost limit, not a principle) |
| Effective field theory | one symbolic example and one estimate | anything beyond |

## Explanatory text only

String theory (page 11). Nothing is computed.

## Planned

A ray-tracing and imaging module. Null geodesic bundles appear in experiment 04, but there is no image-plane code.

## Not implemented, and not claimed

Graviton quantisation, scattering amplitudes, loop calculations, and numerical relativity (3+1 evolution).

**This repository does not solve quantum gravity.**
