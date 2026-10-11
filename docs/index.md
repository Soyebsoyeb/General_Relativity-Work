# general-relativity-from-scratch: documentation

This folder explains the mathematics behind the code and records how each result is checked. It is meant to be read in order, but every page stands on its own.

## What the project does

General relativity describes gravity as the curvature of a four-dimensional spacetime. This repository computes that curvature symbolically, follows particles and light through curved spacetime numerically, and tests every step against an independent calculation.

The central quantities (Christoffel symbols, Riemann and Einstein tensors, geodesics, linearised waves, cosmological expansion) are computed by code in this repository. No existing general-relativity package is used. SymPy provides symbolic algebra, and NumPy and SciPy provide numerical integration.

## Where to start

| If you want to know | Read |
|---|---|
| Which sign and index conventions are used | [Conventions](conventions.md) |
| What is finished, what is partial, what is missing | [Implementation status](status.md) |
| What numbers the experiments produced | [Validation summary](validation.md) |
| How to rerun everything | [Reproducibility](reproducibility.md) |
| Which module does what | [Module overview](api.md) |

## The derivations

Each derivation page follows the same plan: the setting and assumptions, the derivation in numbered steps, the physical meaning, how the code implements it, how the result is checked, and the known limits.

1. [Connection, covariant derivative and curvature](derivations/01_connection_and_curvature.md)
2. [The Einstein-Hilbert action and the field equations](derivations/02_einstein_hilbert.md)
3. [The Schwarzschild solution, orbits and observables](derivations/03_schwarzschild.md)
4. [Charged and rotating black holes](derivations/04_charged_and_rotating.md)
5. [The Newtonian limit](derivations/05_newtonian_limit.md)
6. [Linearised gravity and gravitational waves](derivations/06_linearized_gravity_and_gravitational_waves.md)
7. [FLRW cosmology](derivations/07_flrw_cosmology.md)
8. [Tetrads and Cartan's structure equations](derivations/08_tetrads_and_cartan.md)
9. [A scalar field on a curved background](derivations/09_quantum_fields_on_curved_spacetime.md)
10. [Perturbative gravity and effective field theory](derivations/10_perturbative_gravity_and_eft.md)
11. [A note on string theory](derivations/11_string_theory_bridge.md)

## A note on honesty

Three kinds of statement appear in these pages, and they are kept apart:

- **Exact:** shown by symbolic simplification to the literal value zero, or by exact equality of expressions.
- **Numerical evidence:** checked at sample points or measured from a numerical solution, with the tolerance stated.
- **Truncated series:** valid only to the stated order in a small parameter.

The project does not solve quantum gravity, and nowhere claims to.
