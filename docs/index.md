# general-relativity-from-scratch

A computational framework in which the central objects of general relativity are *derived, implemented, tested and applied*: metric, Christoffel symbols, covariant derivative, Riemann/Ricci/Weyl/Einstein tensors, Killing vectors, geodesics, geodesic deviation, linearised gravity, FLRW cosmology, tetrads, and a scalar field on a black-hole background.

The central calculations are the project's own (Python, SymPy, NumPy, SciPy). No third-party GR package is used or required. SciPy supplies ODE solving and quadrature; SymPy supplies symbolic algebra and simplification.

Start with:

* [Conventions](conventions.md)
* [Implementation status](status.md): what is complete, partial or only planned
* [Validation summary](validation.md): measured numbers from the experiments
* Derivations 1 to 11 in the navigation

Install and run:

```bash
pip install -e ".[dev]"
pytest                      # fast suite
pytest -m slow              # expensive symbolic checks
python experiments/run_all.py
```
