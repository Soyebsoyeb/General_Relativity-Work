# 1. Levi-Civita connection, covariant derivative and curvature

## Definitions and assumptions

A smooth manifold with a metric $g_{\mu\nu}$ (non-degenerate, fixed signature). A connection $\nabla$ is *metric compatible* if $\nabla_\lambda g_{\mu\nu}=0$ and *torsion free* if $\Gamma^\rho{}_{\mu\nu}=\Gamma^\rho{}_{\nu\mu}$ in a coordinate basis. Index conventions are in [conventions](../conventions.md).

## Derivation of the Christoffel formula

Metric compatibility reads
$$\partial_\lambda g_{\mu\nu}=\Gamma^\sigma{}_{\lambda\mu}g_{\sigma\nu}+\Gamma^\sigma{}_{\lambda\nu}g_{\mu\sigma}.$$
Write $\Gamma_{\sigma\lambda\mu}=g_{\sigma\rho}\Gamma^\rho{}_{\lambda\mu}$, so $\partial_\lambda g_{\mu\nu}=\Gamma_{\nu\lambda\mu}+\Gamma_{\mu\lambda\nu}$. Permuting $(\lambda,\mu,\nu)$ cyclically,
$$\partial_\mu g_{\nu\lambda}=\Gamma_{\lambda\mu\nu}+\Gamma_{\nu\mu\lambda},\qquad
\partial_\nu g_{\lambda\mu}=\Gamma_{\mu\nu\lambda}+\Gamma_{\lambda\nu\mu}.$$
Take the first plus the second minus the third. With $\Gamma_{\sigma\mu\nu}=\Gamma_{\sigma\nu\mu}$ everything cancels except $2\Gamma_{\nu\lambda\mu}$:
$$\Gamma^\rho{}_{\mu\nu}=\tfrac12 g^{\rho\sigma}\bigl(\partial_\mu g_{\nu\sigma}+\partial_\nu g_{\mu\sigma}-\partial_\sigma g_{\mu\nu}\bigr).$$
The two assumptions therefore determine $\Gamma$ uniquely.

## Why $\Gamma$ is not a tensor

Under $x\to x'$,
$$\Gamma'^\alpha{}_{\beta\gamma}=\frac{\partial x'^\alpha}{\partial x^\rho}\frac{\partial x^\mu}{\partial x'^\beta}\frac{\partial x^\nu}{\partial x'^\gamma}\Gamma^\rho{}_{\mu\nu}+\frac{\partial x'^\alpha}{\partial x^\rho}\frac{\partial^2 x^\rho}{\partial x'^\beta\partial x'^\gamma}.$$
The inhomogeneous second term is what lets the connection vanish in Cartesian coordinates but not in polar ones. The partial derivative $\partial_\mu V^\nu$ also transforms inhomogeneously, and the two inhomogeneities cancel in
$$\nabla_\mu V^\nu=\partial_\mu V^\nu+\Gamma^\nu{}_{\mu\lambda}V^\lambda,\qquad \nabla_\mu\omega_\nu=\partial_\mu\omega_\nu-\Gamma^\lambda{}_{\mu\nu}\omega_\lambda,$$
so $\nabla V$ is a tensor. Geometrically, the components $V^\nu$ change along a curve partly because the coordinate basis vectors change; $\Gamma$ subtracts that part. A constant Cartesian vector has non-zero partial derivatives in polar coordinates but zero covariant derivative.

## Curvature

Commuting covariant derivatives on a vector,
$$[\nabla_\mu,\nabla_\nu]V^\rho=R^\rho{}_{\sigma\mu\nu}V^\sigma\quad(\text{zero torsion}),$$
which expands to the declared formula for $R^\rho{}_{\sigma\mu\nu}$. Consequences used as tests:

* $R_{\rho\sigma\mu\nu}=-R_{\sigma\rho\mu\nu}=-R_{\rho\sigma\nu\mu}=R_{\mu\nu\rho\sigma}$,
* algebraic Bianchi $R_{\rho[\sigma\mu\nu]}=0$,
* differential Bianchi $\nabla_{[\lambda}R_{\rho\sigma]\mu\nu}=0$, and its contraction $\nabla_\mu G^{\mu\nu}=0$,
* $R_{\mu\nu}=R_{\nu\mu}$.

### Weyl tensor and dimension dependence
$$C_{abcd}=R_{abcd}-\tfrac{2}{n-2}\bigl(g_{a[c}R_{d]b}-g_{b[c}R_{d]a}\bigr)+\tfrac{2}{(n-1)(n-2)}R\,g_{a[c}g_{d]b}.$$
It is trace-free. For $n=3$ it vanishes identically (Riemann is fixed by Ricci). In $n=4$, $K=C^2+2R_{ab}R^{ab}-\tfrac13R^2$. FLRW is conformally flat, so $C=0$. The Gauss–Bonnet density $R^2-4R_{ab}R^{ab}+K$ is topological in $n=4$.

## Implementation
`geometry/connection.py`, `covariant_derivative.py`, `riemann.py`, `ricci.py`, `weyl.py`, `einstein.py`, `bianchi.py`, `invariants.py`. The connection is a plain array, deliberately not a `Tensor`. `covariant_derivative` handles any rank and any `u`/`d` mixture.

## Validation
`tests/test_connection.py`, `test_covariant_derivative.py`, `test_curvature.py`: Cartesian/polar/spherical Minkowski; metric compatibility and zero torsion on Schwarzschild; the inhomogeneous law reproducing polar symbols from Cartesian ones; Leibniz rule; contraction commuting with $\nabla$; all symmetries; both Bianchi identities; Weyl properties; and negative controls (a wrong connection fails compatibility; a wrongly built Riemann fails the cyclic identity).

## Limitations
Tensors store explicit components, with no abstract-index canonicalisation. Simplification is by SymPy; `check_zero` reports when only numeric evidence was obtained.
