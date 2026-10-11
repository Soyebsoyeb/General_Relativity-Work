# Conventions

These are fixed once and used everywhere. `src/gr_from_scratch/core/conventions.py` repeats them and the tests enforce them (e.g. the unit sphere has $R=+2$).

| Item | Choice |
|---|---|
| Signature | $(-,+,+,+)$ |
| Christoffel symbols | $\Gamma^\rho{}_{\mu\nu}=\tfrac12 g^{\rho\sigma}(\partial_\mu g_{\nu\sigma}+\partial_\nu g_{\mu\sigma}-\partial_\sigma g_{\mu\nu})$ |
| Riemann | $R^\rho{}_{\sigma\mu\nu}=\partial_\mu\Gamma^\rho{}_{\nu\sigma}-\partial_\nu\Gamma^\rho{}_{\mu\sigma}+\Gamma^\rho{}_{\mu\lambda}\Gamma^\lambda{}_{\nu\sigma}-\Gamma^\rho{}_{\nu\lambda}\Gamma^\lambda{}_{\mu\sigma}$ |
| Ricci | $R_{\sigma\nu}=R^\rho{}_{\sigma\rho\nu}$ |
| Scalar curvature | $R=g^{\sigma\nu}R_{\sigma\nu}$ |
| Einstein tensor | $G_{\mu\nu}=R_{\mu\nu}-\tfrac12 R g_{\mu\nu}$ |
| Field equations | $G_{\mu\nu}+\Lambda g_{\mu\nu}=8\pi G c^{-4}T_{\mu\nu}$ |
| Covariant derivative | derivative index appended **last**: $(\nabla T)_{ab;c}$ is `array[a, b, c]` |
| Variance string | `"u"` upper, `"d"` lower; Riemann is `"uddd"` |
| Units | $G=c=1$ unless stated; lengths in units of $M$ |
| Electromagnetism | Gaussian: $T_{\mu\nu}=\frac1{4\pi}(F_{\mu\alpha}F_\nu{}^\alpha-\frac14 g_{\mu\nu}F^2)$ |
| Orientation | $\varepsilon^{0123}=-1/\sqrt{-g}$ (only used in the Pontryagin density) |

With these conventions a sphere has positive curvature, $[\nabla_\mu,\nabla_\nu]V^\rho=R^\rho{}_{\sigma\mu\nu}V^\sigma$, and the Einstein–Hilbert action is $\frac{1}{16\pi}\int(R-2\Lambda)\sqrt{-g}\,d^4x$.

## Exact versus approximate

* *Exact symbolic statement* — verified by simplification to the literal $0$, or by exact equality of expressions.
* *Numerical evidence* — verified at sample points (`core.simplification.check_zero` reports which of `exact`, `numeric`, `nonzero` applied). Kerr vacuum is of this kind.
* *Truncated series* — carry an explicit order $\epsilon^n$ (`perturbation_series`).
* *Numerical integration* — carries a tolerance; every experiment records it.
