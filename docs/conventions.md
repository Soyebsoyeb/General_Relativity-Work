# Conventions

Every formula in the code and in these pages uses the choices below. They are fixed once, repeated in `src/gr_from_scratch/core/conventions.py`, and enforced by tests (for example, a unit sphere must come out with positive curvature $R=2$).

## Metric and units

The metric has signature $(-,+,+,+)$: time has a negative square, space positive. Unless a page says otherwise, units are chosen so that $G=c=1$, and lengths are measured in units of the mass $M$.

## Connection and curvature

The Christoffel symbols are

$$
\Gamma^{\rho}{}_{\mu\nu}=\frac{1}{2}\,g^{\rho\sigma}\left(\partial_\mu g_{\nu\sigma}+\partial_\nu g_{\mu\sigma}-\partial_\sigma g_{\mu\nu}\right).
$$

The Riemann tensor is

$$
R^{\rho}{}_{\sigma\mu\nu}=\partial_\mu\Gamma^{\rho}{}_{\nu\sigma}-\partial_\nu\Gamma^{\rho}{}_{\mu\sigma}+\Gamma^{\rho}{}_{\mu\lambda}\Gamma^{\lambda}{}_{\nu\sigma}-\Gamma^{\rho}{}_{\nu\lambda}\Gamma^{\lambda}{}_{\mu\sigma}.
$$

Contracting gives the Ricci tensor and scalar curvature,

$$
R_{\sigma\nu}=R^{\rho}{}_{\sigma\rho\nu},\qquad R=g^{\sigma\nu}R_{\sigma\nu},
$$

and the Einstein tensor is

$$
G_{\mu\nu}=R_{\mu\nu}-\frac{1}{2}\,R\,g_{\mu\nu}.
$$

With these signs, commuting two covariant derivatives of a vector gives

$$
\left[\nabla_\mu,\nabla_\nu\right]V^{\rho}=R^{\rho}{}_{\sigma\mu\nu}\,V^{\sigma},
$$

and a sphere has positive curvature.

## Field equations

$$
G_{\mu\nu}+\Lambda\,g_{\mu\nu}=\frac{8\pi G}{c^{4}}\,T_{\mu\nu}.
$$

For electromagnetism the code uses Gaussian units, so that

$$
T_{\mu\nu}=\frac{1}{4\pi}\left(F_{\mu\alpha}F_{\nu}{}^{\alpha}-\frac{1}{4}\,g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\right).
$$

## How indices are stored in the code

| Item | Choice |
|---|---|
| Upper index | the letter `u` in a variance string |
| Lower index | the letter `d` in a variance string |
| Riemann tensor | stored as $R^{\rho}{}_{\sigma\mu\nu}$, variance string `uddd` |
| Covariant derivative | the new derivative index is appended last, so the array entry `[a, b, c]` of $\nabla T$ means $\nabla_c T_{ab}$ |
| Orientation | $\varepsilon^{0123}=-1/\sqrt{-g}$; used only in the Pontryagin density |

## Three kinds of result

| Kind | Meaning |
|---|---|
| Exact | SymPy reduced the expression to the literal 0, or two expressions are exactly equal |
| Numerical evidence | the expression was evaluated at sample points, or a numerical solution was compared with a reference; the tolerance is stated |
| Truncated series | an expansion in a small parameter, correct only to the order stated |

The helper `check_zero` in `core/simplification.py` reports which of the first two applied. Kerr's vacuum property is of the second kind.
