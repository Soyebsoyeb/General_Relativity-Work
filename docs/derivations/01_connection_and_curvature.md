# 1. Connection, covariant derivative and curvature

**In one paragraph.** On a curved space you cannot simply subtract vectors at different points, because their components depend on the coordinates. The connection tells you how to differentiate in a way that does not depend on that choice. Once a metric is given, requiring the connection to respect it fixes the connection uniquely. The failure of two successive derivatives to commute is the curvature.

## 1. Setting and assumptions

We have a smooth manifold with a metric $g_{\mu\nu}$ that is symmetric and invertible. We look for a derivative operator $\nabla$ with two properties.

- **Metric compatibility:** $\nabla_\lambda g_{\mu\nu}=0$. Lengths and angles are not changed by differentiating.
- **No torsion:** $\Gamma^{\rho}{}_{\mu\nu}=\Gamma^{\rho}{}_{\nu\mu}$ in a coordinate basis.

Sign and index conventions are in [Conventions](../conventions.md).

## 2. Deriving the Christoffel symbols

Lower the first index and write $\Gamma_{\sigma\lambda\mu}=g_{\sigma\rho}\Gamma^{\rho}{}_{\lambda\mu}$. Metric compatibility says

$$
\partial_\lambda g_{\mu\nu}=\Gamma_{\nu\lambda\mu}+\Gamma_{\mu\lambda\nu}.
$$

Write the same equation with the indices $(\lambda,\mu,\nu)$ permuted cyclically:

$$
\partial_\mu g_{\nu\lambda}=\Gamma_{\lambda\mu\nu}+\Gamma_{\nu\mu\lambda},
\qquad
\partial_\nu g_{\lambda\mu}=\Gamma_{\mu\nu\lambda}+\Gamma_{\lambda\nu\mu}.
$$

Add the first two equations and subtract the third. Because $\Gamma_{\sigma\mu\nu}=\Gamma_{\sigma\nu\mu}$, every term cancels except two copies of $\Gamma_{\nu\lambda\mu}$. Raising the index again gives

$$
\Gamma^{\rho}{}_{\mu\nu}=\frac{1}{2}\,g^{\rho\sigma}\left(\partial_\mu g_{\nu\sigma}+\partial_\nu g_{\mu\sigma}-\partial_\sigma g_{\mu\nu}\right).
$$

The two assumptions determine the connection completely.

## 3. The connection is not a tensor

If we change coordinates from $x$ to $x'$, the symbols change by

$$
\Gamma'^{\alpha}{}_{\beta\gamma}=\frac{\partial x'^{\alpha}}{\partial x^{\rho}}\,\frac{\partial x^{\mu}}{\partial x'^{\beta}}\,\frac{\partial x^{\nu}}{\partial x'^{\gamma}}\,\Gamma^{\rho}{}_{\mu\nu}+\frac{\partial x'^{\alpha}}{\partial x^{\rho}}\,\frac{\partial^{2}x^{\rho}}{\partial x'^{\beta}\,\partial x'^{\gamma}}.
$$

The first term is the usual tensor law. The second term, with second derivatives, is what makes $\Gamma$ different from a tensor. It is the reason the symbols can vanish in one coordinate system and not in another.

**Worked example: the flat plane.** In Cartesian coordinates all symbols vanish. In polar coordinates $(r,\varphi)$ the metric is $dr^{2}+r^{2}d\varphi^{2}$ and the only non-zero symbols are

$$
\Gamma^{r}{}_{\varphi\varphi}=-r,\qquad \Gamma^{\varphi}{}_{r\varphi}=\Gamma^{\varphi}{}_{\varphi r}=\frac{1}{r}.
$$

The plane did not change; only the description did. The code reproduces these polar symbols from the zero Cartesian symbols using the law above, which a purely tensorial law could not do (it would map zero to zero).

## 4. The covariant derivative

For a vector and a covector,

$$
\nabla_\mu V^{\nu}=\partial_\mu V^{\nu}+\Gamma^{\nu}{}_{\mu\lambda}V^{\lambda},
\qquad
\nabla_\mu \omega_{\nu}=\partial_\mu \omega_{\nu}-\Gamma^{\lambda}{}_{\mu\nu}\,\omega_{\lambda}.
$$

Each index of a general tensor gets one correction term: plus for an upper index, minus for a lower index. The extra terms exactly cancel the inhomogeneous part of the transformation law, so $\nabla V$ is a genuine tensor even though $\partial V$ is not.

**Why this makes sense.** The components of a vector change along a curve partly because the vector changes and partly because the coordinate basis vectors change. The connection removes the second effect.

**Worked example.** The constant Cartesian vector $\partial_x$ has polar components $\left(\cos\varphi,\,-\sin\varphi/r\right)$. Its partial derivatives are non-zero, but its covariant derivative is zero, as it must be for a vector that is constant. The test suite checks this.

## 5. Curvature

Differentiate a vector twice and take the antisymmetric part:

$$
\left[\nabla_\mu,\nabla_\nu\right]V^{\rho}=R^{\rho}{}_{\sigma\mu\nu}\,V^{\sigma}.
$$

Expanding $\nabla_\mu\nabla_\nu V^{\rho}$, the terms containing derivatives of $V$ cancel in the antisymmetric combination (this uses zero torsion), and what remains is multiplication by

$$
R^{\rho}{}_{\sigma\mu\nu}=\partial_\mu\Gamma^{\rho}{}_{\nu\sigma}-\partial_\nu\Gamma^{\rho}{}_{\mu\sigma}+\Gamma^{\rho}{}_{\mu\lambda}\Gamma^{\lambda}{}_{\nu\sigma}-\Gamma^{\rho}{}_{\nu\lambda}\Gamma^{\lambda}{}_{\mu\sigma}.
$$

A space is flat exactly when this tensor vanishes. Because it is a tensor, it vanishes in one coordinate system only if it vanishes in all of them. That is why flat space gives zero curvature in polar coordinates even though $\Gamma$ does not vanish there.

**Worked example: the unit sphere.** For $d\theta^{2}+\sin^{2}\theta\,d\varphi^{2}$ the code finds the Ricci scalar $R=2$, which is positive as it should be for the convention used.

### Symmetries and identities

The Riemann tensor with all indices down satisfies

$$
R_{\rho\sigma\mu\nu}=-R_{\sigma\rho\mu\nu}=-R_{\rho\sigma\nu\mu}=R_{\mu\nu\rho\sigma},
\qquad
R_{\rho[\sigma\mu\nu]}=0,
$$

and the differential Bianchi identity and its contraction,

$$
\nabla_{[\lambda}R_{\rho\sigma]\mu\nu}=0,
\qquad
\nabla_\mu G^{\mu\nu}=0.
$$

The Ricci tensor is symmetric. None of these is used to build the code; each is a test that the code's output must pass.

## 6. The Weyl tensor and dimension

The part of the curvature not fixed by the Ricci tensor is the Weyl tensor. In $n$ dimensions,

$$
C_{abcd}=R_{abcd}-\frac{2}{n-2}\left(g_{a[c}R_{d]b}-g_{b[c}R_{d]a}\right)+\frac{2}{(n-1)(n-2)}\,R\,g_{a[c}g_{d]b}.
$$

Facts that the tests confirm:

- It is trace-free.
- In three dimensions it vanishes identically, because the Riemann tensor is then determined by the Ricci tensor.
- In two dimensions it is not defined, and the code raises an error.
- In four dimensions the Kretschmann scalar splits as $K=C_{abcd}C^{abcd}+2R_{ab}R^{ab}-\frac{1}{3}R^{2}$.
- FLRW spacetimes are conformally flat, so their Weyl tensor is zero.
- The Gauss-Bonnet combination $R^{2}-4R_{ab}R^{ab}+K$ is a total derivative in four dimensions.

## 7. Implementation

| Task | Module |
|---|---|
| Christoffel symbols, torsion, metric-compatibility residual | `geometry/connection.py` |
| Covariant derivative of any rank, divergence | `geometry/covariant_derivative.py` |
| Riemann tensor | `geometry/riemann.py` |
| Ricci tensor and scalar | `geometry/ricci.py` |
| Einstein tensor | `geometry/einstein.py` |
| Weyl tensor | `geometry/weyl.py` |
| Symmetry and Bianchi residuals | `geometry/bianchi.py` |
| Curvature invariants | `geometry/invariants.py` |

## 8. How it is checked

The tests are in `tests/test_connection.py`, `tests/test_covariant_derivative.py` and `tests/test_curvature.py`. They check:

- flat space in Cartesian, polar, spherical and Minkowski-spherical coordinates;
- the known polar and spherical Christoffel symbols;
- zero torsion and metric compatibility on Schwarzschild;
- that the inhomogeneous transformation law reproduces the polar symbols from the Cartesian ones;
- the Leibniz rule, and that contraction commutes with differentiation;
- all Riemann symmetries and both Bianchi identities;
- the Weyl properties listed above;
- **negative controls:** a deliberately wrong connection must fail the metric-compatibility test, and a deliberately wrong Riemann tensor must fail the cyclic identity.

## 9. Limits

Tensors are stored as explicit component arrays. There is no abstract index algebra and no automatic use of index symmetries. Simplification is done by SymPy, which cannot always prove that an expression is zero. When that happens the code falls back to sampling and reports that the evidence is numerical.
