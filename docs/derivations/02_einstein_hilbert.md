# 2. Einstein–Hilbert action and the field equations

## Setting
$$S=\frac{c^3}{16\pi G}\int(R-2\Lambda)\sqrt{-g}\,d^4x+S_{\rm m},\qquad T_{\mu\nu}=-\frac{2}{\sqrt{-g}}\frac{\delta S_{\rm m}}{\delta g^{\mu\nu}}.$$
Variations are $g_{\mu\nu}\to g_{\mu\nu}+\epsilon h_{\mu\nu}$ with $h$ compactly supported unless a boundary term is discussed.

## Three elementary variations
**Determinant.** $\delta\ln\det g=g^{\mu\nu}\delta g_{\mu\nu}$, hence $\delta\sqrt{-g}=\tfrac12\sqrt{-g}\,g^{\mu\nu}h_{\mu\nu}$.

**Inverse.** From $g^{\mu\alpha}g_{\alpha\nu}=\delta^\mu_\nu$: $\delta g^{\mu\nu}=-g^{\mu\alpha}g^{\nu\beta}h_{\alpha\beta}$.

**Ricci scalar.** The difference of two connections is a tensor, $\delta\Gamma^\rho{}_{\mu\nu}=\tfrac12 g^{\rho\sigma}(\nabla_\mu h_{\nu\sigma}+\nabla_\nu h_{\mu\sigma}-\nabla_\sigma h_{\mu\nu})$. Varying the Riemann formula gives the Palatini identity $\delta R^\rho{}_{\sigma\mu\nu}=\nabla_\mu\delta\Gamma^\rho{}_{\nu\sigma}-\nabla_\nu\delta\Gamma^\rho{}_{\mu\sigma}$; contracting,
$$\delta R=-R^{\mu\nu}h_{\mu\nu}+\nabla_\mu\nabla_\nu h^{\mu\nu}-\Box h,\qquad h=g^{\mu\nu}h_{\mu\nu}.$$

## Boundary term
$\nabla_\mu\nabla_\nu h^{\mu\nu}-\Box h=\nabla_\mu V^\mu$ with $V^\mu=\nabla_\nu h^{\mu\nu}-\nabla^\mu h$, and $\sqrt{-g}\,\nabla_\mu V^\mu=\partial_\mu(\sqrt{-g}V^\mu)$. This integrates to a surface term that vanishes for compactly supported $h$. On a manifold with boundary it does **not** vanish; it is cancelled by the Gibbons–Hawking–York term $\frac1{8\pi}\oint K\sqrt{|\gamma|}\,d^3x$ (discussed here, **not implemented**).

## Result
$$\delta\bigl[\sqrt{-g}(R-2\Lambda)\bigr]=-\sqrt{-g}\bigl(G^{\mu\nu}+\Lambda g^{\mu\nu}\bigr)h_{\mu\nu}+\partial_\mu(\sqrt{-g}V^\mu).$$
Adding $\delta S_{\rm m}=\tfrac12\int\sqrt{-g}\,T^{\mu\nu}h_{\mu\nu}$ and demanding stationarity,
$$G_{\mu\nu}+\Lambda g_{\mu\nu}=8\pi G c^{-4}T_{\mu\nu}.$$
Diffeomorphism invariance of $S_{\rm m}$ gives $\nabla_\mu T^{\mu\nu}=0$ on shell, consistent with $\nabla_\mu G^{\mu\nu}=0$. In vacuum with $\Lambda=0$: $R_{\mu\nu}=0$.

## Implementation and validation (`field_equations/einstein_hilbert.py`)
Each identity is verified symbolically on a curved three-dimensional background with a non-trivial perturbation $h_{\mu\nu}$ (`tests/test_field_equations.py`): $\delta\sqrt g$, $\delta g^{-1}$, $\delta R$ (with a negative control showing a sign error is detected), the total-divergence property of $V^\mu$, and, in the slow suite, the full identity with $\Lambda$. A reduced action for $ds^2=-A\,dt^2+B\,dr^2+r^2d\Omega^2$ is also varied: the Euler–Lagrange equations equal $\sqrt{-g}\,G^{tt}$ and $-\sqrt{-g}\,G^{rr}$ identically. This relies on the principle of symmetric criticality, valid for the compact group $SO(3)$.

## Limitations
There is no general symbolic variational engine: the identities are checked for supplied backgrounds and perturbations, not proven for arbitrary $g$. Boundary terms beyond the compact-support case are not implemented.
