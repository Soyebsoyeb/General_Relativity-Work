# 2. The Einstein-Hilbert action and the field equations

**In one paragraph.** The Einstein field equations can be obtained by asking which metrics make a single number, the action, stationary. This page works through that variation step by step and records exactly which parts the code verifies.

## 1. Setting

The action is

$$
S=\frac{c^{3}}{16\pi G}\int\left(R-2\Lambda\right)\sqrt{-g}\;d^{4}x+S_{\mathrm{m}},
$$

where $S_{\mathrm{m}}$ is the matter action and $\Lambda$ the cosmological constant. The stress-energy tensor is defined by

$$
T_{\mu\nu}=-\frac{2}{\sqrt{-g}}\,\frac{\delta S_{\mathrm{m}}}{\delta g^{\mu\nu}}.
$$

We vary the metric, $g_{\mu\nu}\to g_{\mu\nu}+\epsilon\,h_{\mu\nu}$, where $h_{\mu\nu}$ is a symmetric tensor that vanishes outside a compact region, and keep terms to first order in $\epsilon$.

## 2. Three elementary variations

**The determinant.** For any matrix, $\delta\ln\det g=g^{\mu\nu}\delta g_{\mu\nu}$. Therefore

$$
\delta\sqrt{-g}=\frac{1}{2}\sqrt{-g}\;g^{\mu\nu}h_{\mu\nu}.
$$

**The inverse metric.** Varying $g^{\mu\alpha}g_{\alpha\nu}=\delta^{\mu}{}_{\nu}$ gives

$$
\delta g^{\mu\nu}=-h^{\mu\nu}.
$$

**The Ricci scalar.** The difference between two connections is a tensor, namely

$$
\delta\Gamma^{\rho}{}_{\mu\nu}=\frac{1}{2}\,g^{\rho\sigma}\left(\nabla_\mu h_{\nu\sigma}+\nabla_\nu h_{\mu\sigma}-\nabla_\sigma h_{\mu\nu}\right).
$$

Varying the Riemann formula then gives the Palatini identity $\delta R^{\rho}{}_{\sigma\mu\nu}=\nabla_\mu\delta\Gamma^{\rho}{}_{\nu\sigma}-\nabla_\nu\delta\Gamma^{\rho}{}_{\mu\sigma}$. Contracting it,

$$
\delta R=-R^{\mu\nu}h_{\mu\nu}+\nabla_\mu\nabla_\nu h^{\mu\nu}-\square h,
\qquad h=g^{\mu\nu}h_{\mu\nu}.
$$

## 3. The boundary term

The last two terms of $\delta R$ combine into a divergence,

$$
\nabla_\mu\nabla_\nu h^{\mu\nu}-\square h=\nabla_\mu V^{\mu},
\qquad
V^{\mu}=\nabla_\nu h^{\mu\nu}-\nabla^{\mu}h.
$$

For any vector, $\sqrt{-g}\,\nabla_\mu V^{\mu}=\partial_\mu\left(\sqrt{-g}\,V^{\mu}\right)$, which integrates to a surface term. Because $h$ vanishes at the boundary, this term drops out.

If the region has a boundary where $h$ does not vanish, the term does not drop out. It is then cancelled by adding the Gibbons-Hawking-York boundary term. That case is described here but is **not implemented**.

## 4. The result

Putting the pieces together,

$$
\delta\left[\sqrt{-g}\left(R-2\Lambda\right)\right]=-\sqrt{-g}\left(G^{\mu\nu}+\Lambda g^{\mu\nu}\right)h_{\mu\nu}+\partial_\mu\left(\sqrt{-g}\,V^{\mu}\right).
$$

The matter action varies as $\delta S_{\mathrm m}=\frac{1}{2}\int\sqrt{-g}\,T^{\mu\nu}h_{\mu\nu}\,d^4x$. Requiring the total variation to vanish for every $h_{\mu\nu}$ gives

$$
G_{\mu\nu}+\Lambda\,g_{\mu\nu}=\frac{8\pi G}{c^{4}}\,T_{\mu\nu}.
$$

Two consequences:

- In vacuum with $\Lambda=0$ this reduces to $R_{\mu\nu}=0$.
- Diffeomorphism invariance of the matter action implies $\nabla_\mu T^{\mu\nu}=0$ when the matter equations hold. This agrees with $\nabla_\mu G^{\mu\nu}=0$.

## 5. How the code checks this

The module is `field_equations/einstein_hilbert.py`. It does not contain a general variational engine. Instead it tests each identity above on explicit examples.

The test background is a three-dimensional metric with non-vanishing Einstein tensor, and $h_{\mu\nu}$ is a non-trivial symmetric perturbation with coordinate-dependent entries. The following are verified symbolically:

1. the variation of $\sqrt{|g|}$;
2. the variation of the inverse metric;
3. the variation of $R$, including the $\nabla\nabla h$ and $\square h$ terms;
4. that $\sqrt{-g}\,\nabla_\mu V^{\mu}$ equals $\partial_\mu(\sqrt{-g}\,V^{\mu})$;
5. (slow suite) the complete identity of Section 4 including $\Lambda$.

A **negative control** changes the sign of the $\square h$ term and confirms the test then fails, so the test is capable of detecting an error.

A second check uses a reduced action. For the metric $-A(r)\,dt^{2}+B(r)\,dr^{2}+r^{2}d\Omega^{2}$, the Euler-Lagrange equations of $\sqrt{-g}\,R$ with respect to $A$ and $B$ equal $\sqrt{-g}\,G^{tt}$ and $-\sqrt{-g}\,G^{rr}$ exactly. Doing this reduction before varying is justified by the principle of symmetric criticality, which holds for the compact rotation group.

## 6. Limits

- No general symbolic variational engine: identities are checked for the supplied backgrounds and perturbations, not proved for every metric.
- Boundary terms beyond the compact-support case are not implemented.
