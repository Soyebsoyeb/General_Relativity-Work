# 11. A note on string theory

**In one paragraph.** This page is a short explanation of how string theory connects to general relativity. It is explanatory text only. **Nothing here is implemented or verified by the code**, and it is kept apart from the classical engine on purpose.

## 1. Strings in a curved target space

A string sweeps out a surface $X^{\mu}(\sigma,\tau)$ in a target manifold with metric $G_{\mu\nu}$. The simplest description is a two-dimensional field theory on the surface, the non-linear sigma model, whose action is

$$
S=\frac{1}{4\pi\alpha'}\int d^{2}\sigma\,\sqrt{h}\;h^{ab}\,\partial_a X^{\mu}\,\partial_b X^{\nu}\,G_{\mu\nu}\left(X\right)+\cdots
$$

Here $h_{ab}$ is a metric on the surface and $\alpha'$ is the square of the string length scale.

## 2. How gravity appears

Quantising a string in flat space produces a massless particle with spin two. It has all the properties of a graviton. It comes together with a scalar (the dilaton) and an antisymmetric tensor field. A coherent state of such particles corresponds to a curved target-space metric $G_{\mu\nu}$.

For the quantum theory on the surface to be consistent, its scale symmetry (Weyl invariance) must survive quantisation. This requires the beta functions to vanish. To lowest order in $\alpha'$,

$$
\beta^{G}_{\mu\nu}=\alpha'\left(R_{\mu\nu}+2\nabla_\mu\nabla_\nu\Phi-\frac{1}{4}H_{\mu\lambda\rho}H_{\nu}{}^{\lambda\rho}\right)+O\left(\alpha'^{2}\right).
$$

In the simplest case, with constant dilaton $\Phi$ and no antisymmetric field strength $H$, the condition is

$$
R_{\mu\nu}=0,
$$

the vacuum Einstein equation. In this picture the Einstein equation is not postulated but arises as a consistency condition on the string.

## 3. Corrections

Higher orders in $\alpha'$ add terms quadratic and higher in the curvature. They are the string-theory counterpart of the effective-field-theory corrections of [page 10](10_perturbative_gravity_and_eft.md). Unlike the generic effective theory, they come with a definite ultraviolet completion and extra fields.

## 4. What this repository does and does not do

- It can check that a given metric is Ricci flat, as for Schwarzschild or Kerr. It can also compute the curvature scalars that enter the correction terms.
- It does **not** solve string theory, derive the beta functions, or compute any string amplitude.
- Three things are different and must not be confused: classical general relativity (this engine), general relativity as the low-energy limit of string theory, and a full quantum theory of gravity.
