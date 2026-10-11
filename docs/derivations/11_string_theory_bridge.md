# 11. String theory: a mathematical bridge (explanatory; no code)

This page is deliberately **separate from the classical GR engine**. Nothing here is implemented or verified by the software; it records what the bridge to quantum gravity looks like and where the repository stops.

## Worldsheet and target space
A string sweeps a surface $X^\mu(\sigma,\tau)$ in a target manifold $(\mathcal M,G_{\mu\nu})$. The bosonic non-linear sigma model is
$$S=\frac1{4\pi\alpha'}\int d^2\sigma\sqrt{h}\,h^{ab}\partial_aX^\mu\partial_bX^\nu G_{\mu\nu}(X)+\ldots$$
with the Polyakov metric $h_{ab}$.

## Emergence of gravity
Quantising the flat-space string gives a massless spin-2 state (a graviton) together with a dilaton and antisymmetric tensor. A coherent state of such states corresponds to a curved background $G_{\mu\nu}$. Consistency of the quantum worldsheet theory (Weyl invariance) requires the beta functions to vanish:
$$\beta^G_{\mu\nu}=\alpha'\Bigl(R_{\mu\nu}+2\nabla_\mu\nabla_\nu\Phi-\tfrac14H_{\mu\lambda\rho}H_\nu{}^{\lambda\rho}\Bigr)+O(\alpha'^2).$$
In the simplest case, constant dilaton and vanishing $H$, the leading-order condition is $R_{\mu\nu}=0$: Ricci-flatness, the vacuum Einstein equation, now arising as a consistency condition.

## Corrections
Higher orders in $\alpha'$ add curvature-squared and higher terms (for the bosonic string, $R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$-type corrections at $O(\alpha')$ beyond the leading term). These are the string analogue of the EFT corrections in [page 10](10_perturbative_gravity_and_eft.md), and unlike classical GR they come with a specific UV completion and extra fields.

## What is and is not claimed
* The repository can check that a *given* metric is Ricci flat (Schwarzschild, Kerr) and can compute the curvature invariants entering the correction terms.
* It does **not** solve string theory, derive the beta functions, or compute string amplitudes.
* Distinction: classical GR (this engine) $\ne$ GR as the low-energy limit of string theory $\ne$ a quantum theory of gravity.
