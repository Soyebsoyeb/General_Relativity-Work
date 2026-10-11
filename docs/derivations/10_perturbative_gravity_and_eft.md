# 10. Perturbative gravity and effective field theory

## Status: partially implemented, symbolic only
There is **no graviton quantisation, propagator, loop calculation or scattering amplitude** in this repository.

## Expansion of the action
Writing $g=\eta+\epsilon h$,
$$\sqrt{-g}R=\epsilon\,(\text{total derivative})+\epsilon^2\,\mathcal L^{(2)}[h]+O(\epsilon^3).$$
The $O(\epsilon)$ term is a divergence (the tadpole vanishes about a solution). The code computes $\mathcal L^{(2)}$ with truncated series and verifies, with ten independent component functions $h_{\mu\nu}(t,x)$ in four dimensions:

1. $\delta\mathcal L^{(2)}/\delta h_{\mu\nu}=-G^{(1)\mu\nu}[h]$ (factor 2 off-diagonal for the symmetric pair);
2. $\mathcal L^{(2)}$ and $\tfrac12\mathcal L_{\rm FP}$ have identical Euler–Lagrange equations, where
$$\mathcal L_{\rm FP}=-\tfrac12\partial_\lambda h_{\mu\nu}\partial^\lambda h^{\mu\nu}+\partial_\mu h^{\mu\nu}\partial^\lambda h_{\lambda\nu}-\partial_\mu h^{\mu\nu}\partial_\nu h+\tfrac12\partial_\lambda h\partial^\lambda h,$$
so $\mathcal L^{(2)}=\tfrac12\mathcal L_{\rm FP}$ up to a total derivative (the Fierz–Pauli structure).

The quadratic action inherits the gauge invariance $\delta h=\partial\xi+\partial\xi$; propagating modes are the two transverse-traceless polarisations. Higher-order interaction vertices come from $\mathcal L^{(3)},\mathcal L^{(4)},\dots$ (not computed).

## Dimensional analysis
With $\kappa^2=32\pi G$ and $h\to\kappa h$, interactions scale as $\kappa^n\partial^2h^n$, so the loop expansion parameter is $G E^2/\hbar c^5\sim(E/M_{\rm P})^2$: GR is non-renormalisable, and is treated as an effective field theory valid for $E\ll M_{\rm P}$.

## Higher-curvature corrections
Add $\alpha R^2+\beta R_{\mu\nu}R^{\mu\nu}+\gamma R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}+\delta\,\mathrm{Riem}^3$ with coefficients suppressed by powers of $\ell_P^2$. **Concrete symbolic example** (slow test; reduced static-spherical action): evaluating the Euler–Lagrange equations on Schwarzschild, the $R^2$, $R_{\mu\nu}^2$, $R_{\mu\nu\rho\sigma}^2$ and Gauss–Bonnet contributions all vanish identically, while the cubic term does not. So Schwarzschild is unchanged by all quadratic curvature terms in 4D, and vacuum corrections first appear at cubic order. This is an exact statement about the reduced equations, not about the sizes of the coefficients.

**Validity.** The expansion parameter is $\ell_P^2\sqrt K$ ($\sqrt K=4\sqrt3\,GM/c^2r^3$ for Schwarzschild). `curvature_expansion_parameter` returns values of about $10^{-76}$ at the horizon of a solar-mass hole (and $0.87$ at one Planck mass), i.e. order unity at the Planck mass, where the derivative expansion breaks down. A low-energy expansion is not a fundamental ultraviolet completion.

## Limitations
Quadratic action only for $h$ depending on two coordinates (cost, not principle); no gauge fixing/ghosts; no matter coupling.
