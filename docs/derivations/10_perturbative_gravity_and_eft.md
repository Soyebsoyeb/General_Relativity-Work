# 10. Perturbative gravity and effective field theory

**In one paragraph.** One route towards quantum gravity treats the gravitational field as a small ripple on a fixed background and expands the action in powers of that ripple. This page shows what that expansion gives at lowest order, and what happens when extra curvature terms are added to the action. It does **not** quantise gravity.

## 1. Status

Partly implemented, and only symbolically. There is no graviton quantisation, no propagator, no loop calculation and no scattering amplitude anywhere in this repository.

## 2. Expanding the action

Write $g_{\mu\nu}=\eta_{\mu\nu}+\epsilon\,h_{\mu\nu}$. The expansion of the Einstein-Hilbert Lagrangian has the form

$$
\sqrt{-g}\,R=\epsilon\,\left(\text{total derivative}\right)+\epsilon^{2}\,\mathcal{L}^{(2)}\left[h\right]+O\left(\epsilon^{3}\right).
$$

The linear term is a total derivative, so it does not affect the equations of motion. The quadratic part $\mathcal{L}^{(2)}$ is the free field theory of the ripple.

The code computes $\mathcal{L}^{(2)}$ symbolically, for ten independent component functions $h_{\mu\nu}(t,x)$ in four dimensions, and checks two things.

**First check.** Varying $\mathcal{L}^{(2)}$ with respect to $h_{\mu\nu}$ gives minus the linearised Einstein tensor,

$$
\frac{\delta\mathcal{L}^{(2)}}{\delta h_{\mu\nu}}=-G^{(1)\mu\nu}\left[h\right]
$$

(with a factor of 2 for the off-diagonal components, which each stand for two equal entries $h_{\mu\nu}=h_{\nu\mu}$).

**Second check.** The same equations follow from the Fierz-Pauli Lagrangian,

$$
\mathcal{L}_{\mathrm{FP}}=-\frac{1}{2}\,\partial_\lambda h_{\mu\nu}\,\partial^{\lambda}h^{\mu\nu}+\partial_\mu h^{\mu\nu}\,\partial^{\lambda}h_{\lambda\nu}-\partial_\mu h^{\mu\nu}\,\partial_\nu h+\frac{1}{2}\,\partial_\lambda h\,\partial^{\lambda}h .
$$

The test shows that $\mathcal{L}^{(2)}$ and $\frac{1}{2}\mathcal{L}_{\mathrm{FP}}$ have identical Euler-Lagrange equations, so they differ only by a total derivative.

The quadratic action has the gauge symmetry $h_{\mu\nu}\to h_{\mu\nu}+\partial_\mu\xi_\nu+\partial_\nu\xi_\mu$, and its propagating modes are the two polarisations described on [page 6](06_linearized_gravity_and_gravitational_waves.md). Higher-order interaction terms come from $\mathcal{L}^{(3)}$, $\mathcal{L}^{(4)}$ and so on; these are not computed.

## 3. Dimensional analysis

Gravity has a coupling with dimensions, $\kappa^{2}\sim 32\pi G$. After rescaling $h\to\kappa\,h$, each interaction vertex carries extra powers of $\kappa$ and extra derivatives. The expansion parameter for loop effects is therefore

$$
\frac{G E^{2}}{\hbar c^{5}}\sim\left(\frac{E}{M_{\mathrm{P}}}\right)^{2},
$$

where $E$ is the energy scale and $M_{\mathrm{P}}$ is the Planck mass. The theory is non-renormalisable, which means it can be used as an **effective field theory**, reliable only for $E\ll M_{\mathrm{P}}$.

## 4. Higher-curvature corrections

An effective theory includes every term allowed by symmetry, suppressed by powers of the Planck length $\ell_{P}$. At the next orders the action contains

$$
R^{2},\qquad R_{\mu\nu}R^{\mu\nu},\qquad R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma},\qquad \text{cubic terms in the Riemann tensor}.
$$

**A concrete, exact example.** For static spherical metrics, the code derives the field equations of a reduced Lagrangian containing each of these terms, and evaluates them on the Schwarzschild solution. The result:

- the $R^{2}$ and $R_{\mu\nu}R^{\mu\nu}$ terms contribute nothing;
- the $R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$ term and the Gauss-Bonnet combination contribute nothing;
- a cubic term in the Riemann tensor **does** contribute.

So in four dimensions, Schwarzschild remains an exact solution when any quadratic curvature terms are added, and vacuum corrections first appear at cubic order. This is a statement about the structure of the equations. It says nothing about how large any coefficient is. The check takes about 20 seconds and is in the slow test suite.

## 5. When is the expansion trustworthy?

The natural expansion parameter is $\ell_{P}^{2}\sqrt{K}$, where $K$ is the Kretschmann scalar. For Schwarzschild, $\sqrt{K}=4\sqrt{3}\,GM/\left(c^{2}r^{3}\right)$. The code evaluates this at the horizon:

| Black hole | Expansion parameter at the horizon |
|---|---|
| One solar mass | about $10^{-76}$ |
| $4.3\times10^{6}$ solar masses | about $6\times10^{-90}$ |
| One Planck mass | about $0.87$ |

For astrophysical black holes the effective theory is extremely accurate near the horizon. For a Planck-mass object the expansion parameter is of order one, and the derivative expansion breaks down. A low-energy expansion is not a fundamental ultraviolet completion.

## 6. Limits

The quadratic-action tests restrict $h_{\mu\nu}$ to depend on two coordinates, for cost reasons and not as a matter of principle. There is no gauge fixing and no ghost sector, and no coupling to matter.
