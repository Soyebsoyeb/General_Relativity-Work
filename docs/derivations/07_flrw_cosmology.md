# 7. FLRW cosmology

**In one paragraph.** On large scales the universe is nearly the same everywhere and in every direction. Imposing that symmetry on Einstein's equations gives the Friedmann equations, which describe how the scale of the universe changes in time. This page derives them from the code's Einstein tensor and compares numerical expansion histories with exact solutions.

## 1. Setting

The metric of a homogeneous and isotropic universe is

$$
ds^{2}=-dt^{2}+a(t)^{2}\left[\frac{dr^{2}}{1-kr^{2}}+r^{2}\,d\Omega^{2}\right],
$$

where $a(t)$ is the scale factor and $k=0,\pm1$ describes the spatial curvature. The matter is a perfect fluid with energy density $\rho$ and pressure $p$, at rest in these coordinates:

$$
T_{\mu\nu}=\left(\rho+p\right)u_\mu u_\nu+p\,g_{\mu\nu},\qquad u=\partial_t .
$$

## 2. Deriving the Friedmann equations

The code computes the Einstein tensor of this metric. The relevant components are

$$
G^{t}{}_{t}=-\frac{3\left(\dot a^{2}+k\right)}{a^{2}},
\qquad
G^{r}{}_{r}=-\frac{2a\ddot a+\dot a^{2}+k}{a^{2}}.
$$

The test also confirms that $G^{r}{}_{r}=G^{\theta}{}_{\theta}=G^{\varphi}{}_{\varphi}$, which is the isotropy of the universe.

For the fluid, $T^{t}{}_{t}=-\rho$ and $T^{r}{}_{r}=p$. The field equation $G^{\mu}{}_{\nu}+\Lambda\delta^{\mu}{}_{\nu}=8\pi T^{\mu}{}_{\nu}$ then gives the two Friedmann equations:

$$
\left(\frac{\dot a}{a}\right)^{2}+\frac{k}{a^{2}}=\frac{8\pi}{3}\,\rho+\frac{\Lambda}{3},
$$

$$
\frac{\ddot a}{a}=-\frac{4\pi}{3}\left(\rho+3p\right)+\frac{\Lambda}{3}.
$$

## 3. Conservation of energy comes for free

A third equation, $\dot\rho+3H\left(\rho+p\right)=0$ with $H=\dot a/a$, was **not** imposed. The test solves the two Einstein equations for $\rho$ and $p$ and then shows symbolically that this equation is satisfied automatically. It follows from the contracted Bianchi identity $\nabla_\mu G^{\mu\nu}=0$.

## 4. Numerical solution

Use units where today's Hubble rate is 1 and today's scale factor is 1. For dust, radiation, curvature and a cosmological constant, define

$$
E^{2}(a)=\Omega_{r}\,a^{-4}+\Omega_{m}\,a^{-3}+\Omega_{k}\,a^{-2}+\Omega_{\Lambda},
\qquad
\Omega_{k}=1-\Omega_{r}-\Omega_{m}-\Omega_{\Lambda}.
$$

The code integrates the **second-order** equation for $\ddot a$, which also works through a recollapse, and monitors the first-order constraint $\left(\dot a/a\right)^{2}-E^{2}(a)$ as a measure of accuracy.

## 5. Exact solutions used for comparison

| Model | Exact solution |
|---|---|
| Flat, matter only | $a=\left(\frac{3t}{2}\right)^{2/3}$ |
| Flat, radiation only | $a=\left(2t\right)^{1/2}$ |
| de Sitter ($\Omega_{\Lambda}=1$) | $a=e^{t}$ |
| Flat, matter plus $\Lambda$ | $a=\left(\frac{\Omega_{m}}{\Omega_{\Lambda}}\right)^{1/3}\sinh^{2/3}\left(\frac{3}{2}\sqrt{\Omega_{\Lambda}}\,t\right)$ |
| Closed, matter only ($\Omega_m>1$) | parametric: $a=\frac{a_{\max}}{2}\left(1-\cos\eta\right)$ and $t=\frac{\Omega_{m}}{2\lvert\Omega_{k}\rvert^{3/2}}\left(\eta-\sin\eta\right)$, with $a_{\max}=\Omega_m/\lvert\Omega_k\rvert$ |

The redshift is $1+z=1/a$.

## 6. How it is checked

- The symbolic checks (Sections 2 and 3) are in `tests/test_flrw.py`.
- Each numerical run starts from the exact solution's value and slope at an initial time, then is compared with the exact solution later. All five cases agree to better than about $10^{-9}$ in relative error, with a Friedmann-constraint residual below $10^{-7}$ (experiment 08).
- A **negative control** integrates matter dynamics but compares with the radiation solution; the error is larger than 1%, so the comparison can detect a wrong model.

Note what this does and does not test. The numerical comparison tests the integrator and the form of the equations. That the equations follow from $G=8\pi T$ is tested separately by the symbolic checks.

## 7. Limits

Background expansion only. There are no perturbations, no structure formation and no Boltzmann physics.
