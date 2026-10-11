# 3. The Schwarzschild solution, orbits and observables

**In one paragraph.** Schwarzschild spacetime is the exterior of a non-rotating, uncharged spherical mass. It is the main end-to-end test of the whole code: the metric is derived from the vacuum equations, its curvature is computed, particles and light are followed through it, and each measurable effect is computed in two independent ways.

## 1. Assumptions

The spacetime is static, spherically symmetric, and a vacuum solution. Staticity is **assumed**. That it follows from spherical symmetry is Birkhoff's theorem, which is not proved here.

## 2. Deriving the metric

Start from the most general static spherical form

$$
ds^{2}=-A(r)\,dt^{2}+B(r)\,dr^{2}+r^{2}\,d\Omega^{2},\qquad d\Omega^{2}=d\theta^{2}+\sin^{2}\theta\,d\varphi^{2}.
$$

The code computes the Einstein tensor for this metric. Two of its components are

$$
G^{t}{}_{t}=-\frac{1}{r^{2}}\,\frac{d}{dr}\left[r\left(1-\frac{1}{B}\right)\right],
\qquad
G^{r}{}_{r}-G^{t}{}_{t}=\frac{1}{rB}\left(\frac{A'}{A}+\frac{B'}{B}\right).
$$

Both forms are checked symbolically. Now solve the vacuum equations in three steps.

1. $G^{t}{}_{t}=0$ means $r\left(1-1/B\right)$ is constant, so $1/B=1+C/r$ for a constant $C$.
2. $G^{r}{}_{r}-G^{t}{}_{t}=0$ means $AB$ is constant. Rescaling the time coordinate sets this constant to 1, so $A=1/B$.
3. Far away the metric must reduce to the Newtonian form $g_{tt}\to-(1+2\Phi)$ with $\Phi=-M/r$. This fixes $C=-2M$.

The result is

$$
ds^{2}=-\left(1-\frac{2M}{r}\right)dt^{2}+\left(1-\frac{2M}{r}\right)^{-1}dr^{2}+r^{2}\,d\Omega^{2}.
$$

The test suite confirms that all ten components of the Einstein tensor vanish for this metric, while the Riemann tensor does **not**: vacuum does not mean flat. The Kretschmann scalar is

$$
K=R_{abcd}R^{abcd}=\frac{48\,M^{2}}{r^{6}}.
$$

## 3. The horizon and the singularity

The metric components misbehave at $r=2M$, but this is a feature of the coordinates. Introduce $v=t+r_{*}$ with $r_{*}=r+2M\ln\left(r/2M-1\right)$. Then

$$
ds^{2}=-\left(1-\frac{2M}{r}\right)dv^{2}+2\,dv\,dr+r^{2}\,d\Omega^{2},
$$

whose determinant $-r^{4}\sin^{2}\theta$ is finite and non-zero at $r=2M$. The code obtains this form by transforming the metric, finds it Ricci flat, and finds the same $K$.

At $r=0$ the scalar $K$ diverges, and no coordinate change can remove that. It is a genuine curvature singularity. A word of caution: a single scalar is not a universal detector of singularities. Plane-wave spacetimes, for example, have all polynomial curvature scalars equal to zero.

## 4. Geodesics and conserved quantities

The vector fields $\partial_t$ and $\partial_\varphi$ are Killing vectors, giving two conserved quantities along any geodesic with tangent $u$:

$$
E=\left(1-\frac{2M}{r}\right)\frac{dt}{d\lambda},\qquad L=r^{2}\sin^{2}\theta\,\frac{d\varphi}{d\lambda}.
$$

The code verifies conservation from the definition. For a vector $\xi$ and geodesic tangent $u$,

$$
\frac{d}{d\lambda}\left(g_{\mu\nu}\xi^{\mu}u^{\nu}\right)=\frac{1}{2}\left(\mathcal{L}_{\xi}\,g\right)_{\mu\nu}u^{\mu}u^{\nu},
$$

which vanishes for every geodesic exactly when $\xi$ is Killing. The test also shows that a non-Killing field (for example $\partial_r$) is **not** conserved.

For motion in the equatorial plane, normalisation $g_{\mu\nu}u^{\mu}u^{\nu}=-\kappa$ (with $\kappa=1$ for massive particles and $\kappa=0$ for light) gives

$$
\left(\frac{dr}{d\lambda}\right)^{2}=E^{2}-V_{\mathrm{eff}}(r),
\qquad
V_{\mathrm{eff}}=\left(1-\frac{2M}{r}\right)\left(\kappa+\frac{L^{2}}{r^{2}}\right).
$$

### Circular orbits

- A circular orbit sits where $V_{\mathrm{eff}}'=0$. For massive particles this gives $L^{2}=Mr^{2}/(r-3M)$.
- It is stable if $V_{\mathrm{eff}}''>0$.
- The innermost stable circular orbit (ISCO) is where $V_{\mathrm{eff}}'=V_{\mathrm{eff}}''=0$. Solving symbolically gives $r=6M$ and $E=\sqrt{8/9}$. These are **derived** by the code, not typed in.
- For light ($\kappa=0$), the extremum of $\left(1-2M/r\right)/r^{2}$ is at $r=3M$, the photon sphere. The critical impact parameter is $b_{c}=3\sqrt{3}\,M$.

## 5. Observables

All of these use the orbit equation in the variable $u=1/r$,

$$
\left(\frac{du}{d\varphi}\right)^{2}=\frac{E^{2}-1}{L^{2}}+\frac{2M}{L^{2}}\,u-u^{2}+2M\,u^{3}.
$$

For light the left side equals $b^{-2}-u^{2}+2Mu^{3}$, where $b=L/E$ is the impact parameter.

### Periapsis advance

A bound orbit oscillates between two turning points $u_a$ (far) and $u_p$ (near). The cubic on the right has three roots; the third is $u_3=\frac{1}{2M}-u_a-u_p$. The substitution $u=\frac{u_a+u_p}{2}+\frac{u_p-u_a}{2}\sin\psi$ removes the endpoint singularities and gives

$$
\Delta\varphi=2\int_{-\pi/2}^{\pi/2}\frac{d\psi}{\sqrt{2M\left(u_3-u\right)}},
$$

so that the advance of the periapsis per orbit is $\Delta\varphi-2\pi$. For an orbit with semi-latus rectum $p\gg M$ it tends to

$$
\frac{6\pi M}{p}.
$$

### Deflection of light

$$
\alpha(b)=2\int_{0}^{u_{0}}\frac{du}{\sqrt{b^{-2}-u^{2}+2Mu^{3}}}-\pi
=\frac{4M}{b}+\frac{15\pi}{4}\,\frac{M^{2}}{b^{2}}+O\left(b^{-3}\right).
$$

Close to the critical value the deflection diverges logarithmically:

$$
\alpha\to-\ln\left(\frac{b}{b_{c}}-1\right)+\ln\left[216\left(7-4\sqrt{3}\right)\right]-\pi .
$$

### Gravitational redshift

For a static emitter at $r_e$ and a static receiver at $r_o$,

$$
\frac{\nu_{o}}{\nu_{e}}=\sqrt{\frac{1-2M/r_{e}}{1-2M/r_{o}}},
$$

which for weak fields is approximately $1-\Delta\Phi$.

### Shapiro time delay

For a signal passing at closest approach $r_0$, the extra round-trip time compared with flat space is approximately $4M\left[1+\ln\left(4r_{1}r_{2}/r_{0}^{2}\right)\right]$. The measured remainder scales as $M^{2}/r_{0}$, as expected for a next-order term.

### Tidal forces

In a static orthonormal frame, with $E_{ij}=R^{i}{}_{0j0}$, the tidal tensor is

$$
E_{ij}=\frac{M}{r^{3}}\,\mathrm{diag}\left(-2,\;1,\;1\right).
$$

It is trace-free, as required in vacuum. Neighbouring free-falling particles are stretched along the radial direction ($\ddot\xi=+2M\xi/r^{3}$) and squeezed across it.

## 6. Implementation

| Task | Module |
|---|---|
| Metric, Eddington-Finkelstein form, tortoise coordinate | `spacetimes/schwarzschild.py` |
| Derivation from the ansatz | `field_equations/vacuum_equations.py` |
| Effective potential, ISCO, photon sphere | `dynamics/effective_potential.py` |
| Geodesics derived from a Lagrangian | `dynamics/variational_principles.py` |
| Conservation laws | `dynamics/conserved_quantities.py` |
| Numerical geodesics | `dynamics/geodesics.py` |
| Deflection, precession, redshift, Shapiro delay | `dynamics/observables.py` |
| Tidal tensor | `dynamics/geodesic_deviation.py` |

## 7. How it is checked

The central idea is **two independent routes**.

- The advance and the deflection are each computed by direct quadrature of the exact orbit equation, and by integrating the geodesic equation built from this package's own Christoffel symbols. They agree to roughly $10^{-9}$ to $10^{-11}$ in relative terms (experiments 03 and 04).
- The geodesic equation is itself derived two ways: from the connection, and from the Euler-Lagrange equations of $\frac{1}{2}g_{\mu\nu}\dot x^{\mu}\dot x^{\nu}$. They agree exactly.
- The weak-field formulas show the expected size of their remainders: the residual after $4M/b$ falls as $b^{-2}$, and after the second-order term as $b^{-3}$. The periapsis advance differs from $6\pi M/p$ by a relative amount proportional to $1/p$.
- The strong-field limit is tested against the logarithmic formula and its constant.

## 8. Limits

Only test particles are treated: no spin, no back-reaction on the black hole, no radiation. Birkhoff's theorem is assumed.
