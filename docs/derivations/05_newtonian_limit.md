# 5. The Newtonian limit

**In one paragraph.** For weak fields and slow motion, general relativity must reduce to Newtonian gravity. This page shows how the Newtonian equations come out of the metric and the Einstein equation, and measures how big the relativistic corrections are.

## 1. Assumptions

The gravitational potential is weak, $|\Phi|\ll 1$. Matter moves slowly compared with light, $v\ll 1$. The source is static. Units are $c=1$.

## 2. The weak-field metric

Write

$$
g_{00}=-\left(1+2\Phi\right),\qquad g_{ij}=\left(1-2\Phi\right)\delta_{ij}.
$$

## 3. Equation of motion

For a slow particle $u^{\mu}\approx\left(1,\mathbf{v}\right)$, so the dominant term in the geodesic equation is $\Gamma^{i}{}_{00}$:

$$
\Gamma^{i}{}_{00}=-\frac{1}{2}\,g^{ij}\,\partial_j g_{00}=\partial_i\Phi .
$$

The geodesic equation then reduces to

$$
\frac{d^{2}x^{i}}{dt^{2}}=-\partial_i\Phi .
$$

## 4. Field equation

Treat the metric as flat space plus a small perturbation $h_{\mu\nu}$ with $h_{00}=h_{ii}=-2\Phi$. The linearised Einstein tensor has the time-time component

$$
G^{(1)}_{00}=2\,\nabla^{2}\Phi .
$$

The code checks this symbolically for an arbitrary function $\Phi(x,y,z)$. For dust, $T_{00}\approx\rho$, so $G_{00}=8\pi\,T_{00}$ becomes

$$
\nabla^{2}\Phi=4\pi G\rho,
$$

which is Poisson's equation.

## 5. How large are the corrections?

Three exact results for Schwarzschild show the size of the correction, which is of relative order $M/r$.

- **Static particle.** A particle momentarily at rest has coordinate acceleration $d^{2}r/dt^{2}=-Mf/r^{2}$ with $f=1-2M/r$. This is Newton's value with a relative correction of exactly $2M/r$.
- **Circular orbits.** The period in coordinate time is exactly $2\pi r^{3/2}/\sqrt{M}$, which is Kepler's third law. The period in the particle's own time is shorter by the factor $\sqrt{1-3M/r}$, a relative correction of about $-3M/(2r)$.
- **Precession.** The periapsis advances by about $6\pi M/p$ per orbit. Newtonian orbits are closed ellipses and do not precess.

## 6. Where Newton's description fails

At the 1% level, the proper period of a circular orbit needs $r\gtrsim 150M$ (experiment 06). Even at large radius the precession accumulates linearly with the number of orbits, so a long Newtonian prediction drifts away from the relativistic orbit.

## 7. How it is checked

- The Poisson-equation component test and the coordinate-acceleration test are in `tests/test_quantum_and_numerics.py`.
- The circular-orbit period is measured by integrating the geodesic through $\Delta\varphi=2\pi$ and compared with $2\pi r^{3/2}\sqrt{1-3M/r}$. They agree to about $10^{-14}$ (experiment 06).
- The log-log slope of the relative period correction against $r$ is $-1$.

## 8. Limits

Test-mass limit only. There is no post-Newtonian two-body dynamics.
