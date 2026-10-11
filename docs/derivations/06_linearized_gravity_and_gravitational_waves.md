# 6. Linearised gravity and gravitational waves

**In one paragraph.** When gravity is weak, spacetime is flat space plus a small ripple. The Einstein equations then become linear, and their solutions include waves that travel at the speed of light. This page derives the linearised curvature, explains why most components of the ripple are unphysical, and follows free particles through a passing wave. It ends with a second-order calculation that shows where the linear picture stops.

## 1. Setting

The metric is $g_{\mu\nu}=\eta_{\mu\nu}+h_{\mu\nu}$ with $\left|h_{\mu\nu}\right|\ll1$. Indices on $h$ are raised and lowered with the flat metric $\eta_{\mu\nu}$, and only terms linear in $h$ are kept.

## 2. Linearised curvature

Substituting into the Riemann formula and dropping products of $\Gamma$'s (which are second order), the Riemann tensor with all indices down is

$$
R_{abcd}=\frac{1}{2}\left(\partial_b\partial_c h_{ad}+\partial_a\partial_d h_{bc}-\partial_a\partial_c h_{bd}-\partial_b\partial_d h_{ac}\right).
$$

The Ricci tensor, scalar curvature and Einstein tensor follow by contraction. The code implements these formulas separately from the full curvature engine, and a test confirms that they match the part of the full engine's result that is first order in a small parameter $\epsilon$ for $g=\eta+\epsilon h$.

## 3. Gauge freedom

A small change of coordinates $x^{\mu}\to x^{\mu}+\xi^{\mu}$ changes the perturbation by

$$
h_{\mu\nu}\;\to\;h_{\mu\nu}+\partial_\mu\xi_\nu+\partial_\nu\xi_\mu .
$$

The linearised Riemann tensor does not change. The test confirms this symbolically with completely arbitrary functions for both $h_{\mu\nu}$ and $\xi_\mu$. So the ten components of $h_{\mu\nu}$ contain redundancy that is not physical.

**Counting the physical degrees of freedom.**

1. Define the trace-reversed perturbation $\bar h_{\mu\nu}=h_{\mu\nu}-\frac{1}{2}\eta_{\mu\nu}h$.
2. The Lorenz gauge condition $\partial^{\mu}\bar h_{\mu\nu}=0$ removes 4 components. In this gauge the field equation becomes $G^{(1)}_{\mu\nu}=-\frac{1}{2}\square\bar h_{\mu\nu}$ (verified by the tests).
3. A residual freedom with $\square\xi_\mu=0$ removes 4 more.
4. Ten minus eight leaves **two** physical polarisations.

## 4. Plane waves

In the transverse-traceless gauge, a wave moving along the $z$ axis has only the components

$$
h_{xx}=-h_{yy}=h_{+},\qquad h_{xy}=h_{\times},
$$

which depend on $\omega t-kz$. It solves the vacuum equation only if $k=\omega$, that is, if the wave moves at the speed of light. The test confirms that for $k\neq\omega$ at least one component of the linearised Einstein tensor is non-zero.

## 5. Effect on free particles

In this gauge, particles that start at rest stay at fixed coordinates, because $\Gamma^{i}{}_{00}=0$. What changes is the measured distance between them:

$$
L=\sqrt{g_{ij}\,\Delta^{i}\Delta^{j}} .
$$

A second description uses the geodesic-deviation equation for the separation vector $\xi^{i}$:

$$
\ddot\xi^{\,i}=\frac{1}{2}\,\ddot h^{\mathrm{TT}}_{ij}\,\xi^{\,j}.
$$

The two descriptions should agree to first order in the amplitude $A$ and differ at order $A^{2}$. The measured difference scales with amplitude with log-log slope $2.00$ (experiment 07). The distortion is also trace-free, so a ring of particles changes its shape but its area changes only at second order.

## 6. Beyond linear order

The module `perturbation_series.py` expands the connection, curvature, $\sqrt{-g}$ and Einstein tensor in powers of $\epsilon$, with the truncation order stated explicitly.

**The example.** Take the metric

$$
ds^{2}=-dt^{2}+dz^{2}+\left(1+\epsilon f\right)dx^{2}+\left(1-\epsilon f\right)dy^{2},\qquad f=\sin\left(t-z\right).
$$

Results:

- The first-order Einstein tensor vanishes, as it must for a vacuum wave.
- The second-order Einstein tensor does **not** vanish: for example $G^{(2)}_{tt}=\frac{1}{4}\left(3\cos\left(2t-2z\right)-1\right)$.
- Both statements agree with the exact Einstein tensor of this metric, Taylor-expanded in $\epsilon$.

So a gravitational wave acts as a source of curvature at second order. This is the origin of the effective stress-energy of gravitational waves.

**The error scaling.** The leftover error of the second-order truncation of $G_{tt}$ scales as $\epsilon^{4}$, not $\epsilon^{3}$. The reason: changing $\epsilon\to-\epsilon$ just swaps the $x$ and $y$ directions, so $G_{tt}$ is an even function of $\epsilon$ and has no odd terms. The measured log-log slope is $4.00$.

**Important.** This is a truncated series. It is not a nonlinear evolution of Einstein's equations.

## 7. Implementation and checks

| Task | Module | Tests |
|---|---|---|
| Linearised curvature, gauge, Lorenz gauge, plane waves | `field_equations/linearized_gravity.py` | `tests/test_linearized_gravity.py` |
| Truncated series | `field_equations/perturbation_series.py` | same file |
| Test particles in a wave | `dynamics/gravitational_waves.py` | `tests/test_quantum_and_numerics.py` |

## 8. Limits

The background is flat space only. There is no linearisation about a curved background, no radiation reaction, and no source multipole formulas.
