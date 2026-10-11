# 9. A scalar field on a curved background

**In one paragraph.** A wave equation can be solved on a fixed curved spacetime. Here it is solved on the Schwarzschild background, where it reduces to a one-dimensional scattering problem. The result is a set of transmission coefficients that tell how much of an incoming wave gets through the black hole's potential barrier.

## 1. What this is and is not

It is easy to overclaim in this area, so the distinctions are stated first.

| Setting | In this repository? |
|---|---|
| A classical field on a curved spacetime | **Yes.** This page. |
| A quantum field on a fixed classical background | Described only. No mode sums, no Bogoliubov coefficients, no renormalisation. |
| Gravity itself quantised perturbatively | Not implemented. See [page 10](10_perturbative_gravity_and_eft.md). |
| A complete theory of quantum gravity | Does not exist here and is not claimed. |

The scattering numbers below are classical transmission coefficients. They do not by themselves say anything about Hawking radiation.

## 2. The equation

In natural units a scalar field of mass $m$ coupled to the scalar curvature with strength $\xi$ obeys

$$
\left(\square-m^{2}-\xi R\right)\phi=0,
\qquad
\square\phi=\frac{1}{\sqrt{-g}}\,\partial_\mu\left(\sqrt{-g}\,g^{\mu\nu}\partial_\nu\phi\right).
$$

On Schwarzschild the Ricci scalar is $R=0$, so the coupling $\xi$ has no effect.

## 3. Reduction to a radial problem

Separate the field as

$$
\phi=e^{-i\omega t}\,P(\theta)\,\frac{\psi(r)}{r},
\qquad
P''+\cot\theta\,P'=-\ell(\ell+1)\,P .
$$

Introduce the tortoise coordinate $r_{*}$ by $dr_{*}=dr/f$ with $f=1-2M/r$. The wave equation becomes

$$
\frac{d^{2}\psi}{dr_{*}^{2}}+\left(\omega^{2}-V\right)\psi=0,
\qquad
V=f\left(\frac{\ell(\ell+1)}{r^{2}}+\frac{2M}{r^{3}}+m^{2}\right).
$$

The code verifies this reduction symbolically: it applies the full wave operator to the separated ansatz and shows that the difference from the claimed radial equation is exactly zero.

## 4. The scattering problem

The potential $V$ is a barrier that vanishes at both ends ($r_{*}\to\pm\infty$), so the solution behaves as plane waves there.

- **At the horizon.** Nothing can come back out, so the wave is purely ingoing: $\psi=e^{-i\omega r_{*}}$ with unit amplitude. The integration starts at $r_{*}=-80M$, where the neglected correction is of order $e^{-40}$.
- **At large radius.** The solution is written as an incoming plus an outgoing wave, $\psi=A_{\mathrm{in}}\,u^{-}+A_{\mathrm{out}}\,u^{+}$. The functions $u^{\pm}$ are the Riccati-Hankel solutions of the flat-space centrifugal problem, which are exact solutions when the barrier is only $\ell(\ell+1)/r^{2}$. Using them in place of plain exponentials removes an error that would otherwise be of order $\ell(\ell+1)/\left(\omega^{2}r_{*}^{2}\right)$.

The transmission and reflection coefficients are

$$
T=\frac{1}{\left|A_{\mathrm{in}}\right|^{2}},
\qquad
R=\frac{\left|A_{\mathrm{out}}\right|^{2}}{\left|A_{\mathrm{in}}\right|^{2}} .
$$

Because the Wronskian of $u^{+}$ and $u^{-}$ is constant, energy conservation requires $T+R=1$. The code reports $1-T-R$ as a check. (The more obvious form $\left|A_{\mathrm{in}}\right|^{2}-\left|A_{\mathrm{out}}\right|^{2}-1$ is avoided because it loses all precision through rounding when $T$ is very small.)

## 5. Analytic checks

1. **Low frequency, $\ell=0$.** The absorption cross section $\sigma=\pi T/\omega^{2}$ tends to the area of the horizon, $16\pi M^{2}$. Equivalently, $T/\left(16M^{2}\omega^{2}\right)\to1$ as $\omega\to0$.
2. **Low frequency, general $\ell$.** The transmission scales as $\omega^{2\ell+2}$. The test confirms a slope near 4 for $\ell=1$.
3. **High frequency.** The transmission tends to 1.
4. **Barrier ordering.** At fixed $\omega$, larger $\ell$ means a higher barrier and smaller transmission.

## 6. Results

- The flux residual $1-T-R$ stays below about $5\times10^{-9}$ over the whole grid of experiment 10 ($\ell\le3$, $0.05\le M\omega\le1$).
- Changing the outer boundary radius from $3000M$ to $6000M$ changes a transmission coefficient by about $10^{-11}$.
- $T_0/\left(16M^{2}\omega^{2}\right)$ equals $1.034$ at $M\omega=0.005$ and $1.071$ at $0.01$, approaching 1 from above as $\omega\to0$.

## 7. Implementation

| Task | Module |
|---|---|
| Wave operator, symbolic reduction, scattering solver | `quantum/klein_gordon.py` |
| Horizon and asymptotic boundary data | `numerical/boundary_conditions.py` |
| Tests | `tests/test_quantum_and_numerics.py` |

## 8. Limits

There is no claim about Hawking emission spectra. Only massless fields on Schwarzschild are run. The matching at a finite outer radius ignores corrections of order $M\ln r_{*}/r_{*}^{3}$.
