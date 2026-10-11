# 9. Scalar fields on a fixed curved background

## Classification (important)
| Setting | In this repository? |
|---|---|
| Classical field on curved spacetime | **yes** (the scattering numbers below) |
| Quantum field on a fixed classical background | stated only; no mode sums, no Bogoliubov coefficients, no renormalisation |
| Perturbatively quantised gravity | not implemented (see [page 10](10_perturbative_gravity_and_eft.md)) |
| A complete quantum theory of gravity | does not exist here and is not claimed |

## Equation
$(\Box-m^2-\xi R)\phi=0$, $\Box\phi=\frac1{\sqrt{-g}}\partial_\mu(\sqrt{-g}g^{\mu\nu}\partial_\nu\phi)$, natural units. On Schwarzschild $R=0$ so $\xi$ drops out.

## Reduction (verified symbolically)
Let $\phi=e^{-i\omega t}P(\theta)\psi(r)/r$ with $P''+\cot\theta P'=-\ell(\ell+1)P$. Then $\sqrt{-g}\,\Box\phi$ separates and, with $dr_*=dr/f$,
$$\frac{d^2\psi}{dr_*^2}+(\omega^2-V)\psi=0,\qquad V=f\Bigl(\frac{\ell(\ell+1)}{r^2}+\frac{2M}{r^3}+m^2\Bigr).$$

## Scattering problem
* Horizon ($r_*\to-\infty$): $\psi=e^{-i\omega r_*}$ (purely ingoing, unit amplitude). Starting at $r_*=-80M$ the neglected corrections are $O(f)\sim e^{-40}$.
* Infinity: $\psi=A_{\rm in}u^-+A_{\rm out}u^+$ with $u^\pm$ the Riccati–Hankel solutions of the flat centrifugal problem, which removes the $O(\ell(\ell+1)/\omega^2r_*^2)$ error of plane-wave matching.
* Transmission $T=1/|A_{\rm in}|^2$, reflection $R=|A_{\rm out}|^2/|A_{\rm in}|^2$. The Wronskian of $u^\pm$ is constant, so the numerical check is $1-T-R\approx0$ (a relative form: the absolute form $|A_{\rm in}|^2-|A_{\rm out}|^2-1$ is dominated by round-off when $T\ll1$).

## Analytic checks
* $\omega\to0$, $\ell=0$: $\sigma_{\rm abs}=\pi T/\omega^2\to16\pi M^2$ (horizon area), i.e. $T/(16M^2\omega^2)\to1$.
* $\omega\to0$: $T_\ell\propto\omega^{2\ell+2}$ (slope 4 tested for $\ell=1$).
* $\omega\to\infty$: $T\to1$.

## Validation
`tests/test_quantum_and_numerics.py`: symbolic reduction residual zero; amplitude extraction on exact flat solutions; flux conservation; insensitivity to $r_*^{\max}$; the three limits above. Experiment 10 records the measured values.

## Limitations
No claim about Hawking emission spectra. The finite-$r_*^{\max}$ matching ignores $M\ln r_*/r_*^3$ corrections. Only massless $\xi$-independent Schwarzschild examples are run.
