# 4. Reissner–Nordström and Kerr

For every extension the distinction is explicit: a metric that is **entered** versus a solution that is **verified**.

## Reissner–Nordström (entered, then verified exactly)
Entered: $f=1-2M/r+Q^2/r^2$, $ds^2=-f\,dt^2+f^{-1}dr^2+r^2d\Omega^2$, and $A_\mu=(-Q/r,0,0,0)$. Gaussian units with $G=c=1$:
$$G_{\mu\nu}=8\pi T_{\mu\nu},\quad T_{\mu\nu}=\tfrac1{4\pi}\bigl(F_{\mu\alpha}F_\nu{}^\alpha-\tfrac14g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\bigr),\quad \tfrac1{\sqrt{-g}}\partial_\mu(\sqrt{-g}F^{\mu\nu})=0.$$
All three are verified symbolically (exact). Negative controls: doubling the $Q^2/r^2$ coefficient, or doubling the potential, makes the Einstein equation fail. Horizons $r_\pm=M\pm\sqrt{M^2-Q^2}$ solve $f=0$.

## Kerr (entered; vacuum character verified numerically, not symbolically)
Entered in Boyer–Lindquist form with $\Sigma=r^2+a^2\cos^2\theta$, $\Delta=r^2-2Mr+a^2$, together with its closed-form inverse (checked against $g$ by the constructor).

| Check | Status |
|---|---|
| $R_{\mu\nu}=0$ at sample points ($a/M=0.6$) | numeric, residual below $10^{-25}$ at 30-digit evaluation |
| $K=R_{abcd}R^{abcd}$ against the independently known closed form | numeric, relative agreement about $10^{-12}$ or better |
| $\partial_t,\partial_\phi$ Killing | numeric at sample points |
| horizons $\Delta=0$; ergosurface $g_{tt}=0$ at $r=M+\sqrt{M^2-a^2\cos^2\theta}$ | exact / numeric |
| frame dragging $\omega=-g_{t\phi}/g_{\phi\phi}\to2J/r^3$ | numeric |

**Equatorial circular orbits.** The literature formulas for $E,L,\Omega$ and $r_{\rm ISCO}$ (Bardeen–Press–Teukolsky 1972) are *inputs*. They are validated by (i) the radial acceleration $\Gamma^r{}_{\alpha\beta}u^\alpha u^\beta$ from this package's Christoffel symbols vanishing to round-off, (ii) $E,L$ recomputed from the metric agreeing, and (iii) the radial function $W(r;E,L)=-1-g^{tt}E^2+2g^{t\phi}EL-g^{\phi\phi}L^2$ satisfying $W=W'=0$ with $W''=0$ at the quoted ISCO (experiment 09). ISCO values: $6M$ ($a=0$), $\to M$ prograde and $9M$ retrograde as $a\to M$.

## Limitations
No symbolic Kerr vacuum proof; no off-equatorial orbit analysis; no Kerr–Newman; no tetrad (Carter frame) calculation.
