# 4. Charged and rotating black holes

**In one paragraph.** Reissner-Nordstrom (charged) and Kerr (rotating) black holes extend Schwarzschild. For both, the metric is **entered** into the code rather than derived. The pages below say exactly what is then **verified**, and how strong that verification is.

## 1. Reissner-Nordstrom

### What is entered

$$
ds^{2}=-f\,dt^{2}+f^{-1}dr^{2}+r^{2}\,d\Omega^{2},\qquad f=1-\frac{2M}{r}+\frac{Q^{2}}{r^{2}},
$$

with the electromagnetic potential $A_{\mu}=\left(-Q/r,\,0,\,0,\,0\right)$.

### What is verified (exactly)

Units are $G=c=1$ with Gaussian electromagnetic units. The code checks symbolically that

1. the Einstein equation $G_{\mu\nu}=8\pi\,T_{\mu\nu}$ holds, with $T_{\mu\nu}$ built from the field strength;
2. the source-free Maxwell equation holds, as shown below.

$$
\frac{1}{\sqrt{-g}}\,\partial_\mu\left(\sqrt{-g}\,F^{\mu\nu}\right)=0 .
$$

Both are exact statements. Two **negative controls** show that the checks can fail: doubling the $Q^{2}/r^{2}$ coefficient in $f$, or doubling the potential, makes the Einstein equation fail.

The horizons are the roots of $f=0$, namely $r_{\pm}=M\pm\sqrt{M^{2}-Q^{2}}$.

## 2. Kerr

### What is entered

The Boyer-Lindquist form, with

$$
\Sigma=r^{2}+a^{2}\cos^{2}\theta,\qquad \Delta=r^{2}-2Mr+a^{2},
$$

together with its closed-form inverse metric. The constructor multiplies the metric by the supplied inverse and refuses to continue unless the result is the identity.

### What is verified, and how strongly

| Check | Kind of evidence |
|---|---|
| Ricci tensor vanishes, at sample points ($a/M=0.6$) | numerical: residual below $10^{-25}$ at 30-digit precision |
| Kretschmann scalar equals an independently known closed form | numerical: relative agreement about $10^{-12}$ or better |
| $\partial_t$ and $\partial_\varphi$ are Killing | numerical, at sample points |
| Horizons at $\Delta=0$ | exact |
| The ergosurface is where $g_{tt}=0$, at $r=M+\sqrt{M^{2}-a^{2}\cos^{2}\theta}$ | exact / numerical |
| Frame dragging $\omega=-g_{t\varphi}/g_{\varphi\varphi}$ tends to $2J/r^{3}$ far away | numerical |

The key point: **no symbolic proof** that Kerr solves the vacuum equations is attempted. The evidence is numerical.

### Equatorial circular orbits and the ISCO

The formulas for the energy, angular momentum and angular velocity of circular equatorial orbits, and for the ISCO radius, are taken from the literature (Bardeen, Press and Teukolsky, 1972). They are treated as **inputs** and tested three ways. First, the radial acceleration $\Gamma^{r}{}_{\alpha\beta}u^{\alpha}u^{\beta}$, computed from this package's own Christoffel symbols, vanishes to round-off for the literature orbit.

Second, the energy and angular momentum recomputed from the metric agree with the literature formulas.

Third, define the radial function

$$
W(r;E,L)=-1-g^{tt}E^{2}+2g^{t\varphi}EL-g^{\varphi\varphi}L^{2}
$$

from the engine's inverse metric. It satisfies $W=W'=0$ on a circular orbit, and $W''=0$ at the quoted ISCO radius (experiment 09).

The ISCO lies at $6M$ for $a=0$. As $a\to M$ it moves in to $M$ for prograde orbits and out to $9M$ for retrograde ones. A further test perturbs a prograde circular orbit radially: outside the ISCO it stays bounded, inside it plunges.

## 3. Limits

No symbolic Kerr vacuum proof, no off-equatorial orbit analysis, no Kerr-Newman solution, and no Carter-frame tetrad calculation.
