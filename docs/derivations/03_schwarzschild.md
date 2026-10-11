# 3. Schwarzschild spacetime

## Assumptions
Static, spherically symmetric, asymptotically flat vacuum. Birkhoff's theorem (staticity follows from spherical symmetry) is **not** proved here; staticity is assumed.

## Derivation from the ansatz
Take $ds^2=-A(r)dt^2+B(r)dr^2+r^2d\Omega^2$. The package computes
$$G^t{}_t=-\frac1{r^2}\frac{d}{dr}\Bigl[r\Bigl(1-\frac1B\Bigr)\Bigr],\qquad G^r{}_r-G^t{}_t=\frac1{rB}\Bigl(\frac{A'}A+\frac{B'}B\Bigr)$$
(both forms are checked symbolically in `derive_static_spherical_vacuum`). $G^t{}_t=0$ integrates to $1/B=1+C/r$. The second relation gives $AB=\text{const}$, set to 1 by rescaling $t$. Matching $g_{tt}\to-(1+2\Phi)$ with $\Phi=-M/r$ fixes $C=-2M$:
$$ds^2=-\Bigl(1-\frac{2M}{r}\Bigr)dt^2+\Bigl(1-\frac{2M}{r}\Bigr)^{-1}dr^2+r^2d\Omega^2.$$
All components of $G_{\mu\nu}$ then vanish while the Riemann tensor does not: $K=R_{abcd}R^{abcd}=48M^2/r^6$ (symbolic contraction in the tests).

## Horizon and singularity
$g_{tt}=0$ and $g_{rr}\to\infty$ at $r=2M$ are coordinate effects. With $v=t+r_*$, $r_*=r+2M\ln(r/2M-1)$,
$$ds^2=-\Bigl(1-\frac{2M}{r}\Bigr)dv^2+2\,dv\,dr+r^2d\Omega^2,$$
whose determinant $-r^4\sin^2\theta$ is regular at $r=2M$ (obtained with `Metric.transform`, tested to be Ricci flat with the same $K$). $K\to\infty$ as $r\to0$ cannot be removed by a coordinate change, so this is a curvature singularity. A single scalar is not a universal singularity detector (pp-waves have all polynomial invariants equal to zero).

## Geodesics and constants of motion
$\partial_t$ and $\partial_\phi$ are Killing, giving $E=(1-2M/r)\dot t$ and $L=r^2\sin^2\theta\,\dot\phi$. The conservation law is verified from the definition, $\tfrac{d}{d\lambda}(g_{\mu\nu}\xi^\mu u^\nu)=\tfrac12(\mathcal L_\xi g)_{\mu\nu}u^\mu u^\nu$, and a non-Killing field is shown not to be conserved. For equatorial motion, $g(u,u)=-\kappa$ ($\kappa=1$ timelike, $0$ null) yields
$$\dot r^2=E^2-V_{\rm eff},\qquad V_{\rm eff}=\Bigl(1-\frac{2M}{r}\Bigr)\Bigl(\kappa+\frac{L^2}{r^2}\Bigr).$$

* Circular orbits: $V'=0\Rightarrow L^2=Mr^2/(r-3M)$; stable iff $V''>0$.
* ISCO: $V'=V''=0\Rightarrow r=6M$, $E=\sqrt{8/9}$ (derived symbolically, not input).
* Photon sphere ($\kappa=0$): $\frac{d}{dr}[(1-2M/r)/r^2]=0\Rightarrow r=3M$; critical impact parameter $b_c=3\sqrt3M$.

## Observables
With $u=1/r$: $(du/d\phi)^2=\frac{E^2-1}{L^2}+\frac{2M}{L^2}u-u^2+2Mu^3$.

* **Periapsis advance.** With turning points $u_a,u_p$ and third root $u_3=\frac1{2M}-u_a-u_p$, the substitution $u=\frac{u_a+u_p}2+\frac{u_p-u_a}2\sin\psi$ gives $\Delta\phi=2\int_{-\pi/2}^{\pi/2}d\psi/\sqrt{2M(u_3-u)}$ and advance $\Delta\phi-2\pi\to6\pi M/p$ for $p\gg M$.
* **Deflection.** $\alpha=2\int_0^{u_0}du/\sqrt{b^{-2}-u^2+2Mu^3}-\pi=\frac{4M}b+\frac{15\pi}4\frac{M^2}{b^2}+O(b^{-3})$. Near $b_c$: $\alpha\to-\ln(b/b_c-1)+\ln[216(7-4\sqrt3)]-\pi$.
* **Redshift.** Static emitter and receiver: $\nu_o/\nu_e=\sqrt{f(r_e)/f(r_o)}\approx1-\Delta\Phi$.
* **Shapiro delay.** Round-trip excess $\approx4M[1+\ln(4r_1r_2/r_0^2)]$ with a measured remainder $\propto M^2/r_0$.
* **Tides.** In the static orthonormal frame $E_{ij}=R^i{}_{0j0}=\mathrm{diag}(-2M,M,M)/r^3$: trace-free, and $\ddot\xi^{\hat r}=+2M\xi^{\hat r}/r^3$.

## Implementation
`spacetimes/schwarzschild.py`, `dynamics/effective_potential.py`, `conserved_quantities.py`, `variational_principles.py`, `observables.py`, `geodesic_deviation.py`.

## Validation
Each observable is computed twice, by exact quadrature and by integrating the geodesic equation built from the package's own Christoffel symbols. They agree to roughly $10^{-9}$–$10^{-11}$ (experiments 03–04, `tests/test_geodesics.py`). Weak-field formulas show their expected $O(M/b)$ and $O(M/p)$ remainders.

## Limitations
Test-particle motion only: no spin, back-reaction or radiation.
