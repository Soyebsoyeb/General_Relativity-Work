# 6. Linearised gravity, gravitational waves and higher order

## Assumptions
$g_{\mu\nu}=\eta_{\mu\nu}+h_{\mu\nu}$, $|h|\ll1$; indices moved with $\eta$.

## Linearised curvature
To first order, in the package convention
$$R_{abcd}=\tfrac12\bigl(\partial_b\partial_ch_{ad}+\partial_a\partial_dh_{bc}-\partial_a\partial_ch_{bd}-\partial_b\partial_dh_{ac}\bigr),$$
and $R_{ab}$, $R$, $G_{ab}$ follow by contraction. These are implemented independently of the full engine and tested against the $O(\epsilon)$ coefficient of the full engine for $g=\eta+\epsilon h$.

## Gauge
An infinitesimal diffeomorphism $x\to x+\xi$ changes $h_{\mu\nu}\to h_{\mu\nu}+\partial_\mu\xi_\nu+\partial_\nu\xi_\mu$ and leaves $R^{(1)}_{abcd}$ invariant (verified with arbitrary functions $h_{\mu\nu}$, $\xi_\mu$). So $h$ has 10 components but not 10 physical degrees of freedom. Counting: the Lorenz gauge $\partial^\mu\bar h_{\mu\nu}=0$ (with $\bar h_{\mu\nu}=h_{\mu\nu}-\frac12\eta_{\mu\nu}h$) removes 4 components and leaves a residual gauge freedom with $\Box\xi_\mu=0$ that removes 4 more, leaving **2** physical polarisations. In Lorenz gauge $G^{(1)}_{\mu\nu}=-\tfrac12\Box\bar h_{\mu\nu}$ (tested), so vacuum waves satisfy $\Box\bar h=0$.

## Plane waves
In TT gauge ($h_{0\mu}=0$, $h^i{}_i=0$, $\partial^jh_{ij}=0$) a wave along $z$ has $h_{xx}=-h_{yy}=h_+$, $h_{xy}=h_\times$, functions of $\omega t-kz$. The vacuum equation is satisfied iff $k=\omega$ (tested: with $k\ne\omega$ some $G^{(1)}$ component is non-zero).

## Effect on free particles
Particles at rest in TT coordinates stay at rest ($\Gamma^i{}_{00}=0$), and their proper separation is $L=\sqrt{g_{ij}\Delta^i\Delta^j}$. In the proper detector frame the geodesic-deviation equation is $\ddot\xi^i=\tfrac12\ddot h^{TT}_{ij}\xi^j$. The two descriptions agree to $O(A^2)$ (measured slope 2.00, experiment 07). The distortion is trace-free: ring area changes only at $O(h^2)$.

## Beyond linear order
`perturbation_series.perturbative_geometry` expands $\Gamma$, $R_{\mu\nu}$, $R$, $\sqrt{-g}$ and $G_{\mu\nu}$ for $g=\eta+\epsilon h$ with explicit truncation order. For the plane wave $f=A\sin\omega(t-z)$ in $ds^2=-dt^2+dz^2+(1+\epsilon f)dx^2+(1-\epsilon f)dy^2$: $G^{(1)}=0$ but $G^{(2)}_{tt}=\tfrac14(3\cos2(t-z)-1)$ for $A=\omega=1$, matching the Taylor coefficients of the exact Einstein tensor. The wave therefore sources curvature at second order (the origin of the Isaacson effective stress-energy). Because $\epsilon\to-\epsilon$ is the relabelling $x\leftrightarrow y$, $G_{tt}$ is even in $\epsilon$, and the truncation error of the $O(\epsilon^2)$ result scales as $\epsilon^4$ (measured slope 4.0).

**This is a truncated series, not a nonlinear evolution.**

## Implementation / validation
`field_equations/linearized_gravity.py`, `perturbation_series.py`, `dynamics/gravitational_waves.py`; `tests/test_linearized_gravity.py`; experiments 07 and 11.

## Limitations
Flat background only; no linearisation about curved backgrounds; no radiation reaction, no source multipole formulas.
