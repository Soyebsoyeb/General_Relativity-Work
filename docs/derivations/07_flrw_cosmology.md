# 7. FLRW cosmology

## Assumptions
Homogeneous, isotropic: $ds^2=-dt^2+a(t)^2\bigl[\frac{dr^2}{1-kr^2}+r^2d\Omega^2\bigr]$, perfect fluid $T_{\mu\nu}=(\rho+p)u_\mu u_\nu+pg_{\mu\nu}$ with $u=\partial_t$.

## Derivation
The engine gives (checked symbolically, including isotropy $G^r{}_r=G^\theta{}_\theta=G^\phi{}_\phi$)
$$G^t{}_t=-\frac{3(\dot a^2+k)}{a^2},\qquad G^r{}_r=-\frac{2a\ddot a+\dot a^2+k}{a^2}.$$
$G+\Lambda g=8\pi T$ with $T^t{}_t=-\rho$, $T^r{}_r=p$ gives the Friedmann equations
$$\Bigl(\frac{\dot a}a\Bigr)^2+\frac k{a^2}=\frac{8\pi}3\rho+\frac\Lambda3,\qquad \frac{\ddot a}a=-\frac{4\pi}3(\rho+3p)+\frac\Lambda3.$$
The continuity equation $\dot\rho+3H(\rho+p)=0$ was not imposed: it follows from the contracted Bianchi identity (tested by substituting the $\rho,p$ solved from the two equations).

## Numerics
With $H_0=1$, $a(t_0)=1$: $E^2(a)=\Omega_ra^{-4}+\Omega_ma^{-3}+\Omega_ka^{-2}+\Omega_\Lambda$, $\Omega_k=1-\Omega_r-\Omega_m-\Omega_\Lambda$. The second-order form is integrated (it handles recollapse) and the first-order constraint $(\dot a/a)^2-E^2(a)$ is monitored.

Analytic references: matter $a=(3t/2)^{2/3}$; radiation $a=(2t)^{1/2}$; de Sitter $a=e^t$; flat $\Lambda$CDM $a=(\Omega_m/\Omega_\Lambda)^{1/3}\sinh^{2/3}(\tfrac32\sqrt{\Omega_\Lambda}t)$; closed dust $a=\tfrac{a_{\max}}2(1-\cos\eta)$, $H_0t=\frac{\Omega_m}{2|\Omega_k|^{3/2}}(\eta-\sin\eta)$.

## Validation
Experiment 08 / `tests/test_flrw.py`: relative error against analytic solutions at integrator tolerance for all five cases; a negative control (integrating matter dynamics but comparing with the radiation solution) fails by more than 1%. Redshift is $1+z=1/a$.

## Limitations
Background only. The matching of initial data uses the analytic solution, so these runs test the integrator and the ODE forms, while the derivation of those ODEs from $G=8\pi T$ is tested symbolically. No perturbations, no Boltzmann physics.
