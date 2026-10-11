# 5. Newtonian limit

## Assumptions
Weak field $|\Phi|\ll1$, slow motion $v\ll1$, static source, $c=1$ then restored by dimensional analysis.

## Derivation
Write $g_{00}=-(1+2\Phi)$, $g_{ij}=(1-2\Phi)\delta_{ij}$. For a slow particle $u^\mu\approx(1,\mathbf v)$, so the dominant term is $\Gamma^i{}_{00}=-\tfrac12g^{ij}\partial_jg_{00}=\partial_i\Phi$ (static field, $O(\Phi)$), and the geodesic equation reduces to
$$\frac{d^2x^i}{dt^2}=-\partial_i\Phi.$$
The linearised Einstein tensor for this $h_{\mu\nu}$ has $G^{(1)}_{00}=2\nabla^2\Phi$ (verified symbolically with an arbitrary function $\Phi(x,y,z)$). With $T_{00}\approx\rho$ the 00 equation is $\nabla^2\Phi=4\pi G\rho$.

## Order of the approximation
Corrections are first order in $M/r$ and $v^2$. Exact Schwarzschild statements used in experiments:

* A particle momentarily at rest has $d^2r/dt^2=-Mf/r^2$ with $f=1-2M/r$, i.e. Newton's value with relative correction exactly $2M/r$.
* Circular orbits: coordinate period $2\pi r^{3/2}/\sqrt M$ (Kepler's third law holds exactly in Schwarzschild time); proper period $\times\sqrt{1-3M/r}$, so $\tau/\tau_N-1\approx-3M/2r$.
* Periapsis advance $\approx6\pi M/p$ per orbit with relative corrections $O(M/p)$.

## Where Newton fails
At the 1% level the period needs $r\gtrsim150M$ (experiment 06); the periapsis advance accumulates linearly with orbit number, so even at large $p$ the angular error grows secularly.

## Implementation / validation
`tests/test_quantum_and_numerics.py` (Poisson equation from $G^{(1)}_{00}$; coordinate acceleration $-Mf/r^2$), `tests/test_geodesics.py`, experiment 06.

## Limitations
Test-mass limit only; no post-Newtonian two-body dynamics.
