import nbformat as nbf
from _nb import code, md, notebook

NB = {}

NB["01_metric_and_coordinate_transformations"] = [
    md(r"""# 1. Metric construction and coordinate transformations
A metric $g_{\mu\nu}$ assigns lengths to tangent vectors. Its components depend on the coordinates; the *geometry* does not. This notebook builds the flat plane in Cartesian and polar coordinates and checks that the tensor transformation law connects them.

For $x^\mu=x^\mu(x'^a)$, a tensor with one upper and one lower index transforms as
$$T'^a{}_b=\frac{\partial x'^a}{\partial x^\mu}\frac{\partial x^\nu}{\partial x'^b}T^\mu{}_\nu .$$"""),
    code("""
import sympy as sp
from gr_from_scratch.core.metric import Metric
from gr_from_scratch.core.tensor import Tensor
from gr_from_scratch.core.coordinates import transform_tensor
from gr_from_scratch.spacetimes import minkowski as mk

x, y = sp.symbols('x y', real=True)
cart = Metric([x, y], sp.eye(2), name='plane-cartesian')
polar_from_transform = cart.transform([mk.r, mk.phi], mk.polar_to_cartesian())
print(polar_from_transform.g)
print('same as the hand-written polar metric:', polar_from_transform.g == mk.euclidean_polar().g)
"""),
    md(
        r"""The rotation generator $\xi=-y\,\partial_x+x\,\partial_y$ should become $\partial_\phi$. A vector transforms with $\partial x'/\partial x$ (the inverse Jacobian), a covector with $\partial x/\partial x'$."""
    ),
    code("""
polar = mk.euclidean_polar()
xi = Tensor(cart, [-y, x], 'u')
print('xi in polar components:', [sp.simplify(c) for c in transform_tensor(xi, polar, mk.polar_to_cartesian()).c])

# lowering and raising with the polar metric
g = Tensor.metric_tensor(polar)
print(g.raise_index(0).raise_index(1).c)            # inverse metric
print('signature (neg, pos) at r=1:', polar.signature_at({mk.r: 1}))
"""),
    md(
        r"""**Spacetime check.** The Lorentzian metric $-dt^2+dr^2+r^2d\Omega^2$ has signature $(1,3)$ and determinant $-r^4\sin^2\theta$."""
    ),
    code("""
ms = mk.minkowski_spherical()
print(ms.det, ms.signature_at({mk.r: 2, mk.theta: 1}), ms.is_lorentzian_at({mk.r: 2, mk.theta: 1}))
"""),
]

NB["02_christoffel_from_metric_compatibility"] = [
    md(r"""# 2. Christoffel symbols from metric compatibility
Requiring $\nabla_\lambda g_{\mu\nu}=0$ and zero torsion fixes
$$\Gamma^\rho{}_{\mu\nu}=\tfrac12g^{\rho\sigma}(\partial_\mu g_{\nu\sigma}+\partial_\nu g_{\mu\sigma}-\partial_\sigma g_{\mu\nu}).$$
We compute it for progressively harder metrics, verify both defining properties, and show that $\Gamma$ is **not** a tensor."""),
    code("""
import numpy as np, sympy as sp
from gr_from_scratch.geometry.connection import christoffel, torsion, metric_compatibility_residual
from gr_from_scratch.core.simplification import simplify
from gr_from_scratch.spacetimes import minkowski as mk
from gr_from_scratch.spacetimes.schwarzschild import schwarzschild

def zero(a): return all(simplify(e) == 0 for e in np.asarray(a, dtype=object).ravel())

for name, m in [('Minkowski (Cartesian)', mk.minkowski()), ('plane, polar', mk.euclidean_polar()),
                ('3-space, spherical', mk.euclidean_spherical()), ('Schwarzschild', schwarzschild())]:
    c = christoffel(m)
    print(f'{name:24s} nonzero symbols: {len(c.nonzero()):2d}   torsion=0: {zero(torsion(c))}   nabla g = 0: {zero(metric_compatibility_residual(c))}')
"""),
    md(
        r"""**Non-tensorial character.** In Cartesian coordinates the plane has $\Gamma=0$. Transforming with the *inhomogeneous* law reproduces the polar symbols; the homogeneous tensor law would map zero to zero."""
    ),
    code("""
from gr_from_scratch.core.metric import Metric
from gr_from_scratch.core.coordinates import transform_connection
x, y = sp.symbols('x y', real=True)
cart = Metric([x, y], sp.eye(2))
polar = mk.euclidean_polar()
via = transform_connection(christoffel(cart).gamma, cart.coords, polar.coords, mk.polar_to_cartesian(), simplify)
print('inhomogeneous law reproduces polar symbols:', zero(via - christoffel(polar).gamma))
print(christoffel(polar).nonzero())
"""),
]

NB["03_riemann_curvature_and_flat_space_validation"] = [
    md(r"""# 3. Riemann curvature and flat-space validation
With $R^\rho{}_{\sigma\mu\nu}=\partial_\mu\Gamma^\rho{}_{\nu\sigma}-\partial_\nu\Gamma^\rho{}_{\mu\sigma}+\Gamma^\rho{}_{\mu\lambda}\Gamma^\lambda{}_{\nu\sigma}-\Gamma^\rho{}_{\nu\lambda}\Gamma^\lambda{}_{\mu\sigma}$, flat space must give zero in **every** coordinate system even though $\Gamma\neq0$ in curvilinear ones. A sphere must give $R=2$ for unit radius. Failing either would reveal a sign error."""),
    code("""
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry
from gr_from_scratch.spacetimes import minkowski as mk
import sympy as sp

for name, m in [('Cartesian Minkowski', mk.minkowski()), ('polar plane', mk.euclidean_polar()),
                ('spherical 3-space', mk.euclidean_spherical()), ('Minkowski in spherical', mk.minkowski_spherical())]:
    g = compute_geometry(m)
    print(f'{name:24s} Riemann = 0: {g.riemann.is_zero()}')

s = compute_geometry(mk.unit_sphere())
print('unit sphere  R =', s.scalar, '   R^theta_{phi theta phi} =', sp.simplify(s.riemann.c[0,1,0,1]))
"""),
    md(
        r"""**Symmetries and Bianchi identities** on a curved example (the unit sphere has the full tensor structure in 2D; Schwarzschild is the 4D case)."""
    ),
    code("""
import numpy as np
from gr_from_scratch.geometry import bianchi
from gr_from_scratch.core.simplification import simplify
from gr_from_scratch.spacetimes.schwarzschild import schwarzschild
g = compute_geometry(schwarzschild())
z = lambda a: all(simplify(e) == 0 for e in np.asarray(a, dtype=object).ravel())
for k, v in bianchi.pair_and_antisymmetry_residuals(g.riemann).items(): print(k, z(v))
print('algebraic Bianchi:', z(bianchi.algebraic_bianchi_residual(g.riemann)))
print('differential Bianchi:', z(bianchi.differential_bianchi_residual(g.riemann, g.conn)))
"""),
]

NB["04_schwarzschild_vacuum_verification"] = [
    md(r"""# 4. Schwarzschild vacuum verification
**Derivation.** For $ds^2=-A(r)dt^2+B(r)dr^2+r^2d\Omega^2$ the code finds $G^t{}_t=-\frac1{r^2}\frac{d}{dr}[r(1-1/B)]$. Setting it to zero gives $1/B=1+C/r$; the $rr$ equation gives $AB=$ const. Newtonian matching gives $C=-2M$.

The ansatz assumes staticity and spherical symmetry (Birkhoff's theorem is not proved here)."""),
    code("""
import sympy as sp
from gr_from_scratch.field_equations.vacuum_equations import derive_static_spherical_vacuum, compute_geometry
d = derive_static_spherical_vacuum()
print('G^t_t =', d['G^t_t'])
print('form check (should be 0):', d['check_Gtt_form'], d['check_Grr_minus_Gtt_form'])
print(d['u=1/B solution'])
"""),
    md(
        r"""**Verification.** The entered Schwarzschild metric satisfies all ten vacuum equations, yet its Riemann tensor is non-zero: vacuum is not flat."""
    ),
    code("""
from gr_from_scratch.spacetimes.schwarzschild import schwarzschild
g = compute_geometry(schwarzschild())
print('Einstein tensor zero:', g.einstein.is_zero())
print('Ricci scalar:', g.scalar)
print('Riemann zero:', g.riemann.is_zero())
for idx, val in list(g.riemann.lower_index(0).nonzero_components().items())[:4]:
    print(idx, val)
"""),
]

NB["05_curvature_invariants_and_singularities"] = [
    md(r"""# 5. Curvature invariants and singularities
$K=R_{abcd}R^{abcd}$ is a coordinate-independent scalar. For Schwarzschild $K=48M^2/r^6$: finite at $r=2M$ (the horizon is only a coordinate singularity of these coordinates) and divergent at $r=0$ (a true curvature singularity). Horizon-regular coordinates make the first point explicit.

Caveat: a single invariant is not a universal singularity detector (e.g. pp-waves have all polynomial invariants zero)."""),
    code("""
import sympy as sp, numpy as np, matplotlib.pyplot as plt
from gr_from_scratch.spacetimes.schwarzschild import schwarzschild, eddington_finkelstein, M, r, theta
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry
from gr_from_scratch.geometry import invariants as inv

gs = compute_geometry(schwarzschild())
ge = compute_geometry(eddington_finkelstein())
print('K (Schwarzschild coords):', inv.kretschmann(gs.riemann))
print('K (Eddington-Finkelstein):', inv.kretschmann(ge.riemann))
print('EF metric:', ge.metric.g)
print('EF determinant (regular at r=2M):', ge.metric.det)
print('Euler density == K in vacuum:', sp.simplify(inv.euler_density(gs.riemann) - inv.kretschmann(gs.riemann)) == 0)
print('Pontryagin density:', inv.pontryagin_density(gs.riemann))
"""),
    code("""
rr = np.logspace(np.log10(0.3), 1.5, 200)
fig, ax = plt.subplots(figsize=(4.2, 3))
ax.loglog(rr, 48 / rr**6)
ax.axvline(2, color='0.5', lw=0.7)
ax.set(xlabel='r / M', ylabel=r'$K M^4$', title='Kretschmann scalar (horizon at r = 2M)')
plt.show()
"""),
]

NB["06_geodesic_integration_and_orbit_stability"] = [
    md(r"""# 6. Geodesic integration and orbit stability
Geodesics follow $\ddot x^\mu=-\Gamma^\mu{}_{\alpha\beta}\dot x^\alpha\dot x^\beta$. We integrate with the package's own Christoffel symbols and monitor (i) the normalisation $g(u,u)=-1$ and (ii) the Killing constants $E=-g(\partial_t,u)$, $L=g(\partial_\phi,u)$.

The effective potential $V=(1-2M/r)(1+L^2/r^2)$ gives circular orbits at $V'=0$, stable when $V''>0$; the marginal case defines the ISCO at $6M$."""),
    code("""
import numpy as np, matplotlib.pyplot as plt
from gr_from_scratch.dynamics import observables as ob
from gr_from_scratch.dynamics.effective_potential import circular_orbit_analysis
from gr_from_scratch.spacetimes.schwarzschild import schwarzschild

print('ISCO, photon sphere (derived symbolically):', circular_orbit_analysis(schwarzschild(), 1)['isco'], circular_orbit_analysis(schwarzschild(), 0)['photon_sphere'])

gs = ob.schwarzschild_system()
E, L, up, ua = ob._orbit_constants(20.0, 0.5, 1.0)
rp = 1/up; f = 1 - 2/rp
tr = gs.integrate(np.array([0, rp, np.pi/2, 0.]), np.array([E/f, 0, 0, L/rp**2]), 600, rtol=1e-12, atol=1e-13, max_step=1.0)
Et, Lp = gs.conserved_function([1,0,0,0]), gs.conserved_function([0,0,0,1])
e = np.array([Et(x,u) for x,u in zip(tr.x, tr.u)]); l = np.array([Lp(x,u) for x,u in zip(tr.x, tr.u)])
print('max |g(u,u)+1| =', abs(tr.constraint).max())
print('drift in E, L  =', abs(e-e[0]).max(), abs(l-l[0]).max())
fig, ax = plt.subplots(figsize=(3.6, 3.6))
ax.plot(tr.x[:,1]*np.cos(tr.x[:,3]), tr.x[:,1]*np.sin(tr.x[:,3]), lw=0.7)
ax.add_patch(plt.Circle((0,0), 2, color='k')); ax.set_aspect('equal'); ax.set(xlabel='x/M', ylabel='y/M')
plt.show()
"""),
    md(
        r"""**Stability near the ISCO.** Give a circular orbit a small radial kick. Outside $6M$ it oscillates boundedly; inside, it plunges."""
    ),
    code("""
def run(r0, T=3000):
    ut = 1/np.sqrt(1 - 3/r0)
    x0 = np.array([0, r0, np.pi/2, 0.]); u0 = np.array([ut, 1e-3, 0, ut/r0**1.5])
    ev = lambda lam, y: y[1] - 2.2
    ev.terminal = True
    return gs.integrate(x0, u0, T, rtol=1e-10, atol=1e-12, events=[ev]).x[:,1]
for r0 in (7.5, 6.5, 5.5, 4.5):
    rr = run(r0)
    print(f'r0={r0}: r stays within [{rr.min():.2f}, {rr.max():.2f}] ->', 'bounded (stable)' if rr.min() > 0.9*r0 and rr.max() < 1.1*r0 else 'plunges (unstable)')
"""),
]

for name, cells in NB.items():
    nbf.write(notebook(cells), f"{name}.ipynb")
print("built", len(NB))
