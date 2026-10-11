import nbformat as nbf
from _nb import code, md, notebook

NB = {}

NB["07_light_deflection_and_perihelion_precession"] = [
    md(r"""# 7. Light deflection and perihelion precession
Both observables follow from the orbit equation in $u=1/r$,
$$\Bigl(\frac{du}{d\phi}\Bigr)^2=\frac{E^2-1}{L^2}+\frac{2M}{L^2}u-u^2+2Mu^3 \quad(\text{timelike}),\qquad \Bigl(\frac{du}{d\phi}\Bigr)^2=\frac1{b^2}-u^2+2Mu^3\quad(\text{null}).$$
Each result is computed two independent ways: a **quadrature** of the exact orbit equation and a full **geodesic integration**. Weak-field formulas are compared with both."""),
    code("""
import numpy as np, matplotlib.pyplot as plt
from gr_from_scratch.dynamics import observables as ob

bs = np.array([8, 12, 20, 50, 100, 300, 1000.])
exact = np.array([ob.deflection_quad(b) for b in bs])
geo = np.array([ob.deflection_geodesic(b) for b in bs])
print('b        exact        4M/b     +15pi/4(M/b)^2   geodesic-exact')
for b, e, g in zip(bs, exact, geo):
    print(f'{b:6.0f} {e:12.6e} {ob.deflection_weak_field(b,order=1):12.6e} {ob.deflection_weak_field(b,order=2):12.6e} {g-e: .1e}')
"""),
    md(
        r"""The residual of the first-order formula should fall as $(M/b)^2$ and of the second-order formula as $(M/b)^3$:"""
    ),
    code("""
big = bs > 20
for order in (1, 2):
    res = np.abs(exact - np.array([ob.deflection_weak_field(b, order=order) for b in bs]))
    print('order', order, 'slope of residual vs b:', np.polyfit(np.log(bs[big]), np.log(res[big]), 1)[0])
"""),
    md(
        r"""**Periapsis advance.** For a bound orbit with semi-latus rectum $p$ and eccentricity $e$, the leading result is $6\pi M/p$ per orbit. The exact value follows from a sine substitution that removes the endpoint singularities (see derivation 3)."""
    ),
    code("""
e = 0.3
for p in (20., 100., 1000.):
    q = ob.periapsis_advance_quad(p, e); g, per = ob.periapsis_advance_geodesic(p, e, n_orbits=2)
    print(f'p={p:6.0f}  exact={q:.8f}  geodesic={g:.8f}  6piM/p={ob.periapsis_advance_leading(p):.8f}  ratio exact/leading={q/ob.periapsis_advance_leading(p):.4f}')
"""),
]

NB["08_newtonian_limit_comparisons"] = [
    md(r"""# 8. Newtonian-limit comparisons
For $g_{00}=-(1+2\Phi)$, $g_{ij}=(1-2\Phi)\delta_{ij}$ the linearised Einstein tensor has $G^{(1)}_{00}=2\nabla^2\Phi$, so $G_{00}=8\pi\rho$ is Poisson's equation $\nabla^2\Phi=4\pi\rho$, and $\Gamma^i{}_{00}=\partial_i\Phi$ gives $\ddot x=-\nabla\Phi$."""),
    code("""
import sympy as sp
from gr_from_scratch.field_equations.linearized_gravity import linearized_einstein
t, x, y, z = sp.symbols('t x y z', real=True)
Phi = sp.Function('Phi')(x, y, z)
h = sp.diag(-2*Phi, -2*Phi, -2*Phi, -2*Phi)
G = linearized_einstein(h, [t, x, y, z])
lap = sum(sp.diff(Phi, v, 2) for v in (x, y, z))
print('G00 - 2*laplacian(Phi) =', sp.simplify(G[0,0] - 2*lap))
"""),
    md(
        r"""**Where does the Newtonian description fail?** For circular orbits Kepler's law holds exactly in Schwarzschild time, but the *proper* period differs by $\sqrt{1-3M/r}$, i.e. a relative correction $\simeq 3M/2r$. We measure it by integrating the geodesic through $\Delta\phi=2\pi$."""
    ),
    code("""
import numpy as np, matplotlib.pyplot as plt
from gr_from_scratch.dynamics import observables as ob
from gr_from_scratch.dynamics.geodesics import event_function
gs = ob.schwarzschild_system()
radii = np.logspace(np.log10(8), np.log10(1500), 8)
dev, ex = [], []
for r in radii:
    ut = 1/np.sqrt(1-3/r)
    ev = event_function(3, 2*np.pi, direction=+1, terminal=True)
    tr = gs.integrate(np.array([0,r,np.pi/2,0.]), np.array([ut,0,0,ut/r**1.5]), 3*2*np.pi*r**1.5, rtol=1e-12, atol=1e-13, events=[ev])
    dev.append(tr.sol.t_events[0][0]/(2*np.pi*r**1.5) - 1); ex.append(np.sqrt(1-3/r)-1)
dev = np.array(dev)
print('max |numeric - exact| =', np.abs(dev-np.array(ex)).max())
fig, ax = plt.subplots(figsize=(4.2, 3))
ax.loglog(radii, np.abs(dev), 'o', label='geodesic'); ax.loglog(radii, 1.5/radii, '--', label='3M/2r')
ax.axhline(0.01, color='C1', lw=0.6); ax.set(xlabel='r/M', ylabel=r'$|\\tau/\\tau_N-1|$'); ax.legend(); plt.show()
print('radius where the correction reaches 1% (leading order): r =', 1.5/0.01, 'M')
"""),
]

NB["09_gravitational_wave_perturbations"] = [
    md(r"""# 9. Gravitational-wave perturbations
Write $g=\eta+h$. Under $x\to x+\xi$, $h_{\mu\nu}\to h_{\mu\nu}+\partial_\mu\xi_\nu+\partial_\nu\xi_\mu$. The linearised Riemann tensor is invariant, so most of the ten components of $h$ are gauge. In TT gauge a vacuum wave has two polarisations $h_+,h_\times$, and satisfies the field equations only if $k=\omega$."""),
    code("""
import sympy as sp
from gr_from_scratch.field_equations import linearized_gravity as lg
t, x, y, z = sp.symbols('t x y z', real=True); X = [t, x, y, z]
h = sp.Matrix(4, 4, lambda a, b: sp.Function(f'h{min(a,b)}{max(a,b)}')(*X))
xi = [sp.Function(f'xi{i}')(*X) for i in range(4)]
R1 = lg.linearized_riemann(h, X); R2 = lg.linearized_riemann(lg.gauge_transform(h, xi, X), X)
print('linearised Riemann gauge invariant:', all(sp.simplify(a-b)==0 for a,b in zip(R1.ravel(), R2.ravel())))

w, k, A, B = sp.symbols('omega k A B', positive=True)
hw = lg.plane_wave_tt(X, A, B, w, k=k)
print('TT conditions:', lg.tt_conditions(hw, X))
G = lg.linearized_einstein(hw, X)
print('vacuum if k=omega:', all(sp.simplify(e.subs(k, w))==0 for e in G), '| with k != omega, nonzero:', any(sp.simplify(e)!=0 for e in G))
"""),
    md(
        r"""**Test particles.** The linearised geodesic-deviation equation $\ddot\xi^i=\tfrac12\ddot h_{ij}\xi^j$ is compared with the exact proper separation of particles held at fixed TT coordinates. They should differ at $O(A^2)$."""
    ),
    code("""
import numpy as np, matplotlib.pyplot as plt
from gr_from_scratch.dynamics.gravitational_waves import integrate_deviation, proper_separation
amps = np.array([1e-3, 3e-3, 1e-2, 3e-2, 1e-1]); diffs = []
for A in amps:
    ts, xi_ = integrate_deviation((1,0), A, 0.0, 1.0, 4*np.pi, n=801)
    diffs.append(np.abs(np.linalg.norm(xi_, axis=1) - proper_separation(ts, (1,0), A, 0.0, 1.0)).max())
print('slope of difference vs amplitude:', np.polyfit(np.log(amps), np.log(diffs), 1)[0])

th = np.linspace(0, 2*np.pi, 25)[:-1]
fig, axs = plt.subplots(1, 4, figsize=(8, 2.2))
for ax, ph in zip(axs, [np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]):
    pts = np.array([integrate_deviation((np.cos(a), np.sin(a)), 0.2, 0.0, 1.0, ph, n=3)[1][-1] for a in th])
    ax.plot(*np.vstack([pts, pts[:1]]).T, '-o', ms=2); ax.set_aspect('equal'); ax.set_title(f'wt={ph/np.pi:.1f} pi', fontsize=8)
plt.show()
"""),
]

NB["10_flrw_cosmology"] = [
    md(r"""# 10. FLRW cosmology
For $ds^2=-dt^2+a^2[dr^2/(1-kr^2)+r^2d\Omega^2]$ the engine gives $G^t{}_t=-3(\dot a^2+k)/a^2$ and $G^r{}_r=-(2a\ddot a+\dot a^2+k)/a^2$. With a perfect fluid and $\Lambda$ these are the Friedmann equations. The continuity equation is *not* imposed: it follows from the contracted Bianchi identity."""),
    code("""
import sympy as sp
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry, field_equation_residual
from gr_from_scratch.field_equations import stress_energy as se
from gr_from_scratch.spacetimes import flrw as fl
g = compute_geometry(fl.flrw())
a = fl.a(fl.t); ad = sp.diff(a, fl.t); add = sp.diff(a, fl.t, 2)
Gm = g.einstein.raise_index(0).c
print('G^t_t =', sp.simplify(Gm[0,0])); print('G^r_r =', sp.simplify(Gm[1,1]))
rho = 3*(ad**2 + fl.k)/(8*sp.pi*a**2); p = -(2*a*add + ad**2 + fl.k)/(8*sp.pi*a**2)
print('Einstein eq with this fluid holds:', field_equation_residual(g, se.perfect_fluid(g.metric, rho, p, [1,0,0,0])).is_zero())
print('continuity residual (not imposed):', sp.simplify(sp.diff(rho, fl.t) + 3*ad/a*(rho+p)))
"""),
    md(
        r"""**Numerics versus analytic solutions** ($H_0=1$). Initial data are taken from the analytic solution, so the comparison tests the integrator and the ODE, while the derivation of the ODE from $G=8\pi T$ is the symbolic check above."""
    ),
    code("""
import numpy as np, matplotlib.pyplot as plt
from gr_from_scratch.numerical import friedmann as fr
def run(c, afn, t0, t1):
    h=1e-6; a0=afn(t0); ad0=(afn(t0+h)-afn(t0-h))/(2*h)
    t,a,_,con = fr.integrate(c, a0, ad0, t1-t0); return t+t0, a, afn(t+t0), con
cases = {'matter': run(fr.Cosmology(om_m=1), fr.a_matter_flat, .2, 3), 'radiation': run(fr.Cosmology(om_r=1), fr.a_radiation_flat, .05, 1.5),
         'de Sitter': run(fr.Cosmology(om_l=1), fr.a_de_sitter, 0, 3), 'LCDM': run(fr.Cosmology(om_m=.3, om_l=.7), lambda t: fr.a_lcdm_flat(t,.3), .1, 2.5)}
fig, ax = plt.subplots(1, 2, figsize=(7, 2.6))
for n, (t, a, ref, con) in cases.items():
    ax[0].plot(t, a, label=n); ax[1].semilogy(t, np.abs(a-ref)/ref+1e-16)
    print(f'{n:10s} max rel. error {np.max(np.abs(a-ref)/ref):.1e}   max constraint {np.abs(con).max():.1e}')
ax[0].legend(); ax[0].set(xlabel='H0 t', ylabel='a'); ax[1].set(xlabel='H0 t', ylabel='relative error'); plt.show()
"""),
]

NB["11_tetrad_based_curvature_calculations"] = [
    md(r"""# 11. Tetrad-based curvature
With an orthonormal coframe $e^a$, Cartan's equations $de^a+\omega^a{}_b\wedge e^b=0$ and $\Omega^a{}_b=d\omega^a{}_b+\omega^a{}_c\wedge\omega^c{}_b$ give the curvature **without Christoffel symbols**. The connection is solved algebraically: $\omega_{abc}=\tfrac12(C_{abc}+C_{bca}-C_{cab})$ where $de^a=\tfrac12C^a{}_{bc}e^b\wedge e^c$. We compare with the coordinate calculation."""),
    code("""
import sympy as sp, numpy as np
from gr_from_scratch.geometry.tetrads import Coframe
from gr_from_scratch.geometry.differential_forms import Form
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry
from gr_from_scratch.spacetimes.schwarzschild import schwarzschild, f, r, theta
from gr_from_scratch.dynamics.geodesic_deviation import tidal_tensor

m = schwarzschild(); ff = f()
cf = Coframe(m, sp.diag(sp.sqrt(ff), 1/sp.sqrt(ff), r, r*sp.sin(theta)))
z = lambda a: all(m.simplifier(e) == 0 for e in np.asarray(a, dtype=object).ravel())
print('first structure equation satisfied:', z(cf.first_structure_residual()))
Rf = cf.riemann_frame(); Rc = cf.project_coordinate_riemann(compute_geometry(m).riemann)
print('frame Riemann == coordinate Riemann projected:', z(Rf - Rc))
print('R^{t}_{r t r} (frame) =', Rf[0,1,0,1], '   R^{theta}_{phi theta phi} =', Rf[2,3,2,3])
print('tidal tensor E_ij =', sp.Matrix(tidal_tensor(cf).tolist()).applyfunc(sp.simplify))
"""),
    md(r"""The tidal tensor $E_{ij}=R^i{}_{0j0}=\mathrm{diag}(-2M,M,M)/r^3$ is trace-free as required in vacuum; $\ddot\xi^i=-E_{ij}\xi^j$ stretches radially and compresses transversally.

**Exterior calculus check.** $d^2=0$ on arbitrary forms:"""),
    code("""
x, y, w = sp.symbols('x y w', real=True)
a = Form.one_form([x, y, w], [x*y*w, sp.sin(x)*w**2, y + x**3])
print('d(d a) = 0:', a.d().d().is_zero())
"""),
]

NB["12_quantum_fields_on_curved_backgrounds"] = [
    md(r"""# 12. A scalar field on a curved background
This is a **classical** field on a **fixed** background: it is neither quantised gravity nor (here) a quantised field. The equation is $(\Box-m^2-\xi R)\phi=0$. On Schwarzschild, with $\phi=e^{-i\omega t}P(\theta)\psi(r)/r$ and $P''+\cot\theta P'=-\ell(\ell+1)P$, it reduces to
$$\frac{d^2\psi}{dr_*^2}+(\omega^2-V)\psi=0,\qquad V=\Bigl(1-\frac{2M}r\Bigr)\Bigl(\frac{\ell(\ell+1)}{r^2}+\frac{2M}{r^3}+m^2\Bigr).$$"""),
    code("""
from gr_from_scratch.quantum import klein_gordon as kg
print('residual of the claimed reduction (0 means verified):', kg.schwarzschild_radial_reduction())
"""),
    md(
        r"""**Boundary conditions.** Ingoing at the horizon ($\psi=e^{-i\omega r_*}$), and at large $r_*$ a decomposition into incoming and outgoing Riccati–Hankel waves. Transmission is $T=1/|A_{\rm in}|^2$ and $1-T-R$ monitors flux conservation."""
    ),
    code("""
import numpy as np, matplotlib.pyplot as plt
ws = np.linspace(0.05, 1.0, 14)
fig, ax = plt.subplots(figsize=(4.4, 3))
for l in range(3):
    res = [kg.scattering(l, float(w), rstar_max=max(3000., 80/w)) for w in ws]
    ax.plot(ws, [s.transmission for s in res], label=f'l={l}')
    print(f'l={l}: max |1-T-R| = {max(abs(s.flux_residual) for s in res):.1e}')
ax.set(xlabel=r'$M\\omega$', ylabel='$T_l$'); ax.legend(); plt.show()
"""),
    md(
        r"""**Analytic limit.** As $\omega\to0$, $T_0\to16M^2\omega^2$ so that $\sigma_{\rm abs}=\pi T_0/\omega^2\to16\pi M^2$, the horizon area. The ratio should approach 1:"""
    ),
    code("""
for w in (0.05, 0.02, 0.01, 0.005):
    s = kg.scattering(0, w, rstar_max=max(3000., 120/w)); print(f'M w = {w:6.3f}   T0/(16 M^2 w^2) = {s.transmission/(16*w*w):.4f}')
"""),
]

for name, cells in NB.items():
    nbf.write(notebook(cells), f"{name}.ipynb")
print("built", len(NB))
