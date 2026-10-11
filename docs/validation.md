# Validation summary

This page lists measured results from the twelve experiments in `experiments/`. The numbers are copied from `results/*.json`, which the experiment scripts write, so nothing here is typed in by hand.

Each row compares the code against something independent: an exact quadrature, a closed-form expression, a known analytic limit, or a theoretical convergence rate. The last column says what the value should be.

| Experiment | Quantity | What was measured | Result | Expected |
|---|---|---|---|---|
| 01 | Kretschmann scalar of Schwarzschild | Slope of log K against log r | -6.0000 | $-6$, since $K=48M^{2}/r^{6}$ |
| 01 | Kretschmann scalar of Kerr | Full Riemann contraction against a closed form; largest relative difference | 4.1e-16 | rounding error |
| 02 | ISCO and photon sphere of Schwarzschild | Radii found by solving the effective-potential conditions symbolically | 6*M and 3*M | $6M$ and $3M$ |
| 03 | Periapsis advance | Geodesic integration against exact quadrature; largest relative difference | 4.9e-10 | integrator tolerance |
| 03 | Periapsis advance | Slope of the relative deviation from $6\pi M/p$ against $p$ | -1.05 | $-1$ (relative error of order $M/p$) |
| 04 | Light deflection | Geodesic integration against exact quadrature; largest relative difference | 9.4e-10 | finite start radius |
| 04 | Light deflection | Slopes of the residual after the first and second order terms | -2.02 and -3.02 | $-2$ and $-3$ |
| 05 | Fixed-step Runge-Kutta | Observed order of convergence | 3.97 | 4 |
| 06 | Circular-orbit proper period | Numerical against exact; largest absolute difference | 3.5e-14 | integrator tolerance |
| 06 | Circular-orbit proper period | Slope of the relative correction against $r$ | -1.02 | $-1$ |
| 07 | Gravitational-wave test particles | Slope of the difference (linear equation against exact separation) against amplitude | 2.00 | 2 |
| 08 | FLRW expansion, five models | Worst relative error against the exact solution | 3.3e-10 | integrator tolerance |
| 08 | FLRW expansion, five models | Largest Friedmann-constraint residual | 2.7e-08 | integrator tolerance |
| 09 | Kerr ISCO | Literature radius against the root of the engine's effective function; largest difference | 2.8e-14 | root-finder tolerance |
| 09 | Kerr circular orbit | Radial acceleration from the engine's Christoffel symbols | 2.1e-17 | rounding error |
| 10 | Scalar scattering | Largest value of $1-T-R$ over the grid | 4.4e-09 | integrator tolerance |
| 10 | Scalar scattering | $T_0/\left(16M^{2}\omega^{2}\right)$ at $M\omega=0.005$ | 1.0338 | tends to 1 |
| 11 | Second-order plane-wave curvature | Slope of the truncation error against $\epsilon$ | 4.00 | 4, because $G_{tt}$ is even in $\epsilon$ |

Further notes:

- Experiment 11 also finds that the first-order Einstein tensor of the plane wave vanishes, and that the second-order component is $G^{(2)}_{tt}=\frac{1}{4}\left(3\cos\left(2t-2z\right)-1\right)$ for unit amplitude and frequency.
- Experiment 12 finds the tidal tensor $E_{ij}=\frac{M}{r^{3}}\,\mathrm{diag}\left(-2,1,1\right)$ from Cartan's equations, with zero trace.
- The complete table, in a machine-readable layout, is `results/SUMMARY.md`.

## How to read the tolerances

A result labelled "integrator tolerance" means the difference is limited by the relative tolerance requested from the numerical integrator, typically $10^{-10}$ to $10^{-12}$, and not by any physical approximation. A result labelled "rounding error" is at the level of double-precision arithmetic.
