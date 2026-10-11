# Validation summary

The table is copied from `results/SUMMARY.md`, which is generated from the experiment outputs. Each row has an independent check: a quadrature versus a geodesic integration, a closed form versus the engine, or a known analytic limit.


| Exp. | Quantity | Measured | Reference / expectation |
|---|---|---|---|
| 01 | Schwarzschild K(r) log-log slope | -6.000000 | -6 (K = 48 M^2/r^6) |
| 01 | Kerr K: engine vs closed form, max rel. err. | 4.05e-16 | round-off |
| 02 | ISCO / photon sphere (symbolic) | 6*M / 3*M | 6M / 3M |
| 03 | Periapsis advance: geodesic vs quadrature, max rel. diff. | 4.91e-10 | rtol 1e-12 |
| 03 | Slope of abs(exact/(6 pi M/p) - 1) vs p | -1.05 | -1 (O(M/p)) |
| 04 | Deflection: geodesic vs quadrature, max rel. diff. | 9.38e-10 | finite r0 |
| 04 | Slopes of residual after 1st / 2nd order | -2.02 / -3.02 | -2 / -3 |
| 05 | RK4 observed order | 3.97 | 4 |
| 06 | Circular-orbit proper period: numeric vs exact, max abs. diff. | 3.49e-14 | rtol 1e-12 |
| 06 | Slope of period correction vs r | -1.02 | -1 |
| 07 | Linear deviation vs exact separation: slope in amplitude | 2.00 | 2 |
| 08 | matter (flat): max rel. error / constraint | 3.39e-11 / 3.40e-10 | rtol 1e-11 |
| 08 | radiation (flat): max rel. error / constraint | 3.30e-10 / 9.11e-09 | rtol 1e-11 |
| 08 | de Sitter: max rel. error / constraint | 3.08e-11 / 5.35e-11 | rtol 1e-11 |
| 08 | flat $\Lambda$CDM: max rel. error / constraint | 2.45e-11 / 1.07e-09 | rtol 1e-11 |
| 08 | closed dust ($\Omega_m=2$): max rel. error / constraint | 8.62e-11 / 2.66e-08 | rtol 1e-11 |
| 09 | Kerr ISCO: literature vs W''=0, max discrepancy | 2.75e-14 | root-finder |
| 09 | Circular-orbit radial acceleration (engine Christoffels) | 2.08e-17 | round-off |
| 10 | Max abs(1 - T - R) over grid (l<=3) | 4.43e-09 | rtol 1e-10 |
| 10 | T0/(16 w^2) at w=0.005, 0.01, 0.02, 0.05 | 1.0338, 1.0714, 1.1570, 1.4829 | -> 1 as w -> 0 |
| 10 | Sensitivity to r*_max (l=2, w=0.3) | 3.07e-11 |  |
| 11 | G^(1) vanishes; G^(2)_tt | True; 3*cos(2*t - 2*z)/4 - 1/4 |  |
| 11 | Truncation-error slope in eps | 4.00 | 4 (G_tt even in eps) |
| 12 | Tidal tensor E_ij (Cartan frame) | Matrix([[-2*M/r**3, 0, 0], [0, M/r**3, 0], [0, 0, M/r**3]]) | diag(-2M, M, M)/r^3 |
