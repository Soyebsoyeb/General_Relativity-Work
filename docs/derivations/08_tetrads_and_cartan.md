# 8. Differential forms, tetrads and Cartan's equations

## Definitions
Differential forms are stored as $\{I:\text{coefficient}\}$ with $I$ a sorted index tuple; `Form.wedge` and `Form.d` implement $\wedge$ and exterior derivative ($d^2=0$ and the graded Leibniz rule are tested). An orthonormal coframe $e^a=e^a{}_\mu dx^\mu$ satisfies $g_{\mu\nu}=\eta_{ab}e^a{}_\mu e^b{}_\nu$; $E_a{}^\mu$ is its inverse.

## Cartan structure equations
Zero torsion: $de^a+\omega^a{}_b\wedge e^b=0$, with $\omega_{ab}=-\omega_{ba}$ (frame metric compatibility). Curvature: $\Omega^a{}_b=d\omega^a{}_b+\omega^a{}_c\wedge\omega^c{}_b=\tfrac12R^a{}_{bcd}e^c\wedge e^d$.

**Solution without Christoffel symbols.** Write $de^a=\tfrac12C^a{}_{bc}e^b\wedge e^c$ with $C^a{}_{bc}=E_b{}^\mu E_c{}^\nu(\partial_\mu e^a{}_\nu-\partial_\nu e^a{}_\mu)$ and $\omega^a{}_b=\omega^a{}_{bc}e^c$. Then $C_{abc}=\omega_{abc}-\omega_{acb}$. Using $\omega_{abc}=-\omega_{bac}$,
$$C_{abc}+C_{bca}-C_{cab}=2\omega_{abc}\quad\Longrightarrow\quad\omega_{abc}=\tfrac12(C_{abc}+C_{bca}-C_{cab}).$$

## Equivalence with the coordinate formulation
$R^a{}_{bcd}$ obtained from $\Omega^a{}_b$ must equal $e^a{}_\rho E_b{}^\sigma E_c{}^\mu E_d{}^\nu R^\rho{}_{\sigma\mu\nu}$. This is verified exactly for Schwarzschild (static frame), FLRW with spatial curvature $k$ (symbolic $a(t)$), and the unit sphere ($K=+1$). The first structure equation and the antisymmetry of $\omega_{ab}$ are also verified.

## Why it matters
The frame formulation makes local Lorentz symmetry $e^a\to\Lambda^a{}_b(x)e^b$ explicit: $\omega$ transforms as a $SO(3,1)$ gauge connection, and $\Omega$ is its field strength. Spinors couple to gravity only through the tetrad and spin connection, $\nabla_\mu\psi=\partial_\mu\psi+\tfrac14\omega_{\mu ab}\gamma^{ab}\psi$. First-order (Palatini/Cartan) gravity takes $(e,\omega)$ as independent variables. These points are explanatory only: **spinors and first-order actions are not implemented**.

## Limitations
Frames are supplied by the user (no automatic Gram–Schmidt); the Lorentz-gauge dependence of $\omega$ is not explored; the tidal tensor of experiment 12 uses the static frame.
