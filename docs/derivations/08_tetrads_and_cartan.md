# 8. Tetrads and Cartan's structure equations

**In one paragraph.** Instead of coordinates, one can describe geometry with a set of four orthonormal basis vectors at each point, called a tetrad. In this language, curvature comes from two short equations written with differential forms. The result is an independent route to the same curvature, so agreement with the coordinate calculation is a strong check of both.

## 1. Setting

An orthonormal coframe is a set of four one-forms $e^{a}=e^{a}{}_{\mu}\,dx^{\mu}$ such that

$$
g_{\mu\nu}=\eta_{ab}\,e^{a}{}_{\mu}\,e^{b}{}_{\nu}.
$$

Frame indices $a,b,\dots$ are raised and lowered with the constant matrix $\eta_{ab}$. The matching frame vectors $E_a{}^{\mu}$ are the matrix inverse of $e^{a}{}_{\mu}$.

Differential forms are stored by the code as dictionaries of coefficients on sorted index tuples. The wedge product and the exterior derivative are implemented, and the tests check $d^{2}=0$ and the graded Leibniz rule.

## 2. The structure equations

For a connection without torsion, the **first structure equation** is

$$
de^{a}+\omega^{a}{}_{b}\wedge e^{b}=0,\qquad \omega_{ab}=-\omega_{ba}.
$$

The antisymmetry in $a,b$ expresses compatibility with the frame metric. The **second structure equation** defines the curvature two-form,

$$
\Omega^{a}{}_{b}=d\omega^{a}{}_{b}+\omega^{a}{}_{c}\wedge\omega^{c}{}_{b}=\frac{1}{2}\,R^{a}{}_{bcd}\;e^{c}\wedge e^{d}.
$$

## 3. Solving for the connection without Christoffel symbols

Expand the exterior derivative of the coframe in the frame basis,

$$
de^{a}=\frac{1}{2}\,C^{a}{}_{bc}\;e^{b}\wedge e^{c},
\qquad
C^{a}{}_{bc}=E_b{}^{\mu}E_c{}^{\nu}\left(\partial_\mu e^{a}{}_{\nu}-\partial_\nu e^{a}{}_{\mu}\right).
$$

Write the spin connection as $\omega^{a}{}_{b}=\omega^{a}{}_{bc}\,e^{c}$. The first structure equation then says

$$
C_{abc}=\omega_{abc}-\omega_{acb}.
$$

Write this three times with cyclically permuted indices and use $\omega_{abc}=-\omega_{bac}$. Adding the first two and subtracting the third leaves $2\omega_{abc}$, so

$$
\omega_{abc}=\frac{1}{2}\left(C_{abc}+C_{bca}-C_{cab}\right).
$$

This is purely algebraic in the derivatives of the frame. At no point are Christoffel symbols used.

## 4. Comparing with the coordinate calculation

The frame components of the Riemann tensor from the second structure equation must equal the coordinate Riemann tensor projected onto the frame:

$$
R^{a}{}_{bcd}=e^{a}{}_{\rho}\,E_b{}^{\sigma}\,E_c{}^{\mu}\,E_d{}^{\nu}\,R^{\rho}{}_{\sigma\mu\nu}.
$$

This is verified **exactly** for three cases:

- Schwarzschild, with the static frame;
- FLRW with spatial curvature $k$ and an unspecified function $a(t)$;
- the unit sphere, where the frame curvature $R^{\hat\theta}{}_{\hat\varphi\hat\theta\hat\varphi}$ equals $+1$.

The first structure equation and the antisymmetry of $\omega_{ab}$ are also verified.

The tidal tensor of page 3 is obtained from this frame calculation, with no coordinate Riemann tensor involved.

## 5. Why this formulation matters

These points are explanation only. The code does **not** implement them.

- The frame can be rotated by a Lorentz transformation that varies from point to point, $e^{a}\to\Lambda^{a}{}_{b}(x)\,e^{b}$. The spin connection then transforms like a gauge field, and $\Omega$ is its field strength. This is local Lorentz symmetry.
- Spinors, such as electrons, couple to gravity only through the tetrad and the spin connection, with $\nabla_\mu\psi=\partial_\mu\psi+\frac{1}{4}\omega_{\mu ab}\gamma^{ab}\psi$.
- First-order gravity treats $e$ and $\omega$ as independent fields.

## 6. Implementation

| Task | Module |
|---|---|
| Forms, wedge product, exterior derivative | `geometry/differential_forms.py` |
| Coframe, spin connection, curvature two-form, comparison | `geometry/tetrads.py` |
| Tidal tensor | `dynamics/geodesic_deviation.py` |
| Tests | `tests/test_tetrads.py` |

## 7. Limits

Frames are supplied by the user (there is no automatic orthonormalisation). The dependence of the spin connection on the choice of Lorentz frame is not explored. There are no spinors and no first-order action.
