"""Conventions used throughout the package.

Signature
    Lorentzian metrics use (-, +, +, +).

Curvature
    R^rho_{sigma mu nu} = d_mu Gamma^rho_{nu sigma} - d_nu Gamma^rho_{mu sigma}
                          + Gamma^rho_{mu lam} Gamma^lam_{nu sigma}
                          - Gamma^rho_{nu lam} Gamma^lam_{mu sigma}
    Ricci:  R_{sigma nu} = R^rho_{sigma rho nu}
    Scalar: R = g^{sigma nu} R_{sigma nu}
    With these choices a sphere has positive Ricci scalar and the Einstein
    equations read G_{mu nu} + Lambda g_{mu nu} = 8 pi G T_{mu nu} / c^4.

Index placement
    Tensors store a variance string such as ``"uddd"``: ``u`` is a
    contravariant (upper) slot, ``d`` a covariant (lower) slot.
    The Riemann tensor is stored as R^rho_{sigma mu nu}, variance ``"uddd"``.
    The index added by covariant differentiation is appended last, so
    ``(nabla T)_{a b ; c}`` has the derivative index in the final slot.

Units
    Unless a function says otherwise, G = c = 1 and lengths are measured in
    units of the mass parameter M.
"""

from __future__ import annotations

SIGNATURE_LORENTZIAN = (-1, 1, 1, 1)
UP, DOWN = "u", "d"
