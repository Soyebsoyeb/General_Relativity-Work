"""Friedmann dynamics.

With H0 = 1 and a(t0) = 1 normalisation,
    (adot/a)^2 = Om_r a^-4 + Om_m a^-3 + Om_k a^-2 + Om_L,     Om_k = 1 - Om_r - Om_m - Om_L,
    addot / a  = -(1/2) [ 2 Om_r a^-4 + Om_m a^-3 - 2 Om_L ].
Both follow from G_{mn} + Lambda g_{mn} = 8 pi T_{mn} for FLRW with dust,
radiation and a cosmological constant; the symbolic derivation is in
``tests/test_flrw.py`` and ``docs/derivations/10_flrw.md``. The second-order
form is integrated so that recollapse is handled; the first-order Friedmann
constraint is monitored as a residual.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Cosmology:
    om_m: float = 0.0
    om_r: float = 0.0
    om_l: float = 0.0

    @property
    def om_k(self) -> float:
        return 1.0 - self.om_m - self.om_r - self.om_l

    def e2(self, a):
        return self.om_r / a**4 + self.om_m / a**3 + self.om_k / a**2 + self.om_l

    def addot(self, a):
        return -0.5 * a * (2 * self.om_r / a**4 + self.om_m / a**3 - 2 * self.om_l)


def integrate(cosmo: Cosmology, a0: float, adot0: float, t_end: float, n: int = 801, rtol: float = 1e-11):
    ts = np.linspace(0, t_end, n)

    def rhs(_t, y):
        return [y[1], cosmo.addot(y[0])]

    sol = solve_ivp(rhs, (0, t_end), [a0, adot0], t_eval=ts, rtol=rtol, atol=1e-14, method="DOP853")
    a, ad = sol.y
    constraint = (ad / a) ** 2 - cosmo.e2(a)
    return sol.t, a, ad, constraint


# analytic references ---------------------------------------------------------------
def a_matter_flat(t):
    return (1.5 * t) ** (2 / 3)


def a_radiation_flat(t):
    return (2 * t) ** 0.5


def a_de_sitter(t, a0=1.0):
    return a0 * np.exp(t)


def a_lcdm_flat(t, om_m):
    om_l = 1 - om_m
    return (om_m / om_l) ** (1 / 3) * np.sinh(1.5 * np.sqrt(om_l) * t) ** (2 / 3)


def closed_dust(eta, om_m):
    """Parametric cycloid for Om_m > 1, Lambda = 0: returns (t, a)."""
    ok = abs(1 - om_m)
    a_max = om_m / ok
    return om_m / (2 * ok**1.5) * (eta - np.sin(eta)), 0.5 * a_max * (1 - np.cos(eta))
