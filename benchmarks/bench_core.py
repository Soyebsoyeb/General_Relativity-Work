"""Timing of the main code paths. Numbers depend on the machine; they are written to
results/benchmarks.json together with the Python and package versions."""

import json
import platform
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import scipy
import sympy as sp

from gr_from_scratch.dynamics import observables as ob
from gr_from_scratch.field_equations.vacuum_equations import compute_geometry
from gr_from_scratch.geometry.invariants import kretschmann
from gr_from_scratch.spacetimes import flrw, minkowski
from gr_from_scratch.spacetimes.schwarzschild import schwarzschild


def timed(fn, repeat=1):
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        out = fn()
        best = min(best, time.perf_counter() - t0)
    return best, out


results = {}
results["schwarzschild_full_geometry_s"], g = timed(lambda: compute_geometry(schwarzschild()))
results["schwarzschild_kretschmann_s"], _ = timed(lambda: kretschmann(g.riemann))
results["flrw_full_geometry_s"], _ = timed(lambda: compute_geometry(flrw.flrw()))
results["minkowski_geometry_s"], _ = timed(lambda: compute_geometry(minkowski.minkowski()))
t, gs = timed(ob.schwarzschild_system)
results["geodesic_system_build_s"] = t
E, L, up, ua = ob._orbit_constants(20.0, 0.5, 1.0)
rp = 1 / up
f = 1 - 2 / rp
x0, u0 = np.array([0, rp, np.pi / 2, 0.0]), np.array([E / f, 0, 0, L / rp**2])
tt, tr = timed(lambda: gs.integrate(x0, u0, 5000.0, rtol=1e-10, atol=1e-12), repeat=3)
results["geodesic_5000M_rtol1e-10_s"] = tt
results["geodesic_rhs_evaluations"] = int(tr.sol.nfev)
results["geodesic_rhs_per_second"] = tr.sol.nfev / tt
meta = {
    "python": sys.version.split()[0],
    "platform": platform.platform(),
    "numpy": np.__version__,
    "scipy": scipy.__version__,
    "sympy": sp.__version__,
}
out = {"meta": meta, "timings": results}
(Path(__file__).resolve().parents[1] / "results" / "benchmarks.json").write_text(json.dumps(out, indent=2) + "\n")
for k, v in results.items():
    print(f"{k:36s} {v:.4g}")
