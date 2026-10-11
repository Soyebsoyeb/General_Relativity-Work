# Reproducibility

* Environment: Python $\ge3.10$; dependencies pinned loosely in `pyproject.toml`. Results in this repository were produced with Python 3.12, NumPy 2.4, SciPy 1.17, SymPy 1.14, Matplotlib 3.10.
* Every experiment in `experiments/` is a script with a header giving the research question, equations, assumptions, parameters, numerical method, tolerances, interpretation and limitations. Run `python experiments/run_all.py`; figures go to `figures/`, numbers to `results/*.json`, and `results/SUMMARY.md` is regenerated.
* No random numbers are used in library code. `check_zero` sampling uses a fixed seed.
* Default tests are separated from expensive checks: `pytest` versus `pytest -m slow`.
* Figures can be regenerated from source; only small JSON files are stored.
