# Reproducibility

## Environment

Python 3.10 or newer. The results in the repository were produced with Python 3.12, NumPy 2.4, SciPy 1.17, SymPy 1.14 and Matplotlib 3.10. Install with

```bash
pip install -e ".[dev]"
```

## Tests

```bash
pytest              # fast suite
pytest -m slow      # expensive symbolic checks
```

The fast suite takes about two minutes on one core. The slow tests are separated because they take tens of seconds each.

## Experiments

Each script in `experiments/` begins with a header stating its research question, equations, assumptions, parameters, numerical method, tolerances, interpretation and limitations. To rerun all of them:

```bash
python experiments/run_all.py
```

Figures are written to `figures/`, measured numbers to `results/*.json`, and `results/SUMMARY.md` is regenerated from them. Only small JSON files are stored.

## Notebooks

```bash
python notebooks/execute_all.py
```

executes all twelve notebooks and stores their outputs.

## Determinism

No random numbers are used in the library. The sampling inside `check_zero` uses a fixed seed. Numerical results can differ in the last digits between machines and library versions; the tolerances quoted in the pages are far above that level.
