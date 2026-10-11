"""Run every experiment in order, then regenerate results/SUMMARY.md."""

import subprocess
import sys
from pathlib import Path

here = Path(__file__).parent
for script in sorted(here.glob("exp[0-9][0-9]_*.py")):
    print(f"== {script.name}", flush=True)
    subprocess.run([sys.executable, str(script)], check=True, cwd=here)
subprocess.run([sys.executable, str(here / "summarize.py")], check=True, cwd=here)
