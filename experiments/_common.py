"""Shared helpers for the experiment scripts."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

FIG = ROOT / "figures"
RES = ROOT / "results"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)


def save_fig(fig, name: str) -> Path:
    path = FIG / f"{name}.png"
    fig.savefig(path)
    return path


def _encode(o):
    if hasattr(o, "tolist"):
        return o.tolist()
    return float(o)


def save_results(name: str, data: dict) -> Path:
    path = RES / f"{name}.json"
    path.write_text(json.dumps(data, indent=2, default=_encode) + "\n")
    return path
