"""Execute every notebook in place so outputs are stored. Run from the notebooks/ directory."""

import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

root = Path(__file__).resolve().parent
sys.path.insert(0, str(root.parent / "src"))
failed = []
prefixes = sys.argv[1:]  # optional: only notebooks whose name starts with one of these
for path in sorted(root.glob("[0-9][0-9]_*.ipynb")):
    if prefixes and not any(path.name.startswith(p) for p in prefixes):
        continue
    nb = nbformat.read(path, as_version=4)
    try:
        NotebookClient(nb, timeout=900, kernel_name="python3", resources={"metadata": {"path": str(root)}}).execute()
        nbformat.write(nb, path)
        print("ok  ", path.name)
    except Exception as exc:  # noqa: BLE001
        failed.append(path.name)
        print("FAIL", path.name, str(exc)[-400:])
sys.exit(1 if failed else 0)
