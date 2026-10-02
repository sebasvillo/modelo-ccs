"""Load legacy notebook cells as plain Python namespaces (read-only, notebook untouched).

Used by equivalence tests: the refactored package must reproduce the legacy functions
bit-for-bit until a fix is recorded in docs/errata.md.
"""

import contextlib
import io
import json
import os
import warnings
from pathlib import Path

LEGACY_NOTEBOOK = Path(__file__).resolve().parents[2] / "legacy" / "M9_Integracion.ipynb"


def cell_source(index: int) -> str:
    cells = json.loads(LEGACY_NOTEBOOK.read_text(encoding="utf-8"))["cells"]
    return "".join(cells[index]["source"])


def load_cell1() -> dict:
    """Execute cell 1 definitions (everything except the final `run_case()` call)."""
    os.environ.setdefault("MPLBACKEND", "Agg")
    src = cell_source(1)
    run_line = "best_design, all_designs = run_case()"
    assert src.rstrip().endswith(run_line)
    ns: dict = {"__name__": "legacy_cell1"}
    exec(compile(src.rstrip()[: -len(run_line)], "legacy_cell1", "exec"), ns)
    return ns


def run_legacy_case(ns: dict, **overrides) -> tuple[dict, list]:
    """Call the legacy run_case() with some module-level inputs overridden."""
    import matplotlib.pyplot as plt

    ns.update(overrides)
    with (
        contextlib.redirect_stdout(io.StringIO()),
        warnings.catch_warnings(),
    ):
        warnings.simplefilter("ignore")
        best, candidates = ns["run_case"]()
    plt.close("all")
    return best, candidates
