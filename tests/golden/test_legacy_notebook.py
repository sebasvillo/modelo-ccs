"""Golden master: re-execute the frozen legacy notebook and compare against its saved outputs.

This is a refactoring safety net only. The legacy numbers are known to be wrong
(docs/audit-2026-09.md); later phases move them on purpose, through docs/errata.md.

The notebook runs from a temporary copy, in a temporary directory (so its savefig calls
do not touch the repo), with the Agg backend, top to bottom (the saved execution order,
In[91]..In[101]). Only stdout is compared: matplotlib warnings go to stderr.
"""

import json
import re
import shutil
import time
from pathlib import Path

import nbformat
import pytest
from nbclient import NotebookClient

LEGACY_NOTEBOOK = Path(__file__).resolve().parents[2] / "legacy" / "M9_Integracion.ipynb"

pytestmark = pytest.mark.slow

# Cell indices are positions in nb.cells (markdown cells included).
# Cell 2 reads stale kernel globals (`best`, `design`, `case`) through globals(); the saved
# run (In[92]) picked up a stale `case` with P_total_W = 66 977.8 kW instead of the
# 66 959.4 kW printed by cell 1 in the same run, so its saved outputs cannot be reproduced
# from a clean kernel (audit 6.7). It is checked against clean-kernel values instead.
CELL_STALE_STATE = 2
CELL_SAVED_ERROR = 6

# (cell, label, regex, value as printed). Hand-written second check of the saved outputs.
SAVED_VALUES = [
    (1, "packed height [m]", r"Packed height needed: ([\d.]+) m", "14.86"),
    (1, "pressure drop [Pa]", r"Pressure drop: ([\d.]+) Pa", "1919"),
    (1, "blower power [kW]", r"Gas-side blowers: ([\d.]+) kW", "151.8"),
    (1, "pump power [kW]", r"Pump power: ([\d.]+) kW", "48.2"),
    (1, "KGa [1/s]", r"Volumetric KGa: ([\d.]+) 1/s", "323.4097"),
    (1, "Hatta number", r"Hatta number: ([\d.]+)", "31.3"),
    (1, "cell area [m2]", r"Cell area required = ([\d.]+) m2", "23479.60"),
    (1, "cell power [kW]", r"Cell power = ([\d.]+) kW", "66759.4"),
    (1, "cell energy [kWh/t]", r"Cell energy intensity = ([\d.]+) kWh/tCO2", "1007.3"),
    (10, "coarse LCOC [USD/t]", r"LCOC:\s+\$([\d.]+) / tCO2", "125.82"),
    (10, "coarse diameter [m]", r"Diameter:\s+([\d.]+) m", "5.84"),
    (10, "coarse height [m]", r"Height:\s+([\d.]+) m", "19.78"),
    (12, "refined LCOC [USD/t]", r"LCOC:\s+\$([\d.]+) / tCO2", "125.79"),
]

# Cell 2 from a clean kernel (saved .ipynb says LCOC 123.20, OPEX 64,911,181, 535,822,105 kWh).
CLEAN_KERNEL_CELL2 = [
    ("LCOC [USD/t]", r"LCOC:\s+\$([\d.]+) / tCO2 captured", "123.17"),
    ("annual OPEX [USD/yr]", r"Annual OPEX:\s+\$([\d,]+) / yr", "64,893,539"),
    ("annual electricity [kWh/yr]", r"Annual electricity:\s+([\d,]+) kWh / yr", "535,675,091"),
]


def _stdout(cell) -> str:
    return "".join(
        "".join(o["text"])
        for o in cell.get("outputs", [])
        if o["output_type"] == "stream" and o["name"] == "stdout"
    )


def _value(text: str, pattern: str) -> str:
    match = re.search(pattern, text)
    assert match, f"pattern {pattern!r} not found"
    return match.group(1)


@pytest.fixture(scope="module")
def saved():
    return json.loads(LEGACY_NOTEBOOK.read_text(encoding="utf-8"))["cells"]


@pytest.fixture(scope="module")
def executed(tmp_path_factory):
    workdir = tmp_path_factory.mktemp("legacy_run")
    copy = workdir / LEGACY_NOTEBOOK.name
    shutil.copyfile(LEGACY_NOTEBOOK, copy)
    nb = nbformat.read(copy, as_version=4)
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("MPLBACKEND", "Agg")
        client = NotebookClient(
            nb,
            timeout=1800,
            kernel_name="python3",
            allow_errors=True,  # cell 6 raised KeyError in the saved run too
            resources={"metadata": {"path": str(workdir)}},
        )
        start = time.perf_counter()
        client.execute()
        elapsed = time.perf_counter() - start
    print(f"\nlegacy notebook executed in {elapsed:.1f} s")
    return nb.cells


def test_cell_layout_matches(saved, executed):
    assert [c["cell_type"] for c in executed] == [c["cell_type"] for c in saved]


def test_stdout_identical_except_stale_state_cell(saved, executed):
    for i, (old, new) in enumerate(zip(saved, executed, strict=True)):
        if old["cell_type"] != "code" or i == CELL_STALE_STATE:
            continue
        assert _stdout(new) == _stdout(old), f"cell {i} stdout differs from saved output"


def test_saved_error_is_reproduced(saved, executed):
    def errors(cell):
        return [
            (o["ename"], o["evalue"])
            for o in cell.get("outputs", [])
            if o["output_type"] == "error"
        ]

    assert errors(saved[CELL_SAVED_ERROR]) == [("KeyError", "'T_gas_abs_C'")]
    assert errors(executed[CELL_SAVED_ERROR]) == errors(saved[CELL_SAVED_ERROR])
    others = [i for i, c in enumerate(executed) if i != CELL_SAVED_ERROR and errors(c)]
    assert others == []


@pytest.mark.parametrize(
    ("cell", "label", "pattern", "expected"),
    SAVED_VALUES,
    ids=[f"cell{c}-{label}" for c, label, _, _ in SAVED_VALUES],
)
def test_key_values(saved, executed, cell, label, pattern, expected):
    assert _value(_stdout(saved[cell]), pattern) == expected, "hand-written value vs .ipynb"
    assert _value(_stdout(executed[cell]), pattern) == expected


@pytest.mark.parametrize(
    ("label", "pattern", "expected"),
    CLEAN_KERNEL_CELL2,
    ids=[label for label, _, _ in CLEAN_KERNEL_CELL2],
)
def test_cell2_clean_kernel_values(executed, label, pattern, expected):
    assert _value(_stdout(executed[CELL_STALE_STATE]), pattern) == expected


@pytest.mark.xfail(
    strict=True,
    reason="audit 6.7: saved cell 2 (LCOC 123.20) used stale kernel globals from before "
    "In[91]; a clean kernel gives 123.17",
)
def test_cell2_saved_output_reproduces(saved, executed):
    assert _stdout(executed[CELL_STALE_STATE]) == _stdout(saved[CELL_STALE_STATE])


def test_figures_stay_out_of_repo(executed):
    repo = LEGACY_NOTEBOOK.parents[1]
    for folder in (repo, repo / "legacy", repo / "tests"):
        assert not any(folder.glob("*.png")), f"legacy savefig wrote into {folder}"
