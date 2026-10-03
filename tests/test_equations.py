"""The equation catalog is complete, consistent and exported (CLAUDE.md rule 9)."""

import importlib
import inspect
import json
import re
from pathlib import Path

import pytest

from ccs_predesign import equations
from ccs_predesign.equations import EQUATIONS, REFERENCES, SYMBOLS

EXPORTED = Path(__file__).resolve().parents[1] / "equations.json"
PHYSICS_MODULES = [
    "gas",
    "solvent",
    "kinetics",
    "packing",
    "absorber",
    "cell",
    "tea",
    "optimize",
    "pipeline",
]
# Public helpers that are not equations (unit conversion, bookkeeping, orchestration).
NOT_EQUATIONS = {
    "gas:c_to_k",
    "gas:mean_molar_mass_kg_mol",
    "packing:reynolds",
    "packing:schmidt",
    "solvent:naoh_ions_kmol_m3",
    "tea:safe_div",
    "tea:evaluate_design_point",
    "absorber:design_absorber_by_grid",
    "optimize:evaluate_grid",
    "optimize:coarse_grids",
    "optimize:fine_grids",
    "optimize:optimize_two_stage",
}


def _function(ref: str):
    module, name = ref.split(":")
    return getattr(importlib.import_module(f"ccs_predesign.{module}"), name)


@pytest.mark.parametrize("eq", EQUATIONS.values(), ids=EQUATIONS.keys())
def test_equation_entry_is_consistent(eq):
    assert re.fullmatch(r"[a-z]+\.[a-z0-9_]+", eq.id)
    assert eq.latex.count("{") == eq.latex.count("}")
    assert eq.meaning_es and eq.meaning_en and eq.title_es and eq.title_en
    assert all(s in SYMBOLS for s in eq.symbols), [s for s in eq.symbols if s not in SYMBOLS]
    assert eq.reference is None or eq.reference in REFERENCES
    if eq.provenance == "literature":
        assert eq.reference is not None
    for ref in eq.implemented_by:
        doc = inspect.getdoc(_function(ref)) or ""
        assert eq.id in doc, f"{ref} must cite {eq.id} in its docstring"


def test_every_physics_function_cites_an_equation():
    cited = {ref for eq in EQUATIONS.values() for ref in eq.implemented_by}
    missing = []
    for module in PHYSICS_MODULES:
        mod = importlib.import_module(f"ccs_predesign.{module}")
        for name, fn in inspect.getmembers(mod, inspect.isfunction):
            ref = f"{module}:{name}"
            if (
                fn.__module__ == mod.__name__
                and not name.startswith("_")
                and ref not in NOT_EQUATIONS
            ):
                if ref not in cited:
                    missing.append(ref)
    assert missing == []


def test_exported_catalog_is_up_to_date():
    """Regenerate with `uv run python -m ccs_predesign.equations equations.json`."""
    assert json.loads(EXPORTED.read_text(encoding="utf-8")) == json.loads(
        json.dumps(equations.catalog(), ensure_ascii=False)
    )


def test_errata_ids_exist():
    errata = (Path(__file__).resolve().parents[1] / "docs" / "errata.md").read_text(
        encoding="utf-8"
    )
    for eq in EQUATIONS.values():
        for e_id in eq.errata:
            assert f"| {e_id} |" in errata, (eq.id, e_id)
