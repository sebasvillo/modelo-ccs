"""Snapshot of the model's headline outputs, tied to the latest errata entry.

Regenerate only together with a new docs/errata.md entry:

    uv run python tests/golden/snapshot.py E-00N
"""

import json
import re
import sys
from pathlib import Path

from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.models import AbsorberSpec, CaseInput, CellSpec, FlueGas, Solvent, TEASpec
from ccs_predesign.optimize import optimize_two_stage
from ccs_predesign.pipeline import run_case
from ccs_predesign.tea import design_costs, evaluate_design_point

ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = Path(__file__).with_name("snapshots") / "model_outputs.json"
ERRATA = ROOT / "docs" / "errata.md"

# Bench of thesis Ch. 3 as simulated in legacy cell 17 (audit §8: 15.4 mm column).
BENCH_Q_Nm3_h = (3100.0 / 1e6) * 60.0 * (273.15 / 298.15)
BENCH_Y_CO2 = 500.0 / 3100.0
BENCH_LG = 383.0 / 3100.0
BENCH_RUNS = [(H, C) for H in (0.15, 0.35) for C in (0.10, 0.80, 1.50)]


def latest_errata_id() -> str:
    text = re.sub(r"<!--.*?-->", "", ERRATA.read_text(encoding="utf-8"), flags=re.S)
    ids = re.findall(r"^\|\s*(E-\d{3})\s*\|", text, re.M)
    ids = [i for i in ids if i != "E-000"]
    return ids[-1] if ids else "none"


def _absorber(r) -> dict:
    return {
        "LG_vol": r.LG_vol,
        "height_m": r.height_m,
        "D_col_m": r.D_col_m,
        "capture_achieved": r.capture_achieved,
        "reached_target": r.reached_target,
        "KGa_1_s": r.KGa_1_s,
        "kG_m_s": r.kG_m_s,
        "kL_m_s": r.kL_m_s,
        "Ha": r.Ha,
        "E": r.E,
        "deltaP_Pa": r.deltaP_Pa,
        "blower_power_W": r.blower_power_W,
        "pump_power_W": r.pump_power_W,
        "CO2_captured_mol_s": r.CO2_captured_mol_s,
    }


def _point(p) -> dict:
    keys = ["NaOH_M", "LG_vol", "LCOC_usd_t", "capex_total_usd", "opex_total_usd_y"]
    return {k: getattr(p, k) for k in keys} | {
        "height_m": p.absorber.height_m,
        "D_col_m": p.absorber.D_col_m,
        "capture_achieved": p.absorber.capture_achieved,
    }


def compute_snapshot() -> dict:
    base = run_case(CaseInput())
    best = base.absorber.best
    costs = design_costs(
        best.D_col_m,
        best.height_m,
        best.blower_power_W,
        best.pump_power_W,
        best.CO2_captured_mol_s,
        best.CO2_out_mol_s,
        base.cell,
        CaseInput().gas.P_bar,
        TEASpec(),
    )
    coarse, fine = optimize_two_stage(CaseInput())
    structured = CaseInput(absorber=AbsorberSpec(packing_name="structured_250Y"))
    s_coarse, s_fine = optimize_two_stage(structured)
    bench = []
    for H, C in BENCH_RUNS:
        r = simulate_absorber(
            FlueGas(Q_dry_Nm3_h=BENCH_Q_Nm3_h, y_CO2_dry=BENCH_Y_CO2, y_O2_dry=0.0, T_C=25.0),
            Solvent(NaOH_M=C),
            AbsorberSpec(capture_target=0.999999, max_height_m=H),
            BENCH_LG,
            D_col_fixed_m=0.0154,
        )
        bench.append({"H_bed_m": H, "NaOH_M": C} | _absorber(r))
    default_case = CaseInput()
    j_opt = evaluate_design_point(
        default_case.gas,
        Solvent(NaOH_M=2.5),
        default_case.absorber,
        CellSpec(optimize_j=True, j_max_mA_cm2=1000.0),
        TEASpec(),
        0.0093,
    )
    return {
        "current_density_optimum": {
            "j_mA_cm2": j_opt.cell.j_A_m2 / 10.0,
            "V_cell_V": j_opt.cell.V_cell_V,
            "A_cell_m2": j_opt.cell.A_cell_m2,
            "LCOC_usd_t": j_opt.LCOC_usd_t,
        },
        "base_case": _absorber(best)
        | {
            "feasible": base.absorber.feasible,
            "A_cell_m2": base.cell.A_cell_m2,
            "P_cell_W": base.cell.P_cell_W,
            "P_total_W": base.P_total_W,
            "E_cell_kWh_t": base.E_cell_kWh_t,
            "LCOC_usd_t": costs["LCOC_usd_t"],
            "capex_total_usd": costs["capex_total_usd"],
        },
        "optimisation_coarse": _point(coarse.best) | {"feasible": coarse.feasible},
        "optimisation_fine": _point(fine.best) | {"feasible": fine.feasible},
        "optimisation_structured_250Y": _point(s_fine.best)
        | {"feasible": s_fine.feasible, "n_feasible_coarse": s_coarse.n_feasible},
        "bench": bench,
    }


def main(errata_id: str) -> None:
    if errata_id != latest_errata_id():
        sys.exit(f"{errata_id} is not the latest errata id ({latest_errata_id()})")
    data = {"errata_id": errata_id, "outputs": compute_snapshot()}
    SNAPSHOT.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {SNAPSHOT.relative_to(ROOT)} for {errata_id}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "none")
