"""Refactor safety net for legacy cells 2, 8, 10, 12 and 17 (TEA, optimiser, bench).

The legacy cells run in notebook order in one namespace, so the notebook's state bugs are
reproduced too (cell 2 sets the gas flow to 300 000 Nm3/h for cells 8–12, audit §7).
"""

import numpy as np
import pytest
from legacy_cells import CELL_FIELDS, assert_absorber_equal, exec_cell, load_notebook_state

from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.models import AbsorberSpec, CaseInput, FlueGas, Solvent, TEASpec
from ccs_predesign.optimize import choose_best, coarse_grids, evaluate_grid, fine_grids
from ccs_predesign.pipeline import run_case
from ccs_predesign.tea import evaluate_design_point, tea_summary

pytestmark = pytest.mark.slow

# The optimisation cells ran with the flow left behind by cell 2 (audit §7).
CASE_300K = CaseInput(gas=FlueGas(Q_dry_Nm3_h=300_000.0))

POINT_FIELDS = {
    "NaOH_conc_M": "NaOH_M",
    "LG_vol": "LG_vol",
    "P_total_W": "P_total_W",
    "CO2_captured_tpy": "CO2_captured_t_y",
    "CO2_stack_tpy": "CO2_stack_t_y",
    "CO2_product_tpy": "CO2_product_t_y",
    "E_cell_kWh_tCO2": "E_cell_kWh_t",
    "E_total_kWh_tCO2": "E_total_kWh_t",
    "capex_col_usd": "capex_column_usd",
    "capex_total_usd": "capex_total_usd",
    "capex_cell_usd": "capex_cell_usd",
    "capex_blower_usd": "capex_blower_usd",
    "capex_pump_usd": "capex_pump_usd",
    "annualized_capex_usd_per_year": "annualized_capex_usd_y",
    "opex_fixed_usd_per_year": "opex_fixed_usd_y",
    "opex_electricity_usd_per_year": "opex_electricity_usd_y",
    "opex_solvent_usd_per_year": "opex_solvent_usd_y",
    "opex_water_chem_usd_per_year": "opex_water_chem_usd_y",
    "opex_total_usd_per_year": "opex_total_usd_y",
    "LCOC_usd_per_tCO2": "LCOC_usd_t",
    "indirect_kgCO2e_per_year": "indirect_kgCO2e_y",
    "indirect_tCO2e_per_year": "indirect_tCO2e_y",
    "indirect_kgCO2e_per_tCO2": "indirect_kgCO2e_t",
}


def assert_point_equal(old, new) -> None:
    for old_key, new_key in POINT_FIELDS.items():
        assert getattr(new, new_key) == old[old_key], old_key
    assert_absorber_equal(old, new.absorber, skip=("total_hydraulic_power_W",))
    for old_key, new_key in CELL_FIELDS.items():
        assert getattr(new.cell, new_key) == old[old_key], old_key


@pytest.fixture(scope="module")
def ns():
    return load_notebook_state()


def test_notebook_state_has_the_300k_flow(ns):
    assert ns["Qg_Nm3h"] == 300_000.0
    assert ns["H_packed_m"] == 25.0  # cell 2 alias list misses "height_m" (audit §7)


def test_cell2_tea_summary(ns):
    base = run_case(CaseInput())
    best = base.absorber.best
    s = tea_summary(
        D_m=best.D_col_m,
        H_packed_m=25.0,
        P_bar=1.01,
        CO2_captured_mol_s=best.CO2_captured_mol_s,
        P_total_W=base.P_total_W,
        tea=TEASpec(),
    )
    assert ns["D_col_m"] == best.D_col_m
    assert s.annual_captured_t_y == ns["annual_captured_tCO2"]
    assert s.annual_electricity_kWh_y == ns["annual_electricity_kWh"]
    assert s.column.installed_usd == ns["capex"]["installed_cost_usd"]
    assert s.capex_total_usd == ns["total_installed_capex_usd"]
    assert s.opex_total_usd_y == ns["annual_opex_usd"]
    assert s.annualized_capex_usd_y == ns["annualized_capex_usd_per_year"]
    assert s.LCOC_usd_t == ns["lcoc_usd_per_tCO2"]
    assert s.indirect_tCO2e_y == ns["annual_indirect_emissions_tCO2e"]
    assert s.indirect_kgCO2e_t == ns["indirect_emissions_intensity_kg_per_tCO2"]
    assert f"{s.LCOC_usd_t:,.2f}" == "123.17"  # clean kernel; the saved 123.20 is stale state


@pytest.mark.parametrize(("NaOH_M", "LG_vol"), [(1.5, 0.005), (0.25, 0.07), (4.0, 0.004)])
def test_design_point(ns, NaOH_M, LG_vol):
    old = ns["evaluate_design_point"](NaOH_conc_M=NaOH_M, LG_vol=LG_vol)
    new = evaluate_design_point(
        CASE_300K.gas,
        Solvent(NaOH_M=NaOH_M),
        CASE_300K.absorber,
        CASE_300K.cell,
        CASE_300K.tea,
        LG_vol,
    )
    assert_point_equal(old, new)


@pytest.fixture(scope="module")
def optimisation(ns):
    exec_cell(ns, 10)
    exec_cell(ns, 12)
    coarse = choose_best(evaluate_grid(CASE_300K, *coarse_grids()))
    fine = choose_best(
        evaluate_grid(CASE_300K, *fine_grids(coarse.best.NaOH_M, coarse.best.LG_vol))
    )
    return ns, coarse, fine


@pytest.mark.parametrize("stage", ["coarse", "fine"])
def test_optimisation_stage(optimisation, stage):
    ns, coarse, fine = optimisation
    df, ranked, best_old = (
        (ns["coarse_ok"], ns["coarse_ranked"], ns["best_coarse"])
        if stage == "coarse"
        else (ns["fine_ok"], ns["fine_ranked"], ns["best_fine"])
    )
    new = coarse if stage == "coarse" else fine
    assert "error" not in df.columns
    assert len(df) == len(new.points)
    for (_, row), point in zip(df.iterrows(), new.points, strict=True):
        assert_point_equal(row, point)
    assert np.array_equal(ranked["score"].to_numpy(), np.array(new.scores))
    assert_point_equal(best_old, new.best)


def test_optimisation_printed_values(optimisation):
    _, coarse, fine = optimisation
    c, f = coarse.best, fine.best
    assert f"{c.NaOH_M:.3f}" == "4.000"
    assert f"{c.LG_vol:.4f}" == "0.0040"
    assert f"{c.LCOC_usd_t:.2f}" == "125.82"
    assert f"{c.absorber.D_col_m:.2f}" == "5.84"
    assert f"{c.absorber.height_m:.2f}" == "19.78"
    assert f"{c.capex_total_usd:,.0f}" == "24,484,198"
    assert f"{f.NaOH_M:.3f}" == "4.367"
    assert f"{f.LG_vol:.4f}" == "0.0024"
    assert f"{f.LCOC_usd_t:.2f}" == "125.79"
    assert f"{f.absorber.height_m:.2f}" == "19.88"


@pytest.mark.parametrize("H_real_m", [0.15, 0.35])
@pytest.mark.parametrize("NaOH_M", [0.10, 0.80, 1.50])
def test_bench_cell17(ns, H_real_m, NaOH_M):
    """Bench runs: fixed 15.4 mm column at its physical height (audit §8)."""
    exec_cell(ns, 17)
    kwargs = dict(
        Qg_Nm3h=ns["Q_total_Nm3h"],
        yCO2_dry_in=ns["yCO2_banco"],
        yO2_dry_in=0.0,
        T_gas_abs_C=25.0,
        P_abs_bar=1.01,
        capture_target=0.999999,
        NaOH_conc_M=NaOH_M,
        LG_vol=ns["LG_banco"],
        packing_name="pall_ring_25mm",
        H_eq=0.04,
        max_height_m=H_real_m,
        dz=0.02,
        D_col_m=0.0154,
    )
    old = ns["simulate_absorber_fixed_diameter"](**kwargs)
    new = simulate_absorber(
        FlueGas(
            Q_dry_Nm3_h=kwargs["Qg_Nm3h"], y_CO2_dry=kwargs["yCO2_dry_in"], y_O2_dry=0.0, T_C=25.0
        ),
        Solvent(NaOH_M=NaOH_M),
        AbsorberSpec(capture_target=0.999999, max_height_m=H_real_m),
        kwargs["LG_vol"],
        D_col_fixed_m=0.0154,
    )
    assert_absorber_equal(old, new, skip=("total_hydraulic_power_W",))
    assert new.flood_fraction_actual == old["flood_fraction_actual"]
    assert new.feasible_hydraulically == old["feasible_hydraulically"]
