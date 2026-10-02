"""Refactor safety net: ccs_predesign must equal legacy cell 1 bit-for-bit.

Each case overrides legacy module-level inputs, runs the legacy run_case() and the
package's run_case() with the equivalent CaseInput, and compares every output with ==.
"""

import numpy as np
import pytest
from legacy_cells import (
    CELL_FIELDS,
    assert_absorber_equal,
    load_cell1,
    run_legacy_case,
)

from ccs_predesign.models import AbsorberSpec, CaseInput, CellSpec, FlueGas, Solvent
from ccs_predesign.pipeline import run_case

LEGACY_DEFAULTS = {
    "Qg_Nm3h": 150000.0,
    "yCO2_dry_in": 0.25,
    "yO2_dry_in": 0.03,
    "T_gas_abs_C": 45.0,
    "P_abs_bar": 1.01,
    "capture_target": 0.90,
    "NaOH_conc_M": 1.5,
    "H_eq": 0.04,
    "liquid_resistance_factor": 0.001,
    "k2_CO2_OH_M": 8.5e3,
    "packing_name": "pall_ring_25mm",
    "max_column_height_m": 30.0,
    "design_flood_fraction": 0.60,
    "j_cell_mAcm2": 200.0,
    "cell_voltage_override_V": None,
    "tna_mode": "empirical",
    "tna_const": 0.90,
    "CO2_release_eff": 0.98,
    "pump_eff": 0.70,
    "blower_eff": 0.68,
    "dz_m": 0.02,
    "L_over_G_grid": np.linspace(0.005, 0.060, 28),
}

CASES = {
    "base": {},
    "naoh_3M": {"NaOH_conc_M": 3.0},
    "big_dilute_gas": {"Qg_Nm3h": 300000.0, "yCO2_dry_in": 0.15, "yO2_dry_in": 0.05},
    "structured_250Y": {"packing_name": "structured_250Y"},
    "raschig_capture80": {"packing_name": "raschig_metal_25mm", "capture_target": 0.80},
    "random_hc_flood70": {"packing_name": "random_high_capacity", "design_flood_fraction": 0.7},
    "cell_const_override": {
        "tna_mode": "const",
        "tna_const": 0.85,
        "cell_voltage_override_V": 2.0,
        "j_cell_mAcm2": 100.0,
        "CO2_release_eff": 0.95,
    },
    "cold_pressurised": {"T_gas_abs_C": 30.0, "P_abs_bar": 1.2, "k2_CO2_OH_M": 1.0e4},
    "infeasible_liquid_control": {"liquid_resistance_factor": 1.0},
    "steep_equilibrium": {"H_eq": 0.5, "dz_m": 0.05, "L_over_G_grid": np.linspace(0.01, 0.1, 7)},
    "pump_blower_eff": {"pump_eff": 0.5, "blower_eff": 0.8, "max_column_height_m": 10.0},
}

CASE_FIELDS = {
    "P_total_W": "P_total_W",
    "CO2_captured_tpy": "CO2_captured_t_y",
    "CO2_stack_tpy": "CO2_stack_t_y",
    "CO2_product_tpy": "CO2_product_t_y",
    "E_cell_kWh_tCO2": "E_cell_kWh_t",
}


def case_input(g: dict) -> CaseInput:
    """Map legacy module-level inputs to the package's CaseInput."""
    return CaseInput(
        gas=FlueGas(
            Q_dry_Nm3_h=g["Qg_Nm3h"],
            y_CO2_dry=g["yCO2_dry_in"],
            y_O2_dry=g["yO2_dry_in"],
            T_C=g["T_gas_abs_C"],
            P_bar=g["P_abs_bar"],
        ),
        solvent=Solvent(NaOH_M=g["NaOH_conc_M"]),
        absorber=AbsorberSpec(
            capture_target=g["capture_target"],
            packing_name=g["packing_name"],
            H_eq=g["H_eq"],
            liquid_resistance_factor=g["liquid_resistance_factor"],
            k2_ref_L_mol_s=g["k2_CO2_OH_M"],
            max_height_m=g["max_column_height_m"],
            dz_m=g["dz_m"],
            flood_fraction=g["design_flood_fraction"],
            LG_grid_vol=tuple(float(v) for v in g["L_over_G_grid"]),
            pump_eff=g["pump_eff"],
            blower_eff=g["blower_eff"],
        ),
        cell=CellSpec(
            j_mA_cm2=g["j_cell_mAcm2"],
            V_override_V=g["cell_voltage_override_V"],
            tna_mode=g["tna_mode"],
            tna_const=g["tna_const"],
            CO2_release_eff=g["CO2_release_eff"],
        ),
    )


@pytest.fixture(scope="module")
def legacy_ns():
    return load_cell1()


def test_defaults_match_legacy_inputs(legacy_ns):
    legacy_inputs = {k: legacy_ns[k] for k in LEGACY_DEFAULTS}
    for key, value in LEGACY_DEFAULTS.items():
        assert (
            np.array_equal(legacy_inputs[key], value)
            if key == "L_over_G_grid"
            else (legacy_inputs[key] == value)
        ), key
    assert CaseInput() == case_input(LEGACY_DEFAULTS)


@pytest.mark.parametrize("overrides", CASES.values(), ids=CASES.keys())
def test_run_case_matches_legacy(legacy_ns, overrides):
    inputs = {**LEGACY_DEFAULTS, **overrides}
    best, candidates = run_legacy_case(legacy_ns, **inputs)
    result = run_case(case_input(inputs))

    assert result.absorber.feasible == best["feasible"]
    assert_absorber_equal(best, result.absorber.best)
    assert len(result.absorber.candidates) == len(candidates)
    for old, new in zip(candidates, result.absorber.candidates, strict=True):
        assert_absorber_equal(old, new)
    for old_key, new_key in CELL_FIELDS.items():
        assert getattr(result.cell, new_key) == best["regen"][old_key], old_key
    for old_key, new_key in CASE_FIELDS.items():
        assert getattr(result, new_key) == best[old_key], old_key


def test_base_case_printed_values():
    """Second check of the saved cell-1 printout, from the package alone."""
    r = run_case(CaseInput())
    b = r.absorber.best
    assert f"{b.height_m:.2f}" == "14.86"
    assert f"{b.deltaP_Pa:.0f}" == "1919"
    assert f"{b.blower_power_W / 1000:.1f}" == "151.8"
    assert f"{b.pump_power_W / 1000:.1f}" == "48.2"
    assert f"{b.KGa_1_s:.4f}" == "323.4097"
    assert f"{b.Ha:.1f}" == "31.3"
    assert f"{b.D_col_m:.2f}" == "4.19"
    assert f"{r.cell.A_cell_m2:.2f}" == "23479.60"
    assert f"{r.cell.P_cell_W / 1000:.1f}" == "66759.4"
    assert f"{r.P_total_W / 1000:.1f}" == "66959.4"
    assert f"{r.E_cell_kWh_t:.1f}" == "1007.3"
