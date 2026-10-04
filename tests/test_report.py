"""Explained outputs are truthful: re-evaluating each equation with its inputs gives the value."""

import math

import pytest

from ccs_predesign.models import AbsorberSpec, CaseInput, Solvent
from ccs_predesign.pipeline import run_case
from ccs_predesign.report import explain_case


def Y(y):
    """Solute-free mole ratio (E-018)."""
    return y / (1 - y)


FORMULAS = {
    "height_m": lambda v: v["NTU"] * v["G_flux"] / (v["c_tot"] * v["KGa_mean"]),
    "capture_fraction": lambda v: 1 - Y(v["y_top"]) / Y(v["y"]),
    "diameter_m": lambda v: math.sqrt(
        4.0 * v["Q_G"] / (v["f"] * v["v_flood"]) / (v["n_trains"] * math.pi)
    ),
    "KGa_inlet_1_s": lambda v: v["a_e"] / (1 / v["kG"] + v["H_cc"] / (v["E"] * v["kL"])),
    "hatta_inlet": lambda v: math.sqrt(v["k2"] * v["OH"] * v["D_L"]) / v["kL"],
    "k2_inlet_m3_mol_s": lambda v: (
        10 ** (11.895 - 2382 / v["T"] + 0.221 * v["I"] - 0.016 * v["I"] ** 2) / 1000
    ),
    "enhancement_inlet": lambda v: (
        -(v["Ha"] ** 2) / (2 * (v["E_inf"] - 1))
        + math.sqrt(
            v["Ha"] ** 4 / (4 * (v["E_inf"] - 1) ** 2)
            + v["E_inf"] * v["Ha"] ** 2 / (v["E_inf"] - 1)
            + 1
        )
    ),
    "henry_cc_inlet": lambda v: 1 / (v["H_cp"] * v["R"] * v["T"]),
    "OH_rich_mol_m3": lambda v: (
        v["OH"] - v["nu"] * v["G_I"] * (Y(v["y"]) - Y(v["y_top"])) / v["Q_L"]
    ),
    "T_rich_K": lambda v: (
        v["T_L"] + (-v["dH"]) * v["G_I"] * (Y(v["y"]) - Y(v["y_top"])) / (v["m_L"] * v["cp_L"])
    ),
    "pump_power_W": lambda v: v["rho_L"] * v["g"] * (v["H"] + v["H_extra"]) * v["Q_L"] / v["eta"],
    "blower_power_W": lambda v: v["dP"] * v["Q_G"] / v["eta"],
    "cell_current_A": lambda v: v["n_e"] * v["F"] * v["n_dot"] / v["t_Na"],
    "cell_area_m2": lambda v: v["I_cell"] / v["j"],
    "cell_voltage_V": lambda v: 3.02912 * math.log((v["j_mA"] + 767.22373) / 604.92445),
    "transference_number": lambda v: (
        0.70 + 0.20 * (1 - math.exp(-4 * v["C"])) - 0.05 * min(v["j_mA"], 250) / 250
    ),
    "cell_energy_J_mol": lambda v: v["n_e"] * v["F"] * v["V"] / v["t_Na"],
    "CO2_captured_t_y": lambda v: v["n_dot"] * 44.01e-6 * 3600 * v["h_op"],
    "crf_1_y": lambda v: v["i"] * (1 + v["i"]) ** v["n"] / ((1 + v["i"]) ** v["n"] - 1),
    "cell_capex_usd": lambda v: (v["c_stack"] + v["c_bop"]) * v["A_cell"] * (1 + v["f_u"]),
    "LCOC_usd_t": lambda v: (v["CAPEX"] * v["CRF"] + v["OPEX"]) / v["m_CO2"],
    "indirect_kgCO2e_t": lambda v: v["h_op"] * v["W"] * v["EF"] / (1000 * v["m_CO2"]),
}


@pytest.mark.parametrize(
    "inp",
    [
        CaseInput(),
        CaseInput(solvent=Solvent(NaOH_M=3.0, T_in_C=30.0)),
        CaseInput(absorber=AbsorberSpec(max_diameter_m=5.0)),
    ],
    ids=["base", "3M_cold", "trains"],
)
def test_substituted_values_reproduce_each_output(inp):
    outputs = explain_case(inp, run_case(inp))
    assert {o.key for o in outputs} == set(FORMULAS)
    for o in outputs:
        assert FORMULAS[o.key](o.inputs) == pytest.approx(o.value, rel=1e-6), o.key
