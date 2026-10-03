"""Component-level safety net: each ported function equals its legacy twin bit-for-bit.

Unlike the end-to-end equivalence tests, these keep holding while physics fixes land
elsewhere. When a fix changes one of these functions, its test here is replaced by a test
tied to the errata entry (docs/errata.md).
"""

import itertools

import numpy as np
import pandas as pd
import pytest
from legacy_cells import CELL_FIELDS, load_notebook_state

from ccs_predesign import cell, gas, kinetics, packing, solvent, tea
from ccs_predesign.models import CellResult, CellSpec, TEASpec
from ccs_predesign.optimize import LEGACY_WEIGHTS, weighted_scores

T_K = [283.15, 298.15, 318.15, 333.15]
P_PA = [0.9e5, 1.01e5, 1.5e5]
NAOH_M = [0.1, 1.5, 4.0]
PACKING_NAMES = list(packing.PACKINGS)
COMP = {"CO2": 0.2, "O2": 0.03, "N2": 0.69, "H2O": 0.08}


@pytest.fixture(scope="module")
def ns():
    return load_notebook_state()


def legacy_packing(ns, name):
    return ns["PACKINGS"][name]


def test_packing_database(ns):
    assert set(ns["PACKINGS"]) == set(packing.PACKINGS)
    for name, old in ns["PACKINGS"].items():
        new = packing.PACKINGS[name]
        assert new.a_spec_m2_m3 == old["a_spec_m2m3"]
        assert new.void_fraction == old["void_fraction"]
        assert new.dp_eq_m == old["dp_eq_m"]
        assert new.wetting_ref_uL_m_s == old["wetting_ref_uLm_s"]
        assert new.flood_coeff == old["flood_coeff"]


@pytest.mark.parametrize(("T", "P"), list(itertools.product(T_K, P_PA)))
def test_gas_properties(ns, T, P):
    assert gas.sat_vapor_pressure_water_Pa(T) == ns["sat_vapor_pressure_water_Pa"](T)
    assert gas.density_ideal_kg_m3(T, P, COMP) == ns["gas_density_ideal"](T, P, COMP)
    assert gas.viscosity_sutherland_air_Pa_s(T) == ns["gas_viscosity_sutherland_air"](T)
    assert gas.diffusivity_CO2_in_air_m2_s(T, P) == ns["gas_diffusivity_CO2_in_air"](T, P)
    for y_CO2, y_O2 in [(0.25, 0.03), (0.15, 0.05), (0.16, 0.0)]:
        new = gas.wet_composition_from_dry(y_CO2, y_O2, T, P)
        assert new == ns["wet_gas_composition_from_dry"](y_CO2, y_O2, T, P)
    for Q in [0.17, 150_000.0, 300_000.0]:
        assert gas.dry_molar_flow_mol_s(Q) == ns["dry_mole_flows_from_Nm3h"](Q)


@pytest.mark.parametrize(("T", "C"), list(itertools.product(T_K, NAOH_M)))
def test_solvent_and_kinetics(ns, T, C):
    assert solvent.viscosity_Pa_s(T, C) == ns["liquid_viscosity_Pa_s"](T, C)
    assert solvent.diffusivity_CO2_m2_s(T, C) == ns["liquid_diffusivity_CO2_m2s"](T, C)
    assert solvent.density_kg_m3(T, C) == ns["liquid_density_kgm3"](T, C)
    assert kinetics.k2_m3_mol_s(T, 8.5e3) == ns["reaction_rate_constant_m3mol_s"](T)


@pytest.mark.parametrize("Ha", [1e-6, 0.3, 2.0, 31.3, 49.9, 50.0, 400.0])
def test_enhancement_factor(ns, Ha):
    assert kinetics.enhancement_factor(Ha) == ns["enhancement_factor_ha"](Ha)


@pytest.mark.parametrize("name", PACKING_NAMES)
def test_packing_correlations(ns, name):
    old_pk, new_pk = legacy_packing(ns, name), packing.PACKINGS[name]
    for Re, Sc in itertools.product([1.0, 150.0, 2.0e5], [0.7, 1.2, 500.0]):
        assert packing.sherwood_gas(Re, Sc) == ns["sherwood_gas"](Re, Sc)
        assert packing.sherwood_liquid(Re, Sc) == ns["sherwood_liquid"](Re, Sc)
    for uL in [1e-5, 1e-3, 0.05]:
        assert packing.wetting_fraction(uL, new_pk) == ns["wetting_fraction"](uL, old_pk)
    for rho_g, mu_g, uG in [(1.1, 1.9e-5, 1.5), (1.3, 2.0e-5, 53.8)]:
        assert packing.pressure_drop_ergun_Pa_m(rho_g, mu_g, uG, new_pk) == ns[
            "pressure_drop_Ergun_per_m"
        ](rho_g, mu_g, uG, old_pk)
        assert packing.flooding_velocity_m_s(rho_g, 1030.0, new_pk) == ns[
            "flooding_velocity_simple"
        ](rho_g, 1030.0, old_pk)
        for Q in [5e-5, 53.8]:
            assert packing.column_diameter_m(Q, rho_g, 1030.0, new_pk, 0.6) == ns[
                "choose_diameter"
            ](Q, rho_g, 1030.0, old_pk, flood_frac=0.6)


@pytest.mark.parametrize(
    ("name", "T", "C", "E"),
    [
        (n, T, C, E)
        for n in PACKING_NAMES
        for T in (298.15, 318.15)
        for C in (0.1, 1.5)
        for E in (1.0, 31.3)
    ],
)
def test_mass_transfer_coefficients(ns, name, T, C, E):
    """Legacy passed flows over A_ref = 1 m2, i.e. velocities numerically equal to the flows."""
    old = ns["kL_kG_KGa"](T, 1.01e5, COMP, 53.8, 0.27, legacy_packing(ns, name), C, E=E)
    new = packing.mass_transfer_coefficients(
        T, 1.01e5, COMP, 53.8, 0.27, packing.PACKINGS[name], C, 0.001, E=E
    )
    assert new == old


@pytest.mark.parametrize(
    ("CO2_mol_s", "j", "C", "mode", "V"),
    [
        (418.3, 200.0, 1.5, "empirical", None),
        (2.2e-3, 50.0, 0.1, "empirical", None),
        (800.0, 300.0, 4.0, "const", 2.1),
    ],
)
def test_cell(ns, CO2_mol_s, j, C, mode, V):
    j_A_m2 = j * 10.0
    assert cell.tna_empirical(C * 1000.0, j_A_m2) == ns["tna_empirical"](C * 1000.0, j_A_m2)
    assert cell.voltage_from_current_density_V(j_A_m2) == ns["voltage_from_current_density"](j_A_m2)
    old = ns["electrochemical_regeneration"](
        CO2_capture_mol_s=CO2_mol_s,
        j_cell_mAcm2=j,
        NaOH_conc_M=C,
        tna_mode=mode,
        tna_const=0.9,
        cell_voltage_override_V=V,
        CO2_release_eff=0.97,
    )
    new = cell.electrochemical_regeneration(
        CO2_mol_s,
        C,
        CellSpec(j_mA_cm2=j, tna_mode=mode, tna_const=0.9, V_override_V=V, CO2_release_eff=0.97),
        electrons_per_CO2=1.0,  # legacy sizing (E-008 moved the model to 2 for carbonate)
    )
    for old_key, new_key in CELL_FIELDS.items():
        assert getattr(new, new_key) == old[old_key], old_key


@pytest.mark.parametrize(
    ("D", "H", "P_bar"), [(4.19, 14.86, 1.01), (8.0, 25.0, 1.5), (0.15, 0.35, 0.9)]
)
def test_column_capex(ns, D, H, P_bar):
    old = ns["estimate_absorber_capex_geometry"](D_m=D, H_packed_m=H, pressure_bar=P_bar)
    new = tea.column_capex(D, H, P_bar, TEASpec())
    assert new.installed_usd == old["installed_cost_usd"]
    assert new.purchased_usd == old["purchased_cost_usd"]
    assert tea.crf(0.10, 20) == ns["crf"](0.10, 20)


def legacy_cell_result(old: dict) -> CellResult:
    fields = {new: old[legacy] for legacy, new in CELL_FIELDS.items()}
    return CellResult(**fields, electrons_per_CO2=1.0)


@pytest.mark.parametrize(("C", "LG"), [(1.5, 0.005), (0.25, 0.07), (4.0, 0.004), (2.0, 0.03)])
def test_design_costs(ns, C, LG):
    """TEA of cell 8 applied to the legacy absorber/cell outputs of the same point."""
    old = ns["evaluate_design_point"](NaOH_conc_M=C, LG_vol=LG)
    new = tea.design_costs(
        old["D_col_m"],
        old["height_m"],
        old["blower_power_W"],
        old["pump_power_W"],
        old["CO2_captured_mol_s"],
        old["CO2_out_mol_s"],
        legacy_cell_result(old),
        ns["P_abs_bar"],
        TEASpec(),
    )
    pairs = {
        "P_total_W": "P_total_W",
        "CO2_captured_tpy": "CO2_captured_t_y",
        "CO2_stack_tpy": "CO2_stack_t_y",
        "CO2_product_tpy": "CO2_product_t_y",
        "E_cell_kWh_tCO2": "E_cell_kWh_t",
        "E_total_kWh_tCO2": "E_total_kWh_t",
        "capex_total_usd": "capex_total_usd",
        "opex_total_usd_per_year": "opex_total_usd_y",
        "annualized_capex_usd_per_year": "annualized_capex_usd_y",
        "LCOC_usd_per_tCO2": "LCOC_usd_t",
        "indirect_kgCO2e_per_tCO2": "indirect_kgCO2e_t",
    }
    for old_key, new_key in pairs.items():
        assert new[new_key] == old[old_key], old_key


def test_weighted_scores(ns):
    rows = [
        ns["evaluate_design_point"](NaOH_conc_M=C, LG_vol=LG)
        for C in (0.5, 2.0, 4.0)
        for LG in (0.004, 0.02, 0.06)
    ]
    _, ranked = ns["choose_best_design"](pd.DataFrame(rows))
    legacy_keys = {
        "LCOC_usd_t": "LCOC_usd_per_tCO2",
        "capex_total_usd": "capex_total_usd",
        "opex_total_usd_y": "opex_total_usd_per_year",
        "E_total_kWh_t": "E_total_kWh_tCO2",
        "indirect_kgCO2e_t": "indirect_kgCO2e_per_tCO2",
    }
    objectives = {k: np.array([r[v] for r in rows]) for k, v in legacy_keys.items()}
    assert np.array_equal(weighted_scores(objectives, LEGACY_WEIGHTS), ranked["score"].to_numpy())
