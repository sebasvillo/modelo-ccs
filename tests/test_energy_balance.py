"""Adiabatic liquid energy balance (audit §7, errata E-013)."""

import pytest

from ccs_predesign import solvent
from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent


def test_heat_of_absorption_from_formation_enthalpies():
    """NBS (Wagman 1982): CO3-- −677.1, H2O(l) −285.8, CO2(g) −393.5, OH- −230.0 kJ/mol."""
    assert solvent.DH_ABS_CARBONATE_J_mol == pytest.approx(
        (-677.1 - 285.8 - (-393.5 - 2 * 230.0)) * 1e3, rel=1e-12
    )


@pytest.mark.parametrize("LG", [0.006, 0.02, 0.06])
def test_liquid_energy_balance_closes(LG):
    """Would have caught audit §7: the legacy column was isothermal."""
    r = simulate_absorber(FlueGas(), Solvent(NaOH_M=1.5), AbsorberSpec(), LG)
    m_L = solvent.density_kg_m3(r.T_K, 1.5) * r.Ql_m3_s
    heat_W = -solvent.DH_ABS_CARBONATE_J_mol * r.CO2_captured_mol_s
    dT = r.T_liquid_K[0] - r.T_liquid_K[-1]
    assert r.T_liquid_K[-1] == pytest.approx(r.T_K, abs=1e-9)  # lean solvent enters at gas T
    assert dT == pytest.approx(heat_W / (m_L * solvent.CP_SOLUTION_J_kgK), rel=1e-9)
    assert dT > 0


def test_lower_LG_means_hotter_rich_solvent():
    hot = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(), 0.008)
    cool = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(), 0.06)
    assert hot.T_liquid_K[0] > cool.T_liquid_K[0]


def test_solvent_inlet_temperature_is_an_input():
    r = simulate_absorber(FlueGas(T_C=45.0), Solvent(T_in_C=30.0), AbsorberSpec(), 0.02)
    assert r.T_liquid_K[-1] == pytest.approx(303.15, abs=1e-9)
