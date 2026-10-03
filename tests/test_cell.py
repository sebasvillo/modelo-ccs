"""Electrochemical regeneration: Faraday's law with the absorber's chemistry (audit §6, E-008)."""

import pytest

from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.cell import electrochemical_regeneration
from ccs_predesign.constants import FARADAY_C_mol
from ccs_predesign.models import AbsorberSpec, CaseInput, CellSpec, FlueGas, Solvent
from ccs_predesign.pipeline import run_case


@pytest.mark.parametrize("n_e", [1.0, 2.0])
def test_faraday_sizing(n_e):
    c = electrochemical_regeneration(100.0, 1.5, CellSpec(), electrons_per_CO2=n_e)
    assert c.I_cell_A == pytest.approx(n_e * FARADAY_C_mol * 100.0 / c.tna, rel=1e-12)
    assert c.A_cell_m2 == pytest.approx(c.I_cell_A / c.j_A_m2, rel=1e-12)
    assert c.E_J_mol == pytest.approx(n_e * FARADAY_C_mol * c.V_cell_V / c.tna, rel=1e-12)
    assert c.P_cell_W == pytest.approx(c.E_J_mol * 100.0, rel=1e-12)


def test_invalid_electron_count_is_rejected():
    with pytest.raises(ValueError):
        electrochemical_regeneration(100.0, 1.5, CellSpec(), electrons_per_CO2=1.5)


def test_charge_balance_with_absorber():
    """Would have caught audit §6: Na+ moved by the cell must equal OH- consumed in the column."""
    r = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(), 0.03)
    c = electrochemical_regeneration(
        r.CO2_captured_mol_s, 1.5, CellSpec(), electrons_per_CO2=r.OH_per_CO2
    )
    Na_moved_mol_s = c.I_cell_A * c.tna / FARADAY_C_mol
    OH_consumed_mol_s = (r.OH_mol_m3[-1] - r.OH_mol_m3[0]) * r.Ql_m3_s
    assert Na_moved_mol_s == pytest.approx(OH_consumed_mol_s, rel=1e-9)


def test_pipeline_uses_carbonate_route():
    result = run_case(CaseInput())
    assert result.absorber.best.OH_per_CO2 == 2.0
    assert result.cell.electrons_per_CO2 == 2.0


def test_annual_figures_use_tea_operating_hours():
    """Would have caught bug 9: annual tonnages used 8760 h while the TEA used 8000 h."""
    from ccs_predesign.models import TEASpec
    from ccs_predesign.tea import annual_t_from_mol_s

    for hours in (8000.0, 7000.0):
        r = run_case(CaseInput(tea=TEASpec(hours_per_year=hours)))
        b = r.absorber.best
        assert r.CO2_captured_t_y == pytest.approx(annual_t_from_mol_s(b.CO2_captured_mol_s, hours))
        assert r.CO2_captured_t_y == pytest.approx(b.CO2_captured_mol_s * 44.01e-6 * 3600 * hours)
