"""Integrated case: absorber design + cell sizing + annual figures."""

from .absorber import design_absorber_by_grid
from .cell import electrochemical_regeneration
from .models import CaseInput, CaseResult
from .tea import annual_t_from_mol_s


def run_case(inp: CaseInput) -> CaseResult:
    """Annual figures use inp.tea.hours_per_year, the same basis as the TEA (errata E-009)."""
    design = design_absorber_by_grid(inp.gas, inp.solvent, inp.absorber)
    best = design.best
    cell = electrochemical_regeneration(
        best.CO2_captured_mol_s, inp.solvent.NaOH_M, inp.cell, electrons_per_CO2=best.OH_per_CO2
    )
    hours = inp.tea.hours_per_year
    t_per_hour = annual_t_from_mol_s(best.CO2_captured_mol_s, 1.0)
    return CaseResult(
        absorber=design,
        cell=cell,
        P_total_W=best.blower_power_W + best.pump_power_W + cell.P_cell_W,
        CO2_captured_t_y=annual_t_from_mol_s(best.CO2_captured_mol_s, hours),
        CO2_stack_t_y=annual_t_from_mol_s(best.CO2_out_mol_s, hours),
        CO2_product_t_y=annual_t_from_mol_s(cell.CO2_released_mol_s, hours),
        E_cell_kWh_t=(cell.P_cell_W / 1000.0) / max(t_per_hour, 1e-12),
    )
