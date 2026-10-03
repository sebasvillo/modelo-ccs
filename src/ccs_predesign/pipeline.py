"""Integrated case: absorber design + cell sizing + annual figures (legacy cell 1)."""

from .absorber import design_absorber_by_grid
from .cell import electrochemical_regeneration
from .models import CaseInput, CaseResult

# LEGACY(audit §7): continuous operation (8760 h/y) here; the legacy TEA uses 8000 h/y.
SECONDS_PER_YEAR = 365.0 * 24.0 * 3600.0
_MW_CO2_g_mol = 44.01  # same literal as the legacy annual figures


def run_case(inp: CaseInput) -> CaseResult:
    design = design_absorber_by_grid(inp.gas, inp.solvent, inp.absorber)
    best = design.best
    cell = electrochemical_regeneration(
        best.CO2_captured_mol_s, inp.solvent.NaOH_M, inp.cell, electrons_per_CO2=best.OH_per_CO2
    )

    def t_per_year(mol_s: float) -> float:
        return mol_s * _MW_CO2_g_mol * SECONDS_PER_YEAR / 1e6

    t_per_hour = best.CO2_captured_mol_s * _MW_CO2_g_mol * 3600.0 / 1e6
    return CaseResult(
        absorber=design,
        cell=cell,
        P_total_W=best.blower_power_W + best.pump_power_W + cell.P_cell_W,
        CO2_captured_t_y=t_per_year(best.CO2_captured_mol_s),
        CO2_stack_t_y=t_per_year(best.CO2_out_mol_s),
        CO2_product_t_y=t_per_year(cell.CO2_released_mol_s),
        E_cell_kWh_t=(cell.P_cell_W / 1000.0) / max(t_per_hour, 1e-12),
    )
