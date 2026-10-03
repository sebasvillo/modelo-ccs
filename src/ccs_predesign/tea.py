"""Techno-economics: column CAPEX, auxiliaries, OPEX, LCOC and indirect emissions (legacy)."""

import math

from .absorber import simulate_absorber
from .cell import electrochemical_regeneration
from .models import (
    AbsorberSpec,
    CellResult,
    CellSpec,
    ColumnCapex,
    DesignPoint,
    FlueGas,
    Solvent,
    TEASpec,
    TEASummary,
)

# Fallback specific electricity of legacy cell 2 when no power figure is available.
DEFAULT_SPECIFIC_ELECTRICITY_kWh_t = 250.0


def crf(i: float, n: int) -> float:
    """Capital recovery factor."""
    return i * (1 + i) ** n / ((1 + i) ** n - 1)


def safe_div(a: float, b: float, eps: float = 1e-12) -> float:
    return a / max(b, eps)


def annual_t_from_mol_s(CO2_mol_s: float, hours_per_year: float) -> float:
    """CO2 mol/s → t/y at 44.01 g/mol."""
    return CO2_mol_s * 44.01 * 3600.0 * hours_per_year / 1e6


def column_capex(D_m: float, H_packed_m: float, P_bar: float, tea: TEASpec) -> ColumnCapex:
    """Packed-column CAPEX: (shell + packing + internals) × pressure factor × installation."""
    A_cs_m2 = math.pi * D_m**2 / 4.0
    shell_area_m2 = math.pi * D_m * H_packed_m
    packed_volume_m3 = A_cs_m2 * H_packed_m
    pressure_factor = 1.0 + 0.03 * max(P_bar - 1.0, 0.0)
    shell = tea.shell_usd_m2 * shell_area_m2
    packing_ = tea.packing_usd_m3 * packed_volume_m3
    internals = tea.internals_usd_m2 * A_cs_m2
    purchased = pressure_factor * (shell + packing_ + internals)
    return ColumnCapex(
        A_cs_m2=A_cs_m2,
        shell_area_m2=shell_area_m2,
        packed_volume_m3=packed_volume_m3,
        pressure_factor=pressure_factor,
        shell_usd=shell,
        packing_usd=packing_,
        internals_usd=internals,
        purchased_usd=purchased,
        installed_usd=purchased * tea.installation_factor,
    )


def design_costs(
    D_col_m: float,
    height_m: float,
    blower_power_W: float,
    pump_power_W: float,
    CO2_captured_mol_s: float,
    CO2_out_mol_s: float,
    cell: CellResult,
    P_bar: float,
    tea: TEASpec,
) -> dict[str, float]:
    """CAPEX, OPEX, LCOC and indirect emissions of one design (legacy cell 8)."""
    hours = tea.hours_per_year
    P_total_W = blower_power_W + pump_power_W + cell.P_cell_W
    P_total_kW = P_total_W / 1000.0

    captured_t_y = annual_t_from_mol_s(CO2_captured_mol_s, hours)
    t_per_h = safe_div(captured_t_y, hours)

    col = column_capex(D_col_m, height_m, P_bar, tea)
    capex_cell = tea.cell_capex_usd_m2 * cell.A_cell_m2
    capex_blower = tea.blower_capex_usd_kW * (blower_power_W / 1000.0)
    capex_pump = tea.pump_capex_usd_kW * (pump_power_W / 1000.0)
    capex_total = col.installed_usd + capex_cell + capex_blower + capex_pump

    opex_fixed = tea.fixed_om_fraction * capex_total
    opex_elec = hours * P_total_kW * tea.electricity_usd_kWh
    opex_solvent = tea.solvent_makeup_usd_t * captured_t_y
    opex_water = tea.water_chem_usd_t * captured_t_y
    opex_total = opex_fixed + opex_elec + opex_solvent + opex_water

    annualized = capex_total * crf(tea.discount_rate, tea.project_life_y)
    indirect_kg_y = hours * P_total_kW * tea.grid_EF_kgCO2e_kWh

    return {
        "P_total_W": P_total_W,
        "CO2_captured_t_y": captured_t_y,
        "CO2_stack_t_y": annual_t_from_mol_s(CO2_out_mol_s, hours),
        "CO2_product_t_y": annual_t_from_mol_s(cell.CO2_released_mol_s, hours),
        "E_cell_kWh_t": safe_div(cell.P_cell_W / 1000.0, t_per_h),
        "E_total_kWh_t": safe_div(P_total_W / 1000.0, t_per_h),
        "capex_column_usd": col.installed_usd,
        "capex_cell_usd": capex_cell,
        "capex_blower_usd": capex_blower,
        "capex_pump_usd": capex_pump,
        "capex_total_usd": capex_total,
        "annualized_capex_usd_y": annualized,
        "opex_fixed_usd_y": opex_fixed,
        "opex_electricity_usd_y": opex_elec,
        "opex_solvent_usd_y": opex_solvent,
        "opex_water_chem_usd_y": opex_water,
        "opex_total_usd_y": opex_total,
        "LCOC_usd_t": safe_div(annualized + opex_total, captured_t_y),
        "indirect_kgCO2e_y": indirect_kg_y,
        "indirect_tCO2e_y": indirect_kg_y / 1000.0,
        "indirect_kgCO2e_t": safe_div(indirect_kg_y, captured_t_y),
    }


def evaluate_design_point(
    flue: FlueGas,
    solv: Solvent,
    absorber: AbsorberSpec,
    cell_spec: CellSpec,
    tea: TEASpec,
    LG_vol: float,
) -> DesignPoint:
    """Absorber + cell + TEA for one (NaOH, L/G) point (legacy cell 8).

    LEGACY(audit §7): points that miss the capture target are costed like any other.
    """
    res = simulate_absorber(flue, solv, absorber, LG_vol)
    cell = electrochemical_regeneration(res.CO2_captured_mol_s, solv.NaOH_M, cell_spec)
    costs = design_costs(
        res.D_col_m,
        res.height_m,
        res.blower_power_W,
        res.pump_power_W,
        res.CO2_captured_mol_s,
        res.CO2_out_mol_s,
        cell,
        flue.P_bar,
        tea,
    )
    return DesignPoint(NaOH_M=solv.NaOH_M, LG_vol=LG_vol, absorber=res, cell=cell, **costs)


def tea_summary(
    D_m: float,
    H_packed_m: float,
    P_bar: float,
    CO2_captured_mol_s: float,
    P_total_W: float | None,
    tea: TEASpec,
) -> TEASummary:
    """Single-design TEA of legacy cell 2.

    LEGACY(audit §7): uses 44.0 g/mol (cell 8 uses 44.01), ignores cell/blower/pump CAPEX,
    and in the notebook receives H_packed_m = 25 m (default) instead of the model height.
    """
    hours = tea.hours_per_year
    annual_t = CO2_captured_mol_s * 44.0e-3 * 3600.0 * hours / 1000.0
    if P_total_W is None or P_total_W <= 0:
        kWh_y = annual_t * DEFAULT_SPECIFIC_ELECTRICITY_kWh_t
    else:
        kWh_y = (P_total_W / 1000.0) * hours

    col = column_capex(D_m, H_packed_m, P_bar, tea)
    capex_total = col.installed_usd + 0.0 + 0.0 + 0.0  # no auxiliary CAPEX in cell 2

    opex_fixed = tea.fixed_om_fraction * capex_total
    opex_elec = kWh_y * tea.electricity_usd_kWh
    opex_solvent = tea.solvent_makeup_usd_t * annual_t
    opex_water = tea.water_chem_usd_t * annual_t
    opex_total = opex_fixed + opex_elec + opex_solvent + opex_water
    annualized = capex_total * crf(tea.discount_rate, tea.project_life_y)

    return TEASummary(
        annual_captured_t_y=annual_t,
        annual_electricity_kWh_y=kWh_y,
        column=col,
        capex_total_usd=capex_total,
        opex_fixed_usd_y=opex_fixed,
        opex_electricity_usd_y=opex_elec,
        opex_solvent_usd_y=opex_solvent,
        opex_water_chem_usd_y=opex_water,
        opex_total_usd_y=opex_total,
        annualized_capex_usd_y=annualized,
        LCOC_usd_t=(annualized + opex_total) / annual_t,
        indirect_tCO2e_y=(kWh_y * tea.grid_EF_kgCO2e_kWh) / 1000.0,
        indirect_kgCO2e_t=(kWh_y * tea.grid_EF_kgCO2e_kWh) / annual_t,
    )
