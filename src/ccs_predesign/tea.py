"""Techno-economics: column CAPEX, auxiliaries, OPEX, LCOC and indirect emissions (legacy)."""

import math

from scipy import optimize

from .absorber import simulate_absorber
from .cell import electrochemical_regeneration
from .constants import FARADAY_C_mol
from .models import (
    AbsorberSpec,
    CellResult,
    CellSpec,
    ColumnCapex,
    DesignPoint,
    FlueGas,
    Solvent,
    TEASpec,
)


def crf(i: float, n: int) -> float:
    """Capital recovery factor.

    Equations: tea.crf.
    """
    return i * (1 + i) ** n / ((1 + i) ** n - 1)


def safe_div(a: float, b: float, eps: float = 1e-12) -> float:
    return a / max(b, eps)


def annual_t_from_mol_s(CO2_mol_s: float, hours_per_year: float) -> float:
    """CO2 mol/s → t/y at 44.01 g/mol.

    Equations: tea.annual_tonnage.
    """
    return CO2_mol_s * 44.01 * 3600.0 * hours_per_year / 1e6


def column_capex(D_m: float, H_packed_m: float, P_bar: float, tea: TEASpec) -> ColumnCapex:
    """Packed-column CAPEX: (shell + packing + internals) × pressure factor × installation.

    Equations: tea.column_capex.
    """
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
    """CAPEX, OPEX, LCOC and indirect emissions of one design (legacy cell 8).

    Equations: tea.lcoc, tea.indirect_emissions, tea.cell_capex, tea.cell_opex.
    """
    hours = tea.hours_per_year
    P_total_W = blower_power_W + pump_power_W + cell.P_cell_W
    P_total_kW = P_total_W / 1000.0

    captured_t_y = annual_t_from_mol_s(CO2_captured_mol_s, hours)
    t_per_h = safe_div(captured_t_y, hours)

    col = column_capex(D_col_m, height_m, P_bar, tea)
    # Cell: stack + balance of plant per m2, plus the uninstalled-cost factor (errata E-015)
    stack_usd = tea.cell_stack_usd_m2 * cell.A_cell_m2 * (1.0 + tea.cell_uninstalled_factor)
    bop_usd = tea.cell_bop_usd_m2 * cell.A_cell_m2 * (1.0 + tea.cell_uninstalled_factor)
    capex_cell = stack_usd + bop_usd
    capex_blower = tea.blower_capex_usd_kW * (blower_power_W / 1000.0)
    capex_pump = tea.pump_capex_usd_kW * (pump_power_W / 1000.0)
    capex_total = col.installed_usd + capex_cell + capex_blower + capex_pump

    # Fixed O&M on the absorber side; the cell has its own O&M and stack replacement (Zhang)
    opex_fixed = tea.fixed_om_fraction * (capex_total - capex_cell)
    opex_elec = hours * P_total_kW * tea.electricity_usd_kWh
    opex_solvent = tea.solvent_makeup_usd_t * captured_t_y
    opex_water = tea.water_chem_usd_t * captured_t_y
    opex_stack = tea.stack_replacement_fraction * stack_usd / tea.stack_replacement_interval_y
    opex_cell_om = tea.cell_om_fraction * stack_usd
    h2_kg_y = cell.I_cell_A / (2.0 * FARADAY_C_mol) * 2.016e-3 * 3600.0 * hours
    opex_h2 = tea.h2_loss_fraction * h2_kg_y * tea.h2_price_usd_kg
    opex_total = (
        opex_fixed + opex_elec + opex_solvent + opex_water + opex_stack + opex_cell_om + opex_h2
    )

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
        "opex_stack_replacement_usd_y": opex_stack,
        "opex_cell_om_usd_y": opex_cell_om,
        "opex_h2_loss_usd_y": opex_h2,
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

    Points that miss the capture target are costed too; optimize.choose_best excludes them.
    """
    res = simulate_absorber(flue, solv, absorber, LG_vol)

    def cost_at(spec: CellSpec):
        cell = electrochemical_regeneration(
            res.CO2_captured_mol_s, solv.NaOH_M, spec, electrons_per_CO2=res.OH_per_CO2
        )
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
        return cell, costs

    warnings: list[str] = []
    if cell_spec.optimize_j:
        # bug 9: trade cell CAPEX (area ∝ 1/j) against electricity (V and t_Na depend on j)
        lo, hi = cell_spec.j_min_mA_cm2, cell_spec.j_max_mA_cm2
        best = optimize.minimize_scalar(
            lambda j: cost_at(cell_spec.model_copy(update={"j_mA_cm2": j}))[1]["LCOC_usd_t"],
            bounds=(lo, hi),
            method="bounded",
            options={"xatol": 1e-3},
        )
        cell_spec = cell_spec.model_copy(update={"j_mA_cm2": float(best.x)})
        if min(best.x - lo, hi - best.x) < 1e-2 * (hi - lo):
            warnings.append(
                f"optimal current density {best.x:.0f} mA/cm2 sits at the search bound "
                f"[{lo:g}, {hi:g}]; the j–V fit range of the thesis cell is not documented"
            )
    cell, costs = cost_at(cell_spec)
    return DesignPoint(
        NaOH_M=solv.NaOH_M,
        LG_vol=LG_vol,
        absorber=res,
        cell=cell,
        warnings=tuple(warnings),
        **costs,
    )
