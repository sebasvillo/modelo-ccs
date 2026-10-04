"""Explained outputs: each headline number with its equation id and substituted values.

This is what the API returns so the website can answer "where does this number come from?"
(CLAUDE.md rule 9). Values are taken from the model results; the few derived ones (E_inf,
k2 and H_cp at the liquid inlet) reuse the same package functions as the model.
"""

from importlib.metadata import version

from . import ERRATA_ID, kinetics, solvent
from .constants import FARADAY_C_mol, G_m_s2, R_J_molK
from .models import (
    CaseInput,
    CaseResponse,
    CaseResult,
    CaseSummary,
    CostShare,
    ExplainedValue,
    Profiles,
    SensitivityPoint,
)
from .tea import crf, design_costs

# Electricity prices for the sensitivity curve [USD/kWh]: from very cheap renewables to
# expensive grid power. Colombia's 2024 solar auction cleared near 0.018 USD/kWh.
ELECTRICITY_PRICES_USD_KWH = (0.0, 0.0182, 0.03, 0.05, 0.075, 0.10, 0.12, 0.15, 0.20)


def explain_case(inp: CaseInput, result: CaseResult) -> tuple[ExplainedValue, ...]:
    a = result.absorber.best
    c = result.cell
    k = result.costs
    spec, tea = inp.absorber, inp.tea
    C = inp.solvent.NaOH_M
    OH_in = a.OH_mol_m3[-1]
    T_in = a.T_liquid_K[-1]
    y_in, y_out = a.yCO2_wet_in, a.yCO2_out
    rho_l_in = solvent.density_kg_m3(T_in, C)
    m_L = rho_l_in * a.Ql_m3_s
    k2_in = kinetics.k2_pohorecki_moniuk_m3_mol_s(T_in, C)
    E_inf_in = kinetics.enhancement_infinite(
        solvent.diffusivity_NaOH_m2_s(T_in, C),
        OH_in,
        a.D_l_m2_s,
        y_out * a.c_tot_mol_m3 / a.H_cc_CO2,
    )
    CRF = crf(tea.discount_rate, tea.project_life_y)
    j_mA = c.j_A_m2 / 10.0

    def ev(key, value, unit, eq, **inputs) -> ExplainedValue:
        return ExplainedValue(key=key, value=value, unit=unit, equation_id=eq, inputs=inputs)

    return (
        ev(
            "height_m",
            a.height_m,
            "m",
            "absorber.height",
            NTU=a.NTU,
            G_flux=a.G_flux_mol_m2_s,
            c_tot=a.c_tot_mol_m3,
            KGa_mean=a.KGa_mean_1_s,
        ),
        ev(
            "capture_fraction",
            a.capture_achieved,
            "-",
            "absorber.capture",
            y=y_in,
            y_top=y_out,
        ),
        ev(
            "diameter_m",
            a.D_col_m,
            "m",
            "packing.column_diameter",
            n_trains=float(a.n_trains),
            Q_G=a.Qg_actual_m3_s,
            f=spec.flood_fraction,
            v_flood=a.v_flood_m_s,
        ),
        ev(
            "KGa_inlet_1_s",
            a.KGa_1_s,
            "1/s",
            "absorber.two_film",
            kG=a.kG_m_s,
            kL=a.kL_m_s,
            E=a.E,
            H_cc=a.H_cc_CO2,
            a_e=a.a_eff_m2_m3,
        ),
        ev(
            "hatta_inlet",
            a.Ha,
            "-",
            "kinetics.hatta",
            k2=k2_in,
            OH=OH_in,
            D_L=a.D_l_m2_s,
            kL=a.kL_m_s,
        ),
        ev("k2_inlet_m3_mol_s", k2_in, "m³/(mol·s)", "kinetics.k2", T=T_in, I=C),
        ev("enhancement_inlet", a.E, "-", "kinetics.decoursey", Ha=a.Ha, E_inf=E_inf_in),
        ev(
            "henry_cc_inlet",
            a.H_cc_CO2,
            "-",
            "solvent.henry_cc",
            H_cp=1.0 / (a.H_cc_CO2 * R_J_molK * T_in),
            R=R_J_molK,
            T=T_in,
        ),
        ev(
            "OH_rich_mol_m3",
            a.OH_mol_m3[0],
            "mol/m³",
            "absorber.hydroxide_balance",
            OH=OH_in,
            nu=a.OH_per_CO2,
            G_I=a.G_inert_mol_s,
            y=y_in,
            y_top=y_out,
            Q_L=a.Ql_m3_s,
        ),
        ev(
            "T_rich_K",
            a.T_liquid_K[0],
            "K",
            "absorber.energy_balance",
            T_L=T_in,
            dH=solvent.DH_ABS_CARBONATE_J_mol,
            G_I=a.G_inert_mol_s,
            y=y_in,
            y_top=y_out,
            m_L=m_L,
            cp_L=solvent.CP_SOLUTION_J_kgK,
        ),
        ev(
            "pump_power_W",
            a.pump_power_W,
            "W",
            "absorber.pump_power",
            rho_L=a.rho_l_kg_m3,
            g=G_m_s2,
            H=a.height_m,
            H_extra=spec.pump_extra_head_m,
            Q_L=a.Ql_m3_s,
            eta=spec.pump_eff,
        ),
        ev(
            "blower_power_W",
            a.blower_power_W,
            "W",
            "absorber.blower_power",
            dP=a.deltaP_Pa,
            Q_G=a.Qg_actual_m3_s,
            eta=spec.blower_eff,
        ),
        ev(
            "cell_current_A",
            c.I_cell_A,
            "A",
            "cell.faraday",
            n_e=c.electrons_per_CO2,
            F=FARADAY_C_mol,
            n_dot=a.CO2_captured_mol_s,
            t_Na=c.tna,
        ),
        ev("cell_area_m2", c.A_cell_m2, "m²", "cell.faraday", I_cell=c.I_cell_A, j=c.j_A_m2),
        ev("cell_voltage_V", c.V_cell_V, "V", "cell.jv_curve", j_mA=j_mA),
        ev("transference_number", c.tna, "-", "cell.transference", C=C, j_mA=j_mA),
        ev(
            "cell_energy_J_mol",
            c.E_J_mol,
            "J/mol",
            "cell.energy",
            n_e=c.electrons_per_CO2,
            F=FARADAY_C_mol,
            V=c.V_cell_V,
            t_Na=c.tna,
        ),
        ev(
            "CO2_captured_t_y",
            k.CO2_captured_t_y,
            "t/año",
            "tea.annual_tonnage",
            n_dot=a.CO2_captured_mol_s,
            h_op=tea.hours_per_year,
        ),
        ev("crf_1_y", CRF, "1/año", "tea.crf", i=tea.discount_rate, n=float(tea.project_life_y)),
        ev(
            "cell_capex_usd",
            k.capex_cell_usd,
            "USD",
            "tea.cell_capex",
            c_stack=tea.cell_stack_usd_m2,
            c_bop=tea.cell_bop_usd_m2,
            A_cell=c.A_cell_m2,
            f_u=tea.cell_uninstalled_factor,
        ),
        ev(
            "LCOC_usd_t",
            k.LCOC_usd_t,
            "USD/t",
            "tea.lcoc",
            CAPEX=k.capex_total_usd,
            CRF=CRF,
            OPEX=k.opex_total_usd_y,
            m_CO2=k.CO2_captured_t_y,
        ),
        ev(
            "indirect_kgCO2e_t",
            k.indirect_kgCO2e_t,
            "kgCO₂e/t",
            "tea.indirect_emissions",
            h_op=tea.hours_per_year,
            W=k.P_total_W,
            EF=tea.grid_EF_kgCO2e_kWh,
            m_CO2=k.CO2_captured_t_y,
        ),
    )


def lcoc_breakdown(inp: CaseInput, result: CaseResult) -> tuple[CostShare, ...]:
    """LCOC split by cost item [USD/t]; the shares add up to the LCOC (Equations: tea.lcoc)."""
    k, tea = result.costs, inp.tea
    CRF = crf(tea.discount_rate, tea.project_life_y)
    t = k.CO2_captured_t_y
    items = [
        ("capex_cell", k.capex_cell_usd * CRF),
        ("capex_absorber", (k.capex_total_usd - k.capex_cell_usd) * CRF),
        ("electricity", k.opex_electricity_usd_y),
        ("stack_replacement", k.opex_stack_replacement_usd_y),
        ("cell_om", k.opex_cell_om_usd_y),
        ("fixed_om", k.opex_fixed_usd_y),
        ("chemicals", k.opex_solvent_usd_y + k.opex_water_chem_usd_y + k.opex_h2_loss_usd_y),
    ]
    return tuple(CostShare(key=key, usd_t=usd_y / t) for key, usd_y in items)


def lcoc_vs_electricity(inp: CaseInput, result: CaseResult) -> tuple[SensitivityPoint, ...]:
    """LCOC of the same design at other electricity prices (costing only; Equations: tea.lcoc)."""
    a = result.absorber.best
    points = []
    for price in ELECTRICITY_PRICES_USD_KWH:
        tea = inp.tea.model_copy(update={"electricity_usd_kWh": price})
        costs = design_costs(
            a.D_col_m, a.height_m, a.blower_power_W, a.pump_power_W, a.CO2_captured_mol_s,
            a.CO2_out_mol_s, result.cell, inp.gas.P_bar, tea, a.n_trains,
        )  # fmt: skip
        points.append(SensitivityPoint(electricity_usd_kWh=price, LCOC_usd_t=costs["LCOC_usd_t"]))
    return tuple(points)


def case_response(inp: CaseInput, result: CaseResult) -> CaseResponse:
    a = result.absorber.best
    k = result.costs
    warnings = a.warnings + (
        () if result.absorber.feasible else ("no feasible design in the L/G scan",)
    )
    return CaseResponse(
        model_version=version("ccs-predesign"),
        errata_id=ERRATA_ID,
        feasible=result.absorber.feasible,
        warnings=warnings,
        summary=CaseSummary(
            packing_name=a.packing_name,
            NaOH_M=inp.solvent.NaOH_M,
            LG_vol=a.LG_vol,
            n_trains=a.n_trains,
            D_col_m=a.D_col_m,
            height_m=a.height_m,
            capture_fraction=a.capture_achieved,
            CO2_captured_t_y=k.CO2_captured_t_y,
            indirect_tCO2e_y=k.indirect_tCO2e_y,
            net_captured_t_y=k.CO2_captured_t_y - k.indirect_tCO2e_y,
            E_total_kWh_t=k.E_total_kWh_t,
            A_cell_m2=result.cell.A_cell_m2,
            P_total_MW=k.P_total_W / 1e6,
            T_rich_C=a.T_liquid_K[0] - 273.15,
            LCOC_usd_t=k.LCOC_usd_t,
            capex_total_usd=k.capex_total_usd,
            opex_total_usd_y=k.opex_total_usd_y,
        ),
        lcoc_breakdown=lcoc_breakdown(inp, result),
        lcoc_vs_electricity=lcoc_vs_electricity(inp, result),
        outputs=explain_case(inp, result),
        profiles=Profiles(
            z_m=a.z_m,
            yCO2=a.yCO2,
            OH_mol_m3=a.OH_mol_m3,
            T_liquid_K=a.T_liquid_K,
            KGa_1_s=a.KGa_profile_1_s,
        ),
    )
