"""Explained outputs: each headline number with its equation id and substituted values.

This is what the API returns so the website can answer "where does this number come from?"
(CLAUDE.md rule 9). Values are taken from the model results; the few derived ones (E_inf,
k2 and H_cp at the liquid inlet) reuse the same package functions as the model.
"""

from importlib.metadata import version

from . import ERRATA_ID, kinetics, solvent
from .constants import FARADAY_C_mol, G_m_s2, R_J_molK
from .models import CaseInput, CaseResponse, CaseResult, ExplainedValue, Profiles
from .tea import crf


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
            G=a.G_mol_s,
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
            G=a.G_mol_s,
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


def case_response(inp: CaseInput, result: CaseResult) -> CaseResponse:
    a = result.absorber.best
    warnings = a.warnings + (
        () if result.absorber.feasible else ("no feasible design in the L/G scan",)
    )
    return CaseResponse(
        model_version=version("ccs-predesign"),
        errata_id=ERRATA_ID,
        feasible=result.absorber.feasible,
        warnings=warnings,
        outputs=explain_case(inp, result),
        profiles=Profiles(
            z_m=a.z_m,
            yCO2=a.yCO2,
            OH_mol_m3=a.OH_mol_m3,
            T_liquid_K=a.T_liquid_K,
            KGa_1_s=a.KGa_profile_1_s,
        ),
    )
