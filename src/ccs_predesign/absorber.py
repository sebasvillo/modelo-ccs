"""Packed absorber: column profile for one L/G and the L/G design scan (legacy model)."""

import logging
import math

from . import gas, kinetics, packing, solvent
from .constants import G_m_s2, MW_kg_mol, R_J_molK
from .models import AbsorberDesign, AbsorberResult, AbsorberSpec, FlueGas, Solvent

logger = logging.getLogger(__name__)


def ntu_increment(
    KGa_1_s: float, c_tot_mol_m3: float, G_flux_mol_m2_s: float, dz_m: float
) -> float:
    """Gas-phase transfer units in a slice: dNTU = K_G·a·c_tot·dz / G''  [dimensionless].

    From the gas balance G''·dy/dz = −K_G·a·c_tot·(y − y*), with K_G in m/s, a in 1/m,
    c_tot = P/(R·T) in mol/m3 and G'' the gas molar flux in mol/(m2·s). Constant G'' is assumed
    (dilute-gas approximation). Provenance: theory (film model, e.g. Seader & Henley,
    Separation Process Principles, ch. 6). Fixes audit §1 (errata E-001).
    """
    return KGa_1_s * c_tot_mol_m3 * dz_m / G_flux_mol_m2_s


def simulate_absorber(
    flue: FlueGas,
    solv: Solvent,
    spec: AbsorberSpec,
    LG_vol: float,
    D_col_fixed_m: float | None = None,
) -> AbsorberResult:
    """Integrate the column for a volumetric L/G (m3 liquid / m3 actual gas).

    The diameter comes from spec.flood_fraction of the flooding velocity, unless
    D_col_fixed_m is given (legacy cell 4 variant, used for the bench column).

    Ha, E and KGa follow the local OH- concentration and absorption stops at the
    stoichiometric OH- capacity (audit §4, errata E-004).
    LEGACY(audit §5): lean solvent (x = 0) imposed at the gas inlet, i.e. co-current.
    """
    pk = packing.PACKINGS[spec.packing_name]
    T_K = gas.c_to_k(flue.T_C)
    P_Pa = flue.P_bar * 1e5

    comp_wet = gas.wet_composition_from_dry(flue.y_CO2_dry, flue.y_O2_dry, T_K, P_Pa)
    y_in = comp_wet["CO2"]
    y_H2O = comp_wet["H2O"]
    y_out_target = y_in * (1.0 - spec.capture_target)

    n_dry = gas.dry_molar_flow_mol_s(flue.Q_dry_Nm3_h)
    n_wet = n_dry / max(1.0 - y_H2O, 1e-12)
    Qg_m3_s = n_wet * R_J_molK * T_K / P_Pa
    rho_g = gas.density_ideal_kg_m3(T_K, P_Pa, comp_wet)

    Ql_m3_s = LG_vol * Qg_m3_s
    rho_l = solvent.density_kg_m3(T_K, solv.NaOH_M)
    if D_col_fixed_m is None:
        D_col, v_flood, v_oper = packing.column_diameter_m(
            Qg_m3_s, rho_g, rho_l, pk, flood_frac=spec.flood_fraction
        )
    else:
        D_col = float(D_col_fixed_m)
        v_flood = packing.flooding_velocity_m_s(rho_g, rho_l, pk)
        v_oper = Qg_m3_s / max(math.pi * D_col**2 / 4.0, 1e-12)
    A_col = math.pi * D_col**2 / 4.0
    uG = Qg_m3_s / A_col
    uL = Ql_m3_s / A_col

    k2 = kinetics.k2_m3_mol_s(T_K, spec.k2_ref_L_mol_s)
    OH_mol_m3 = solv.NaOH_M * 1000.0
    k1_pseudo = k2 * OH_mol_m3

    H_cc = solvent.henry_cc_CO2(T_K, solvent.naoh_ions_kmol_m3(solv.NaOH_M))

    def coefficients(E: float) -> dict[str, float]:
        return packing.mass_transfer_coefficients(
            T_K, P_Pa, comp_wet, uG, uL, pk, solv.NaOH_M,
            H_cc, E=E,
        )  # fmt: skip

    props0 = coefficients(1.0)
    Ha = kinetics.hatta_number(k1_pseudo, props0["Dl"], props0["kL"])
    E = kinetics.enhancement_factor(Ha)
    props = coefficients(E)
    KGa = props["KGa"]  # at the liquid inlet (fresh solvent)

    L_mol_s = rho_l * Ql_m3_s / MW_kg_mol["water"]  # LEGACY: water-proxy molar mass
    G_mol_s = n_wet

    c_tot = P_Pa / (R_J_molK * T_K)
    G_flux = G_mol_s / A_col
    logger.debug("KGa_in=%.3f 1/s | G''=%.3e mol/m2/s", KGa, G_flux)

    def local_KGa(OH: float) -> tuple[float, float, float]:
        """Ha, E and KGa with the local hydroxide concentration (audit §4, errata E-004)."""
        Ha_loc = kinetics.hatta_number(k2 * OH, props0["Dl"], props0["kL"])
        E_loc = kinetics.enhancement_factor(Ha_loc)
        KG = 1.0 / (1.0 / max(props0["kG"], 1e-12) + H_cc / max(props0["kL"] * E_loc, 1e-12))
        return Ha_loc, E_loc, KG * props0["a_eff"]

    # Carbonate route CO2 + 2 OH- -> CO3-- + H2O: each mol of CO2 consumes 2 mol of OH-.
    OH_PER_CO2 = 2.0
    OH_in_mol_s = OH_mol_m3 * Ql_m3_s
    capacity_mol_s = OH_in_mol_s / OH_PER_CO2

    z_vals, y_vals, x_vals = [0.0], [y_in], [0.0]
    ystar_vals, rate_vals = [spec.H_eq * 0.0], [0.0]
    OH_vals, Ha_vals, E_vals, KGa_vals = [OH_mol_m3], [Ha], [E], [KGa]
    y, x, z, OH, ntu_total = y_in, 0.0, 0.0, OH_mol_m3, 0.0
    reached = exhausted = False
    KGa_loc = KGa
    while z < spec.max_height_m:
        # explicit step: chemistry evaluated at the start of the slice (node values)
        step = ntu_increment(KGa_loc, c_tot, G_flux, spec.dz_m)
        y_star = spec.H_eq * x
        driving0 = max(y - y_star, 0.0)
        y_new = max(y_star + (y - y_star) * math.exp(-step), 0.0)
        dn_mol_s = G_mol_s * (y - y_new)
        available_mol_s = OH * Ql_m3_s / OH_PER_CO2
        if dn_mol_s >= available_mol_s:  # stoichiometric limit: solvent exhausted
            dn_mol_s = available_mol_s
            y_new = y - dn_mol_s / G_mol_s
            exhausted = True
        OH = 0.0 if exhausted else OH - OH_PER_CO2 * dn_mol_s / Ql_m3_s
        x = x + dn_mol_s / max(L_mol_s, 1e-12)
        y = y_new
        z += spec.dz_m
        ntu_total += step

        z_vals.append(z)
        y_vals.append(y)
        x_vals.append(x)
        ystar_vals.append(spec.H_eq * x)
        rate_vals.append(KGa_loc * driving0)
        Ha_loc, E_loc, KGa_loc = local_KGa(OH)
        OH_vals.append(OH)
        Ha_vals.append(Ha_loc)
        E_vals.append(E_loc)
        KGa_vals.append(KGa_loc)

        if y <= y_out_target:
            reached = True
            break
        if exhausted:
            break

    warnings = []
    if exhausted:
        warnings.append(
            f"solvent exhausted at z = {z:.2f} m: OH- capacity {capacity_mol_s:.4g} mol/s CO2 "
            "(carbonate route, 2 OH- per CO2); the slow bicarbonate route is not modelled"
        )
    if not reached and not exhausted:
        warnings.append(
            f"capture target not reached within max_height_m = {spec.max_height_m:g} m "
            f"(capture {(y_in - y) / y_in:.1%})"
        )
    if capacity_mol_s < G_mol_s * (y_in - y_out_target):
        warnings.append(
            f"stoichiometric OH- capacity ({capacity_mol_s:.4g} mol/s CO2) is below the capture "
            f"target ({G_mol_s * (y_in - y_out_target):.4g} mol/s): raise L/G or NaOH"
        )

    height_m = z_vals[-1]
    y_out = y_vals[-1]

    dpdz = packing.pressure_drop_ergun_Pa_m(props["rho_g"], props["mu_g"], uG, pk)
    deltaP = dpdz * height_m
    blower_W = deltaP * Qg_m3_s / max(spec.blower_eff, 1e-12)
    pump_head_m = 5.0 + 0.5 * height_m  # LEGACY(audit §10): should be >= H + losses
    pump_W = rho_l * G_m_s2 * pump_head_m * Ql_m3_s / max(spec.pump_eff, 1e-12)

    CO2_in = G_mol_s * y_in
    CO2_out = G_mol_s * y_out

    return AbsorberResult(
        LG_vol=LG_vol,
        packing_name=spec.packing_name,
        T_K=T_K,
        P_Pa=P_Pa,
        comp_wet_in=comp_wet,
        yCO2_wet_in=y_in,
        yCO2_out=y_out,
        yCO2_out_target=y_out_target,
        capture_achieved=(y_in - y_out) / max(y_in, 1e-12),
        reached_target=reached,
        Qg_actual_m3_s=Qg_m3_s,
        Ql_m3_s=Ql_m3_s,
        n_dry_mol_s=n_dry,
        n_wet_mol_s=n_wet,
        G_mol_s=G_mol_s,
        G_flux_mol_m2_s=G_flux,
        c_tot_mol_m3=c_tot,
        NTU=ntu_total,
        L_mol_s=L_mol_s,
        rho_g_kg_m3=rho_g,
        rho_l_kg_m3=rho_l,
        mu_g_Pa_s=props["mu_g"],
        mu_l_Pa_s=props["mu_l"],
        kG_m_s=props["kG"],
        kL_m_s=props["kL"],
        KGa_1_s=KGa,
        a_eff_m2_m3=props["a_eff"],
        wetting_fraction=props["wet"],
        Ha=Ha,
        E=E,
        H_cc_CO2=H_cc,
        k1_pseudo_1_s=k1_pseudo,
        D_g_m2_s=props["Dg"],
        D_l_m2_s=props["Dl"],
        ReG=props["ReG"],
        ReL=props["ReL"],
        ScG=props["ScG"],
        ScL=props["ScL"],
        uG_m_s=uG,
        uL_m_s=uL,
        D_col_m=D_col,
        v_flood_m_s=v_flood,
        v_oper_m_s=v_oper,
        flood_fraction_actual=uG / max(v_flood, 1e-12),
        feasible_hydraulically=uG <= spec.flood_fraction * v_flood,
        height_m=height_m,
        dpdz_Pa_m=dpdz,
        deltaP_Pa=deltaP,
        blower_power_W=blower_W,
        pump_power_W=pump_W,
        hydraulic_power_W=blower_W + pump_W,
        CO2_in_mol_s=CO2_in,
        CO2_out_mol_s=CO2_out,
        CO2_captured_mol_s=max(CO2_in - CO2_out, 0.0),
        stack_CO2_ppmv=y_out * 1e6,
        z_m=tuple(z_vals),
        yCO2=tuple(y_vals),
        x_loading=tuple(x_vals),
        y_star=tuple(ystar_vals),
        rate_indicator=tuple(rate_vals),
        OH_mol_m3=tuple(OH_vals),
        Ha_profile=tuple(Ha_vals),
        E_profile=tuple(E_vals),
        KGa_profile_1_s=tuple(KGa_vals),
        OH_capacity_CO2_mol_s=capacity_mol_s,
        solvent_exhausted=exhausted,
        warnings=tuple(warnings),
    )


def design_absorber_by_grid(flue: FlueGas, solv: Solvent, spec: AbsorberSpec) -> AbsorberDesign:
    """Scan spec.LG_grid_vol and pick a feasible design by spec.objective.

    If no point reaches the target, the best capture is returned with feasible=False.
    """
    candidates = tuple(simulate_absorber(flue, solv, spec, float(lg)) for lg in spec.LG_grid_vol)
    feasible = [r for r in candidates if r.reached_target and r.height_m <= spec.max_height_m]
    if not feasible:
        best = max(candidates, key=lambda r: (r.capture_achieved, -r.height_m))
        return AbsorberDesign(best=best, feasible=False, candidates=candidates)
    if spec.objective == "min_height":
        best = min(feasible, key=lambda r: r.height_m)
    else:
        best = min(feasible, key=lambda r: r.hydraulic_power_W)
    return AbsorberDesign(best=best, feasible=True, candidates=candidates)
