"""Packed absorber: column profile for one L/G and the L/G design scan (legacy model)."""

import logging
import math

from scipy import optimize

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

    Counter-current (audit §5, errata E-007): gas enters at z = 0, fresh solvent at the top.
    Design mode integrates up from the rich end until the target is met; if it is not met
    within max_height_m, the column is rated at that height by shooting on the gas outlet.
    Ha, E, KGa and y* follow the local liquid composition (audit §4, errata E-004/E-006).
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

    # Carbonate route CO2 + 2 OH- -> CO3-- + H2O: each mol of CO2 consumes 2 mol of OH-.
    OH_PER_CO2 = 2.0
    capacity_mol_s = OH_mol_m3 * Ql_m3_s / OH_PER_CO2
    kG, kL, a_eff, Dl = props0["kG"], props0["kL"], props0["a_eff"], props0["Dl"]
    ntu_per_KGa_m = c_tot / G_flux  # dNTU = KGa · c_tot · dz / G''
    K1, K2, Kw = solvent.carbonate_constants(T_K)
    ystar_coeff = H_cc * 1000.0 * Kw**2 / (K1 * K2 * c_tot)  # y* = coeff·[CO3]_M/[OH]_M²

    def node(y: float, y_top: float) -> tuple[float, float, float, float, float]:
        """Local liquid state and rates at a height where the gas mole fraction is y.

        Counter-current balance over the section above: the liquid has absorbed
        G·(y − y_top) since entering fresh at the top (audit §5, errata E-007).
        """
        OH = OH_mol_m3 - OH_PER_CO2 * G_mol_s * (y - y_top) / Ql_m3_s
        if OH <= 0.0:  # no hydroxide left: no reaction and no driving force
            KG = 1.0 / (1.0 / kG + H_cc / kL)
            return OH, 0.0, 1.0, KG * a_eff, y
        Ha_loc = kinetics.hatta_number(k2 * OH, Dl, kL)
        E_loc = kinetics.enhancement_factor(Ha_loc)
        KG = 1.0 / (1.0 / max(kG, 1e-12) + H_cc / max(kL * E_loc, 1e-12))
        CO3_M = (OH_mol_m3 - OH) / OH_PER_CO2 / 1000.0
        y_star = min(ystar_coeff * CO3_M / (OH / 1000.0) ** 2, y)
        return OH, Ha_loc, E_loc, KG * a_eff, y_star

    def pinch_gap(y: float, y_top: float) -> float:
        """y*(liquid in contact with gas y) − y; it crosses zero at a rich-end pinch."""
        OH = OH_mol_m3 - OH_PER_CO2 * G_mol_s * (y - y_top) / Ql_m3_s
        CO3_M = (OH_mol_m3 - OH) / OH_PER_CO2 / 1000.0
        return ystar_coeff * CO3_M / (OH / 1000.0) ** 2 - y

    def march(y_top: float, z_end: float, y_stop: float | None = None, keep: bool = False):
        """Integrate upward from the gas inlet (z = 0) with explicit exponential steps.

        With y_stop (design mode) the last slice is cut exactly where y = y_stop (audit §11);
        otherwise the march ends exactly at z_end (rating mode).
        """
        y, z, ntu, reached = y_in, 0.0, 0.0, False
        prof: dict[str, list[float]] = {
            k: [] for k in ("z", "y", "OH", "Ha", "E", "KGa", "ys", "r")
        }
        while z < z_end - 1e-12:
            OH, Ha_loc, E_loc, KGa_loc, y_star = node(y, y_top)
            if keep:
                for k, v in zip(
                    prof,
                    (z, y, OH, Ha_loc, E_loc, KGa_loc, y_star, KGa_loc * (y - y_star)),
                    strict=True,
                ):
                    prof[k].append(v)
            k_m = KGa_loc * ntu_per_KGa_m  # transfer units per metre
            h = min(spec.dz_m, z_end - z)
            y_new = y_star + (y - y_star) * math.exp(-k_m * h)
            if y_stop is not None and y_new <= y_stop:
                h = math.log((y - y_star) / (y_stop - y_star)) / k_m
                y_new, reached = y_stop, True
            y, z, ntu = y_new, z + h, ntu + k_m * h
            if reached:
                break
        if keep:
            OH, Ha_loc, E_loc, KGa_loc, y_star = node(y, y_top)
            for k, v in zip(
                prof,
                (z, y, OH, Ha_loc, E_loc, KGa_loc, y_star, KGa_loc * (y - y_star)),
                strict=True,
            ):
                prof[k].append(v)
        return y, z, ntu, reached, prof

    def march_down(y_top: float, height: float, keep: bool = False):
        """Rating mode: integrate down from the lean end (fresh solvent, gas outlet y_top).

        Going down, y_(z−h) = y* + (y − y*)·exp(+k·h). Well-conditioned near a rich-end pinch,
        which an upward march would start inside of.
        """
        y, z, ntu = y_top, height, 0.0
        prof: dict[str, list[float]] = {
            k: [] for k in ("z", "y", "OH", "Ha", "E", "KGa", "ys", "r")
        }
        while z > 1e-12:
            OH, Ha_loc, E_loc, KGa_loc, y_star = node(y, y_top)
            if keep:
                for k, v in zip(
                    prof,
                    (z, y, OH, Ha_loc, E_loc, KGa_loc, y_star, KGa_loc * (y - y_star)),
                    strict=True,
                ):
                    prof[k].append(v)
            k_m = KGa_loc * ntu_per_KGa_m
            h = min(spec.dz_m, z)
            y_new = y_star + (y - y_star) * math.exp(k_m * h)
            y_dry = y_top + OH_mol_m3 * Ql_m3_s / (OH_PER_CO2 * G_mol_s) * (1.0 - 1e-12)
            y_hi = min(y_new, y_dry)  # gas level at which the liquid would hold no OH- left
            if pinch_gap(y, y_top) >= 0.0:  # already sitting on the pinch: no driving force
                y_new = y
            elif y_hi < y_new or pinch_gap(y_hi, y_top) > 0.0:
                # the explicit step crossed the equilibrium pinch: stop exactly on it
                y_new = optimize.brentq(pinch_gap, y, y_hi, args=(y_top,), xtol=1e-15, rtol=1e-13)
            y = y_new
            z, ntu = z - h, ntu + k_m * h
        if keep:
            OH, Ha_loc, E_loc, KGa_loc, y_star = node(y, y_top)
            for k, v in zip(
                prof,
                (0.0, y, OH, Ha_loc, E_loc, KGa_loc, y_star, KGa_loc * (y - y_star)),
                strict=True,
            ):
                prof[k].append(v)
            prof = {k: v[::-1] for k, v in prof.items()}
        return y, ntu, prof

    OH_bottom_design = OH_mol_m3 - OH_PER_CO2 * G_mol_s * (y_in - y_out_target) / Ql_m3_s
    reached = False
    if OH_bottom_design > 0.0 and node(y_in, y_out_target)[4] < y_in:
        y_end, z_end, ntu_total, reached, prof = march(
            y_out_target, spec.max_height_m, y_stop=y_out_target, keep=True
        )
    if reached:
        y_out = y_out_target
    else:
        # rating at the maximum height: shoot on the gas outlet so that the gas reaching the
        # bottom equals the feed (the solvent enters fresh at the top by construction)
        y_lo = max(y_in - capacity_mol_s / G_mol_s, 0.0)
        H = spec.max_height_m

        def residual(y_top: float) -> float:
            return march_down(y_top, H)[0] - y_in

        y_out = (
            y_lo
            if residual(y_lo) >= 0.0
            else optimize.brentq(residual, y_lo, y_in, xtol=1e-15, rtol=1e-13)
        )
        _, ntu_total, prof = march_down(y_out, H, keep=True)
        z_end = H

    z_vals, y_vals = prof["z"], prof["y"]
    OH_vals = prof["OH"]
    x_vals = [G_mol_s * (yv - y_out) / max(L_mol_s, 1e-12) for yv in y_vals]
    exhausted = OH_vals[0] <= 0.01 * OH_mol_m3

    warnings = []
    if exhausted:
        warnings.append(
            f"rich solvent leaves with {max(OH_vals[0], 0.0) / OH_mol_m3:.1%} of the inlet OH- "
            "(carbonate route, 2 OH- per CO2); the slow bicarbonate route is not modelled"
        )
    if not reached:
        warnings.append(
            f"capture target not reached within max_height_m = {spec.max_height_m:g} m "
            f"(capture {(y_in - y_out) / y_in:.1%})"
        )
    if capacity_mol_s < G_mol_s * (y_in - y_out_target):
        warnings.append(
            f"stoichiometric OH- capacity ({capacity_mol_s:.4g} mol/s CO2) is below the capture "
            f"target ({G_mol_s * (y_in - y_out_target):.4g} mol/s): raise L/G or NaOH"
        )

    height_m = z_end

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
        y_star=tuple(prof["ys"]),
        rate_indicator=tuple(prof["r"]),
        OH_mol_m3=tuple(OH_vals),
        Ha_profile=tuple(prof["Ha"]),
        E_profile=tuple(prof["E"]),
        KGa_profile_1_s=tuple(prof["KGa"]),
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
