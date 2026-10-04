"""NaOH solution properties: density, viscosity, heat capacity, diffusivities, CO2 solubility."""

import math

from .constants import MW_kg_mol, R_J_molK

# Laliberté, M. (2009). J. Chem. Eng. Data 54, 1725–1760, Table 1, row NaOH: apparent-density
# coefficients c0–c4 (Laliberté & Cooper 2004 model), solute-viscosity coefficients v1–v6
# (Laliberté 2007 model) and apparent heat-capacity coefficients a1–a6.
LALIBERTE_NAOH_C = (
    319.020509469838,
    528.592358475315,
    -0.102197896602724,
    0.000350420706415566,
    765.970470238438,
)
LALIBERTE_NAOH_V = (
    448.457566713375,
    0.00871452408983102,
    -431.97212334697,
    0.0160144202049452,
    104.011738670148,
    4.64493684488816,
)
LALIBERTE_NAOH_A = (
    -0.922780764834469,
    -0.0412353462450485,
    1.8722524604359,
    -5.94223565147303,
    3.13007617842649,
    0.141040805508813,
)
# Fitted range of each correlation (Laliberté 2009, Table 1): t_min °C, t_max °C, max w.
LALIBERTE_NAOH_VALID = {
    "density": (4.0, 120.0, 0.5029),
    "viscosity": (12.5, 70.0, 0.56),
    "heat_capacity": (4.0, 120.0, 0.3035),
}
# Liquid water specific heat at 0.1 MPa, 20–70 °C (IAPWS-95: 4178–4190 J/(kg·K), ±0.2 %).
CP_WATER_J_kgK = 4182.0


def _water_density_kg_m3(t_C: float) -> float:
    """Pure water density, Kell (1975) form used by Laliberté & Cooper (2004)."""
    num = (
        (((-2.8054253e-10 * t_C + 1.0556302e-7) * t_C - 4.6170461e-5) * t_C - 0.0079870401) * t_C
        + 16.945176
    ) * t_C + 999.83952
    return num / (1.0 + 0.01687985 * t_C)


def _water_viscosity_Pa_s(t_C: float) -> float:
    """Pure water viscosity, Laliberté (2007) eq. 9 (mPa·s converted to Pa·s)."""
    return (t_C + 246.0) / ((0.05594 * t_C + 5.2842) * t_C + 137.37) * 1e-3


def _density_at_mass_fraction(t_C: float, w: float) -> float:
    c0, c1, c2, c3, c4 = LALIBERTE_NAOH_C
    rho_app = (c0 * w + c1) * math.exp(1e-6 * (t_C + c4) ** 2) / (w + c2 + c3 * t_C)
    return 1.0 / ((1.0 - w) / _water_density_kg_m3(t_C) + w / rho_app)


def mass_fraction_NaOH(T_K: float, C_NaOH_M: float) -> float:
    """NaOH mass fraction at molarity C and temperature T: w = C·M_NaOH/ρ(T, w), iterated.

    Equations: solvent.density.
    """
    t_C = T_K - 273.15
    c_kg_m3 = 1000.0 * C_NaOH_M * MW_kg_mol["NaOH"]
    w = c_kg_m3 / 1000.0
    for _ in range(50):
        w_new = c_kg_m3 / _density_at_mass_fraction(t_C, w)
        if abs(w_new - w) < 1e-13:
            return w_new
        w = w_new
    return w


def density_kg_m3(T_K: float, C_NaOH_M: float) -> float:
    """NaOH solution density, Laliberté & Cooper (2004) model with Laliberté (2009) coefficients.

    Equations: solvent.density.
    """
    return _density_at_mass_fraction(T_K - 273.15, mass_fraction_NaOH(T_K, C_NaOH_M))


def viscosity_Pa_s(T_K: float, C_NaOH_M: float) -> float:
    """NaOH solution viscosity, Laliberté (2007) model with Laliberté (2009) coefficients.

    ln μ = w_w·ln μ_w + w·ln μ_s, μ_s = exp[(v1·w^v2 + v3)/(v4·t + 1)]/(v5·w^v6 + 1), t in °C.

    Equations: solvent.viscosity.
    """
    t_C = T_K - 273.15
    w = mass_fraction_NaOH(T_K, C_NaOH_M)
    mu_w = _water_viscosity_Pa_s(t_C)
    if w <= 0.0:
        return mu_w
    v1, v2, v3, v4, v5, v6 = LALIBERTE_NAOH_V
    mu_s = math.exp((v1 * w**v2 + v3) / (v4 * t_C + 1.0)) / (v5 * w**v6 + 1.0) * 1e-3
    return math.exp((1.0 - w) * math.log(mu_w) + w * math.log(mu_s))


def heat_capacity_J_kgK(T_K: float, C_NaOH_M: float) -> float:
    """NaOH solution specific heat, Laliberté (2009): c_p = w_w·c_p,w + w·c_p,app.

    c_p,app = a1·exp(α) + a5·(1 − w_w)^a6, α = a2·t + a3·exp(0.01·t) + a4·(1 − w_w), kJ/(kg·K).

    Equations: solvent.heat_capacity.
    """
    t_C = T_K - 273.15
    w = mass_fraction_NaOH(T_K, C_NaOH_M)
    a1, a2, a3, a4, a5, a6 = LALIBERTE_NAOH_A
    alpha = a2 * t_C + a3 * math.exp(0.01 * t_C) + a4 * w
    cp_app = (a1 * math.exp(alpha) + a5 * w**a6) * 1000.0
    return (1.0 - w) * CP_WATER_J_kgK + w * cp_app


# CO2 diffusivity in water, Versteeg & van Swaaij (1988), J. Chem. Eng. Data 33, 29–34.
D_CO2_WATER_PREEXP_m2_s = 2.35e-6
D_CO2_WATER_E_K = 2119.0
# Modified Stokes–Einstein relation of the same paper: D·μ^0.8 = constant.
D_CO2_VISCOSITY_EXPONENT = 0.8


def diffusivity_CO2_m2_s(T_K: float, C_NaOH_M: float) -> float:
    """CO2 diffusivity in the alkaline solution: water value corrected with D·μ^0.8 = const.

    Equations: solvent.diffusivity_co2.
    """
    D_w = D_CO2_WATER_PREEXP_m2_s * math.exp(-D_CO2_WATER_E_K / T_K)
    mu_ratio = _water_viscosity_Pa_s(T_K - 273.15) / viscosity_Pa_s(T_K, C_NaOH_M)
    return D_w * mu_ratio**D_CO2_VISCOSITY_EXPONENT


def property_range_warnings(T_K: float, C_NaOH_M: float) -> list[str]:
    """Warnings when the Laliberté correlations are used outside their fitted range."""
    t_C = T_K - 273.15
    w = mass_fraction_NaOH(T_K, C_NaOH_M)
    out = []
    for name, (t_lo, t_hi, w_hi) in LALIBERTE_NAOH_VALID.items():
        if not (t_lo <= t_C <= t_hi) or w > w_hi:
            out.append(
                f"NaOH {name} (Laliberté 2009) extrapolated: {t_C:.1f} °C, w = {w:.3f} "
                f"outside {t_lo:g}–{t_hi:g} °C, w ≤ {w_hi:g}"
            )
    return out


# CO2 solubility in pure water (Sander 2015, Atmos. Chem. Phys. 15, 4399, recommended value).
HCP_CO2_WATER_298_mol_m3Pa = 3.3e-4
HCP_CO2_dlnH_d1T_K = 2400.0

# Weisenberger & Schumpe (1996), AIChE J. 42, 298: ion-specific h_i [m3/kmol] and the
# gas-specific parameter of CO2, h_G = h_G0 + h_T·(T − 298.15 K).
SCHUMPE_H_ION_m3_kmol = {
    "Na+": 0.1143,
    "K+": 0.0922,
    "OH-": 0.0839,
    "HCO3-": 0.0967,
    "CO3--": 0.1423,
}
SCHUMPE_HG0_CO2_m3_kmol = -0.0172
SCHUMPE_HT_CO2_m3_kmolK = -0.338e-3


def henry_cp_CO2_water_mol_m3Pa(T_K: float) -> float:
    """CO2 solubility in water H_cp = c_L/p [mol/(m3·Pa)], van 't Hoff form around 298.15 K.

    Provenance: literature (Sander 2015). Validity: about 273–353 K.

    Equations: solvent.henry_water.
    """
    return HCP_CO2_WATER_298_mol_m3Pa * math.exp(HCP_CO2_dlnH_d1T_K * (1.0 / T_K - 1.0 / 298.15))


def salting_out_log10(T_K: float, ions_kmol_m3: dict[str, float]) -> float:
    """Sechenov term log10(H_cp,water / H_cp,solution) = Σ (h_i + h_G)·c_i.

    Provenance: literature (Weisenberger & Schumpe 1996). Validity: 273–363 K, ionic strength
    up to about 5 kmol/m3.

    Equations: solvent.salting_out.
    """
    h_G = SCHUMPE_HG0_CO2_m3_kmol + SCHUMPE_HT_CO2_m3_kmolK * (T_K - 298.15)
    return sum((SCHUMPE_H_ION_m3_kmol[ion] + h_G) * c for ion, c in ions_kmol_m3.items())


def henry_cc_CO2(T_K: float, ions_kmol_m3: dict[str, float]) -> float:
    """Dimensionless Henry constant H_cc = c_G/c_L of CO2 in the electrolyte solution.

    H_cc = 1/(H_cp·R·T), with H_cp corrected for salting-out (audit §3, errata E-003).

    Equations: solvent.henry_cc.
    """
    H_cp = henry_cp_CO2_water_mol_m3Pa(T_K) * 10.0 ** (-salting_out_log10(T_K, ions_kmol_m3))
    return 1.0 / (H_cp * R_J_molK * T_K)


def naoh_ions_kmol_m3(C_NaOH_M: float) -> dict[str, float]:
    """Fresh NaOH solution: [Na+] = [OH-] = C (mol/L = kmol/m3).

    LEGACY(audit §4): carbonate formed along the column is not accounted for yet.
    """
    return {"Na+": C_NaOH_M, "OH-": C_NaOH_M}


def carbonate_constants(T_K: float) -> tuple[float, float, float]:
    """Ideal-solution constants K1, K2 [mol/L] of carbonic acid and Kw [mol2/L2] of water.

    K1, K2: Plummer & Busenberg (1982), Geochim. Cosmochim. Acta 46, 1011 (0–90 °C).
    Kw: Harned & Owen form, log Kw = −4470.99/T + 6.0875 − 0.01706·T (0–60 °C).
    Provenance: literature. At 25 °C: pK1 6.352, pK2 10.329, pKw 13.995.

    Equations: solvent.carbonate_constants.
    """
    lgT = math.log10(T_K)
    log_K1 = -356.3094 - 0.06091964 * T_K + 21834.37 / T_K + 126.8339 * lgT - 1684915.0 / T_K**2
    log_K2 = -107.8871 - 0.03252849 * T_K + 5151.79 / T_K + 38.92561 * lgT - 563713.9 / T_K**2
    log_Kw = -4470.99 / T_K + 6.0875 - 0.01706 * T_K
    return 10.0**log_K1, 10.0**log_K2, 10.0**log_Kw


def free_CO2_equilibrium_mol_m3(T_K: float, OH_mol_m3: float, CO3_mol_m3: float) -> float:
    """Dissolved CO2 in equilibrium with OH-/CO3-- : [CO2] = Kw²·[CO3--]/(K1·K2·[OH-]²).

    Ideal solution (activity coefficients = 1). Infinite when no hydroxide is left.

    Equations: solvent.free_co2.
    """
    if OH_mol_m3 <= 0.0:
        return math.inf
    K1, K2, Kw = carbonate_constants(T_K)
    OH_M, CO3_M = OH_mol_m3 / 1000.0, CO3_mol_m3 / 1000.0
    return 1000.0 * Kw**2 * CO3_M / (K1 * K2 * OH_M**2)


# Limiting (infinite-dilution) diffusivity of NaOH at 25 °C, CRC Handbook of Chemistry and Physics.
D_NAOH_298_m2_s = 2.13e-9


def diffusivity_NaOH_m2_s(T_K: float, C_NaOH_M: float) -> float:
    """NaOH (OH- with its counter-ion) diffusivity, Stokes–Einstein scaling D·μ/T = const.

    Provenance: literature value + theory scaling, using this module's viscosity correlation.

    Equations: solvent.naoh_diffusivity.
    """
    return (
        D_NAOH_298_m2_s
        * (T_K / 298.15)
        * viscosity_Pa_s(298.15, 0.0)
        / viscosity_Pa_s(T_K, C_NaOH_M)
    )


# Heat of absorption, CO2(g) + 2 OH-(aq) -> CO3--(aq) + H2O(l), from standard enthalpies of
# formation (Wagman et al. 1982, NBS tables): −677.1 − 285.8 − (−393.5 − 2·230.0) kJ/mol.
DH_ABS_CARBONATE_J_mol = -109.4e3


def surface_tension_N_m(T_K: float, C_NaOH_M: float) -> float:
    """Surface tension of the NaOH solution.

    Water: IAPWS (2014) σ = 0.2358·τ^1.256·(1 − 0.625·τ), τ = 1 − T/647.096 (literature).
    NaOH raises σ by about 1.8 mN/m per mol/L (own closure, typical of 1:1 hydroxide salts).

    Equations: solvent.surface_tension.
    """
    tau = 1.0 - T_K / 647.096
    return 0.2358 * tau**1.256 * (1.0 - 0.625 * tau) + 1.8e-3 * C_NaOH_M
