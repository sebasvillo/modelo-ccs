"""NaOH solution properties: legacy simplified correlations and CO2 solubility."""

import math

from .constants import R_J_molK


def viscosity_Pa_s(T_K: float, C_NaOH_M: float) -> float:
    """Water viscosity (Vogel-type) times a linear NaOH correction."""
    T_C = T_K - 273.15
    mu_water = 2.414e-5 * 10 ** (247.8 / (T_C + 133.15))
    return mu_water * (1.0 + 0.08 * C_NaOH_M)


def diffusivity_CO2_m2_s(T_K: float, C_NaOH_M: float) -> float:
    """CO2 diffusivity in the alkaline solution; decreases with NaOH concentration."""
    D0 = 1.9e-9 * (T_K / 298.15) ** 1.25
    return D0 / (1.0 + 0.30 * C_NaOH_M)


def density_kg_m3(T_K: float, C_NaOH_M: float) -> float:
    """Approximate density of NaOH solution near 1–2 M."""
    return 1000.0 + 24.0 * C_NaOH_M - 0.30 * (T_K - 298.15)


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
    """
    return HCP_CO2_WATER_298_mol_m3Pa * math.exp(HCP_CO2_dlnH_d1T_K * (1.0 / T_K - 1.0 / 298.15))


def salting_out_log10(T_K: float, ions_kmol_m3: dict[str, float]) -> float:
    """Sechenov term log10(H_cp,water / H_cp,solution) = Σ (h_i + h_G)·c_i.

    Provenance: literature (Weisenberger & Schumpe 1996). Validity: 273–363 K, ionic strength
    up to about 5 kmol/m3.
    """
    h_G = SCHUMPE_HG0_CO2_m3_kmol + SCHUMPE_HT_CO2_m3_kmolK * (T_K - 298.15)
    return sum((SCHUMPE_H_ION_m3_kmol[ion] + h_G) * c for ion, c in ions_kmol_m3.items())


def henry_cc_CO2(T_K: float, ions_kmol_m3: dict[str, float]) -> float:
    """Dimensionless Henry constant H_cc = c_G/c_L of CO2 in the electrolyte solution.

    H_cc = 1/(H_cp·R·T), with H_cp corrected for salting-out (audit §3, errata E-003).
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
    """
    lgT = math.log10(T_K)
    log_K1 = -356.3094 - 0.06091964 * T_K + 21834.37 / T_K + 126.8339 * lgT - 1684915.0 / T_K**2
    log_K2 = -107.8871 - 0.03252849 * T_K + 5151.79 / T_K + 38.92561 * lgT - 563713.9 / T_K**2
    log_Kw = -4470.99 / T_K + 6.0875 - 0.01706 * T_K
    return 10.0**log_K1, 10.0**log_K2, 10.0**log_Kw


def free_CO2_equilibrium_mol_m3(T_K: float, OH_mol_m3: float, CO3_mol_m3: float) -> float:
    """Dissolved CO2 in equilibrium with OH-/CO3-- : [CO2] = Kw²·[CO3--]/(K1·K2·[OH-]²).

    Ideal solution (activity coefficients = 1). Infinite when no hydroxide is left.
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

# Specific heat of 1–3 M NaOH solution. Own closure: typical value, ±5 % over that range.
CP_SOLUTION_J_kgK = 3900.0
