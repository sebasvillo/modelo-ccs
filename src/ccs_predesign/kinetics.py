"""CO2 + OH- kinetics: rate constant, Hatta number and enhancement factor."""

import math

# Pohorecki & Moniuk (1988), Chem. Eng. Sci. 43, 1677: CO2 + OH- in electrolyte solutions.
PM_VALID_T_K = (291.0, 314.0)


def k2_pohorecki_moniuk_m3_mol_s(T_K: float, ionic_strength_kmol_m3: float) -> float:
    """Second-order rate constant of CO2 + OH- in NaOH solution [m3/(mol·s)].

    log10 k2_inf = 11.895 − 2382/T  [L/(mol·s)];  log10(k2/k2_inf) = 0.221·I − 0.016·I² (NaOH).
    Provenance: literature (Pohorecki & Moniuk 1988). Validity: 291–314 K (warned outside).
    Replaces the legacy 8.5e3·(T/298.15)^1.5 L/(mol·s) (errata E-012).

    Equations: kinetics.k2.
    """
    ion = ionic_strength_kmol_m3
    log_k2_L = 11.895 - 2382.0 / T_K + 0.221 * ion - 0.016 * ion**2
    return 10.0**log_k2_L / 1000.0


def hatta_number(k1_pseudo_1_s: float, D_l_m2_s: float, kL_m_s: float) -> float:
    """Equations: kinetics.hatta."""
    return math.sqrt(max(k1_pseudo_1_s * D_l_m2_s, 1e-24)) / max(kL_m_s, 1e-12)


def enhancement_factor(Ha: float) -> float:
    """Pseudo-first-order enhancement E = Ha/tanh(Ha), E = Ha for Ha >= 50.

    LEGACY(audit §9): no instantaneous-reaction limit E_inf.

    Equations: kinetics.enhancement_pfo.
    """
    Ha = max(Ha, 1e-12)
    if Ha < 50.0:
        return Ha / math.tanh(Ha)
    return Ha


def enhancement_infinite(
    D_OH_m2_s: float, OH_mol_m3: float, D_CO2_m2_s: float, CO2_i_mol_m3: float, nu: float = 2.0
) -> float:
    """Instantaneous-reaction limit E_inf = 1 + D_OH·[OH-]/(ν·D_CO2·[CO2]_i), ν = 2 for CO2 + 2 OH-.

    Provenance: theory (film theory; Danckwerts 1970, Gas-Liquid Reactions).

    Equations: kinetics.e_infinite.
    """
    return 1.0 + D_OH_m2_s * OH_mol_m3 / (nu * D_CO2_m2_s * max(CO2_i_mol_m3, 1e-30))


def enhancement_factor_decoursey(Ha: float, E_inf: float) -> float:
    """Explicit enhancement factor bounded by E_inf (DeCoursey 1974, Chem. Eng. Sci. 29, 1867).

    E = −Ha²/(2(E_inf−1)) + sqrt(Ha⁴/(4(E_inf−1)²) + E_inf·Ha²/(E_inf−1) + 1).
    Limits: sqrt(1 + Ha²) for E_inf → ∞ (fast pseudo-first order), E_inf for Ha → ∞.
    Provenance: literature (surface-renewal theory). Fixes audit §9, errata E-011.

    Equations: kinetics.decoursey.
    """
    if E_inf <= 1.0 + 1e-12:
        return 1.0
    a = Ha**2 / (E_inf - 1.0)
    return -a / 2.0 + math.sqrt(a**2 / 4.0 + E_inf * a + 1.0)
