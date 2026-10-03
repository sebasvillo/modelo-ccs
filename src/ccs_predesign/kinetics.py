"""CO2 + OH- kinetics: rate constant, Hatta number and enhancement factor."""

import math


def k2_m3_mol_s(T_K: float, k2_ref_L_mol_s: float) -> float:
    """Second-order rate constant [m3/mol/s], scaled as (T/298.15)^1.5 from k2_ref [L/mol/s]."""
    k2_L = k2_ref_L_mol_s * (T_K / 298.15) ** 1.5
    return k2_L / 1000.0


def hatta_number(k1_pseudo_1_s: float, D_l_m2_s: float, kL_m_s: float) -> float:
    return math.sqrt(max(k1_pseudo_1_s * D_l_m2_s, 1e-24)) / max(kL_m_s, 1e-12)


def enhancement_factor(Ha: float) -> float:
    """Pseudo-first-order enhancement E = Ha/tanh(Ha), E = Ha for Ha >= 50.

    LEGACY(audit §9): no instantaneous-reaction limit E_inf.
    """
    Ha = max(Ha, 1e-12)
    if Ha < 50.0:
        return Ha / math.tanh(Ha)
    return Ha
