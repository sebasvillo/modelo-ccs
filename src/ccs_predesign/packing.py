"""Packing database plus mass-transfer and hydraulic correlations (legacy forms)."""

import math

import numpy as np

from . import gas, solvent
from .constants import G_m_s2
from .models import Packing

PACKINGS: dict[str, Packing] = {
    "raschig_metal_25mm": Packing(
        a_spec_m2_m3=200.0,
        void_fraction=0.94,
        dp_eq_m=0.025,
        wetting_ref_uL_m_s=0.003,
        flood_coeff=0.42,
    ),
    "pall_ring_25mm": Packing(
        a_spec_m2_m3=210.0,
        void_fraction=0.92,
        dp_eq_m=0.025,
        wetting_ref_uL_m_s=0.003,
        flood_coeff=0.45,
    ),
    "structured_250Y": Packing(
        a_spec_m2_m3=250.0,
        void_fraction=0.97,
        dp_eq_m=0.008,
        wetting_ref_uL_m_s=0.0015,
        flood_coeff=0.60,
    ),
    "random_high_capacity": Packing(
        a_spec_m2_m3=160.0,
        void_fraction=0.95,
        dp_eq_m=0.035,
        wetting_ref_uL_m_s=0.004,
        flood_coeff=0.40,
    ),
}


def reynolds(rho_kg_m3: float, u_m_s: float, d_m: float, mu_Pa_s: float) -> float:
    return rho_kg_m3 * u_m_s * d_m / max(mu_Pa_s, 1e-12)


def schmidt(mu_Pa_s: float, rho_kg_m3: float, D_m2_s: float) -> float:
    return mu_Pa_s / max(rho_kg_m3 * D_m2_s, 1e-18)


def sherwood_gas(Re: float, Sc: float) -> float:
    """LEGACY(audit §9): Ranz–Marshall/Wakao particle form, cited in the thesis as Onda."""
    return 2.0 + 1.10 * (Re**0.70) * (Sc ** (1.0 / 3.0))


def sherwood_liquid(Re: float, Sc: float) -> float:
    """LEGACY(audit §9): same particle form as sherwood_gas, exponent 0.60."""
    return 2.0 + 1.10 * (Re**0.60) * (Sc ** (1.0 / 3.0))


def wetting_fraction(uL_m_s: float, packing: Packing) -> float:
    """LEGACY(audit §9): exponential own closure, not Onda's wetting correlation."""
    u_ref = packing.wetting_ref_uL_m_s
    return float(np.clip(0.20 + 0.80 * (1.0 - np.exp(-uL_m_s / max(u_ref, 1e-9))), 0.20, 1.0))


def liquid_film_kL(
    T_K: float, C_NaOH_M: float, uL_m_s: float, packing: Packing
) -> tuple[float, float]:
    """Liquid-film coefficient kL [m/s] and CO2 diffusivity D_l [m2/s] at the liquid temperature.

    Same correlation as mass_transfer_coefficients, evaluated on its own so the absorber can
    follow the local liquid temperature (errata E-013).
    """
    rho_l = solvent.density_kg_m3(T_K, C_NaOH_M)
    mu_l = solvent.viscosity_Pa_s(T_K, C_NaOH_M)
    Dl = solvent.diffusivity_CO2_m2_s(T_K, C_NaOH_M)
    d_h = 4.0 * packing.void_fraction / max(packing.a_spec_m2_m3, 1e-12)
    ReL = reynolds(rho_l, max(uL_m_s, 1e-12), d_h, mu_l)
    ScL = schmidt(mu_l, rho_l, Dl)
    return sherwood_liquid(ReL, ScL) * Dl / d_h, Dl


def mass_transfer_coefficients(
    T_K: float,
    P_Pa: float,
    comp_gas: dict[str, float],
    uG_m_s: float,
    uL_m_s: float,
    packing: Packing,
    C_NaOH_M: float,
    henry_cc: float,
    E: float = 1.0,
) -> dict[str, float]:
    """Film model: kG, kL and the overall volumetric coefficient KGa [1/s].

    uG_m_s and uL_m_s are superficial velocities in the actual column (audit §2, E-002).
    Overall gas-side coefficient 1/K_G = 1/k_G + H_cc/(E·k_L), with H_cc = c_G/c_L the
    dimensionless Henry constant (two-film theory; audit §3, errata E-003).
    """
    rho_g = gas.density_ideal_kg_m3(T_K, P_Pa, comp_gas)
    mu_g = gas.viscosity_sutherland_air_Pa_s(T_K)
    rho_l = solvent.density_kg_m3(T_K, C_NaOH_M)
    mu_l = solvent.viscosity_Pa_s(T_K, C_NaOH_M)
    Dg = gas.diffusivity_CO2_in_air_m2_s(T_K, P_Pa)
    Dl = solvent.diffusivity_CO2_m2_s(T_K, C_NaOH_M)

    eps = packing.void_fraction
    a_spec = packing.a_spec_m2_m3
    d_h = 4.0 * eps / max(a_spec, 1e-12)

    uG = uG_m_s
    uL = max(uL_m_s, 1e-12)

    ReG = reynolds(rho_g, uG, d_h, mu_g)
    ReL = reynolds(rho_l, uL, d_h, mu_l)
    ScG = schmidt(mu_g, rho_g, Dg)
    ScL = schmidt(mu_l, rho_l, Dl)
    ShG = sherwood_gas(ReG, ScG)
    ShL = sherwood_liquid(ReL, ScL)
    kG = ShG * Dg / d_h
    kL = ShL * Dl / d_h

    wet = wetting_fraction(uL, packing)
    a_eff = a_spec * wet

    KG = 1.0 / (1.0 / max(kG, 1e-12) + henry_cc / max(kL * E, 1e-12))
    KGa = KG * a_eff

    return {
        "rho_g": rho_g,
        "mu_g": mu_g,
        "rho_l": rho_l,
        "mu_l": mu_l,
        "Dg": Dg,
        "Dl": Dl,
        "d_h": d_h,
        "uG": uG,
        "uL": uL,
        "ReG": ReG,
        "ReL": ReL,
        "ScG": ScG,
        "ScL": ScL,
        "ShG": ShG,
        "ShL": ShL,
        "kG": kG,
        "kL": kL,
        "a_eff": a_eff,
        "wet": wet,
        "KG": KG,
        "KGa": KGa,
    }


def pressure_drop_ergun_Pa_m(rho_g: float, mu_g: float, uG_m_s: float, packing: Packing) -> float:
    """LEGACY(audit §9): Ergun is for particle beds, also applied to structured packings."""
    eps = packing.void_fraction
    dp = packing.dp_eq_m
    term1 = 150.0 * (1.0 - eps) ** 2 * mu_g * uG_m_s / max((eps**3) * dp**2, 1e-18)
    term2 = 1.75 * (1.0 - eps) * rho_g * uG_m_s**2 / max((eps**3) * dp, 1e-18)
    return term1 + term2


def flooding_velocity_m_s(rho_g: float, rho_l: float, packing: Packing) -> float:
    eps = packing.void_fraction
    dp = packing.dp_eq_m
    C = packing.flood_coeff
    return (
        C * math.sqrt(max((rho_l - rho_g), 1e-12) * G_m_s2 * dp / max(rho_g, 1e-12)) * (eps**0.15)
    )


def column_diameter_m(
    Qg_m3_s: float, rho_g: float, rho_l: float, packing: Packing, flood_frac: float
) -> tuple[float, float, float]:
    """Diameter from a fraction of the flooding velocity. Returns (D, v_flood, v_oper)."""
    v_flood = flooding_velocity_m_s(rho_g, rho_l, packing)
    v_oper = max(0.05, flood_frac * v_flood)
    A = Qg_m3_s / v_oper
    D = math.sqrt(4.0 * A / math.pi)
    return D, v_flood, v_oper
