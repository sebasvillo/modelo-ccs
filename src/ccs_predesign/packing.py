"""Packings: Billet & Schultes (1999) database, hydraulics and mass transfer.

Billet & Schultes, Trans IChemE 77A (1999) 498–504, is the default method: it is the most reliable
general model according to the most recent review in the bibliography (Flagiello et al. 2021,
ChemEngineering 5, 43). Onda et al. (1968) is kept as an alternative mass-transfer model for
comparison (errata E-016, E-017).
"""

import math
from typing import Literal

from scipy import optimize

from . import gas, solvent
from .constants import G_m_s2
from .models import Packing

# Critical surface tension of packing materials for Onda (1968) [N/m].
SIGMA_C_N_m = {"metal": 0.075, "ceramic": 0.061, "plastic": 0.033}

# Billet & Schultes (1999) Table 2a (dumped) and 2b (regular): a [m2/m3], ε, C_S, C_Fl, C_P0,
# C_L, C_V. Only packings with the full set of constants are listed.
PACKINGS: dict[str, Packing] = {
    "pall_ring_metal_25mm": Packing(
        kind="random", material="metal", nominal_size_m=0.025, a_spec_m2_m3=223.5,
        void_fraction=0.954, C_S=2.627, C_Fl=2.083, C_P0=0.957, C_L=1.440, C_V=0.336,
    ),
    "pall_ring_metal_50mm": Packing(
        kind="random", material="metal", nominal_size_m=0.050, a_spec_m2_m3=112.6,
        void_fraction=0.951, C_S=2.725, C_Fl=1.580, C_P0=0.763, C_L=1.192, C_V=0.410,
    ),
    "raschig_super_ring_metal_no1": Packing(
        kind="random", material="metal", nominal_size_m=0.025, a_spec_m2_m3=160.0,
        void_fraction=0.980, C_S=3.491, C_Fl=2.200, C_P0=0.500, C_L=1.290, C_V=0.440,
    ),
    "hiflow_ring_metal_50mm": Packing(
        kind="random", material="metal", nominal_size_m=0.050, a_spec_m2_m3=92.3,
        void_fraction=0.977, C_S=2.702, C_Fl=1.626, C_P0=0.421, C_L=1.168, C_V=0.408,
    ),
    "raschig_ring_ceramic_25mm": Packing(
        kind="random", material="ceramic", nominal_size_m=0.025, a_spec_m2_m3=190.0,
        void_fraction=0.680, C_S=2.454, C_Fl=1.899, C_P0=1.329, C_L=1.361, C_V=0.412,
    ),
    "ralu_pak_metal_yc250": Packing(
        kind="structured", material="metal", nominal_size_m=None, a_spec_m2_m3=250.0,
        void_fraction=0.945, C_S=3.178, C_Fl=2.558, C_P0=0.191, C_L=1.334, C_V=0.385,
    ),
}  # fmt: skip

# Ranges of the Billet & Schultes database (their Table 1).
BS_MASS_TRANSFER_F_MAX = 2.77  # Pa^0.5
BS_MASS_TRANSFER_UL_RANGE_m3_m2h = (0.256, 118.0)


def hydraulic_diameter_m(packing: Packing) -> float:
    """d_h = 4ε/a. Equations: packing.bs_holdup."""
    return 4.0 * packing.void_fraction / packing.a_spec_m2_m3


def bs_holdup_below_loading(uL_m_s: float, rho_l: float, mu_l: float, packing: Packing) -> float:
    """Liquid holdup below the loading point, h_L = (12·η_L·u_L·a²/(g·ρ_L))^(1/3) (BS eq. 11).

    Equations: packing.bs_holdup.
    """
    return (12.0 * mu_l * uL_m_s * packing.a_spec_m2_m3**2 / (G_m_s2 * rho_l)) ** (1.0 / 3.0)


def _flow_parameter(LV_mass: float, rho_g: float, rho_l: float) -> float:
    return LV_mass * math.sqrt(rho_g / rho_l)


def bs_flooding(
    LV_mass: float, rho_g: float, mu_g: float, rho_l: float, mu_l: float, packing: Packing
) -> tuple[float, float]:
    """Gas velocity and liquid holdup at the flooding point (BS eqs. 13 and 36–40).

    u_V,Fl = √2·√(g/ψ_Fl)·(ε − h_Fl)^1.5/ε^0.5·√(h_Fl/a)·√(ρ_L/ρ_V), with ψ_Fl from the flow
    parameter and h_Fl from h³(3h − ε) = 6/g·a²·ε·(η_L/ρ_L)·(L/V)·(ρ_V/ρ_L)·u_V,Fl.

    Equations: packing.bs_flooding.
    """
    eps, a = packing.void_fraction, packing.a_spec_m2_m3
    X = _flow_parameter(LV_mass, rho_g, rho_l)
    if X <= 0.4:
        n_fl, C_fl = -0.194, packing.C_Fl
    else:
        n_fl, C_fl = -0.708, 0.6244 * packing.C_Fl * (mu_l / mu_g) ** 0.1028
    psi_fl = G_m_s2 / C_fl**2 * (X * (mu_l / mu_g) ** 0.2) ** (-2.0 * n_fl)
    R = 6.0 / G_m_s2 * a**2 * eps * (mu_l / rho_l) * LV_mass * (rho_g / rho_l)

    def holdup(u: float) -> float:
        rhs = R * u
        if rhs >= 2.0 * eps**4:
            return eps
        return optimize.brentq(lambda h: h**3 * (3.0 * h - eps) - rhs, eps / 3.0, eps)

    def velocity(h: float) -> float:
        return (
            math.sqrt(2.0 * G_m_s2 / psi_fl)
            * (eps - h) ** 1.5
            / math.sqrt(eps)
            * math.sqrt(h / a)
            * math.sqrt(rho_l / rho_g)
        )

    u_hi = velocity(eps / 3.0)
    u_fl = optimize.brentq(lambda u: u - velocity(holdup(u)), 1e-9 * u_hi, u_hi, xtol=1e-12)
    return u_fl, holdup(u_fl)


def bs_loading_velocity(
    LV_mass: float, rho_g: float, mu_g: float, rho_l: float, mu_l: float, packing: Packing
) -> float:
    """Gas velocity at the loading point (BS eqs. 31–35). Equations: packing.bs_loading."""
    eps, a = packing.void_fraction, packing.a_spec_m2_m3
    X = _flow_parameter(LV_mass, rho_g, rho_l)
    if X <= 0.4:
        n_s, C_s = -0.326, packing.C_S
    else:
        n_s, C_s = -0.723, 0.695 * packing.C_S * (mu_l / mu_g) ** 0.1588
    psi_s = G_m_s2 / C_s**2 * (X * (mu_l / mu_g) ** 0.4) ** (-2.0 * n_s)

    def velocity(u_vs: float) -> float:
        u_ls = rho_g / rho_l * LV_mass * u_vs
        film = 12.0 / G_m_s2 * mu_l / rho_l * u_ls
        return (
            math.sqrt(G_m_s2 / psi_s)
            * (eps / a ** (1.0 / 6.0) - math.sqrt(a) * film ** (1.0 / 3.0))
            * film ** (1.0 / 6.0)
            * math.sqrt(rho_l / rho_g)
        )

    return optimize.brentq(lambda u: u - velocity(u), 1e-9, 50.0, xtol=1e-12)


def bs_holdup(uG_m_s: float, u_flood_m_s: float, h_L_S: float, h_L_Fl: float) -> float:
    """Holdup between loading and flooding, h_L = h_S + (h_Fl − h_S)·(u_V/u_V,Fl)^13 (BS eq. 12).

    Below the loading point the second term is negligible. Equations: packing.bs_holdup.
    """
    return h_L_S + (h_L_Fl - h_L_S) * (uG_m_s / u_flood_m_s) ** 13


def bs_pressure_drop_Pa_m(
    uG_m_s: float,
    uL_m_s: float,
    rho_g: float,
    mu_g: float,
    h_L: float,
    h_L_S: float,
    D_col_m: float,
    packing: Packing,
) -> float:
    """Irrigated pressure drop Δp/H = ψ_L·a/(ε − h_L)³·F_V²/2·1/K (BS eqs. 23–30).

    Equations: packing.bs_pressure_drop.
    """
    eps, a = packing.void_fraction, packing.a_spec_m2_m3
    d_p = 6.0 * (1.0 - eps) / a
    K = 1.0 / (1.0 + 2.0 / 3.0 / (1.0 - eps) * d_p / D_col_m)
    Re_V = uG_m_s * d_p / ((1.0 - eps) * mu_g / rho_g) * K
    Fr_L = uL_m_s**2 * a / G_m_s2
    C1 = 13300.0 / a**1.5
    psi_L = (
        packing.C_P0
        * (64.0 / Re_V + 1.8 / Re_V**0.08)
        * ((eps - h_L) / eps) ** 1.5
        * (h_L / h_L_S) ** 0.3
        * math.exp(C1 * math.sqrt(Fr_L))
    )
    F_V = uG_m_s * math.sqrt(rho_g)
    return psi_L * a / (eps - h_L) ** 3 * F_V**2 / 2.0 / K


def column_diameter_m(Qg_m3_s: float, u_flood_m_s: float, flood_frac: float) -> float:
    """D from operating at a fraction of the flooding velocity.

    Equations: packing.column_diameter.
    """
    return math.sqrt(4.0 * Qg_m3_s / (flood_frac * u_flood_m_s) / math.pi)


def bs_wetted_area_m2_m3(
    uL_m_s: float, rho_l: float, mu_l: float, sigma_l: float, packing: Packing
) -> float:
    """a_Ph/a = 1.5 (a·d_h)^−0.5 Re_L^−0.2 We_L^0.75 Fr_L^−0.45 below loading (BS eq. 9).

    σ_L is floored at 0.03 N/m as Billet & Schultes recommend. Equations: packing.bs_area.
    """
    a = packing.a_spec_m2_m3
    d_h = hydraulic_diameter_m(packing)
    sigma = max(sigma_l, 0.03)
    Re = uL_m_s * d_h * rho_l / mu_l
    We = uL_m_s**2 * rho_l * d_h / sigma
    Fr = uL_m_s**2 / (G_m_s2 * d_h)
    return a * 1.5 * (a * d_h) ** -0.5 * Re**-0.2 * We**0.75 * Fr**-0.45


def bs_kL_m_s(uL_m_s: float, h_L: float, D_l: float, packing: Packing) -> float:
    """β_L = C_L·12^(1/6)·(ū_L·D_L/d_h)^(1/2), ū_L = u_L/h_L (BS eqs. 6–7).

    Equations: packing.bs_kl.
    """
    d_h = hydraulic_diameter_m(packing)
    return packing.C_L * 12.0 ** (1.0 / 6.0) * math.sqrt(uL_m_s / h_L * D_l / d_h)


def bs_kG_m_s(
    uG_m_s: float, h_L: float, rho_g: float, mu_g: float, D_g: float, packing: Packing
) -> float:
    """β_V = C_V·(ε − h_L)^−½·(a/d_h)^½·D_V·(u_V/(a·ν_V))^¾·(ν_V/D_V)^⅓ (BS eq. 8).

    Equations: packing.bs_kg.
    """
    a = packing.a_spec_m2_m3
    nu_g = mu_g / rho_g
    return (
        packing.C_V
        / math.sqrt(packing.void_fraction - h_L)
        * math.sqrt(a / hydraulic_diameter_m(packing))
        * D_g
        * (uG_m_s / (a * nu_g)) ** 0.75
        * (nu_g / D_g) ** (1.0 / 3.0)
    )


def onda_wetted_area_m2_m3(
    L_kg_m2s: float, rho_l: float, mu_l: float, sigma_l: float, packing: Packing
) -> float:
    """Wetted area a_w = a·{1 − exp[−1.45 (σc/σ)^0.75 Re^0.1 Fr^−0.05 We^0.2]} (Onda 1968).

    Re = L/(a·μ), Fr = L²·a/(ρ²·g), We = L²/(ρ·σ·a), L the liquid mass flux [kg/(m2·s)].
    Provenance: literature (Onda, Takeuchi & Okumoto 1968, J. Chem. Eng. Japan 1, 56).

    Equations: packing.onda_wetted_area.
    """
    a = packing.a_spec_m2_m3
    sigma_c = SIGMA_C_N_m[packing.material]
    Re = L_kg_m2s / (a * mu_l)
    Fr = L_kg_m2s**2 * a / (rho_l**2 * G_m_s2)
    We = L_kg_m2s**2 / (rho_l * sigma_l * a)
    x = 1.45 * (sigma_c / sigma_l) ** 0.75 * Re**0.1 * Fr**-0.05 * We**0.2
    return a * (1.0 - math.exp(-x))


def onda_kL_m_s(
    L_kg_m2s: float, a_w: float, rho_l: float, mu_l: float, D_l: float, packing: Packing
) -> float:
    """kL·(ρ/(μ·g))^(1/3) = 0.0051 (L/(a_w·μ))^(2/3) Sc^(−1/2) (a·dp)^0.4 (Onda 1968).

    Equations: packing.onda_kl.
    """
    Sc = mu_l / (rho_l * D_l)
    group = 0.0051 * (L_kg_m2s / (a_w * mu_l)) ** (2.0 / 3.0) * Sc**-0.5
    ad = packing.a_spec_m2_m3 * packing.nominal_size_m
    return group * ad**0.4 * (mu_l * G_m_s2 / rho_l) ** (1 / 3)


def onda_kG_m_s(G_kg_m2s: float, rho_g: float, mu_g: float, D_g: float, packing: Packing) -> float:
    """kG/(a·D_G) = C (G/(a·μ))^0.7 Sc^(1/3) (a·dp)^−2, C = 5.23 (2.00 below 15 mm) (Onda 1968).

    kG in concentration units [m/s] (kG,p·R·T of the original).

    Equations: packing.onda_kg.
    """
    a = packing.a_spec_m2_m3
    dp = packing.nominal_size_m
    C = 5.23 if dp >= 0.015 else 2.00
    Sc = mu_g / (rho_g * D_g)
    return C * a * D_g * (G_kg_m2s / (a * mu_g)) ** 0.7 * Sc ** (1 / 3) * (a * dp) ** -2.0


MassTransferModel = Literal["billet_schultes", "onda"]


def gas_side_kG_m_s(
    T_K: float,
    P_Pa: float,
    comp_gas: dict[str, float],
    uG_m_s: float,
    h_L: float,
    packing: Packing,
    model: MassTransferModel = "billet_schultes",
) -> float:
    """Gas-film coefficient at the gas conditions.

    Equations: packing.bs_kg, packing.onda_kg.
    """
    rho_g = gas.density_ideal_kg_m3(T_K, P_Pa, comp_gas)
    mu_g = gas.viscosity_sutherland_air_Pa_s(T_K)
    Dg = gas.diffusivity_CO2_in_air_m2_s(T_K, P_Pa)
    if model == "onda":
        return onda_kG_m_s(rho_g * uG_m_s, rho_g, mu_g, Dg, packing)
    return bs_kG_m_s(uG_m_s, h_L, rho_g, mu_g, Dg, packing)


def liquid_side(
    T_K: float,
    C_NaOH_M: float,
    uL_m_s: float,
    packing: Packing,
    model: MassTransferModel = "billet_schultes",
) -> tuple[float, float, float, float]:
    """kL [m/s], CO2 diffusivity D_l [m2/s], interfacial area a_e [m2/m3] and holdup h_L at T.

    Equations: packing.bs_kl, packing.bs_area, packing.bs_holdup, packing.onda_kl,
    packing.onda_wetted_area.
    """
    rho_l = solvent.density_kg_m3(T_K, C_NaOH_M)
    mu_l = solvent.viscosity_Pa_s(T_K, C_NaOH_M)
    D_l = solvent.diffusivity_CO2_m2_s(T_K, C_NaOH_M)
    sigma = solvent.surface_tension_N_m(T_K, C_NaOH_M)
    h_L = bs_holdup_below_loading(uL_m_s, rho_l, mu_l, packing)
    if model == "onda":
        L = rho_l * uL_m_s
        a_w = onda_wetted_area_m2_m3(L, rho_l, mu_l, sigma, packing)
        return onda_kL_m_s(L, a_w, rho_l, mu_l, D_l, packing), D_l, a_w, h_L
    a_e = bs_wetted_area_m2_m3(uL_m_s, rho_l, mu_l, sigma, packing)
    return bs_kL_m_s(uL_m_s, h_L, D_l, packing), D_l, a_e, h_L
