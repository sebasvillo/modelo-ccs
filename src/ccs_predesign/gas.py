"""Flue-gas conditioning and gas-phase properties."""

from .constants import MOL_PER_NM3, MW_kg_mol, P_STD_Pa, R_J_molK


def c_to_k(T_C: float) -> float:
    return T_C + 273.15


def sat_vapor_pressure_water_Pa(T_K: float) -> float:
    """Antoine equation for water (roughly 1–100 °C), coefficients in mmHg.

    Equations: gas.saturation_pressure.
    """
    T_C = T_K - 273.15
    A, B, C = 8.07131, 1730.63, 233.426
    P_mmHg = 10 ** (A - B / (T_C + C))
    return P_mmHg * 133.322368


def dry_molar_flow_mol_s(Q_dry_Nm3_h: float) -> float:
    """Total dry molar flow [mol/s] from normal volumetric flow [Nm3/h].

    Equations: gas.molar_flow.
    """
    return Q_dry_Nm3_h * MOL_PER_NM3 / 3600.0


def wet_composition_from_dry(
    y_CO2_dry: float,
    y_O2_dry: float,
    T_K: float,
    P_Pa: float,
    y_H2O_cap: float | None = None,
) -> dict[str, float]:
    """Wet-basis mole fractions, humidified to saturation (capped at 0.25) at T_K, P_Pa.

    Equations: gas.wet_composition.
    """
    y_H2O_sat = min(0.25, sat_vapor_pressure_water_Pa(T_K) / P_Pa)
    y_H2O = y_H2O_sat if y_H2O_cap is None else min(y_H2O_cap, y_H2O_sat)
    scale = 1.0 - y_H2O
    y_CO2_wet = y_CO2_dry * scale
    y_O2_wet = y_O2_dry * scale
    y_N2_wet = max(0.0, 1.0 - y_CO2_wet - y_O2_wet - y_H2O)
    total = y_CO2_wet + y_O2_wet + y_N2_wet + y_H2O
    return {
        "CO2": y_CO2_wet / total,
        "O2": y_O2_wet / total,
        "N2": y_N2_wet / total,
        "H2O": y_H2O / total,
    }


def mean_molar_mass_kg_mol(comp: dict[str, float]) -> float:
    return sum(MW_kg_mol[k] * y for k, y in comp.items())


def density_ideal_kg_m3(T_K: float, P_Pa: float, comp: dict[str, float]) -> float:
    """Equations: gas.ideal_density."""
    return P_Pa * mean_molar_mass_kg_mol(comp) / (R_J_molK * T_K)


def viscosity_sutherland_air_Pa_s(T_K: float) -> float:
    """Flue gas approximated as air (Sutherland's law).

    Equations: gas.viscosity.
    """
    mu0 = 1.716e-5
    T0 = 273.15
    S = 111.0
    return mu0 * ((T_K / T0) ** 1.5) * (T0 + S) / (T_K + S)


def diffusivity_CO2_in_air_m2_s(T_K: float, P_Pa: float) -> float:
    """Rough temperature/pressure scaling of the CO2–air diffusivity (1.6e-5 m2/s at 298 K).

    Equations: gas.diffusivity.
    """
    D_ref = 1.6e-5
    return D_ref * (T_K / 298.15) ** 1.75 * (P_STD_Pa / P_Pa)
