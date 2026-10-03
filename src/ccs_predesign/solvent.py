"""NaOH solution properties (simplified correlations from the legacy model)."""


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
