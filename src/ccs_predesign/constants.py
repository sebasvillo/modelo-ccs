"""Physical constants and molar masses (SI)."""

R_J_molK = 8.31446261815324
FARADAY_C_mol = 96485.3329
P_STD_Pa = 101325.0
T_STD_K = 273.15
G_m_s2 = 9.81

# Molar volume basis used by the legacy model for normal cubic metres.
# LEGACY(audit §7): 1 Nm3 = 44.615 mol (0 °C, 1 atm); the legacy TEA also uses 22.414 L/mol
# elsewhere. Kept until the TEA is unified.
MOL_PER_NM3 = 44.615

MW_kg_mol = {
    "CO2": 44.01e-3,
    "N2": 28.0134e-3,
    "O2": 31.9988e-3,
    "H2O": 18.01528e-3,
    "air": 28.97e-3,
    "water": 18.01528e-3,
    "NaOH": 40.00e-3,
}
