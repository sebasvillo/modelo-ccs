"""Electrochemical solvent regeneration (PSE cell), legacy sizing."""

import math

import numpy as np

from .constants import FARADAY_C_mol
from .models import CellResult, CellSpec


def tna_empirical(C_NaOH_mol_m3: float, j_A_m2: float) -> float:
    """Empirical Na+ transference number vs NaOH concentration and current density."""
    C_M = C_NaOH_mol_m3 / 1000.0
    j_mA_cm2 = j_A_m2 * 1e-1
    t = 0.70 + 0.20 * (1.0 - math.exp(-4.0 * C_M))
    t -= 0.05 * min(j_mA_cm2, 250.0) / 250.0
    return float(np.clip(t, 0.10, 0.98))


def voltage_from_current_density_V(j_A_m2: float) -> float:
    """Inverse of the fitted j–V curve j[mA/cm2] = 604.92445·exp(V/3.02912) − 767.22373."""
    j_mA_cm2 = j_A_m2 * 1e-1
    return 3.02912 * math.log((j_mA_cm2 + 767.22373) / 604.92445)


def energy_per_mol_CO2_J(V_cell_V: float, tna: float, electrons_per_CO2: float) -> float:
    """Electrical energy per mol of CO2 regenerated: n_e·F·V/t_Na [J/mol] (Faraday)."""
    return electrons_per_CO2 * FARADAY_C_mol * V_cell_V / max(tna, 1e-12)


def electrochemical_regeneration(
    CO2_mol_s: float, NaOH_M: float, spec: CellSpec, *, electrons_per_CO2: float
) -> CellResult:
    """Size the cell for the captured CO2 flow.

    Each mol of OH- regenerated needs one Na+ across the membrane, i.e. one Faraday divided by
    the Na+ transference number: I = n_e·F·n_CO2/t_Na. n_e is the OH- consumed per CO2 in the
    absorber: 2 for the carbonate route, 1 for bicarbonate (audit §6, errata E-008).
    """
    if electrons_per_CO2 not in (1.0, 2.0):
        raise ValueError(
            f"electrons_per_CO2 must be 1 (bicarbonate) or 2 (carbonate), got {electrons_per_CO2}"
        )
    j_A_m2 = spec.j_mA_cm2 * 10.0
    C_NaOH_mol_m3 = NaOH_M * 1000.0
    tna = (
        float(spec.tna_const) if spec.tna_mode == "const" else tna_empirical(C_NaOH_mol_m3, j_A_m2)
    )
    V_cell = (
        voltage_from_current_density_V(j_A_m2)
        if spec.V_override_V is None
        else float(spec.V_override_V)
    )
    A_cell_m2 = electrons_per_CO2 * CO2_mol_s * FARADAY_C_mol / max(tna * j_A_m2, 1e-12)
    I_cell_A = j_A_m2 * A_cell_m2
    CO2_released = CO2_mol_s * spec.CO2_release_eff
    return CellResult(
        j_A_m2=j_A_m2,
        tna=tna,
        V_cell_V=V_cell,
        electrons_per_CO2=electrons_per_CO2,
        E_J_mol=energy_per_mol_CO2_J(V_cell, tna, electrons_per_CO2),
        A_cell_m2=A_cell_m2,
        I_cell_A=I_cell_A,
        P_cell_W=I_cell_A * V_cell,
        CO2_released_mol_s=CO2_released,
        CO2_slip_mol_s=CO2_mol_s - CO2_released,
    )
