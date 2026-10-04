"""Validation against independent data from the bibliography (see docs/validation.md)."""

import math

import pytest

from ccs_predesign import cell, gas, kinetics, packing
from ccs_predesign.models import CellSpec

ELLEBRACHT_250Y = {"a_eff_m": 170.0, "kL_m_s": 0.23e-3}  # Energy Environ. Sci. 16 (2023) Table 1


def test_mellapak_area_vs_ellebracht_2023():
    """Hanley & Chen a_e for 250Y within ±35 % of the area fitted to NaOH absorption data.

    Conditions: 74 mm column, 28 SLPM gas (10 % CO2), 2.75 L/min of 0.4 M NaOH, ~22 °C. The
    fitted k_L is not compared: with Ha ≈ 11 the flux depends on E·k_L ≈ √(k2[OH-]D), so k_L is
    poorly identified by the fit (the authors say the area dominates).
    """
    pk = packing.PACKINGS["mellapak_metal_250y"]
    T, P, A = 295.15, 1.01325e5, math.pi * 0.074**2 / 4
    uG, uL = 28e-3 / 60 * (T / 293.15) / A, 2.75e-3 / 60 / A
    comp = {"CO2": 0.1, "N2": 0.9, "O2": 0.0, "H2O": 0.0}
    rho_g, mu_g = gas.density_ideal_kg_m3(T, P, comp), gas.viscosity_sutherland_air_Pa_s(T)
    _, _, a_e, _ = packing.liquid_side(T, 0.4, uL, pk, gas_props=(uG, rho_g, mu_g))
    assert a_e == pytest.approx(ELLEBRACHT_250Y["a_eff_m"], rel=0.35)


def test_pall_ring_flooding_capacity():
    """Billet & Schultes flooding of 25 mm metal Pall rings, air/water at low liquid load:
    F ≈ 3 Pa^0.5, the order of vendor capacity charts (Kister 1992)."""
    u, _ = packing.bs_flooding(
        1.0, 1.2, 1.8e-5, 1000.0, 1.0e-3, packing.PACKINGS["pall_ring_metal_25mm"]
    )
    assert 2.5 < u * math.sqrt(1.2) < 3.5


def test_cell_at_zhang_operating_point():
    """Zhang et al. (2024) SI Note 8: 200 mA/cm2 → V = 1.42 V, t_Na = 0.87 (the thesis fits)."""
    j_A_m2 = 2000.0
    assert cell.voltage_from_current_density_V(j_A_m2) == pytest.approx(1.42, abs=0.01)
    assert cell.tna_empirical(1500.0, j_A_m2) == pytest.approx(0.87, abs=0.015)
    c = cell.electrochemical_regeneration(
        1e6 / 86400 / 44.0, 1.5, CellSpec(), electrons_per_CO2=1.0
    )  # 1 t/day, bicarbonate route as in their example
    assert c.A_cell_m2 == pytest.approx(14.58, rel=0.02)


def test_kinetics_reference_value():
    """Pohorecki & Moniuk (1988): k2∞ ≈ 8.0e3 L/(mol·s) at 25 °C."""
    assert kinetics.k2_pohorecki_moniuk_m3_mol_s(298.15, 0.0) * 1000 == pytest.approx(
        8.05e3, rel=0.01
    )
