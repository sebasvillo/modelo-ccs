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


def test_ellebracht_250Y_capture_is_flagged_outside_hanley_chen_range():
    """Ellebracht et al. (2023) ESI: 74 mm column, 3 × 150 mm of 250Y, v_G 0.105 m/s, v_L 45 m/h,
    10 % CO2, fresh 0.4 M NaOH, 25 °C → initial capture 67 %. The model gives ≈ 40 % because the
    Hanley & Chen gas-film correlation falls below Sh = 2 at this tiny gas load; it must say so.
    """
    from ccs_predesign.absorber import simulate_absorber
    from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent

    A = math.pi * 0.074**2 / 4
    Q_N_h = 0.105 * A * 273.15 / 298.15 * 3600
    r = simulate_absorber(
        FlueGas(Q_dry_Nm3_h=Q_N_h, y_CO2_dry=0.10, y_O2_dry=0.0, T_C=25.0, P_bar=1.01325),
        Solvent(NaOH_M=0.4),
        AbsorberSpec(
            packing_name="mellapak_metal_250y", capture_target=0.999, max_height_m=0.45, dz_m=0.005
        ),
        (45 / 3600) / 0.105,
        D_col_fixed_m=0.074,
    )
    assert 0.3 < r.capture_achieved < 0.67  # documented under-prediction (docs/validation.md)
    assert any("too low for Hanley & Chen" in w for w in r.warnings)
