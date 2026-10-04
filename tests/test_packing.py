"""Billet & Schultes (1999) packing model, with Onda (1968) as alternative (errata E-016, E-017)."""

import math

import pytest

from ccs_predesign import packing, solvent
from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent

PALL = packing.PACKINGS["pall_ring_metal_25mm"]
AIR_WATER = dict(rho_g=1.2, mu_g=1.8e-5, rho_l=1000.0, mu_l=1.0e-3)


def test_database_matches_billet_schultes_table_2a():
    assert (PALL.a_spec_m2_m3, PALL.void_fraction) == (223.5, 0.954)
    assert (PALL.C_S, PALL.C_Fl, PALL.C_P0, PALL.C_L, PALL.C_V) == (
        2.627,
        2.083,
        0.957,
        1.440,
        0.336,
    )


def test_flooding_of_25mm_pall_rings_air_water():
    """Would have caught the legacy flooding closure: F_fl ≈ 3 Pa^0.5 at low liquid load."""
    u_fl, h_fl = packing.bs_flooding(1.0, packing=PALL, **AIR_WATER)
    assert 2.6 < u_fl * math.sqrt(1.2) < 3.4
    assert h_fl == pytest.approx(PALL.void_fraction / 3, rel=0.05)  # BS eq. 13 lower bound
    u_fl_wet, _ = packing.bs_flooding(10.0, packing=PALL, **AIR_WATER)
    assert u_fl_wet < u_fl  # more liquid, earlier flooding


def test_flooding_point_solves_eqs_13_and_36():
    eps, a = PALL.void_fraction, PALL.a_spec_m2_m3
    LV, w = 3.0, AIR_WATER
    u, h = packing.bs_flooding(LV, packing=PALL, **w)
    X = LV * math.sqrt(w["rho_g"] / w["rho_l"])
    psi = 9.81 / PALL.C_Fl**2 * (X * (w["mu_l"] / w["mu_g"]) ** 0.2) ** (2 * 0.194)
    lhs13 = h**3 * (3 * h - eps)
    rhs13 = 6 / 9.81 * a**2 * eps * w["mu_l"] / w["rho_l"] * LV * w["rho_g"] / w["rho_l"] * u
    assert lhs13 == pytest.approx(rhs13, rel=1e-6)
    u36 = math.sqrt(2) * math.sqrt(9.81 / psi) * (eps - h) ** 1.5 / math.sqrt(eps)
    u36 *= math.sqrt(h / a) * math.sqrt(w["rho_l"] / w["rho_g"])
    assert u == pytest.approx(u36, rel=1e-6)


def test_loading_point_is_about_70_percent_of_flooding():
    for LV in (0.5, 2.0, 10.0):
        u_fl, _ = packing.bs_flooding(LV, packing=PALL, **AIR_WATER)
        u_s = packing.bs_loading_velocity(LV, packing=PALL, **AIR_WATER)
        assert 0.55 < u_s / u_fl < 0.8


def test_pressure_drop_in_typical_design_range():
    w = AIR_WATER
    uL = 0.005
    h = packing.bs_holdup_below_loading(uL, w["rho_l"], w["mu_l"], PALL)
    dpdz = packing.bs_pressure_drop_Pa_m(1.5, uL, w["rho_g"], w["mu_g"], h, h, 2.0, PALL)
    dry = packing.bs_pressure_drop_Pa_m(1.5, 1e-9, w["rho_g"], w["mu_g"], 1e-6, 1e-6, 2.0, PALL)
    assert 50 < dpdz < 600  # Pa/m, irrigated 25 mm rings near 50 % of flooding
    assert dpdz > dry


def test_mass_transfer_coefficients_formulas():
    uL, h_L, D_l = 0.01, 0.08, 1.5e-9
    d_h = 4 * PALL.void_fraction / PALL.a_spec_m2_m3
    kL = PALL.C_L * 12 ** (1 / 6) * math.sqrt(uL / h_L * D_l / d_h)
    assert packing.bs_kL_m_s(uL, h_L, D_l, PALL) == pytest.approx(kL, rel=1e-12)
    rho_g, mu_g, Dg, uG = 1.1, 1.9e-5, 1.7e-5, 1.0
    nu = mu_g / rho_g
    kG = PALL.C_V / math.sqrt(PALL.void_fraction - h_L) * math.sqrt(PALL.a_spec_m2_m3 / d_h) * Dg
    kG *= (uG / (PALL.a_spec_m2_m3 * nu)) ** 0.75 * (nu / Dg) ** (1 / 3)
    assert packing.bs_kG_m_s(uG, h_L, rho_g, mu_g, Dg, PALL) == pytest.approx(kG, rel=1e-12)


def test_column_coefficients_in_typical_ranges():
    r = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(), 0.02)
    assert 1e-3 < r.kG_m_s < 5e-2  # CLAUDE.md sanity range
    assert 1e-4 < r.kL_m_s < 1e-3
    assert r.KGa_1_s < 1.0
    assert r.v_oper_m_s == pytest.approx(0.6 * r.v_flood_m_s, rel=1e-9)
    assert 50 < r.dpdz_Pa_m < 600


def test_onda_alternative_and_structured_packing():
    onda = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(mass_transfer_model="onda"), 0.02)
    bs = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(), 0.02)
    assert onda.kL_m_s != bs.kL_m_s and onda.D_col_m == bs.D_col_m
    ralu = simulate_absorber(
        FlueGas(), Solvent(), AbsorberSpec(packing_name="ralu_pak_metal_yc250"), 0.02
    )
    assert ralu.KGa_1_s > 0
    with pytest.raises(ValueError):
        simulate_absorber(
            FlueGas(),
            Solvent(),
            AbsorberSpec(packing_name="ralu_pak_metal_yc250", mass_transfer_model="onda"),
            0.02,
        )


def test_surface_tension_of_water_IAPWS():
    assert solvent.surface_tension_N_m(298.15, 0.0) == pytest.approx(0.07197, abs=2e-5)


MELLAPAK = packing.PACKINGS["mellapak_metal_250y"]


def test_hanley_chen_sheet_metal_formulas():
    """Hanley & Chen (2012) eqs. 70–72, as printed in the original (not the review's version)."""
    d_e = 4 * 0.970 / 250.0
    rho_g, mu_g, Dg, uG = 1.1, 1.9e-5, 1.7e-5, 1.5
    rho_l, mu_l, Dl, uL, sigma = 1030.0, 8e-4, 1.6e-9, 0.01, 0.075
    Re_V, Re_L = d_e * uG * rho_g / mu_g, d_e * uL * rho_l / mu_l
    kG = 0.0084 * Re_V * (mu_g / (rho_g * Dg)) ** (1 / 3) * Dg / d_e
    kL = 0.33 * Re_L * (mu_l / (rho_l * Dl)) ** (1 / 3) * Dl / d_e
    am = 0.539 * Re_V**0.145 * Re_L**-0.153 * (d_e * rho_l * uL**2 / sigma) ** 0.2
    am *= (uL**2 / (9.81 * d_e)) ** -0.2 * (rho_g / rho_l) ** -0.033 * (mu_g / mu_l) ** 0.090
    assert packing.hc_kG_m_s(uG, rho_g, mu_g, Dg, MELLAPAK) == pytest.approx(kG, rel=1e-12)
    assert packing.hc_kL_m_s(uL, rho_l, mu_l, Dl, MELLAPAK) == pytest.approx(kL, rel=1e-12)
    a_m = packing.hc_area_m2_m3(uG, rho_g, mu_g, uL, rho_l, mu_l, sigma, MELLAPAK)
    assert a_m == pytest.approx(250.0 * am, rel=1e-12)


def test_corrugation_angle_exponents_follow_the_original():
    steep = MELLAPAK.model_copy(update={"corrugation_angle_deg": 30.0})
    factor = math.cos(math.radians(30)) / math.cos(math.radians(45))
    args = (1.5, 1.1, 1.9e-5, 1.7e-5)
    ratio_kG = packing.hc_kG_m_s(*args, steep) / packing.hc_kG_m_s(*args, MELLAPAK)
    assert ratio_kG == pytest.approx(factor**-7.15, rel=1e-12)
    area = (1.5, 1.1, 1.9e-5, 0.01, 1030.0, 8e-4, 0.075)
    ratio_a = packing.hc_area_m2_m3(*area, steep) / packing.hc_area_m2_m3(*area, MELLAPAK)
    assert ratio_a == pytest.approx(factor**4.078, rel=1e-12)


def test_mellapak_uses_its_recommended_model():
    """Billet & Schultes has no C_L/C_V for Mellapak 250Y: 'auto' picks Hanley & Chen."""
    r = simulate_absorber(
        FlueGas(), Solvent(), AbsorberSpec(packing_name="mellapak_metal_250y"), 0.02
    )
    assert r.reached_target and r.KGa_1_s > 0
    with pytest.raises(ValueError):
        simulate_absorber(
            FlueGas(),
            Solvent(),
            AbsorberSpec(packing_name="mellapak_metal_250y", mass_transfer_model="billet_schultes"),
            0.02,
        )
    with pytest.raises(ValueError):
        simulate_absorber(
            FlueGas(), Solvent(), AbsorberSpec(mass_transfer_model="hanley_chen"), 0.02
        )
