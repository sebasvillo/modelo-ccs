"""Onda et al. (1968) mass-transfer correlations for random packings (audit §8, errata E-014)."""

import math

import pytest

from ccs_predesign import gas, packing, solvent
from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent

PALL = packing.PACKINGS["pall_ring_25mm"]


def test_surface_tension_of_water_IAPWS():
    assert solvent.surface_tension_N_m(298.15, 0.0) == pytest.approx(0.07197, abs=2e-5)
    assert solvent.surface_tension_N_m(298.15, 1.5) > solvent.surface_tension_N_m(298.15, 0.0)


def test_wetted_area_formula_and_bounds():
    rho, mu, sigma, L = 1030.0, 8e-4, 0.075, 10.0
    a = PALL.a_spec_m2_m3
    x = 1.45 * (0.075 / sigma) ** 0.75 * (L / (a * mu)) ** 0.1
    x *= (L**2 * a / (rho**2 * 9.81)) ** -0.05 * (L**2 / (rho * sigma * a)) ** 0.2
    a_w = packing.onda_wetted_area_m2_m3(L, rho, mu, sigma, PALL)
    assert a_w == pytest.approx(a * (1 - math.exp(-x)), rel=1e-12)
    assert 0 < a_w < a
    assert packing.onda_wetted_area_m2_m3(30.0, rho, mu, sigma, PALL) > a_w


def test_kL_and_kG_formulas():
    rho, mu, D, L, a_w = 1030.0, 8e-4, 1.6e-9, 10.0, 150.0
    a, dp = PALL.a_spec_m2_m3, PALL.dp_eq_m
    kL = 0.0051 * (L / (a_w * mu)) ** (2 / 3) * (mu / (rho * D)) ** -0.5 * (a * dp) ** 0.4
    kL *= (mu * 9.81 / rho) ** (1 / 3)
    assert packing.onda_kL_m_s(L, a_w, rho, mu, D, PALL) == pytest.approx(kL, rel=1e-12)
    rho_g, mu_g, Dg, G = 1.1, 1.9e-5, 1.7e-5, 2.0
    kG = 5.23 * a * Dg * (G / (a * mu_g)) ** 0.7 * (mu_g / (rho_g * Dg)) ** (1 / 3) * (a * dp) ** -2
    assert packing.onda_kG_m_s(G, rho_g, mu_g, Dg, PALL) == pytest.approx(kG, rel=1e-12)


def test_random_packing_coefficients_in_typical_ranges():
    """Would have caught audit §8/§9: the particle Sherwood form gave kG = 0.38 m/s (E-002)."""
    r = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(), 0.011)
    assert 1e-4 < r.kL_m_s < 1e-3
    assert r.kG_m_s < 0.15
    assert 0.3 < r.wetting_fraction < 1.0
    assert r.KGa_1_s < 1.0


def test_onda_kG_uses_gas_inlet_conditions():
    T, P = 318.15, 1.01e5
    comp = gas.wet_composition_from_dry(0.25, 0.03, T, P)
    rho_g = gas.density_ideal_kg_m3(T, P, comp)
    kG = packing.gas_side_kG_m_s(T, P, comp, 2.0, PALL)
    expected = packing.onda_kG_m_s(
        rho_g * 2.0,
        rho_g,
        gas.viscosity_sutherland_air_Pa_s(T),
        gas.diffusivity_CO2_in_air_m2_s(T, P),
        PALL,
    )
    assert kG == expected


def test_structured_and_high_gas_load_are_flagged():
    s = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(packing_name="structured_250Y"), 0.01)
    assert any("structured packing" in w for w in s.warnings)
    r = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(), 0.011)
    assert any("F-factor" in w for w in r.warnings)  # legacy flooding gives uG ≈ 3.9 m/s
