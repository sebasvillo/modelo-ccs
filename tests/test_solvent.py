"""CO2 solubility in NaOH solution (audit §3, errata E-003)."""

import math

import pytest

from ccs_predesign import solvent
from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.constants import R_J_molK
from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent


def test_henry_water_25C_matches_literature():
    """Sander (2015): H_cp = 3.3e-4 mol/(m3·Pa) → c_G/c_L ≈ 1.22 (Ostwald coefficient ≈ 0.82)."""
    H_cc = solvent.henry_cc_CO2(298.15, {})
    assert H_cc == pytest.approx(1.0 / (3.3e-4 * R_J_molK * 298.15), rel=1e-12)
    assert 1.15 < H_cc < 1.30


def test_henry_increases_with_temperature():
    assert solvent.henry_cc_CO2(318.15, {}) / solvent.henry_cc_CO2(298.15, {}) == pytest.approx(
        math.exp(2400.0 * (1 / 298.15 - 1 / 318.15)) * 298.15 / 318.15, rel=1e-12
    )


def test_salting_out_1M_NaOH_25C():
    """Weisenberger & Schumpe: (h_Na + h_G)·c + (h_OH + h_G)·c = 0.1638 at c = 1 kmol/m3."""
    assert solvent.salting_out_log10(298.15, solvent.naoh_ions_kmol_m3(1.0)) == pytest.approx(
        (0.1143 - 0.0172) + (0.0839 - 0.0172), rel=1e-12
    )
    assert solvent.henry_cc_CO2(298.15, solvent.naoh_ions_kmol_m3(1.0)) > solvent.henry_cc_CO2(
        298.15, {}
    )


def test_overall_coefficient_uses_henry_constant():
    r = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(), 0.02)
    KG = 1.0 / (1.0 / r.kG_m_s + r.H_cc_CO2 / (r.E * r.kL_m_s))
    assert r.KGa_1_s == pytest.approx(KG * r.a_eff_m2_m3, rel=1e-12)
    liquid_share = (r.H_cc_CO2 / (r.E * r.kL_m_s)) * KG
    assert liquid_share > 0.9  # fast reaction in NaOH is still liquid-film controlled


def test_absorption_rate_has_interior_maximum_in_NaOH():
    """Would have caught audit §3: with the 0.001 factor KGa rose monotonically with NaOH.

    With the real Henry constant, salting-out, viscosity and diffusivity make the rate peak
    at intermediate concentration (near 2.5 M with the Laliberté properties, errata E-021).
    """
    KGa = {
        C: simulate_absorber(FlueGas(), Solvent(NaOH_M=C), AbsorberSpec(), 0.02).KGa_1_s
        for C in (0.5, 2.5, 5.0)
    }
    assert KGa[2.5] > KGa[0.5]
    assert KGa[2.5] > KGa[5.0]
