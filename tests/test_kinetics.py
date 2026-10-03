"""CO2 + OH- kinetics against Pohorecki & Moniuk (1988) (CLAUDE.md validation target, E-012)."""

import math

import pytest

from ccs_predesign import kinetics
from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.constants import R_J_molK
from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent


def test_infinite_dilution_rate_constant_25C():
    """log10 k2_inf = 11.895 − 2382/T → 8.05e3 L/(mol·s) at 25 °C."""
    k2_L = kinetics.k2_pohorecki_moniuk_m3_mol_s(298.15, 0.0) * 1000.0
    assert k2_L == pytest.approx(10 ** (11.895 - 2382 / 298.15), rel=1e-12)
    assert 7.5e3 < k2_L < 8.5e3


def test_activation_energy():
    """Would have caught the legacy (T/298.15)^1.5 scaling (Ea ≈ 4 kJ/mol instead of ≈ 46)."""
    k_a = kinetics.k2_pohorecki_moniuk_m3_mol_s(291.15, 0.0)
    k_b = kinetics.k2_pohorecki_moniuk_m3_mol_s(313.15, 0.0)
    Ea = R_J_molK * math.log(k_b / k_a) / (1 / 291.15 - 1 / 313.15)
    assert Ea == pytest.approx(2382 * math.log(10) * R_J_molK, rel=1e-9)
    assert 40e3 < Ea < 50e3


def test_ionic_strength_correction_for_NaOH():
    ratio = kinetics.k2_pohorecki_moniuk_m3_mol_s(298.15, 1.5) / (
        kinetics.k2_pohorecki_moniuk_m3_mol_s(298.15, 0.0)
    )
    assert math.log10(ratio) == pytest.approx(0.221 * 1.5 - 0.016 * 1.5**2, rel=1e-12)


def test_extrapolation_outside_validity_is_warned():
    hot = simulate_absorber(FlueGas(T_C=45.0), Solvent(), AbsorberSpec(), 0.03)
    ok = simulate_absorber(FlueGas(T_C=30.0), Solvent(), AbsorberSpec(), 0.03)
    assert any("Pohorecki" in w for w in hot.warnings)
    assert not any("Pohorecki" in w for w in ok.warnings)
