"""Enhancement factor bounded by the instantaneous-reaction limit (audit §9, errata E-011)."""

import math

import pytest

from ccs_predesign import kinetics, solvent
from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent


@pytest.mark.parametrize("Ha", [0.5, 3.0, 34.0, 300.0])
def test_decoursey_limits(Ha):
    assert kinetics.enhancement_factor_decoursey(Ha, 1e12) == pytest.approx(
        math.sqrt(1 + Ha**2), rel=1e-6
    )
    assert kinetics.enhancement_factor_decoursey(1e6, 50.0) == pytest.approx(50.0, rel=1e-3)
    for E_inf in (2.0, 20.0, 200.0):
        E = kinetics.enhancement_factor_decoursey(Ha, E_inf)
        assert 1.0 <= E <= min(E_inf, math.sqrt(1 + Ha**2)) * (1 + 1e-12)


def test_einf_definition():
    E_inf = kinetics.enhancement_infinite(2e-9, 1000.0, 1e-9, 4.0, nu=2.0)
    assert E_inf == pytest.approx(1 + 2e-9 * 1000.0 / (2 * 1e-9 * 4.0), rel=1e-12)


def test_naoh_diffusivity_scaling():
    assert solvent.diffusivity_NaOH_m2_s(298.15, 0.0) == pytest.approx(2.13e-9, rel=1e-12)
    assert solvent.diffusivity_NaOH_m2_s(318.15, 0.0) > solvent.diffusivity_NaOH_m2_s(298.15, 0.0)
    assert solvent.diffusivity_NaOH_m2_s(298.15, 2.0) < solvent.diffusivity_NaOH_m2_s(298.15, 0.0)


@pytest.mark.parametrize(("C", "LG"), [(0.1, 0.1235), (1.5, 0.005), (1.5, 0.06)])
def test_profile_never_exceeds_instantaneous_limit(C, LG):
    """Would have caught audit §9 / rule 5: legacy E = Ha had no E <= E_inf bound."""
    r = simulate_absorber(FlueGas(), Solvent(NaOH_M=C), AbsorberSpec(), LG)
    D_OH = solvent.diffusivity_NaOH_m2_s(r.T_K, C)
    for y, OH, E in zip(r.yCO2, r.OH_mol_m3, r.E_profile, strict=True):
        if OH <= 0:
            continue
        E_inf = kinetics.enhancement_infinite(D_OH, OH, r.D_l_m2_s, y * r.c_tot_mol_m3 / r.H_cc_CO2)
        assert E <= E_inf * (1 + 1e-12)
