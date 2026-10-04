"""NaOH solution properties from the literature (errata E-021).

Would have caught the own-closure properties of the thesis (μ·(1 + 0.08·C), ρ linear in C,
D_CO2/(1 + 0.30·C), constant c_p = 3900 J/(kg·K)), which had no source.
"""

import pytest

from ccs_predesign import solvent
from ccs_predesign.constants import MW_kg_mol


def _molarity(T_K: float, w: float) -> float:
    rho = solvent._density_at_mass_fraction(T_K - 273.15, w)
    return w * rho / MW_kg_mol["NaOH"] / 1000.0


# Independent implementation of Laliberté (2009) in the `thermo` package (Bell, MIT licence),
# evaluated once with its NaOH row: (t °C, w) -> (ρ kg/m3, μ Pa·s, c_p J/(kg·K)).
THERMO_REFERENCE = {
    (25.0, 0.06): (1062.8949882029765, 0.0011406402855163172, 3904.2433478598414),
    (60.0, 0.12): (1112.109073377024, 0.0008478647357289572, 3820.7876009460806),
}


@pytest.mark.parametrize(("t_C", "w"), THERMO_REFERENCE)
def test_laliberte_matches_independent_implementation(t_C, w):
    rho_ref, mu_ref, cp_ref = THERMO_REFERENCE[(t_C, w)]
    T = t_C + 273.15
    C = _molarity(T, w)
    assert solvent.mass_fraction_NaOH(T, C) == pytest.approx(w, rel=1e-10)
    assert solvent.density_kg_m3(T, C) == pytest.approx(rho_ref, rel=1e-10)
    assert solvent.viscosity_Pa_s(T, C) == pytest.approx(mu_ref, rel=1e-10)
    # `thermo` uses an IAPWS fit for water c_p; this model uses a constant 4182 J/(kg·K).
    assert solvent.heat_capacity_J_kgK(T, C) == pytest.approx(cp_ref, rel=2e-3)


def test_density_close_to_handbook_values():
    """CRC Handbook, concentrative properties at 20 °C: 4 % NaOH ≈ 1.043, 10 % ≈ 1.109 g/cm3."""
    T = 293.15
    assert solvent.density_kg_m3(T, _molarity(T, 0.04)) == pytest.approx(1042.8, rel=5e-3)
    assert solvent.density_kg_m3(T, _molarity(T, 0.10)) == pytest.approx(1108.9, rel=5e-3)


def test_pure_water_limits():
    T = 298.15
    assert solvent.density_kg_m3(T, 0.0) == pytest.approx(997.05, abs=0.05)
    assert solvent.viscosity_Pa_s(T, 0.0) == pytest.approx(0.890e-3, rel=2e-3)
    assert solvent.heat_capacity_J_kgK(T, 0.0) == solvent.CP_WATER_J_kgK


def test_co2_diffusivity_versteeg_and_viscosity_correction():
    """Versteeg & van Swaaij (1988): D_w = 2.35e-6·exp(−2119/T), and D·μ^0.8 = const."""
    T = 298.15
    assert solvent.diffusivity_CO2_m2_s(T, 0.0) == pytest.approx(1.92e-9, rel=5e-3)
    for C in (1.0, 2.0):
        ratio = solvent.diffusivity_CO2_m2_s(T, C) / solvent.diffusivity_CO2_m2_s(T, 0.0)
        mu_ratio = solvent.viscosity_Pa_s(T, 0.0) / solvent.viscosity_Pa_s(T, C)
        assert ratio == pytest.approx(mu_ratio**0.8, rel=1e-12)
        assert ratio < 1.0


def test_properties_trend_with_concentration_and_temperature():
    assert solvent.viscosity_Pa_s(298.15, 2.0) > solvent.viscosity_Pa_s(298.15, 1.0)
    assert solvent.viscosity_Pa_s(333.15, 1.5) < solvent.viscosity_Pa_s(298.15, 1.5)
    assert solvent.density_kg_m3(298.15, 2.0) > solvent.density_kg_m3(298.15, 1.0)
    assert solvent.heat_capacity_J_kgK(298.15, 2.0) < solvent.heat_capacity_J_kgK(298.15, 1.0)


def test_out_of_range_use_is_warned():
    assert solvent.property_range_warnings(318.15, 1.5) == []
    hot = solvent.property_range_warnings(353.15, 1.5)
    assert len(hot) == 1 and "viscosity" in hot[0]
