"""CO2 back-pressure over NaOH/Na2CO3 solution (replaces the legacy y* = H_eq·x, errata E-006)."""

import math

import pytest

from ccs_predesign import solvent
from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent


@pytest.mark.parametrize(
    ("T_K", "pK1", "pK2", "pKw"),
    [
        (273.15, 6.579, 10.625, 14.94),
        (298.15, 6.352, 10.329, 13.995),
        (323.15, 6.285, 10.172, 13.26),
    ],
)
def test_carbonate_constants_match_tabulated_values(T_K, pK1, pK2, pKw):
    """Tabulated values: Harned & Davis (1943), Harned & Scholes (1941), Harned & Owen (1958)."""
    K1, K2, Kw = solvent.carbonate_constants(T_K)
    assert -math.log10(K1) == pytest.approx(pK1, abs=0.01)
    assert -math.log10(K2) == pytest.approx(pK2, abs=0.01)
    assert -math.log10(Kw) == pytest.approx(pKw, abs=0.02)


def test_free_CO2_from_speciation_identity():
    """[CO2] from Kw²·[CO3]/(K1·K2·[OH]²) satisfies K1 and K2 with [H+] = Kw/[OH-]."""
    T, OH, CO3 = 318.15, 500.0, 500.0
    K1, K2, Kw = solvent.carbonate_constants(T)
    H = Kw / (OH / 1000.0)
    HCO3 = H * (CO3 / 1000.0) / K2
    CO2 = solvent.free_CO2_equilibrium_mol_m3(T, OH, CO3) / 1000.0
    assert H * HCO3 / CO2 == pytest.approx(K1, rel=1e-12)


def test_back_pressure_negligible_with_excess_hydroxide_and_rising_near_exhaustion():
    T, H_cc, c_tot = 318.15, 3.0, 38.0
    y = [H_cc * solvent.free_CO2_equilibrium_mol_m3(T, OH, (1500 - OH) / 2) / c_tot
         for OH in (1400.0, 750.0, 10.0, 0.1)]  # fmt: skip
    assert y[0] < 1e-9
    assert all(a < b for a, b in zip(y, y[1:], strict=False))
    assert solvent.free_CO2_equilibrium_mol_m3(T, 0.0, 750.0) == math.inf


def test_legacy_equilibrium_slope_is_gone():
    """Would have caught the legacy closure: y* = 0.04·x was an arbitrary tunable slope."""
    assert "H_eq" not in AbsorberSpec.model_fields
    r = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(), 0.03)
    assert max(r.y_star) < 1e-8
