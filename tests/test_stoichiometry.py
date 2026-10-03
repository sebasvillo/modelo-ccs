"""Hydroxide stoichiometry along the column (audit §4, errata E-004)."""

import itertools

import pytest

from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent

GRID = list(itertools.product([0.25, 1.5, 4.0], [0.002, 0.005, 0.02, 0.06]))


@pytest.mark.parametrize(("C", "LG"), GRID)
def test_capture_never_exceeds_stoichiometric_capacity(C, LG):
    """Would have caught audit §4: legacy base case captured 418 mol/s with capacity 201."""
    r = simulate_absorber(FlueGas(), Solvent(NaOH_M=C), AbsorberSpec(), LG)
    capacity = C * 1000.0 * r.Ql_m3_s / 2.0
    assert r.OH_capacity_CO2_mol_s == pytest.approx(capacity, rel=1e-12)
    assert r.CO2_captured_mol_s <= capacity * (1 + 1e-9)
    assert min(r.OH_mol_m3) >= 0.0


@pytest.mark.parametrize(("C", "LG"), GRID)
def test_hydroxide_balance_closes(C, LG):
    r = simulate_absorber(FlueGas(), Solvent(NaOH_M=C), AbsorberSpec(), LG)
    consumed = (r.OH_mol_m3[0] - r.OH_mol_m3[-1]) * r.Ql_m3_s
    assert consumed == pytest.approx(2.0 * r.CO2_captured_mol_s, rel=1e-9)


def test_legacy_base_case_is_flagged_as_exhausted():
    r = simulate_absorber(FlueGas(), Solvent(NaOH_M=1.5), AbsorberSpec(), 0.005)
    assert r.solvent_exhausted
    assert not r.reached_target
    assert any("exhausted" in w for w in r.warnings)
    assert any("below the capture target" in w for w in r.warnings)


def test_local_chemistry_follows_hydroxide_depletion():
    r = simulate_absorber(FlueGas(), Solvent(NaOH_M=1.5), AbsorberSpec(), 0.01)
    assert all(a >= b for a, b in itertools.pairwise(r.OH_mol_m3))
    assert all(a >= b for a, b in itertools.pairwise(r.Ha_profile))
    assert r.KGa_profile_1_s[-1] < r.KGa_profile_1_s[0]
    assert r.KGa_1_s == r.KGa_profile_1_s[0]


def test_unreached_target_is_warned():
    r = simulate_absorber(FlueGas(), Solvent(NaOH_M=1.5), AbsorberSpec(max_height_m=2.0), 0.06)
    assert not r.reached_target and not r.solvent_exhausted
    assert any("not reached" in w for w in r.warnings)
