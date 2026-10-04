"""Parallel trains when one column would exceed the maximum diameter (audit §13)."""

import math

import pytest

from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.models import AbsorberSpec, CaseInput, FlueGas, Solvent
from ccs_predesign.pipeline import run_case


def test_no_limit_means_a_single_column():
    r = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(max_diameter_m=None), 0.02)
    assert r.n_trains == 1


def test_default_limit_is_8_m_after_hossain_2026():
    """Would have caught single columns of any diameter (audit §13, errata E-020)."""
    assert AbsorberSpec().max_diameter_m == 8.0
    r = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(), 0.02)
    assert r.D_col_m <= 8.0 and r.n_trains >= 2


@pytest.mark.parametrize("D_max", [8.0, 5.0, 3.0])
def test_trains_keep_velocities_and_height(D_max):
    one = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(max_diameter_m=None), 0.02)
    split = simulate_absorber(FlueGas(), Solvent(), AbsorberSpec(max_diameter_m=D_max), 0.02)
    assert split.D_col_m <= D_max
    assert split.n_trains == math.ceil((one.D_col_m / D_max) ** 2)
    assert split.uG_m_s == pytest.approx(one.uG_m_s, rel=1e-12)
    assert split.height_m == pytest.approx(one.height_m, rel=1e-6)
    assert split.n_trains * split.D_col_m**2 == pytest.approx(one.D_col_m**2, rel=1e-12)


def test_trains_cost_more_column_steel():
    one = run_case(CaseInput(absorber=AbsorberSpec(max_diameter_m=None)))
    split = run_case(CaseInput(absorber=AbsorberSpec(max_diameter_m=5.0)))
    assert split.absorber.best.n_trains > 1
    assert split.costs.capex_column_usd > one.costs.capex_column_usd
