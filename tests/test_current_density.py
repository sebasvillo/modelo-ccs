"""Current density as a design variable (bug 9: it was fixed at 200 mA/cm2)."""

import pytest
from pydantic import ValidationError

from ccs_predesign.models import CaseInput, CellSpec, Solvent, TEASpec
from ccs_predesign.tea import evaluate_design_point

CASE = CaseInput()
DEFAULT_TEA = TEASpec()


def point(cell: CellSpec, tea: TEASpec = DEFAULT_TEA):
    return evaluate_design_point(CASE.gas, Solvent(NaOH_M=2.5), CASE.absorber, cell, tea, 0.0093)


def test_default_keeps_the_thesis_current_density():
    p = point(CellSpec())
    assert p.cell.j_A_m2 == pytest.approx(2000.0)
    assert p.warnings == ()


def test_interior_optimum_is_a_minimum():
    """Would have caught bug 9: LCOC(j) has an interior optimum the legacy never searched."""
    tea = DEFAULT_TEA  # Zhang et al. (2024) cell costs: interior optimum near 285 mA/cm2
    best = point(CellSpec(optimize_j=True, j_max_mA_cm2=1000.0), tea)
    j_star = best.cell.j_A_m2 / 10.0
    assert 50.0 < j_star < 1000.0 and best.warnings == ()
    for dj in (-10.0, 10.0):
        other = point(CellSpec(j_mA_cm2=j_star + dj), tea)
        assert other.LCOC_usd_t > best.LCOC_usd_t
    assert best.LCOC_usd_t < point(CellSpec(), tea).LCOC_usd_t


def test_optimum_on_a_bound_is_warned():
    cheap = TEASpec(cell_stack_usd_m2=420.0, cell_bop_usd_m2=0.0)  # legacy cell cost
    p = point(CellSpec(optimize_j=True), cheap)  # cheap membrane: lowest j wins
    assert p.cell.j_A_m2 / 10.0 == pytest.approx(50.0, abs=0.5)
    assert any("search bound" in w for w in p.warnings)


def test_bounds_are_validated():
    with pytest.raises(ValidationError):
        CellSpec(j_min_mA_cm2=200.0, j_max_mA_cm2=100.0)
