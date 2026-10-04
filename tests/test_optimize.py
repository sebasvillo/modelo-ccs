"""Design-space search only ranks feasible designs (audit §7, errata E-005)."""

import math

import pytest

from ccs_predesign.models import AbsorberSpec, CaseInput
from ccs_predesign.optimize import choose_best, evaluate_grid

NAOH = [0.5, 1.5, 3.0]
LG = [0.004, 0.02, 0.06]


@pytest.fixture(scope="module")
def structured_points():
    return evaluate_grid(
        CaseInput(absorber=AbsorberSpec(packing_name="ralu_pak_metal_yc250")), NAOH, LG
    )


def test_best_design_reaches_the_target(structured_points):
    """Would have caught audit §7: the legacy score picked designs far below the target."""
    result = choose_best(structured_points)
    assert any(not p.absorber.reached_target for p in structured_points)
    assert result.feasible
    assert result.best.absorber.reached_target
    assert result.n_feasible == sum(p.absorber.reached_target for p in structured_points)


def test_infeasible_points_are_not_scored(structured_points):
    result = choose_best(structured_points)
    for p, s in zip(result.points, result.scores, strict=True):
        assert math.isnan(s) != p.absorber.reached_target


def test_no_feasible_design_is_flagged():
    points = evaluate_grid(CaseInput(absorber=AbsorberSpec(max_height_m=0.5)), NAOH, LG)
    result = choose_best(points)
    assert not result.feasible and result.n_feasible == 0
    assert result.best.absorber.capture_achieved == max(p.absorber.capture_achieved for p in points)
    assert result.warnings and "no design reaches" in result.warnings[0]
