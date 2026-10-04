"""Cell CAPEX/OPEX after Zhang et al. (2024), Supplementary Note 8 (errata E-015)."""

import pytest

from ccs_predesign.constants import FARADAY_C_mol
from ccs_predesign.models import CaseInput, TEASpec
from ccs_predesign.pipeline import run_case


def test_zhang_reference_case_per_m2():
    """SI Note 8: 14.58 m2 at 1.12 × 7722 $/m2 = $126,100 (stack), × 4914 $/m2 = $80,244 (BOP)."""
    tea = TEASpec()
    assert 14.58 * tea.cell_stack_usd_m2 * (1 + tea.cell_uninstalled_factor) == pytest.approx(
        126_100, rel=1e-3
    )
    assert 14.58 * tea.cell_bop_usd_m2 * (1 + tea.cell_uninstalled_factor) == pytest.approx(
        80_244, rel=1e-3
    )


def test_cell_cost_terms():
    """Would have caught the legacy 420 $/m2 cell cost (≈ 34× below the DOE H2A-based estimate)."""
    r = run_case(CaseInput())
    k, c, tea = r.costs, r.cell, TEASpec()
    stack = tea.cell_stack_usd_m2 * c.A_cell_m2 * 1.12
    bop = tea.cell_bop_usd_m2 * c.A_cell_m2 * 1.12
    assert k.capex_cell_usd == pytest.approx(stack + bop, rel=1e-12)
    assert k.opex_stack_replacement_usd_y == pytest.approx(0.30 * stack / 7.0, rel=1e-12)
    assert k.opex_cell_om_usd_y == pytest.approx(0.025 * stack, rel=1e-12)
    h2_kg_y = c.I_cell_A / (2 * FARADAY_C_mol) * 2.016e-3 * 3600 * tea.hours_per_year
    assert k.opex_h2_loss_usd_y == pytest.approx(0.025 * 5.0 * h2_kg_y, rel=1e-12)
    assert k.opex_fixed_usd_y == pytest.approx(
        0.04 * (k.capex_total_usd - k.capex_cell_usd), rel=1e-12
    )
    assert k.capex_cell_usd / c.A_cell_m2 > 10_000  # installed cell cost per m2
