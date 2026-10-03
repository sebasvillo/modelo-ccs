"""Counter-current column: boundary conditions and design/rating consistency (audit §5, E-007)."""

import pytest

from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent


@pytest.mark.parametrize(
    ("max_height_m", "LG"),
    [(60.0, 0.06), (30.0, 0.06), (30.0, 0.005)],
    ids=["design", "rating", "pinched"],
)
def test_boundary_conditions(max_height_m, LG):
    """Would have caught audit §5: legacy put the lean solvent (x = 0) at the gas inlet."""
    r = simulate_absorber(
        FlueGas(), Solvent(NaOH_M=1.5), AbsorberSpec(max_height_m=max_height_m), LG
    )
    assert r.z_m[0] == 0.0 and r.z_m[-1] == pytest.approx(r.height_m, abs=1e-9)
    assert r.yCO2[0] == pytest.approx(r.yCO2_wet_in, rel=1e-9)  # feed gas at the bottom
    assert r.yCO2[-1] == pytest.approx(r.yCO2_out, rel=1e-9)  # treated gas at the top
    assert r.OH_mol_m3[-1] == pytest.approx(1500.0, rel=1e-9)  # fresh solvent at the top
    assert r.x_loading[-1] == pytest.approx(0.0, abs=1e-12)
    assert r.x_loading[0] == max(r.x_loading)  # rich solvent leaves at the bottom


def test_rating_at_design_height_recovers_design_capture():
    flue, solv = FlueGas(), Solvent(NaOH_M=1.5)
    design = simulate_absorber(flue, solv, AbsorberSpec(max_height_m=60.0), 0.06)
    assert design.reached_target
    rating = simulate_absorber(
        flue, solv, AbsorberSpec(max_height_m=design.height_m, capture_target=0.95), 0.06
    )
    assert not rating.reached_target
    assert rating.capture_achieved == pytest.approx(design.capture_achieved, abs=2e-4)


def test_pinched_column_reaches_stoichiometric_capacity():
    """Solvent-limited and tall: the rich end pinches and capture tends to the OH- capacity."""
    r = simulate_absorber(FlueGas(), Solvent(NaOH_M=1.5), AbsorberSpec(), 0.005)
    assert r.CO2_captured_mol_s == pytest.approx(r.OH_capacity_CO2_mol_s, rel=1e-3)
    assert r.CO2_captured_mol_s <= r.OH_capacity_CO2_mol_s
    assert r.y_star[0] == pytest.approx(r.yCO2[0], rel=1e-6)  # equilibrium at the rich end
