"""Gas-phase balance on an inert (CO2-free) basis (errata E-018)."""

import pytest

from ccs_predesign.absorber import simulate_absorber
from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent


@pytest.mark.parametrize("target", [0.5, 0.9, 0.97])
def test_capture_is_defined_on_molar_flows(target):
    """Would have caught the dilute-gas shortcut: y_out = y_in(1 − φ) over-captures at 25 % CO2."""
    r = simulate_absorber(FlueGas(), Solvent(NaOH_M=2.5), AbsorberSpec(capture_target=target), 0.03)
    assert r.reached_target
    assert r.CO2_captured_mol_s / r.CO2_in_mol_s == pytest.approx(target, rel=1e-9)
    assert r.yCO2_out > r.yCO2_wet_in * (1 - target)  # mole fraction falls less than the flow


def test_inert_gas_flow_is_conserved():
    r = simulate_absorber(FlueGas(), Solvent(NaOH_M=2.5), AbsorberSpec(), 0.03)
    assert r.G_inert_mol_s == pytest.approx(r.G_mol_s * (1 - r.yCO2_wet_in), rel=1e-12)
    G_out = r.G_inert_mol_s / (1 - r.yCO2_out)
    assert r.G_mol_s - G_out == pytest.approx(r.CO2_captured_mol_s, rel=1e-9)  # only CO2 leaves
