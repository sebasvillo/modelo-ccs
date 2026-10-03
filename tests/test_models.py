import pytest
from pydantic import ValidationError

from ccs_predesign.models import AbsorberSpec, CaseInput, FlueGas


def test_defaults_are_the_legacy_base_case():
    inp = CaseInput()
    assert inp.gas.Q_dry_Nm3_h == 150_000.0
    assert inp.solvent.NaOH_M == 1.5
    assert inp.absorber.packing_name == "pall_ring_25mm"
    assert len(inp.absorber.LG_grid_vol) == 28


@pytest.mark.parametrize(
    "bad",
    [
        lambda: FlueGas(y_CO2_dry=1.5),
        lambda: FlueGas(Q_dry_Nm3_h=-1.0),
        lambda: AbsorberSpec(capture_target=1.0),
        lambda: AbsorberSpec(unknown_field=1),
    ],
)
def test_invalid_inputs_are_rejected(bad):
    with pytest.raises(ValidationError):
        bad()


def test_inputs_are_immutable():
    with pytest.raises(ValidationError):
        CaseInput().gas.T_C = 20.0
