"""Load legacy notebook cells as plain Python namespaces (read-only, notebook untouched).

Used by equivalence tests: the refactored package must reproduce the legacy functions
bit-for-bit until a fix is recorded in docs/errata.md.
"""

import contextlib
import io
import json
import warnings
from pathlib import Path

import matplotlib
import numpy as np

# Legacy cells call plt.show(); a GUI backend would block the test run.
matplotlib.use("Agg")

LEGACY_NOTEBOOK = Path(__file__).resolve().parents[2] / "legacy" / "M9_Integracion.ipynb"


def cell_source(index: int) -> str:
    cells = json.loads(LEGACY_NOTEBOOK.read_text(encoding="utf-8"))["cells"]
    return "".join(cells[index]["source"])


def load_notebook_state(cells=(1, 2, 4, 6, 8)) -> dict:
    """Execute whole cells in notebook order in one namespace, like the saved kernel run.

    Cell 6 raises KeyError('T_gas_abs_C') in the notebook too; that error is expected here.
    """
    import matplotlib.pyplot as plt

    ns: dict = {"__name__": "legacy_notebook", "display": lambda *a, **k: None}
    for index in cells:
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                exec(compile(cell_source(index), f"legacy_cell{index}", "exec"), ns)
            except KeyError as err:
                if not (index == 6 and err.args == ("T_gas_abs_C",)):
                    raise
        plt.close("all")
    return ns


ABSORBER_FIELDS = {
    "yCO2_wet_in": "yCO2_wet_in",
    "yCO2_out": "yCO2_out",
    "yCO2_out_target": "yCO2_out_target",
    "capture_achieved": "capture_achieved",
    "reached_target": "reached_target",
    "Qg_actual_m3s": "Qg_actual_m3_s",
    "Ql_m3s": "Ql_m3_s",
    "n_dry_mol_s": "n_dry_mol_s",
    "n_wet_mol_s": "n_wet_mol_s",
    "G_molar": "G_mol_s",
    "L_molar": "L_mol_s",
    "rho_g": "rho_g_kg_m3",
    "rho_l": "rho_l_kg_m3",
    "mu_g": "mu_g_Pa_s",
    "mu_l": "mu_l_Pa_s",
    "kG": "kG_m_s",
    "kL": "kL_m_s",
    "KGa": "KGa_1_s",
    "a_eff": "a_eff_m2_m3",
    "wetting_fraction": "wetting_fraction",
    "Ha": "Ha",
    "E": "E",
    "k1_pseudo": "k1_pseudo_1_s",
    "D_g": "D_g_m2_s",
    "D_l": "D_l_m2_s",
    "ReG": "ReG",
    "ReL": "ReL",
    "ScG": "ScG",
    "ScL": "ScL",
    "uG": "uG_m_s",
    "uL": "uL_m_s",
    "D_col_m": "D_col_m",
    "v_flood_m_s": "v_flood_m_s",
    "v_oper_m_s": "v_oper_m_s",
    "height_m": "height_m",
    "dpdz_Pa_m": "dpdz_Pa_m",
    "deltaP_Pa": "deltaP_Pa",
    "blower_power_W": "blower_power_W",
    "pump_power_W": "pump_power_W",
    "total_hydraulic_power_W": "hydraulic_power_W",
    "CO2_in_mol_s": "CO2_in_mol_s",
    "CO2_out_mol_s": "CO2_out_mol_s",
    "CO2_captured_mol_s": "CO2_captured_mol_s",
    "stack_CO2_ppmv": "stack_CO2_ppmv",
    "comp_wet_in": "comp_wet_in",
    "T_K": "T_K",
    "P_Pa": "P_Pa",
}
PROFILE_FIELDS = ["z_m", "yCO2", "x_loading", "y_star", "rate_indicator"]
CELL_FIELDS = {
    "j_A_m2": "j_A_m2",
    "tna": "tna",
    "V_cell": "V_cell_V",
    "E_mol_J": "E_J_mol",
    "A_cell_m2": "A_cell_m2",
    "I_cell_A": "I_cell_A",
    "P_cell_W": "P_cell_W",
    "CO2_released_mol_s": "CO2_released_mol_s",
    "CO2_slip_mol_s": "CO2_slip_mol_s",
}


def assert_absorber_equal(old, new, skip=()) -> None:
    for old_key, new_key in ABSORBER_FIELDS.items():
        if old_key in skip:
            continue
        assert getattr(new, new_key) == old[old_key], old_key
    for key in PROFILE_FIELDS:
        assert np.array_equal(np.asarray(getattr(new, key)), old[key]), key
