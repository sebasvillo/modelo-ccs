"""Pydantic input/output schemas (single source of truth for the model interface).

Defaults reproduce the legacy base case (thesis notebook, cell 1).
"""

from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field


class _Frozen(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


# ----------------------------------------------------------------- inputs


class Packing(_Frozen):
    a_spec_m2_m3: float = Field(gt=0, description="specific surface area")
    void_fraction: float = Field(gt=0, lt=1)
    dp_eq_m: float = Field(gt=0, description="equivalent particle diameter")
    wetting_ref_uL_m_s: float = Field(gt=0, description="reference liquid velocity for wetting")
    flood_coeff: float = Field(gt=0)


class FlueGas(_Frozen):
    Q_dry_Nm3_h: float = Field(150_000.0, gt=0, description="dry flue-gas flow")
    y_CO2_dry: float = Field(0.25, gt=0, lt=1)
    y_O2_dry: float = Field(0.03, ge=0, lt=1)
    T_C: float = Field(45.0, description="absorber inlet gas temperature")
    P_bar: float = Field(1.01, gt=0, description="absorber pressure")


class Solvent(_Frozen):
    NaOH_M: float = Field(1.5, gt=0, description="circulating NaOH concentration [mol/L]")


def _legacy_LG_grid() -> tuple[float, ...]:
    return tuple(float(v) for v in np.linspace(0.005, 0.060, 28))


class AbsorberSpec(_Frozen):
    capture_target: float = Field(0.90, gt=0, lt=1)
    packing_name: str = "pall_ring_25mm"
    H_eq: float = Field(0.04, ge=0, description="LEGACY: equilibrium slope y* = H_eq·x")
    liquid_resistance_factor: float = Field(
        0.001, gt=0, description="LEGACY(audit §3): tuning factor on the liquid resistance"
    )
    k2_ref_L_mol_s: float = Field(8.5e3, gt=0, description="CO2 + OH- rate constant at 25 °C")
    max_height_m: float = Field(30.0, gt=0)
    dz_m: float = Field(0.02, gt=0)
    flood_fraction: float = Field(0.60, gt=0, lt=1)
    LG_grid_vol: tuple[float, ...] = Field(
        default_factory=_legacy_LG_grid, description="L/G scan [m3 liquid / m3 actual gas]"
    )
    objective: Literal["min_total_power", "min_height"] = "min_total_power"
    pump_eff: float = Field(0.70, gt=0, le=1)
    blower_eff: float = Field(0.68, gt=0, le=1)


class CellSpec(_Frozen):
    j_mA_cm2: float = Field(200.0, gt=0, description="cell current density")
    V_override_V: float | None = Field(None, description="if None, from the j–V correlation")
    tna_mode: Literal["empirical", "const"] = "empirical"
    tna_const: float = Field(0.90, gt=0, le=1)
    CO2_release_eff: float = Field(0.98, gt=0, le=1)


class CaseInput(_Frozen):
    gas: FlueGas = FlueGas()
    solvent: Solvent = Solvent()
    absorber: AbsorberSpec = AbsorberSpec()
    cell: CellSpec = CellSpec()


# ---------------------------------------------------------------- outputs


class AbsorberResult(_Frozen):
    LG_vol: float
    packing_name: str
    T_K: float
    P_Pa: float
    comp_wet_in: dict[str, float]
    yCO2_wet_in: float
    yCO2_out: float
    yCO2_out_target: float
    capture_achieved: float
    reached_target: bool
    Qg_actual_m3_s: float
    Ql_m3_s: float
    n_dry_mol_s: float
    n_wet_mol_s: float
    G_mol_s: float = Field(description="total gas molar flow (not a flux)")
    L_mol_s: float = Field(description="liquid molar flow, water-proxy molar mass")
    rho_g_kg_m3: float
    rho_l_kg_m3: float
    mu_g_Pa_s: float
    mu_l_Pa_s: float
    kG_m_s: float
    kL_m_s: float
    KGa_1_s: float
    a_eff_m2_m3: float
    wetting_fraction: float
    Ha: float
    E: float
    k1_pseudo_1_s: float
    D_g_m2_s: float
    D_l_m2_s: float
    ReG: float
    ReL: float
    ScG: float
    ScL: float
    uG_m_s: float
    uL_m_s: float
    D_col_m: float
    v_flood_m_s: float
    v_oper_m_s: float
    height_m: float
    dpdz_Pa_m: float
    deltaP_Pa: float
    blower_power_W: float
    pump_power_W: float
    hydraulic_power_W: float
    CO2_in_mol_s: float
    CO2_out_mol_s: float
    CO2_captured_mol_s: float
    stack_CO2_ppmv: float
    z_m: tuple[float, ...]
    yCO2: tuple[float, ...]
    x_loading: tuple[float, ...]
    y_star: tuple[float, ...]
    rate_indicator: tuple[float, ...]

    @property
    def Qg_actual_m3_h(self) -> float:
        return self.Qg_actual_m3_s * 3600.0

    @property
    def Ql_m3_h(self) -> float:
        return self.Ql_m3_s * 3600.0


class AbsorberDesign(_Frozen):
    best: AbsorberResult
    feasible: bool
    candidates: tuple[AbsorberResult, ...]


class CellResult(_Frozen):
    j_A_m2: float
    tna: float
    V_cell_V: float
    E_J_mol: float
    A_cell_m2: float
    I_cell_A: float
    P_cell_W: float
    CO2_released_mol_s: float
    CO2_slip_mol_s: float


class CaseResult(_Frozen):
    absorber: AbsorberDesign
    cell: CellResult
    P_total_W: float
    CO2_captured_t_y: float
    CO2_stack_t_y: float
    CO2_product_t_y: float
    E_cell_kWh_t: float
