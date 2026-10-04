"""Pydantic input/output schemas (single source of truth for the model interface).

Defaults reproduce the legacy base case (thesis notebook, cell 1).
"""

from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator


class _Frozen(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


# ----------------------------------------------------------------- inputs


class Packing(_Frozen):
    a_spec_m2_m3: float = Field(gt=0, description="specific surface area")
    void_fraction: float = Field(gt=0, lt=1)
    dp_eq_m: float = Field(gt=0, description="equivalent particle (nominal) diameter")
    wetting_ref_uL_m_s: float = Field(gt=0, description="LEGACY wetting closure parameter")
    flood_coeff: float = Field(gt=0)
    kind: Literal["random", "structured"] = "random"
    sigma_c_N_m: float | None = Field(
        None, gt=0, description="critical surface tension of the packing material (Onda)"
    )


class FlueGas(_Frozen):
    Q_dry_Nm3_h: float = Field(150_000.0, gt=0, description="dry flue-gas flow")
    y_CO2_dry: float = Field(0.25, gt=0, lt=1)
    y_O2_dry: float = Field(0.03, ge=0, lt=1)
    T_C: float = Field(45.0, description="absorber inlet gas temperature")
    P_bar: float = Field(1.01, gt=0, description="absorber pressure")


class Solvent(_Frozen):
    NaOH_M: float = Field(1.5, gt=0, description="circulating NaOH concentration [mol/L]")
    T_in_C: float | None = Field(None, description="lean solvent inlet temperature; None = gas T")


def _legacy_LG_grid() -> tuple[float, ...]:
    return tuple(float(v) for v in np.linspace(0.005, 0.060, 28))


class AbsorberSpec(_Frozen):
    capture_target: float = Field(0.90, gt=0, lt=1)
    packing_name: str = "pall_ring_25mm"
    max_height_m: float = Field(30.0, gt=0, le=100.0)
    dz_m: float = Field(0.02, ge=0.005, le=1.0, description="integration step")
    flood_fraction: float = Field(0.60, gt=0, lt=1)
    LG_grid_vol: tuple[float, ...] = Field(
        default_factory=_legacy_LG_grid,
        min_length=1,
        max_length=100,
        description="L/G scan [m3 liquid / m3 actual gas]",
    )
    objective: Literal["min_total_power", "min_height"] = "min_total_power"
    pump_eff: float = Field(0.70, gt=0, le=1)
    pump_extra_head_m: float = Field(
        5.0, ge=0, description="own closure: distributor, freeboard and piping losses above the bed"
    )
    blower_eff: float = Field(0.68, gt=0, le=1)


class CellSpec(_Frozen):
    j_mA_cm2: float = Field(200.0, gt=0, description="cell current density")
    V_override_V: float | None = Field(None, description="if None, from the j–V correlation")
    tna_mode: Literal["empirical", "const"] = "empirical"
    tna_const: float = Field(0.90, gt=0, le=1)
    CO2_release_eff: float = Field(0.98, gt=0, le=1)
    optimize_j: bool = Field(False, description="choose j per design point to minimise the LCOC")
    j_min_mA_cm2: float = Field(50.0, gt=0)
    j_max_mA_cm2: float = Field(250.0, gt=0, description="t_Na correlation caps j at 250")

    @model_validator(mode="after")
    def _j_bounds(self):
        if not self.j_min_mA_cm2 < self.j_max_mA_cm2:
            raise ValueError("j_min_mA_cm2 must be below j_max_mA_cm2")
        return self


class TEASpec(_Frozen):
    """Cost and emission parameters of the legacy TEA (notebook cell 8)."""

    hours_per_year: float = Field(
        8000.0, gt=0, le=8760, description="operating hours per year (capacity factor × 8760)"
    )
    discount_rate: float = Field(0.10, gt=0)
    project_life_y: int = Field(20, gt=0)
    electricity_usd_kWh: float = Field(0.12, ge=0)
    grid_EF_kgCO2e_kWh: float = Field(0.2104, ge=0, description="Colombian grid default")
    fixed_om_fraction: float = Field(0.04, ge=0, description="fixed O&M / installed CAPEX")
    solvent_makeup_usd_t: float = Field(0.75, ge=0)
    water_chem_usd_t: float = Field(0.15, ge=0)
    # Electrochemical cell, Zhang et al. (2024) Supplementary Note 8 (DOE H2A PEM model, 2019).
    cell_stack_usd_m2: float = Field(7722.0, ge=0, description="stack cost per electrode area")
    cell_bop_usd_m2: float = Field(4914.0, ge=0, description="balance of plant per electrode area")
    cell_uninstalled_factor: float = Field(0.12, ge=0, description="installation on top of cost")
    stack_replacement_fraction: float = Field(0.30, ge=0, description="of the stack cost")
    stack_replacement_interval_y: float = Field(7.0, gt=0)
    cell_om_fraction: float = Field(0.025, ge=0, description="labour + maintenance per year")
    h2_loss_fraction: float = Field(0.025, ge=0, description="H2 lost from the HER/HOR loop")
    h2_price_usd_kg: float = Field(5.0, ge=0)
    blower_capex_usd_kW: float = Field(280.0, ge=0)
    pump_capex_usd_kW: float = Field(180.0, ge=0)
    shell_usd_m2: float = Field(2500.0, ge=0)
    packing_usd_m3: float = Field(3500.0, ge=0)
    internals_usd_m2: float = Field(1800.0, ge=0)
    installation_factor: float = Field(1.65, gt=0)


class CaseInput(_Frozen):
    gas: FlueGas = FlueGas()
    solvent: Solvent = Solvent()
    absorber: AbsorberSpec = AbsorberSpec()
    cell: CellSpec = CellSpec()
    tea: TEASpec = TEASpec()


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
    G_flux_mol_m2_s: float = Field(description="gas molar flux G/A_col")
    c_tot_mol_m3: float = Field(description="gas molar concentration P/(R·T)")
    NTU: float = Field(description="gas-phase transfer units over the packed height")
    KGa_mean_1_s: float = Field(description="height-averaged K_G·a_e = NTU·G''/(c_tot·H)")
    L_mol_s: float = Field(description="liquid molar flow, water-proxy molar mass")
    rho_g_kg_m3: float
    rho_l_kg_m3: float
    mu_g_Pa_s: float
    mu_l_Pa_s: float
    kG_m_s: float
    kL_m_s: float
    KGa_1_s: float = Field(description="at the liquid inlet; see KGa_profile_1_s")
    a_eff_m2_m3: float
    wetting_fraction: float
    Ha: float = Field(description="at the liquid inlet (top); see Ha_profile")
    E: float = Field(description="at the liquid inlet (top), DeCoursey with E_inf; see E_profile")
    H_cc_CO2: float = Field(description="dimensionless Henry constant c_G/c_L (salting-out)")
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
    flood_fraction_actual: float = Field(description="uG / v_flood")
    feasible_hydraulically: bool = Field(description="uG <= flood_fraction · v_flood")
    height_m: float
    dpdz_Pa_m: float
    deltaP_Pa: float
    blower_power_W: float
    pump_power_W: float
    pump_head_m: float
    hydraulic_power_W: float
    CO2_in_mol_s: float
    CO2_out_mol_s: float
    CO2_captured_mol_s: float
    stack_CO2_ppmv: float
    z_m: tuple[float, ...]
    yCO2: tuple[float, ...]
    x_loading: tuple[float, ...]
    y_star: tuple[float, ...] = Field(description="CO2 equilibrium over the local liquid")
    rate_indicator: tuple[float, ...]
    OH_mol_m3: tuple[float, ...] = Field(description="local hydroxide concentration")
    T_liquid_K: tuple[float, ...] = Field(description="local liquid temperature (adiabatic)")
    Ha_profile: tuple[float, ...]
    E_profile: tuple[float, ...]
    KGa_profile_1_s: tuple[float, ...]
    OH_capacity_CO2_mol_s: float = Field(description="max CO2 the solvent can take as carbonate")
    OH_per_CO2: float = Field(description="OH- consumed per CO2 absorbed (2: carbonate route)")
    solvent_exhausted: bool
    warnings: tuple[str, ...] = ()

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
    electrons_per_CO2: float = Field(description="Na+ moved per CO2: 2 carbonate, 1 bicarbonate")
    E_J_mol: float
    A_cell_m2: float
    I_cell_A: float
    P_cell_W: float
    CO2_released_mol_s: float
    CO2_slip_mol_s: float


class DesignCosts(_Frozen):
    """Output of tea.design_costs for one design."""

    P_total_W: float
    CO2_captured_t_y: float
    CO2_stack_t_y: float
    CO2_product_t_y: float
    E_cell_kWh_t: float
    E_total_kWh_t: float
    capex_column_usd: float
    capex_cell_usd: float
    capex_blower_usd: float
    capex_pump_usd: float
    capex_total_usd: float
    annualized_capex_usd_y: float
    opex_fixed_usd_y: float
    opex_electricity_usd_y: float
    opex_solvent_usd_y: float
    opex_water_chem_usd_y: float
    opex_stack_replacement_usd_y: float
    opex_cell_om_usd_y: float
    opex_h2_loss_usd_y: float
    opex_total_usd_y: float
    LCOC_usd_t: float
    indirect_kgCO2e_y: float
    indirect_tCO2e_y: float
    indirect_kgCO2e_t: float


class CaseResult(_Frozen):
    absorber: AbsorberDesign
    cell: CellResult
    costs: DesignCosts
    P_total_W: float
    CO2_captured_t_y: float
    CO2_stack_t_y: float
    CO2_product_t_y: float
    E_cell_kWh_t: float


class ExplainedValue(_Frozen):
    """One output with the equation that produced it and the values substituted into it."""

    key: str
    value: float
    unit: str
    equation_id: str
    inputs: dict[str, float] = Field(description="symbol key (equations.SYMBOLS) -> value")


class Profiles(_Frozen):
    z_m: tuple[float, ...]
    yCO2: tuple[float, ...]
    OH_mol_m3: tuple[float, ...]
    T_liquid_K: tuple[float, ...]
    KGa_1_s: tuple[float, ...]


class CaseResponse(_Frozen):
    """What the API returns for one case (CLAUDE.md rule 9)."""

    model_version: str
    errata_id: str
    feasible: bool
    warnings: tuple[str, ...]
    outputs: tuple[ExplainedValue, ...]
    profiles: Profiles


class ColumnCapex(_Frozen):
    A_cs_m2: float
    shell_area_m2: float
    packed_volume_m3: float
    pressure_factor: float
    shell_usd: float
    packing_usd: float
    internals_usd: float
    purchased_usd: float
    installed_usd: float


class DesignPoint(_Frozen):
    """One (NaOH, L/G) point of the design space with its TEA (legacy cell 8)."""

    NaOH_M: float
    LG_vol: float
    absorber: AbsorberResult
    cell: CellResult
    P_total_W: float
    CO2_captured_t_y: float
    CO2_stack_t_y: float
    CO2_product_t_y: float
    E_cell_kWh_t: float
    E_total_kWh_t: float
    capex_column_usd: float
    capex_cell_usd: float
    capex_blower_usd: float
    capex_pump_usd: float
    capex_total_usd: float
    annualized_capex_usd_y: float
    opex_fixed_usd_y: float
    opex_electricity_usd_y: float
    opex_solvent_usd_y: float
    opex_water_chem_usd_y: float
    opex_stack_replacement_usd_y: float
    opex_cell_om_usd_y: float
    opex_h2_loss_usd_y: float
    opex_total_usd_y: float
    LCOC_usd_t: float
    indirect_kgCO2e_y: float
    indirect_tCO2e_y: float
    indirect_kgCO2e_t: float
    warnings: tuple[str, ...] = ()


class OptimizationResult(_Frozen):
    best: DesignPoint
    best_score: float
    points: tuple[DesignPoint, ...]
    scores: tuple[float, ...] = Field(description="NaN for points that miss the capture target")
    feasible: bool
    n_feasible: int
    warnings: tuple[str, ...] = ()
