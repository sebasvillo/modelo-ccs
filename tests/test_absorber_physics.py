"""Physics checks of the absorber integration (audit §1, errata E-001)."""

import itertools
import math

import pint
import pytest

from ccs_predesign.absorber import ntu_increment, simulate_absorber
from ccs_predesign.constants import R_J_molK
from ccs_predesign.models import AbsorberSpec, FlueGas, Solvent

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

# Bench-like case: capture well below 100 %, so exp(-NTU) is well resolved.
BENCH_GAS = FlueGas(Q_dry_Nm3_h=0.1704, y_CO2_dry=0.1613, y_O2_dry=0.0, T_C=25.0)


def test_ntu_increment_is_dimensionless():
    dntu = ntu_increment(
        Q_(300.0, "1/s"), Q_(38.2, "mol/m**3"), Q_(149.0, "mol/(m**2*s)"), Q_(0.02, "m")
    )
    assert dntu.dimensionless


def test_legacy_exponent_was_not_dimensionless():
    """Audit §1: the legacy step KGa/G·dz used the total molar flow G [mol/s]."""
    legacy = Q_(300.0, "1/s") / Q_(2053.0, "mol/s") * Q_(0.02, "m")
    assert not legacy.dimensionless
    assert legacy.units == ureg("m/mol").units


@pytest.mark.parametrize("Z_m", [0.15, 0.25, 0.35])
def test_capture_matches_analytic_ntu(Z_m):
    """y_out/y_in = exp(−Σ K_G·a·c_tot·dz/G''), using the local K_G·a profile.

    Exact for y* = 0; with excess OH- the equilibrium back-pressure is ~1e-10 (E-006).
    """
    r = simulate_absorber(
        BENCH_GAS,
        Solvent(NaOH_M=0.8),
        AbsorberSpec(capture_target=0.999999, max_height_m=Z_m, dz_m=0.01),
        LG_vol=0.1235,
        D_col_fixed_m=0.150,
    )
    c_tot = BENCH_GAS.P_bar * 1e5 / (R_J_molK * (BENCH_GAS.T_C + 273.15))
    G_flux = r.n_wet_mol_s / (math.pi * 0.150**2 / 4.0)
    dz = [b - a for a, b in itertools.pairwise(r.z_m)]
    # local transfer rate includes the shrinking gas flow (1 − y)²/(1 − y_in) (E-018)
    k = [
        K * c_tot / G_flux * (1 - y) ** 2 / (1 - r.yCO2_wet_in)
        for K, y in zip(r.KGa_profile_1_s, r.yCO2, strict=True)
    ]
    ntu_lower = sum(ki * h for ki, h in zip(k[:-1], dz, strict=True))
    ntu_upper = sum(ki * h for ki, h in zip(k[1:], dz, strict=True))
    lo, hi = sorted((ntu_lower, ntu_upper))
    assert not r.solvent_exhausted
    assert max(r.y_star) < 1e-8
    assert lo * (1 - 1e-12) <= r.NTU <= hi * (1 + 1e-12)
    assert r.yCO2_out / r.yCO2_wet_in == pytest.approx(math.exp(-r.NTU), rel=1e-5)


def test_co2_mass_balance_closes():
    r = simulate_absorber(
        BENCH_GAS, Solvent(NaOH_M=0.8), AbsorberSpec(max_height_m=0.35), 0.1235, D_col_fixed_m=0.15
    )
    gas_side = r.CO2_in_mol_s - r.CO2_out_mol_s  # inert-gas basis since E-018
    Y = [y / (1 - y) for y in (r.yCO2_wet_in, r.yCO2_out)]
    assert gas_side == pytest.approx(r.G_inert_mol_s * (Y[0] - Y[1]), rel=1e-12)
    hydroxide_side = (r.OH_mol_m3[-1] - r.OH_mol_m3[0]) * r.Ql_m3_s / 2.0
    assert hydroxide_side == pytest.approx(gas_side, rel=1e-9)
    assert r.x_loading[-1] == pytest.approx(0.0, abs=1e-12)  # lean solvent at the top
    assert r.L_mol_s * r.x_loading[0] == pytest.approx(gas_side, rel=1e-9)


@pytest.mark.parametrize("scale", [0.01, 4.0])
def test_column_scale_up_is_invariant(scale):
    """Audit §2 (E-002): same superficial velocities and flux → same coefficients and height.

    The flooding-based diameter makes A_col ∝ Q, so scaling the gas flow must not change
    kG, kL, KGa or the packed height. With the legacy 1 m2 reference area it did.
    """
    spec, solv = AbsorberSpec(), Solvent()
    base = simulate_absorber(FlueGas(), solv, spec, 0.005)
    scaled = simulate_absorber(FlueGas(Q_dry_Nm3_h=150_000.0 * scale), solv, spec, 0.005)
    assert scaled.uG_m_s == pytest.approx(base.uG_m_s, rel=1e-12)
    for field in ("kG_m_s", "kL_m_s", "KGa_1_s", "Ha", "NTU"):
        assert getattr(scaled, field) == pytest.approx(getattr(base, field), rel=1e-9), field
    assert scaled.height_m == pytest.approx(base.height_m, abs=1e-9)


@pytest.mark.parametrize("H_max", [2.0, 10.0, 30.0])
def test_pump_lifts_solvent_above_the_bed(H_max):
    """Would have caught audit §10: legacy head 5 + 0.5·H is below H for H > 10 m."""
    spec = AbsorberSpec(max_height_m=H_max)
    r = simulate_absorber(FlueGas(), Solvent(), spec, 0.06)
    assert r.pump_head_m >= r.height_m
    assert r.pump_head_m == pytest.approx(r.height_m + spec.pump_extra_head_m, rel=1e-12)
    expected_W = r.rho_l_kg_m3 * 9.81 * r.pump_head_m * r.Ql_m3_s / spec.pump_eff
    assert r.pump_power_W == pytest.approx(expected_W, rel=1e-12)
