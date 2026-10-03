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
    """With y* = 0, y_out/y_in = exp(−Σ K_G·a·c_tot·dz/G''), using the local K_G·a profile."""
    r = simulate_absorber(
        BENCH_GAS,
        Solvent(NaOH_M=0.8),
        AbsorberSpec(capture_target=0.999999, max_height_m=Z_m, H_eq=0.0, dz_m=0.01),
        LG_vol=0.1235,
        D_col_fixed_m=0.150,
    )
    c_tot = BENCH_GAS.P_bar * 1e5 / (R_J_molK * (BENCH_GAS.T_C + 273.15))
    G_flux = r.n_wet_mol_s / (math.pi * 0.150**2 / 4.0)
    dz = [b - a for a, b in itertools.pairwise(r.z_m)]
    ntu = sum(K * c_tot * h / G_flux for K, h in zip(r.KGa_profile_1_s[:-1], dz, strict=True))
    assert not r.solvent_exhausted
    assert r.yCO2_out / r.yCO2_wet_in == pytest.approx(math.exp(-ntu), rel=1e-9)
    assert r.NTU == pytest.approx(ntu, rel=1e-9)


def test_co2_mass_balance_closes():
    r = simulate_absorber(
        BENCH_GAS, Solvent(NaOH_M=0.8), AbsorberSpec(max_height_m=0.35), 0.1235, D_col_fixed_m=0.15
    )
    gas_side = r.G_mol_s * (r.yCO2_wet_in - r.yCO2_out)
    liquid_side = r.L_mol_s * r.x_loading[-1]
    assert liquid_side == pytest.approx(gas_side, rel=1e-9)


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
