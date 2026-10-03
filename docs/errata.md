# Errata: numerical changes versus the thesis

Every change that moves a number reported by the thesis or by the legacy notebook gets one entry here,
plus an updated test. Golden masters are never rebaselined silently, and numbers are never tuned back
toward the thesis values. Audit references point to [`audit-2026-09.md`](audit-2026-09.md).

After adding a row, regenerate the snapshot with `uv run python tests/golden/snapshot.py E-00N`;
`tests/golden/test_snapshot.py` fails if the snapshot is not tied to the latest row.

Old values refer to the legacy notebook (`legacy/M9_Integracion.ipynb`, clean kernel) unless
marked *thesis*, which refers to the thesis document.

| id | equation / module | old value | new value | cause | test | PR |
|----|-------------------|-----------|-----------|-------|------|----|
| E-000 | *(template)* `absorber.py`, Eq. 2.xx | H = … m | H = … m | audit §… : short reason | `tests/…::test_…` | #… |
| E-001 | `absorber.ntu_increment` (integration step) | base case: H = 14.86 m (*thesis* Table 2.2: 16.36 m), capture 90.01 %, ΔP 1919 Pa, blower 151.8 kW, pump 48.2 kW, cell 66 759 kW, LCOC 125.68 USD/t (cell-8 TEA); optimum 4.50 M, L/G 0.0024, H 13.72 m; bench capture 99.8–100 % | base case: H = 0.04 m (2 steps of dz), capture 96.01 % (overshoot of the last step), ΔP 5.2 Pa, blower 0.41 kW, pump 19.5 kW, cell 71 205 kW, LCOC 124.75 USD/t; optimum 1.35 M, L/G 0.001, H 0.04 m; bench capture 4.5–9.9 % | audit §1: step exponent KGa/G·dz [m/mol] replaced by K_G·a·c_tot·dz/G'' (dimensionless). Heights are still wrong because of §2 (A_ref = 1 m²) and §3 (factor 0.001); next entries | `tests/test_absorber_physics.py` | #5 |
| E-002 | `packing.mass_transfer_coefficients` (velocities), `absorber.simulate_absorber` | after E-001: uG for Re = 53.8 m/s, kG 2.35 m/s, kL 1.43e-4 m/s, KGa 323 1/s, Ha 31.3, H 0.04 m, LCOC 124.75 USD/t; optimum 1.35 M, L/G 0.001; bench (15.4 mm column) capture 4.5–9.9 %, KGa 0.08 1/s | uG = 3.91 m/s, kG 0.377 m/s, kL 2.97e-5 m/s, KGa 72.8 1/s, Ha 150, H 0.14 m, capture 92.3 %, LCOC 124.76 USD/t; optimum 2.06 M, L/G 0.001, H 0.18 m, LCOC 124.67 USD/t; bench capture 99.85–99.98 %, KGa 12.0–12.4 1/s | audit §2: Re and wetting used flows over a 1 m² reference area instead of superficial velocities in the actual column. kG is still ~10× typical (particle Sherwood form, §9) and the liquid side is still scaled by 0.001 (§3) | `tests/test_absorber_physics.py::test_column_scale_up_is_invariant` | #6 |
| E-003 | `solvent.henry_cc_CO2`, `packing.mass_transfer_coefficients` (liquid resistance) | after E-002: liquid resistance × 0.001; KGa 72.8 1/s, Ha 150, H 0.14 m at L/G 0.005, capture 92.3 %, feasible; LCOC 124.76 USD/t; optimum 2.06 M, H 0.18 m, capture 91.1 %; bench capture 99.9 % at any NaOH | H_cc = 3.19 (45 °C, 1.5 M); KGa 0.293 1/s, Ha 34; **infeasible**: best capture 89.4 % at the 30 m limit (L/G 0.06); ΔP 3877 Pa, blower 307 kW, pump 932 kW, LCOC 128.25 USD/t; optimiser picks 3.83 M at 63 % capture (infeasible points not filtered, audit §7); bench capture 11–39 %, now rising with NaOH up to ~0.8 M | audit §3: the undocumented factor 0.001 is replaced by the dimensionless Henry constant (Sander 2015) with Weisenberger–Schumpe (1996) salting-out for fresh NaOH. The column is now liquid-film controlled (> 90 % of the resistance) | `tests/test_solvent.py` | #7 |

<!--
Template for a new entry (copy the row, increment the id):
| E-001 | <equation id or module> | <value + units, thesis/legacy> | <value + units> | <cause, audit §> | <tests/path::test_name> | <PR link> |
-->
