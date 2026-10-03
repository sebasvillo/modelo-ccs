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

<!--
Template for a new entry (copy the row, increment the id):
| E-001 | <equation id or module> | <value + units, thesis/legacy> | <value + units> | <cause, audit §> | <tests/path::test_name> | <PR link> |
-->
