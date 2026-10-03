# Session log

Three lines per working session: what was done, what was found, what comes next.

## 2026-10-02 · Session 1 · Phase 0 "freeze and measure"
- Froze the thesis notebook in `legacy/` (hash test), set up uv + Python 3.13 (`uv.lock`), repo CLAUDE.md, audit, errata template, MIT license, and CI (ruff + pytest incl. slow).
- Golden master replays the notebook (~5 s warm, ~22 s first run on a Mac): all cells match byte-for-byte except cell 2, which read stale kernel globals (LCOC 123.20 saved vs 123.17 clean; strict xfail, audit §7).
- Next (session 2): start the refactor of cell 1 into `src/ccs_predesign` as pure functions, golden master kept green.

## 2026-10-02 · Session 2 · Phase 1 refactor, part 1 (cell 1)
- Ported legacy cell 1 into `src/ccs_predesign` (pydantic inputs/outputs, pure functions, no globals/print); legacy bugs kept and tagged `LEGACY(audit §n)`.
- Equivalence test runs the legacy cell-1 functions from the frozen notebook and compares 11 input variations with `==` (bit-for-bit, incl. infeasible path); a 1e-10 perturbation breaks it.
- Next (session 3): port TEA/optimiser/bench (cells 2, 4/8, 10, 12, 17) the same way, unifying the duplicated cells.

## 2026-10-02 · Session 3 · Phase 1 refactor, part 2 (TEA, optimiser, bench)
- Ported cells 2, 8, 10, 12 and the fixed-diameter variant of cell 4 into `tea.py`, `optimize.py` and `absorber.py`; equivalence test replays the notebook state and matches every grid point, score and bench run bit-for-bit.
- New audit facts: cell 2's TEA uses the default H = 25 m (alias list misses `height_m`); cell 4's optimiser is dead code. The legacy helper must force the Agg backend (macOS backend blocks on plt.show).
- Next (session 4): first physics PR, audit §1–2 (dimensionless NTU exponent and real column area), with errata entries.

## 2026-10-03 · Session 4 · First physics fixes (audit §1–§2)
- Safety nets for physics work: component-level legacy equivalence (80 tests) and a headline-output snapshot tied to the latest errata id (#4).
- E-001 (#5): dimensionless NTU step K_G·a·c_tot·dz/G''; E-002 (#6): Re/wetting with real superficial velocities. Base case H 14.86 → 0.04 → 0.14 m, as the audit predicted (compensating errors removed one at a time).
- Next (session 5): §3, replace the 0.001 liquid-resistance factor by the CO2 Henry constant with Weisenberger–Schumpe salting-out.

## 2026-10-03 · Session 5 · Liquid side, stoichiometry, feasible optimisation (audit §3, §4, §7)
- E-003 (#7): Henry constant (Sander 2015) with Weisenberger–Schumpe salting-out replaces the 0.001 factor; column now liquid-film controlled, KGa ≈ 0.3 1/s, rate peaks near 1.5 M NaOH.
- E-004 (#8): local OH⁻ balance, local Ha/E/KGa and a hard stoichiometric limit with warnings; E-005 (#9): optimiser ranks only designs that reach the target.
- Result: pall ring 25 mm cannot reach 90 % in 30 m; structured 250Y optimum 1.79 M, H 25.1 m, LCOC 126.7 USD/t. Next: §5 counter-current BVP and a physical y* (session 6).

## 2026-10-03 · Session 6 · Equilibrium y* and counter-current column (audit §4, §5, §11)
- E-006 (#10): y* from carbonate speciation (Plummer & Busenberg 1982, Harned & Owen) replaces the tunable H_eq slope; y* ≈ 1e-10 with excess OH⁻.
- E-007 (#11): counter-current column; design mode marches up from the rich end to the exact target height, rating mode shoots on the outlet marching down from the lean end with exact pinch handling (scipy brentq).
- Structured 250Y optimum now H 22.5 m, LCOC 126.4 USD/t; pall ring still infeasible at 30 m (88.7 %). Next: §6 cell electrons per CO2 (session 7).

## 2026-10-03 · Session 7 · Cell electrons and TEA consistency (audit §6, §7, §10)
- E-008 (#12): the cell moves 2 Na⁺ per CO2 (carbonate route, taken from the absorber chemistry): 2015 kWh/t, LCOC ≈ 250 USD/t, as the audit estimated.
- E-009 (#13): annual tonnages use the TEA operating hours (8000 h/y); cell-2 TEA removed; 22.414 vs 44.615 was not an error. E-010 (#14): pump head = H + 5 m.
- Next (session 8): §7 energy balance (heat of absorption) and the E_inf limit; then §8 packing correlations.

## 2026-10-03 · Session 8 · Reaction limits, kinetics and energy balance (audit §7, §9)
- E-011 (#15): E ≤ E∞ via DeCoursey (1974); E-012 (#16): Pohorecki & Moniuk (1988) k2 with local ionic strength (×5 at 45 °C, 1.5 M) makes the pall-ring base case feasible again (H 17.3 m).
- E-013 (#17): adiabatic liquid energy balance (ΔH = −109.4 kJ/mol); rich solvent leaves at ~64 °C; liquid-side properties follow T_L.
- Now: pall optimum 3.0 M, H 15.0 m, LCOC 249.5 USD/t (2 F/CO2 dominates). Next: §8 packing correlations (Onda, Billet–Schultes), then current-density optimisation.
