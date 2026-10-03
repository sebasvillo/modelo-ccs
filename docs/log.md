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
