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
