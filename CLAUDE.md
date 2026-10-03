# CLAUDE.md: modelo-ccs (open CCS pre-design tool)

## What this repo is
Open-source pre-design tool for point-source CO2 capture: packed absorption column (NaOH/KOH) + electrochemical solvent regeneration (PSE cell, Zhang et al.) + techno-economics (LCOC). It comes from Sebastián Villota's undergraduate thesis (Universidad de los Andes, 2026; advisors A. González Mancera and R. Sierra Ramírez). Target users: engineers, regulators, students, plants in developing countries. Output is **pre-feasibility (AACE Class 5)**, never investment-grade.

Owner: Sebas (mechanical engineer, not a web dev). Explain web/devops choices in plain Spanish. Code, identifiers, docstrings and repo docs are in English. The future UI is bilingual ES/EN.

## Priority: correct before matching
The goal is a **physically correct** model, not one that reproduces the thesis numbers or the bench experiments.
- The golden master (`tests/golden/`) is only a safety net for the refactor (phases 0–1). It must stay green while code is restructured without changing results.
- After that, numbers change **only** through `docs/errata.md`, and are **never** tuned, fitted or calibrated back toward thesis or experimental values. Report honest deviations.

## Decisions already taken (do not reopen)
- IP is the author's; no institutional clearance needed.
- License: MIT for code (`LICENSE`), CC BY 4.0 for documentation.
- Python 3.13 everywhere (project, CI, local), managed with **uv**. `requires-python = ">=3.13,<3.14"`, versions locked in `uv.lock`.
- Later phases: the API runs in a free container (cold-start wait is acceptable); a separate website repo consumes `equations.json` and calls the API. Never commit website code here, or model code there.

## Layout
```
legacy/M9_Integracion.ipynb   # frozen thesis notebook. NEVER edit (hash-checked by tests). Source of golden-master numbers.
src/ccs_predesign/
  gas.py          # flue-gas conditioning, properties
  solvent.py      # NaOH/KOH properties, speciation/equilibrium
  kinetics.py     # k2, Hatta, enhancement (incl. E_inf)
  packing.py      # packing database + mass-transfer/hydraulic correlations
  absorber.py     # counter-current column solver (BVP)
  cell.py         # electrochemical regeneration
  tea.py          # CAPEX/OPEX/LCOC, cost database with sources
  optimize.py     # constrained optimisation (scipy)
  models.py       # pydantic input/output schemas (single source of truth)
  equations.py    # equation catalog (LaTeX, meaning ES/EN, symbols, provenance, refs) → exported to equations.json
api/              # FastAPI app. Thin: validates, calls ccs_predesign, returns JSON.
tests/            # pytest: unit, golden master (tests/golden/), validation vs literature
docs/audit-2026-09.md  # code audit of the legacy notebook
docs/errata.md    # every numerical change vs the thesis, with the reason
docs/log.md       # 3-line summary per working session
```
Status: cell 1 of the legacy notebook (absorber + L/G scan + cell sizing) lives in `constants.py`, `gas.py`, `solvent.py`, `kinetics.py`, `packing.py`, `absorber.py`, `cell.py`, `models.py` and `pipeline.py` (`run_case(CaseInput)`), ported bit-for-bit in phase 1. Remaining legacy bugs are tagged `LEGACY(audit §n)` in the code; each one is removed only by a physics PR with an errata entry. Fixed so far: §1 (E-001), §2 (E-002), §3 (E-003), §4 stoichiometry and local chemistry (E-004), physical y* from carbonate equilibrium (E-006), counter-current design/rating with exact last slice (E-007, §5 and §11), cell electrons per CO2 from the absorber chemistry (E-008, §6), hours and pump head (E-009, E-010), E ≤ E_inf via DeCoursey (E-011), Pohorecki & Moniuk kinetics with local ionic strength (E-012), adiabatic liquid energy balance (E-013, §7), Onda (1968) mass transfer for random packings (E-014, §8 partial: hydraulics and structured packings still legacy), optimiser feasibility filter (E-005, from bug 9, moved ahead because it made optimisation results meaningless). TEA (`tea.py`: cell 8 `evaluate_design_point`, cell 2 `tea_summary`), optimiser (`optimize.py`: cells 10/12 two-stage grid) and the fixed-diameter bench variant (cell 17) are ported too, ported bit-for-bit in phase 1. Cell 4's optimiser is dead code (only driven by failing cell 6) and was not ported. The equation catalog lives in `equations.py` (48 equations, bilingual, with provenance and validity) and is exported to `equations.json` at the repo root (`uv run python -m ccs_predesign.equations equations.json`); `tests/test_equations.py` fails if a physics function does not cite its equation id or the export is stale.

## Hard rules
1. **Physics lives only in `src/ccs_predesign`.** API and frontend never compute physics.
2. **SI units internally**, explicit suffixes in names (`_m`, `_Pa`, `_mol_s`, `_mol_m2s`). Tests check dimensions with `pint`. Never mix total molar flow (mol/s) with molar flux (mol/m²/s).
3. **No hidden globals, no notebook state.** Everything is a pure function of a pydantic input model. No `globals()` lookups, no `print` in library code (use `logging`).
4. **Every equation carries provenance**: `theory | literature | own_closure`, a reference key and a validity range. Outputs return `warnings` whenever a correlation is used outside its range. This traceability is a core feature (the thesis's Table 2.1), not decoration.
5. **Hard physical constraints are checked, not assumed**: stoichiometric OH⁻ capacity, flooding fraction, max height, E ≤ E_inf, electrons per CO2 (2 for carbonate, 1 for bicarbonate).
6. **Any change that moves a thesis number** gets an entry in `docs/errata.md` (old value, new value, cause, test, PR) and an updated test. Never silently rebaseline golden masters.
7. **Equations are a product feature.** `equations.py` is the single catalog: id, LaTeX, plain-language meaning (ES/EN), symbols with units, provenance, reference and validity range. Every physics function cites its equation id. API responses include, per output, the equation id plus the intermediate values used. Equations are never hand-written downstream.
8. Small PRs: one physics fix per PR, each with a test that would have caught the bug.
9. Before finishing a task: `uv run ruff check` and `uv run pytest -q` (slow tests included) must be green.
10. Never edit `legacy/M9_Integracion.ipynb`. If the golden master breaks because of library versions, fix the environment (pin versions), never the notebook.

## Known bugs inherited from the thesis notebook (fix in this order)
Details in `docs/audit-2026-09.md`.
1. Integration exponent uses KGa/G_molar·dz (not dimensionless) → use KG·a·c_tot/(G/A).
2. `A_ref = 1` in mass-transfer coefficients → use the real column area.
3. `liquid_resistance_factor = 0.001` → replace with the Henry constant (salting-out corrected).
4. No stoichiometric limit, Ha/E not updated as OH⁻ depletes → local chemistry each step.
5. Co-current integration presented as counter-current → BVP solver.
6. Cell assumes 1 e⁻ per CO2 → configurable route (carbonate 2 e⁻ / bicarbonate 1 e⁻).
7. Isothermal → add energy balance (heat of absorption).
8. Correlations mislabelled (Sh "Onda", wetting "Onda", Ergun for structured) → proper Onda / Billet–Schultes / Rocha–Bravo–Fair; TPMS packings from Ellebracht et al. 2023 data.
9. TEA: 8760 vs 8000 h, 22.414 vs 44.615, optimiser ignores infeasible points, current density not optimised.

## Safety nets (tests/golden/)
1. `test_legacy_notebook.py` (slow): replays the frozen notebook with nbclient and compares its stdout with the saved outputs. Always green; it only tests the legacy file. All cells reproduce byte-for-byte except cell 2 (stale kernel state in the saved run: LCOC 123.20 saved vs 123.17 clean; strict xfail, audit §7). Cell 6 raises `KeyError: 'T_gas_abs_C'` in the saved run too.
2. `test_components_equivalence.py`: every ported function equals its legacy twin bit-for-bit (properties, correlations, cell, TEA costing, scores). When a physics PR changes a function, replace its test here by one tied to the errata entry.
3. `test_snapshot.py` + `snapshots/model_outputs.json`: headline outputs (base case, optimisation, bench) at rel. tol. 1e-9, tied to the latest id in `docs/errata.md`. A physics PR adds the errata row, then regenerates with `uv run python tests/golden/snapshot.py E-00N`. Never regenerate without an errata row.
The end-to-end refactor proofs (`test_cell1_equivalence.py`, `test_tea_equivalence.py`, PRs #2–#3) were retired by E-001; they live in git history.

## Validation targets
- Kinetics vs Pohorecki & Moniuk (1988). kLa vs Ellebracht et al. (2023) NaOH/TPMS data. Bench data: `data/lab2_titration.csv` (12 runs, thesis Table 3.1; to be added).
- Sanity ranges: kG ~1e-3–5e-2 m/s, KGa (gas) typically < 1 1/s, cell 1–4 V.

## Commands
```
uv sync                        # install (Python 3.13 + locked deps)
uv run pytest -q               # all tests, including the slow golden master
uv run pytest -q -m "not slow" # fast loop
uv run ruff check              # lint
uv run ruff format             # format
uvicorn api.main:app --reload  # (later) local API at http://127.0.0.1:8000/docs
```

## Working style for Claude
- Use plan mode before any physics change; state the equation, units and reference before coding.
- Prefer boring, well-known libraries (numpy, scipy, pydantic, FastAPI). No new dependency without saying why.
- Sessions are short: finish one closed, tested PR per session and write a 3-line summary in `docs/log.md`.
- Commits small and in Spanish.
