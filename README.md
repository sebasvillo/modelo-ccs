# modelo-ccs

Open-source pre-design tool for point-source CO2 capture:

- packed absorption column with NaOH/KOH,
- electrochemical solvent regeneration (PSE cell, after Zhang et al.),
- techno-economics (levelized cost of capture, LCOC).

It grows out of the undergraduate thesis *Herramientas de Diseño para Captura de Carbono por
Absorción en Sectores Difíciles de Abatir* (Sebastián Villota, Universidad de los Andes, 2026;
advisors A. González Mancera and R. Sierra Ramírez).

> **Status: phase 0 (frozen legacy + golden master).** Results are pre-feasibility estimates
> (AACE Class 5), never investment-grade. The legacy notebook has known issues listed in
> [`docs/audit-2026-09.md`](docs/audit-2026-09.md); its numbers are not reliable and are kept
> only as a refactoring safety net.

## Layout

```
legacy/M9_Integracion.ipynb   frozen thesis notebook (never edited)
src/ccs_predesign/            the model package (being built)
tests/                        pytest suite; tests/golden/ replays the legacy notebook
docs/                         audit, errata, session log
```

## Quick start

Requires [uv](https://docs.astral.sh/uv/) (Python 3.13 is installed by uv).

```bash
uv sync
uv run pytest -q          # includes the slow golden-master test
uv run pytest -q -m "not slow"
uv run ruff check
```

## API

```bash
uv run uvicorn api.main:app --reload   # http://127.0.0.1:8000/docs
```

`POST /v1/case` takes a `CaseInput` (send `{}` for the default case) and returns every headline
number with the id of its equation and the values substituted into it, plus warnings and column
profiles. `GET /v1/equations` serves the equation catalog.

## Equations

Every equation the model uses, with its meaning (ES/EN), symbols and units, provenance
(`theory`, `literature` or `own_closure`), reference and validity range, is listed in
[`src/ccs_predesign/equations.py`](src/ccs_predesign/equations.py) and exported to
[`equations.json`](equations.json).

## License

Code: [MIT](LICENSE). Documentation (`docs/`, README): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
