"""HTTP API for the CCS pre-design model. Thin: validate, call ccs_predesign, return JSON.

Run locally:  uv run uvicorn api.main:app --reload   →  http://127.0.0.1:8000/docs
"""

import os
from importlib.metadata import version

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from ccs_predesign.equations import catalog
from ccs_predesign.models import CaseInput, CaseResponse
from ccs_predesign.pipeline import run_case
from ccs_predesign.report import case_response

app = FastAPI(
    title="modelo-ccs",
    summary="Open pre-design of CO2 capture with NaOH and electrochemical regeneration.",
    description="Pre-feasibility estimates (AACE Class 5), never investment-grade. Every output "
    "carries the id of the equation that produced it and the values substituted into it.",
    version=version("ccs-predesign"),
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=os.environ.get("CCS_CORS_ORIGINS", "*").split(","),
    allow_methods=["GET", "POST"],
    allow_headers=["Content-Type"],
)


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "model_version": version("ccs-predesign")}


@app.get("/v1/equations")
def equations() -> dict:
    """The equation catalog (same content as equations.json)."""
    return catalog()


@app.post("/v1/case", response_model=CaseResponse)
def case(inp: CaseInput) -> CaseResponse:
    """Design the absorber over the L/G scan, size the cell and cost the plant.

    Send `{}` for the default case; any field of CaseInput can be overridden.
    """
    return case_response(inp, run_case(inp))
