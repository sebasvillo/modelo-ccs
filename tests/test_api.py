"""HTTP API: thin wrapper around the package (CLAUDE.md rules 1 and 9)."""

import pytest
from fastapi.testclient import TestClient

from api.main import app
from ccs_predesign.equations import EQUATIONS, SYMBOLS


@pytest.fixture(scope="module")
def client():
    return TestClient(app)


def test_health(client):
    r = client.get("/health")
    assert r.status_code == 200 and r.json()["status"] == "ok"


def test_equations_endpoint_serves_the_catalog(client):
    body = client.get("/v1/equations").json()
    assert {e["id"] for e in body["equations"]} == set(EQUATIONS)


@pytest.fixture(scope="module")
def default_case(client):
    r = client.post("/v1/case", json={})
    assert r.status_code == 200
    return r.json()


def test_every_output_is_explained(default_case):
    assert default_case["feasible"] is True
    keys = {o["key"] for o in default_case["outputs"]}
    assert {"height_m", "diameter_m", "LCOC_usd_t", "cell_area_m2", "T_rich_K"} <= keys
    for o in default_case["outputs"]:
        eq = EQUATIONS[o["equation_id"]]
        assert set(o["inputs"]) <= set(eq.symbols), (o["key"], set(o["inputs"]) - set(eq.symbols))
        assert all(s in SYMBOLS for s in o["inputs"])


def test_profiles_and_warnings(default_case):
    p = default_case["profiles"]
    assert len(p["z_m"]) == len(p["yCO2"]) == len(p["T_liquid_K"]) > 10
    assert any("Pohorecki" in w for w in default_case["warnings"])


def test_inputs_are_validated(client):
    assert client.post("/v1/case", json={"gas": {"y_CO2_dry": 1.5}}).status_code == 422
    assert client.post("/v1/case", json={"absorber": {"dz_m": 1e-6}}).status_code == 422
    assert client.post("/v1/case", json={"unknown": 1}).status_code == 422


def test_overrides_reach_the_model(client):
    body = client.post("/v1/case", json={"solvent": {"NaOH_M": 3.0}}).json()
    oh = next(o for o in body["outputs"] if o["key"] == "OH_rich_mol_m3")
    assert oh["inputs"]["OH"] == pytest.approx(3000.0)


def test_breakdown_and_sensitivity(default_case):
    lcoc = default_case["summary"]["LCOC_usd_t"]
    assert sum(c["usd_t"] for c in default_case["lcoc_breakdown"]) == pytest.approx(lcoc, rel=1e-9)
    curve = default_case["lcoc_vs_electricity"]
    prices = [p["electricity_usd_kWh"] for p in curve]
    values = [p["LCOC_usd_t"] for p in curve]
    assert prices == sorted(prices) and values == sorted(values)  # costlier power, costlier CO2
    at_default = next(p for p in curve if p["electricity_usd_kWh"] == 0.12)["LCOC_usd_t"]
    assert at_default == pytest.approx(lcoc, rel=1e-9)
    s = default_case["summary"]
    assert s["net_captured_t_y"] == pytest.approx(s["CO2_captured_t_y"] - s["indirect_tCO2e_y"])
