"""The model's headline numbers only move together with an errata entry."""

import json
import math

from snapshot import SNAPSHOT, compute_snapshot, latest_errata_id


def _assert_close(expected, actual, path="outputs"):
    if isinstance(expected, dict):
        assert expected.keys() == actual.keys(), path
        for k in expected:
            _assert_close(expected[k], actual[k], f"{path}.{k}")
    elif isinstance(expected, list):
        assert len(expected) == len(actual), path
        for i, (e, a) in enumerate(zip(expected, actual, strict=True)):
            _assert_close(e, a, f"{path}[{i}]")
    elif isinstance(expected, bool) or isinstance(expected, str):
        assert expected == actual, path
    else:
        assert math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-15), (
            f"{path}: {actual!r} != snapshot {expected!r}. A number moved: add an entry to "
            "docs/errata.md and regenerate with `uv run python tests/golden/snapshot.py E-00N`."
        )


def test_snapshot_is_tied_to_latest_errata():
    data = json.loads(SNAPSHOT.read_text(encoding="utf-8"))
    assert data["errata_id"] == latest_errata_id(), (
        "docs/errata.md has a newer entry than the snapshot: regenerate the snapshot"
    )


def test_model_matches_snapshot():
    data = json.loads(SNAPSHOT.read_text(encoding="utf-8"))
    _assert_close(data["outputs"], compute_snapshot())
