"""The legacy thesis notebook is frozen: any byte change must fail the suite."""

import hashlib
from pathlib import Path

LEGACY_NOTEBOOK = Path(__file__).resolve().parents[1] / "legacy" / "M9_Integracion.ipynb"

# SHA-256 of the notebook as committed to the thesis repo
# (git blob sha1 a0070aa612b1165f311b62653de69f07e1d74bd6).
LEGACY_SHA256 = "2fcec85b45ce580505e6c8019eb4cfa1dbf3f745208a006f18106d348597ff37"


def test_legacy_notebook_is_unchanged():
    digest = hashlib.sha256(LEGACY_NOTEBOOK.read_bytes()).hexdigest()
    assert digest == LEGACY_SHA256, (
        "legacy/M9_Integracion.ipynb changed. It is the frozen source of the golden master "
        "and must never be edited (see CLAUDE.md)."
    )
