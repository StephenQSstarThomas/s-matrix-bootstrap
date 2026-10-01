"""Kernel-layer checks for the 2403 reproduction: provenance, anchors, three-way agreement."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SUBMODULE = REPO / "third_party/gauge-theory-bootstrap"


def _provenance_table():
    text = (REPO / "third_party/PROVENANCE.md").read_text()
    rows = re.findall(r"^\|\s*`([^`]+)`\s*\|\s*([0-9a-f]{64})\s*\|", text, flags=re.M)
    assert rows, "PROVENANCE.md has no file/sha256 table"
    return rows


def test_provenance_hashes():
    if not (SUBMODULE / "README.md").exists():
        pytest.skip("submodule not checked out")
    for rel, digest in _provenance_table():
        assert hashlib.sha256((SUBMODULE / rel).read_bytes()).hexdigest() == digest, rel
