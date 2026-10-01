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


def _write_export(tmp, name, values, generator):
    from smatrix_bootstrap.gtb.kernels import ExportWriter
    w = ExportWriter(tmp / generator, {"M": 2, "L": 2, "nu0": -20, "digits": 40, "generator": generator})
    w.write(name, [len(values)], [values])
    w.finish()
    return tmp / generator


def test_compare_flags_a_40_digit_difference(tmp_path):
    from smatrix_bootstrap.gtb.compare import compare_exports
    a = _write_export(tmp_path, "K", ["1234567890123456789012345678901234567890e-40", "0", "-5e-60"], "a")
    b = _write_export(tmp_path, "K", ["1234567890123456789012345678901234567891e-40", "0", "-5e-60"], "b")
    c = _write_export(tmp_path, "K", ["1234567890123456789012345678901234567890e-40", "1e-50", "-5e-60"], "c")
    same = compare_exports(a, a)["K"]
    assert same["max_rel"] == 0 and same["zero_mismatch"] == 0
    near = compare_exports(a, b)["K"]
    assert 5e-40 < near["max_rel"] < 2e-39
    zero = compare_exports(a, c)["K"]
    assert zero["zero_mismatch"] == 1 and zero["max_abs_small"] == 1e-50
