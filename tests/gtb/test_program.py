"""Conic-program layer of the 2403 reproduction: block algebra, eliminations, PMP round trip through SDPB."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

SCRATCH = Path("/playpen1/shiqiu/gtb2403_runs/scratch_tests")


def _docker_ok():
    return shutil.which("docker") is not None and subprocess.run(
        ["docker", "image", "inspect", "bootstrapcollaboration/sdpb:3.1.0"], capture_output=True).returncode == 0


@pytest.mark.skipif(not _docker_ok(), reason="SDPB docker image unavailable")
def test_pmp_roundtrip_known_optimum():
    """max y1 s.t. [[1, y1], [y1, 1]] >= 0 and 0 <= y2 <= 2: optimum y1 = 1 exactly."""
    from smatrix_bootstrap.gtb.program import Block
    from smatrix_bootstrap.gtb import pmpwriter
    n = 3                                   # [y0, y1, y2], y0 pinned to 1
    blocks = [Block("disk", 2, {(0, 0): ["1", "0", "0"], (0, 1): ["0", "1", "0"], (1, 1): ["1", "0", "0"]}),
              Block("y2_lo", 1, {(0, 0): ["0", "0", "1"]}),
              Block("y2_hi", 1, {(0, 0): ["2", "0", "-1"]})]
    work = SCRATCH / "pmp_roundtrip"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    info = pmpwriter.write(blocks, ["0", "1", "0"], work / "pmp.json", n_vars=n)
    assert info["n_blocks"] == 3 and info["n_vars"] == 3
    res = pmpwriter.solve(work, precision=128, nproc=1, timeout=600, duality_gap=1e-30, primal_error=1e-30, dual_error=1e-30)
    assert res["terminateReason"] == "found primal-dual optimal solution"
    y = pmpwriter.read_y(work / "out", n_vars=n)
    assert abs(float(y["text"][0]) - 1.0) < 1e-20
    assert 0.0 <= float(y["text"][1]) <= 2.0
