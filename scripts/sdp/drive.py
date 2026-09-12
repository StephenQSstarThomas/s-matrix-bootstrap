#!/usr/bin/env python
"""Fan a stage of the SDP reproduction out over independent solver processes.

Usage:  python scripts/sdp/drive.py <stage> --root <dir> [--jobs N]

Each worker is a separate ``python -m smatrix_bootstrap.sdp solve`` invocation
writing its own ``report.json``; the radial directions of a boundary sweep are
independent, so they are simply sliced across workers (task section 5:
BLAS threads 4-8 per process, at most 16 processes).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

PY = sys.executable
# The pre-registered density caliber: ||rho||_4 <= 100 * 3775 = 377500, the value
# used both by the authors' 2403 code and by this repository's historical
# Newton mainline (basis.py:233), which reproduced the Fig.3 tip to 0.1%.
BASE = ["-m", "smatrix_bootstrap.sdp", "solve", "--cone-scaling", "rownorm",
        "--B", "377500", "--B-norm", "l4", "--generate"]
EPS_GRID = (6e-3, 4e-3, 2e-3, 1e-3, 6e-4, 2e-4)


def chunks(n: int, k: int):
    step = max(1, (n + k - 1) // k)
    return [(i, min(i + step, n)) for i in range(0, n, step)]


def stage_jobs(stage: str, root: str, args) -> list[tuple[str, list[str]]]:
    out = []
    if stage == "p0":
        out.append(("p0", BASE + ["--points", "tip", "--out", f"{root}/p0"]))
    elif stage == "fig3":
        for a, b in chunks(args.ndir, args.jobs):
            out.append((f"fig3_{a:03d}", BASE + ["--ndir", str(args.ndir),
                                                 "--slice", f"{a}:{b}",
                                                 "--out", f"{root}/fig3/{a:03d}"]))
    elif stage == "fig4":
        for eps in EPS_GRID:
            cals = ["chi-a", "chi-b", "chi-c"] if abs(eps - 2e-3) < 1e-12 else ["chi-b"]
            for cal in cals:
                tag = f"eps{eps:.0e}_{cal}"
                for a, b in chunks(args.ndir, max(1, args.jobs // 2)):
                    out.append((f"fig4_{tag}_{a:03d}",
                                BASE + ["--chiral", "--chi", cal, "--eps-chi", str(eps),
                                        "--ndir", str(args.ndir), "--half",
                                        "--slice", f"{a}:{b}",
                                        "--out", f"{root}/fig4/{tag}/{a:03d}"]))
    elif stage == "fig8":
        for a, b in chunks(args.ndir, args.jobs):
            out.append((f"fig8_{a:03d}",
                        BASE + ["--chiral", "--uv", "--ndir", str(args.ndir), "--half",
                                "--slice", f"{a}:{b}",
                                "--out", f"{root}/fig8/{a:03d}"]))
    elif stage == "points":
        for chi in ("chi-a", "chi-b", "chi-c"):
            for sr in ("SR-a", "SR-b", "SR-c"):
                out.append((f"pt_{chi}_{sr}",
                            BASE + ["--chiral", "--uv", "--chi", chi, "--sr", sr,
                                    "--points", "tip", "ref", "mid",
                                    "--x-tip", str(args.x_tip),
                                    "--out", f"{root}/points/{chi}_{sr}"]))
    elif stage == "fig11":
        base = [x for x in BASE if x not in ("--M", "--L", str(args.M), str(args.L))]
        for M, L in ((50, 8), (50, 10), (50, 12), (45, 10), (60, 10)):
            out.append((f"ml_{M}_{L}",
                        base + ["--chiral", "--uv", "--M", str(M), "--L", str(L),
                                "--points", "tip", "ref", "mid",
                                "--x-tip", str(args.x_tip),
                                "--out", f"{root}/fig11/M{M}L{L}"]))
    else:
        raise SystemExit(f"unknown stage {stage}")
    return out


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("stage")
    p.add_argument("--root", required=True)
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--ndir", type=int, default=24)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--M", type=int, default=50)
    p.add_argument("--L", type=int, default=10)
    p.add_argument("--x-tip", type=float, default=0.0826)
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()

    env = dict(os.environ)
    # Clarabel's default linear-algebra backend (faer) threads with rayon, which
    # ignores OMP_NUM_THREADS; without RAYON_NUM_THREADS each worker grabs the
    # whole machine and the workers thrash each other.
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "RAYON_NUM_THREADS"):
        env[k] = str(a.threads)
    global BASE
    BASE = BASE + ["--M", str(a.M), "--L", str(a.L)]
    jobs = stage_jobs(a.stage, a.root, a)
    print(f"{len(jobs)} worker(s), {a.threads} BLAS threads each")
    if a.dry_run:
        for name, cmd in jobs:
            print(name, " ".join(cmd))
        return 0
    os.makedirs(f"{a.root}/logs", exist_ok=True)
    running = []
    for name, cmd in jobs:
        while len(running) >= a.jobs:
            running = [q for q in running if q[1].poll() is None]
            if len(running) >= a.jobs:
                time.sleep(2.0)          # don't busy-spin a core on poll()
        log = open(f"{a.root}/logs/{name}.log", "w")
        running.append((name, subprocess.Popen([PY] + cmd, env=env,
                                               stdout=log, stderr=subprocess.STDOUT)))
        print("started", name, flush=True)
    bad = 0
    for name, q in running:
        if q.wait() != 0:
            print("FAILED", name)
            bad += 1
    print("done, failures:", bad)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
