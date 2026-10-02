"""Certify (V1) and analyse every finished 2403 leaf under ROOT, once.  Usage:
    nice python scripts/gtb/postprocess.py ROOT [--loop SECONDS]
A leaf is finished when report.json holds an optimal SDPB termination.  Watson leaves and the initial leaves they were
seeded from also get `gtb analyze` (with the L = 24 export for omitted waves when M = 50, L = 20).  Results go to
certificate.json / analysis.json in each leaf; a one-line summary per action goes to ROOT/postprocess.log.
"""
import json
import subprocess
import sys
import time
from pathlib import Path

OPTIMAL = "found primal-dual optimal solution"
RUN = [sys.executable, "-m", "smatrix_bootstrap.run", "gtb"]


def _log(root, msg):
    with open(root / "postprocess.log", "a") as f:
        f.write(time.strftime("%Y-%m-%dT%H:%M:%S ") + msg + "\n")


def _seed_refs(root: Path) -> str:
    """Text of every queued line and every running command: a leaf named there may still seed a hot start."""
    txt = "".join(q.read_text() for q in root.glob("queue*.txt"))
    ps = subprocess.run(["ps", "-eo", "args"], capture_output=True, text=True).stdout
    return txt + ps


def _prune(root: Path, leaf: Path, refs: str):
    """Delete sdp.ck of a certified leaf that no queued or running command can still use as --hot-start seed."""
    ck = leaf / "sdp.ck"
    if ck.is_dir() and (leaf / "certificate.json").exists() and f"--hot-start {leaf} " not in refs + " " \
            and not refs.rstrip().endswith(f"--hot-start {leaf}"):
        subprocess.run(["rm", "-rf", str(ck)], check=True)
        _log(root, f"pruned {ck}")


def once(root: Path):
    seeds = set()
    for rep_path in root.rglob("report.json"):
        rep = json.loads(rep_path.read_text())
        if rep.get("watson_from"):
            seeds.add(str(Path(rep["watson_from"]).resolve()))
    for rep_path in sorted(root.rglob("report.json")):
        leaf = rep_path.parent
        rep = json.loads(rep_path.read_text())
        if (rep.get("sdpb_result") or {}).get("terminateReason") != OPTIMAL:
            continue
        k = rep["kernels"]
        if not (leaf / "certificate.json").exists():
            r = subprocess.run(RUN + ["certify", "--leaf", str(leaf), "--kernels", k], capture_output=True, text=True)
            _log(root, f"certify {leaf} rc={r.returncode} {r.stdout.strip()[-200:]}")
        wants = rep.get("watson_from") or str(leaf.resolve()) in seeds
        if wants and not (leaf / "analysis.json").exists():
            extra = []
            spec = rep["spec"]
            big = Path(k + "_L24")
            if spec["M"] == 50 and spec["L"] == 20 and big.exists():
                extra = ["--omitted-kernels", str(big)]
            r = subprocess.run(RUN + ["analyze", "--leaf", str(leaf), "--kernels", k] + extra,
                               capture_output=True, text=True)
            _log(root, f"analyze {leaf} rc={r.returncode} {r.stdout.strip()[-200:]}")
        if (leaf / "analysis.json").exists() or not wants:
            _prune(root, leaf, _seed_refs(root))


if __name__ == "__main__":
    root = Path(sys.argv[1])
    loop = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[2] == "--loop" else 0
    while True:
        once(root)
        if not loop:
            break
        time.sleep(loop)
