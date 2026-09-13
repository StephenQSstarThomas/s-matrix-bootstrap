"""Checked SDPB 3.1.0 checkpoint copies for the local, single-host launcher.

Cross-input reuse is an initial iterate, never transferred feasibility. Input
and output checkpoints remain separate, and no source file is modified.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import time

import numpy as np

VERSION = "3.1.0"
IMAGE = "bootstrapcollaboration/sdpb:3.1.0"


def _hash(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _json(path):
    return json.loads(Path(path).read_text())


def _canonical(value):
    return json.dumps(value, sort_keys=True, default=lambda v: v.tolist())


def _structure(root, record, precision):
    """Ordered converted matrix/sampling shapes, including variable counts."""
    sdp, meta = root / "sdp", record["pmp"]
    control = _json(sdp / "control.json")
    command = control.get("command", "")
    if not re.search(r"--precision(?:=|\s+)" + str(precision) + r"(?:\s|$)", command):
        raise ValueError("Converted SDP precision does not match the checkpoint")
    infos = _json(sdp / "pmp_info.json")
    count = meta["n_blocks"]
    if control["num_blocks"] != count or len(infos) != count:
        raise ValueError("Converted SDP block count mismatch")
    shapes = []
    for i, info in enumerate(infos):
        shape = _json(sdp / f"block_info_{i}.json")
        if (info.get("index") != i or shape.get("dim", 0) <= 0 or
                shape.get("num_points", 0) <= 0 or info.get("dim") != shape["dim"] or
                len(info.get("samplePoints", [])) != shape["num_points"]):
            raise ValueError(f"Converted SDP matrix shape mismatch at block {i}")
        shapes.append((shape, {k: v for k, v in info.items() if k != "path"}))
    norm = _json(sdp / "normalization.json")["normalization"]
    objective = _json(sdp / "objectives.json")
    if len(norm) != meta["n_vars"] or len(objective["b"]) != meta["n_vars"] - 1:
        raise ValueError("Converted SDP variable layout mismatch")
    return {"blocks": shapes, "normalization": norm, "n_vars": meta["n_vars"]}


def _copy_checked(source, target, expected):
    """Independent inode plus before/after hashes detects torn or changed input."""
    if source.is_symlink() or not source.is_file() or source.stat().st_size == 0:
        raise ValueError(f"Missing or invalid checkpoint file: {source.name}")
    if _hash(source) != expected:
        raise ValueError(f"Source changed before copying: {source.name}")
    shutil.copy2(source, target)
    if source.samefile(target) or _hash(target) != expected or _hash(source) != expected:
        raise ValueError(f"Source changed while copying: {source.name}")


def validate_resume_source(source_report, settings):
    """Preflight a terminal source; only the new wall-clock budget may change."""
    source_report = Path(source_report).resolve()
    source = _json(source_report)
    if (source.get("solver") != "SDPB" or source.get("status") not in
            ("not_accepted", "solver_failed", "readback_failed") or
            source.get("accepted") is not False or "seconds" not in source):
        raise ValueError("Resume requires a terminal, unaccepted SDPB source leaf")
    process = _json(source_report.parent / "sdpb_process.json")
    if (process.get("terminal") is not True or "returncode" not in process or
            process["returncode"] != source.get("sdpb", {}).get("returncode")):
        raise ValueError("Resume source solver process is not terminal")
    pid = process.get("pid")
    if type(pid) is not int or pid <= 0:
        raise ValueError("Resume source process identity is missing")
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        pass
    except PermissionError as exc:
        raise ValueError("Resume source process may still be live") from exc
    else:
        raise ValueError("Resume source process is still live (or PID was reused)")
    if type(settings.timeout) is not int or settings.timeout <= 0:
        raise ValueError("Resume timeout must be a positive integer")
    for key, value in vars(settings).items():
        if key != "timeout" and source["settings"].get(key) != value:
            raise ValueError(f"Same-input resume cannot change setting {key}")
    return source


def prepare_warm_start(source_report, workdir, record, settings, *, resume=False):
    """Stage accepted cross-input reuse or a terminal same-input continuation.

    This deliberately supports only the repository's SDPB 3.1.0 Docker
    launcher on one host. Identical block timings, rank count and granularity
    preserve its deterministic block distribution. General cluster reuse is
    outside this tested scope. SDPB still validates the binary checkpoint.
    """
    started = time.monotonic()
    source_report, dest = Path(source_report).resolve(), Path(workdir).resolve()
    root, source_hash = source_report.parent, _hash(source_report)
    source = _json(source_report)
    if root == dest or root in dest.parents:
        raise ValueError("Warm-start output must be independent of the source run")
    if resume:
        validate_resume_source(source_report, settings)
        if (source["pmp"]["sha256"] != record["pmp"]["sha256"] or
                source["direction"] != record["direction"] or
                source.get("fix_f00") != record.get("fix_f00")):
            raise ValueError("Same-input resume requires identical PMP, objective and section")
    elif (source.get("solver") != "SDPB" or source.get("status") != "numerically_accepted" or
            source.get("accepted") is not True or
            source.get("verification", {}).get("primal_feasible") is not True or
            source.get("convergence", {}).get("solver_optimal") is not True):
        raise ValueError("Warm start requires an accepted source leaf")
    if (os.environ.get("SDPB_NATIVE", "0") != "0" or os.environ.get("SDPB_HOSTFILE") or
            os.environ.get("SDPB_IMAGE", IMAGE) != IMAGE):
        raise ValueError("Warm start is scoped to the local SDPB 3.1.0 Docker launcher")
    if _canonical(source["spec"]) != _canonical(record["spec"]):
        raise ValueError("Warm start cannot change the model contract")
    if (source.get("fix_f00") is None) != (record.get("fix_f00") is None):
        raise ValueError("Warm start cannot change section presence")
    for key in ("precision", "nproc"):
        if source["settings"][key] != getattr(settings, key):
            raise ValueError(f"Checkpoint {key} mismatch")
    for key in ("n_a", "n_vars", "n_blocks"):
        if source["pmp"][key] != record["pmp"][key]:
            raise ValueError(f"Checkpoint {key} layout mismatch")
    for directory, rec in ((root, source), (dest, record)):
        if _hash(directory / "pmp.json") != rec["pmp"]["sha256"]:
            raise ValueError("PMP changed before warm start")
    if source.get("basis", {}).get("sha256") != record.get("basis", {}).get("sha256"):
        raise ValueError("Checkpoint coordinate basis identity mismatch")
    identities = {source_report: source_hash, root / "pmp.json": source["pmp"]["sha256"]}
    if resume:
        identities[root / "sdpb_process.json"] = _hash(root / "sdpb_process.json")
    if source.get("basis"):
        for directory, rec in ((root, source), (dest, record)):
            path = directory / "basis.npy"
            if _hash(path) != rec["basis"]["sha256"]:
                raise ValueError("Coordinate basis changed before warm start")
            if np.load(path, mmap_mode="r", allow_pickle=False).shape[1:] != (rec["pmp"]["n_a"],):
                raise ValueError("Coordinate basis matrix shape mismatch")
        identities[root / "basis.npy"] = source["basis"]["sha256"]
    checkpoint = root / "sdp.ck"
    meta_path = checkpoint / "checkpoint.json"
    meta_hash, meta = _hash(meta_path), _json(meta_path)
    options = meta.get("options", {})
    if meta.get("version") != VERSION:
        raise ValueError("Checkpoint SDPB version mismatch")
    if any(int(options.get(key, -1)) != settings.precision for key in ("precision", "precision_actual")):
        raise ValueError("Checkpoint precision metadata mismatch")
    if int(options.get("procGranularity", 1)) != 1:
        raise ValueError("Checkpoint rank distribution requires procGranularity=1")
    log = (root / "sdpb.log").read_text()
    mpi = re.search(r"MPI processes:\s*(\d+),\s*nodes:\s*(\d+)", log)
    if (not re.search(r"SDPB version:\s*3\.1\.0(?:\s|$)", log) or not mpi or
            tuple(map(int, mpi.groups())) != (settings.nproc, 1)):
        raise ValueError("Source solver version or single-host rank distribution mismatch")
    structure = _structure(root, source, settings.precision)
    if structure != _structure(dest, record, settings.precision):
        raise ValueError("Ordered SDP block/sampling structure mismatch")
    generations = {meta.get("current"), meta.get("backup")}
    if any(type(g) is not int or g < 0 for g in generations):
        raise ValueError("Invalid checkpoint generation metadata")
    files = [meta_path, checkpoint / "block_timings"]
    for generation in sorted(generations):
        expected = {f"checkpoint_{generation}_{rank}" for rank in range(settings.nproc)}
        actual = {p.name for p in checkpoint.glob(f"checkpoint_{generation}_*")}
        if actual != expected:
            raise ValueError(f"Missing or extra rank files in checkpoint generation {generation}")
        files.extend(checkpoint / name for name in sorted(expected))
    timings = [float(v) for v in (checkpoint / "block_timings").read_text().split()]
    if len(timings) != source["pmp"]["n_blocks"] or any(not math.isfinite(t) or t < 0 for t in timings):
        raise ValueError("Checkpoint block_timings shape or values mismatch")
    hashes = {p.name: _hash(p) for p in files}
    if hashes["checkpoint.json"] != meta_hash:
        raise ValueError("Checkpoint generation changed during validation")
    initial = dest / "warm_start_input"
    initial.mkdir(exist_ok=False)
    for path in files:
        _copy_checked(path, initial / path.name, hashes[path.name])
    target_timings = dest / "sdp/block_timings"
    if target_timings.exists():
        raise ValueError("Fresh target SDP unexpectedly contains block_timings")
    _copy_checked(checkpoint / "block_timings", target_timings, hashes["block_timings"])
    if any(_hash(path) != digest for path, digest in identities.items()) or _hash(meta_path) != meta_hash:
        raise ValueError("Source identity or checkpoint generation changed while staging")
    result = {"mode": "same-input-resume" if resume else "accepted-cross-input-warm-start",
              "source_report": str(source_report), "source_report_sha256": source_hash,
              "source_checkpoint": str(checkpoint), "input_directory": initial.name,
              "output_checkpoint_directory": "sdp.ck", "current_generation": meta["current"],
              "backup_generation": meta["backup"], "version": VERSION,
              "precision": settings.precision, "nproc": settings.nproc, "nodes": 1,
              "image": IMAGE, "proc_granularity": 1, "source_file_sha256": hashes,
              "ordered_structure_sha256": hashlib.sha256(_canonical(structure).encode()).hexdigest(),
              "source_pmp_sha256": source["pmp"]["sha256"],
              "target_pmp_sha256": record["pmp"]["sha256"],
              "basis_sha256": source.get("basis", {}).get("sha256"),
              "copy_method": "copy2 with independent inodes and before/after SHA256; no hard links",
              "feasibility_transferred": False, "seconds": time.monotonic() - started,
              "scope": "local single-host SDPB 3.1.0 Docker; matching ranks, granularity and block_timings; no general multi-node validation",
              "binary_scope": "all current/backup rank files preserved and hashed; binary decoding is performed by SDPB"}
    (dest / "warm_start_manifest.json").write_text(json.dumps(result, indent=2))
    return result
