"""SDPB-only execution, immutable round evidence and numerical acceptance.

SDPB's y is the maximization candidate (dualObjective). primalObjective is
the candidate upper support bound. Neither is a rigorous interval certificate.
No failed or partial iterate is used as a relaxation optimum.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import shutil
import time

import numpy as np

from .pmp import Pmp
ROOT = Path(__file__).resolve().parents[3]

@dataclass
class Settings:
    precision: int = 192
    nproc: int = 8
    duality_gap: float = 1e-6
    primal_error: float = 1e-10
    dual_error: float = 1e-10
    initial_scale: float = 1e4
    max_iterations: int = 500
    timeout: int = 7200
    checkpoint_interval: int = 300
    digits: int = 17

def read_out(path):
    out = {}
    if not Path(path).exists():
        return out
    for line in Path(path).read_text().splitlines():
        match = re.match(r'\s*(\w+)\s*=\s*"?([^";]+)"?;', line)
        if match:
            key, value = match.group(1), match.group(2).strip()
            try:
                out[key] = float(value)
            except ValueError:
                out[key] = value
    return out

def convergence_report(data, returncode, objective=None, settings=None):
    cfg = settings or Settings()
    p, d = data.get("primalObjective"), data.get("dualObjective")
    finite = all(isinstance(x, (float, int)) and np.isfinite(x) for x in
                 (p, d, data.get("primalError"), data.get("dualError"),
                  data.get("dualityGap"), objective))
    gap = abs(p-d) / max(1., abs(p)+abs(d)) if finite else None
    match = abs(objective-d) if finite else None
    ok = (returncode == 0 and finite and
          data.get("terminateReason") == "found primal-dual optimal solution" and
          data["primalError"] <= cfg.primal_error and
          data["dualError"] <= cfg.dual_error and
          data["dualityGap"] <= cfg.duality_gap and gap <= cfg.duality_gap and
          match <= 1e-8 * max(1., abs(d)) and p >= d-1e-10*max(1., abs(d)))
    return {"solver_optimal": bool(ok), "independent_dual_certified": False,
            "candidate_upper": p, "candidate_lower": d,
            "absolute_gap": abs(p-d) if finite else None,
            "recomputed_gap": gap, "objective_readback_error": match,
            "bound_scope": "solver numerical residuals; not an interval bound",
            "terminate_reason": data.get("terminateReason")}

def write_json(path, data):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, default=lambda x:
                             x.tolist() if isinstance(x, np.ndarray) else x.item()))
    tmp.replace(path)

def stop_owned_process(proc, cidfile, log):
    """Release only the invocation we created, including an already-exited group."""
    import signal
    if cidfile.exists() and cidfile.read_text().strip():
        try:
            subprocess.run(["docker", "stop", "-t", "10", cidfile.read_text().strip()],
                           stdout=log, stderr=log, timeout=30)
        except (OSError, subprocess.SubprocessError):
            pass
    for sig, wait in ((signal.SIGTERM, 15), (signal.SIGKILL, 15)):
        if proc.poll() is not None:
            return
        try:
            os.killpg(proc.pid, sig)
        except ProcessLookupError:
            pass
        try:
            proc.wait(timeout=wait)
            return
        except subprocess.TimeoutExpired:
            pass

def wait_process(proc, timeout):
    """Turn timeouts and CLI interruption into persisted terminal return codes."""
    import signal
    previous = signal.getsignal(signal.SIGTERM)
    def interrupt(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupt)
    try:
        return proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        return 124
    except KeyboardInterrupt:
        return 130
    finally:
        signal.signal(signal.SIGTERM, previous)

def execute(workdir, stage, args, cfg, mpi=False):
    """Persist full logs; a timeout terminates only this invocation's container."""
    workdir = Path(workdir).resolve()
    cmd = [str(ROOT / "scripts/sdpb/sdpb.sh"), str(workdir)]
    if mpi:
        cmd += ["-n", str(cfg.nproc)]
    cmd += args
    cidfile = workdir / (stage + ".cid")
    env = dict(os.environ, SDPB_CIDFILE=str(cidfile))
    start = time.monotonic()
    with (workdir / (stage + ".log")).open("w") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env,
                                start_new_session=True)
        write_json(workdir / (stage + "_process.json"),
                   {"pid": proc.pid, "command": cmd, "started_unix": time.time()})
        try:
            rc = wait_process(proc, cfg.timeout + 60)
        finally:
            if proc.poll() is None:
                stop_owned_process(proc, cidfile, log)
    write_json(workdir / (stage + "_process.json"),
               {"pid": proc.pid, "command": cmd, "returncode": rc, "terminal": True})
    return {"returncode": rc, "seconds": time.monotonic()-start, "command": cmd}

def run_once(spec, workdir, direction=(1., 0.), fix_f00=None, settings=None, prepared=None,basis_source_report=None,exact_input=False,warm_start=False,resume=False,face=None,functional=None):
    if prepared and basis_source_report:raise ValueError('Prepared constraints and transferred basis are distinct reuse modes')
    if (face is not None or functional is not None) and (exact_input or warm_start):
        raise ValueError('The face diagnostic adds a block and changes the objective; no exact-input or warm-start reuse')
    if warm_start and (not prepared or exact_input):raise ValueError('Warm start requires the same accepted prepared support source')
    if resume and (not prepared or not exact_input or warm_start):raise ValueError('Resume requires a distinct exact-input refinement')
    origin=json.loads(Path(prepared).read_text()) if exact_input and prepared else None
    if exact_input and (origin is None or list(direction)!=origin['direction'] or fix_f00!=origin.get('fix_f00')):
        raise ValueError('Exact-input refinement requires the original objective and section')
    cfg, dest = settings or Settings(), Path(workdir).resolve()
    if resume:
        from .checkpoint import validate_resume_source
        validate_resume_source(prepared,cfg)
    dest.mkdir(parents=True, exist_ok=True)
    if any(dest.iterdir()):
        raise FileExistsError(f"Use a fresh result directory; preserving {dest}")
    rec = {"spec": asdict(spec), "settings": asdict(cfg), "direction": list(direction),
           "fix_f00": fix_f00, "face": face, "functional": functional, "status": "assembling",
           "accepted": False, "certified": False, "solver": "SDPB", "pmp_degree": 0}
    report_path = dest / "report.json"
    write_json(report_path, rec)
    start = time.monotonic()
    w = (Pmp.from_saved(prepared,direction,fix_f00,face,functional) if prepared else
         Pmp(spec, direction=direction, fix_f00=fix_f00, digits=cfg.digits,basis_source_report=basis_source_report,
             face=face, functional=functional))
    rec['face'],rec['functional'] = w.face,w.functional
    if prepared:
        if any(not np.array_equal(getattr(spec,k),getattr(w.spec,k)) for k in asdict(spec)):
            raise ValueError('Prepared constraints do not match the requested model')
        if cfg.digits != w.digits:
            raise ValueError('Prepared constraints retain their original coefficient precision')
        rec['reused_source_report'] = str(Path(prepared).resolve())
    if exact_input:
        shutil.copyfile(Path(prepared).parent/'pmp.json',dest/'pmp.json')
        rec['pmp']=dict(origin['pmp'],path=str(dest/'pmp.json'))
        rec['exact_input_source']={'report':str(Path(prepared).resolve()),'sha256':hashlib.sha256(Path(prepared).read_bytes()).hexdigest(),
                                   'checkpoint_loaded':False,'scope':'new solve of identical input; no source iterate restored'}
    else:rec["pmp"] = w.write(str(dest / "pmp.json"))
    rec["pmp"]["sha256"] = hashlib.file_digest((dest / "pmp.json").open("rb"), "sha256").hexdigest()
    if exact_input and rec['pmp']['sha256']!=origin['pmp']['sha256']:raise ValueError('Exact input changed while copying')
    if w.basis is not None:
        transfer=getattr(w,'basis_source',None)
        if exact_input:
            transfer={'source_basis_path':str(Path(prepared).parent/'basis.npy'),
                      'source_basis_sha256':origin['basis']['sha256'],'source_tolerance':origin['basis']['tolerance']}
        if transfer:
            raw=Path(transfer['source_basis_path']).read_bytes()
            if hashlib.sha256(raw).hexdigest()!=transfer['source_basis_sha256']:raise ValueError('Transferred basis changed during assembly')
            (dest/'basis.npy').write_bytes(raw)
        else:np.save(dest / "basis.npy", w.basis)
        from .assembly import BASIS_PACKING
        rec["basis"] = {"normalized_projection_error": w.ops.basis_projection_error,
                         "sha256": hashlib.sha256((dest / "basis.npy").read_bytes()).hexdigest(),
                         "tolerance": transfer['source_tolerance'] if transfer else spec.basis_tol,
                         "requested_tolerance":spec.basis_tol,"packing":BASIS_PACKING,
                         "construction":'saved coordinates' if transfer else 'prepared constraints' if prepared else 'canonical unscaled effective rows',
                         "source":transfer,
                         "scope": "fixed saved coordinates; fresh constraints and verification" if transfer else
                                  "floating-point SVD truncation, not an exact rank certificate"}
        if exact_input:rec['basis'].update(construction='identical-input refinement',scope='exact saved input coordinates; no new SVD')
    np.save(dest / "disk_mask.npy", w.keep)
    rec["source_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in sorted(Path(__file__).parent.glob("*.py"))}
    if exact_input:
        rec['verification_source_sha256']=rec['source_sha256'];rec['source_sha256']=origin['source_sha256']
    rec["status"] = "pmp2sdp"
    write_json(report_path, rec)
    rec["pmp2sdp"] = execute(dest, "pmp2sdp", ["pmp2sdp", f"--precision={cfg.precision}",
                              "--input=pmp.json", "--output=sdp"], cfg)
    if rec["pmp2sdp"]["returncode"]:
        rec["status"] = "preprocess_failed"
        rec["convergence"] = convergence_report({}, rec["pmp2sdp"]["returncode"])
        write_json(report_path, rec)
        return rec, w, None
    if warm_start or resume:
        from .checkpoint import prepare_warm_start
        try:rec['warm_start']=prepare_warm_start(prepared,dest,rec,cfg,resume=resume)
        except (ValueError,OSError,KeyError,TypeError) as exc:
            rec.update(status='warm_start_rejected',warm_start_error=str(exc))
            write_json(report_path,rec);return rec,w,None
    rec["status"] = "solving"
    write_json(report_path, rec)
    args = ["sdpb", f"--precision={cfg.precision}", "-s", "sdp", "-o", "out",
            f"--maxIterations={cfg.max_iterations}", f"--maxRuntime={cfg.timeout}",
            f"--checkpointInterval={cfg.checkpoint_interval}",
            f"--dualityGapThreshold={cfg.duality_gap}",
            f"--primalErrorThreshold={cfg.primal_error}",
            f"--dualErrorThreshold={cfg.dual_error}",
            f"--initialMatrixScalePrimal={cfg.initial_scale}",
            f"--initialMatrixScaleDual={cfg.initial_scale}"]
    if warm_start or resume:args += ['-i',rec['warm_start']['input_directory'],'-c',rec['warm_start']['output_checkpoint_directory']]
    rec["sdpb"] = execute(dest, "sdpb", args, cfg, mpi=True)
    if resume:rec['exact_input_source'].update(checkpoint_loaded='Loading binary checkpoint from' in (dest/'sdpb.log').read_text(),scope='same PMP and saved coordinates; checkpoint initialization attempted in an independent run; full acceptance repeated')
    rec["sdpb_result"] = read_out(dest / "out/out.txt")
    if rec['sdpb']['returncode']:
        rec['untrusted_partial_output']=rec['sdpb_result'];rec['sdpb_result']={}
        log=(dest/'sdpb.log').read_text()
        rec['solver_failure']='Cholesky(Q) numerical HPD failure' if 'Cholesky(Q)' in log else 'solver process failed'
        from .accuracy import trace_report
        rec['failure_trace']=trace_report(dest/'out/iterations.json',cfg.precision)
    sol = None
    yp = dest / "out/y.txt"
    if not rec['sdpb']['returncode'] and yp.exists() and yp.stat().st_size:
        try:
            sol = w.read_solution(yp)
            rec['y_sha256']=hashlib.sha256(yp.read_bytes()).hexdigest()
            np.savez_compressed(dest / "solution.npz", **{k: v for k, v in sol.items() if k != "y_text"})
            rec["solution_sha256"] = hashlib.sha256((dest / "solution.npz").read_bytes()).hexdigest()
            rec["verification"] = w.verify(sol)
        except (ValueError, OSError, FloatingPointError) as exc:
            rec["readback_error"] = str(exc)
            sol = None
    rec["convergence"] = convergence_report(rec["sdpb_result"],
        rec["sdpb"]["returncode"], rec.get("verification", {}).get("objective_recomputed"), cfg)
    rec["accepted"] = (rec["convergence"]["solver_optimal"] and
                       rec.get("verification", {}).get("primal_feasible", False))
    rec["status"] = "numerically_accepted" if rec["accepted"] else "not_accepted"
    if 'solver_failure' in rec:
        rec['status']='solver_failed';rec['convergence']['terminate_reason']=rec['solver_failure']
    if "readback_error" in rec:
        rec["status"] = "readback_failed"
    rec["seconds"] = time.monotonic()-start
    if sol is not None:
        from .observables import resonance_report, subthreshold_curves
        if w.precise:
            rec["observables"] = rec['verification'].pop('observables')
            rec["subthreshold"] = rec['verification'].pop('subthreshold')
        else:
            rec["observables"] = resonance_report(w.ops, sol["c"])
            rec["subthreshold"] = subthreshold_curves(w.ops, sol["c"], np.linspace(.05, 3.95, 40))
    write_json(report_path, rec)
    print(f"{dest.name}: {rec['status']} objective="
          f"{rec.get('verification', {}).get('objective_recomputed')} "
          f"termination={rec['convergence']['terminate_reason']}", flush=True)
    return rec, w, sol

def support_from_saved(source_report, outdir, point=None, gap=None,*,direction=None,fix_f00=None,warm_start=False,timeout=None,
                       face_margin=None,functional=None):
    """Registered or arbitrary support; reuse base constraints and exact saved basis.

    With ``face_margin`` and ``functional`` the run keeps the source direction and
    section, adds the slab ``d.(f00,f11) >= source value - margin`` and optimises
    the registered node functional over that near-optimal face instead.
    """
    from .spec import ModelSpec
    from .constraints import chiral_reference_point
    if face_margin is not None and functional is None:
        raise ValueError('A face slab needs a functional to optimise over it')
    if functional is not None:
        if point is not None or direction is not None or fix_f00 is not None:
            raise ValueError('A functional run keeps the source direction and section')
        if warm_start:raise ValueError('Warm start is not available for functional runs')
        if face_margin is not None and (isinstance(face_margin,bool) or not isinstance(face_margin,(int,float))
                                        or not np.isfinite(face_margin) or face_margin<=0):
            raise ValueError('face margin must be a positive finite number')
    if point is not None and (direction is not None or fix_f00 is not None):
        raise ValueError('Registered point conflicts with an arbitrary direction or section')
    if direction is None and fix_f00 is not None:raise ValueError('A custom section requires a direction')
    if direction is not None:
        direction=np.asarray(direction,dtype=float)
        if direction.shape!=(2,) or not np.all(np.isfinite(direction)) or not np.any(direction!=0):
            raise ValueError('Support direction must have two finite components and be nonzero')
        direction=tuple(direction.tolist())
    if fix_f00 is not None and not np.isfinite(fix_f00):raise ValueError('Fixed f00 must be finite')
    source = json.loads(Path(source_report).read_text())
    if (source.get('solver')!='SDPB' or source.get('status')!='numerically_accepted' or source.get('accepted') is not True
            or source.get('verification',{}).get('primal_feasible') is not True
            or source.get('convergence',{}).get('solver_optimal') is not True):
        raise ValueError('Supports require a completed numerically accepted source leaf')
    fixed,face=fix_f00,None
    if functional is not None:
        direction=tuple(float(v) for v in source['direction']);fixed=source.get('fix_f00')
        if face_margin is not None:
            face={'direction':list(direction),'value':source['verification']['objective_recomputed'],'margin':float(face_margin)}
    elif direction is None:
        point=point or 'ref';xr=chiral_reference_point()[0]
        if point not in ('tip','ref','mid'):raise ValueError('Unknown registered point')
        if point=='mid' and (source.get('fix_f00') is not None or source['direction']!=[1.,0.]):
            raise ValueError('Mid requires the accepted tip of this same model')
        direction = (1.,0.) if point=='tip' else (0.,1.)
        fixed = None if point=='tip' else xr if point=='ref' else (xr+source['verification']['f00_3'])/2
    cfg = Settings(**source['settings'])
    if gap is not None:cfg.duality_gap=gap
    if timeout is not None:
        if type(timeout) is not int or timeout<=0:raise ValueError('timeout must be a positive integer number of seconds')
        cfg.timeout=timeout
    result,_,_ = run_once(ModelSpec(**source['spec']),outdir,direction,fixed,cfg,prepared=source_report,warm_start=warm_start,
                          face=face,functional=functional)
    return result

def generate(spec, workdir, direction=(1., 0.), fix_f00=None, settings=None,
             start_tol=0., max_rounds=12, add_per_round=120,basis_source_report=None):
    """Nested disk relaxations; default all disks avoids unsupported sparse seeds.

    A failed optimization stops. A larger returned objective is not a proof of
    unboundedness. Every round is retained in its own directory.
    """
    from .assembly import Operators
    ops = Operators(spec.M,spec.L,spec.scattering_prescription,spec.operator_dps)
    mask = np.ones(3 * spec.L * spec.M, dtype=bool) if start_tol == 0 else (
        ops.nu_measured > start_tol * ops.nu_measured.max())
    history, previous = [], None
    root = Path(workdir)
    for rnd in range(max_rounds):
        rec, w, sol = run_once(replace(spec, disk_mask=mask.copy()),
                               root / f"round_{rnd:03d}", direction, fix_f00, settings,basis_source_report=basis_source_report)
        cv = rec.get("convergence", {"solver_optimal": False})
        history.append({"round": rnd, "report": f"round_{rnd:03d}/report.json",
                        "accepted": rec["accepted"], "convergence": cv})
        write_json(root / "generation.json", history)
        if not cv["solver_optimal"]:
            break
        # Intervals for optima of nested relaxations must permit nonincrease.
        if previous is not None and cv["candidate_lower"] > previous + 1e-8*max(1., abs(previous)):
            history[-1]["failure"] = "nested_relaxation_intervals_inconsistent"
            history[-1]["accepted"] = False
            write_json(root / "generation.json", history)
            rec["accepted"] = False
            rec["status"] = "nested_relaxation_intervals_inconsistent"
            write_json(root / f"round_{rnd:03d}/report.json", rec)
            break
        previous = cv["candidate_upper"]
        if rec["accepted"]:
            break
        hre, him = ops.h_re @ sol["c"], ops.h_im @ sol["c"]
        mag2 = hre**2 + him**2
        rel = (mag2-2*him) / np.maximum(mag2, 1e-300)
        bad = ((rel > 1e-6) & (mag2 > 1e-12)) | (mag2-2*him > 1e-8)
        new = np.flatnonzero(bad & ~mask)
        if not len(new):
            break
        new = new[np.argsort(rel[new])[::-1][:add_per_round]]
        mask[new] = True
    return rec
