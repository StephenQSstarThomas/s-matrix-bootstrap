"""Compare Python's 2309 operators with independent 100-digit Mathematica formulas."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

import numpy as np

from . import constraints as C, formfactor as F, grid

MMA_TIMEOUT_SECONDS = 18000
NONZERO_RELATIVE_TOLERANCE = 1e-10
ZERO_ROUNDOFF_EPS_FACTOR = 32

def _compare_values(value, target):
    """Separate relative accuracy from float64 rounding at exact reference zeros.

    The zero budget is 32 machine epsilons at the reference array scale;
    it is a rounding diagnostic, never an absolute floor for nonzero values.
    """
    value, target = np.asarray(value, dtype=float), np.asarray(target, dtype=float)
    if (target.shape != value.shape or not target.size or
            not np.all(np.isfinite(target)) or not np.all(np.isfinite(value))):
        raise ValueError('Mathematica comparison shape/finite-value mismatch')
    zero = target == 0
    absolute = np.abs(value-target)
    relative = float(np.max(absolute[~zero]/np.abs(target[~zero]), initial=0.))
    zero_error = float(np.max(absolute[zero], initial=0.))
    scale = max(1., float(np.max(np.abs(target))))
    budget = ZERO_ROUNDOFF_EPS_FACTOR*np.finfo(np.float64).eps*scale
    return {'max_nonzero_relative_error':relative,
            'nonzero_relative_tolerance':NONZERO_RELATIVE_TOLERANCE,
            'max_absolute_error_at_exact_zeros':zero_error,
            'exact_zero_absolute_budget':budget, 'reference_scale':scale,
            'exact_zero_count':int(zero.sum()), 'nonzero_count':int((~zero).sum()),
            'max_absolute_error':float(np.max(absolute)),
            'passed':bool(relative < NONZERO_RELATIVE_TOLERANCE and zero_error <= budget)}

def _comparison_fields(comparison):
    return {'comparison':comparison,
            'relative_errors':{key:row['max_nonzero_relative_error'] for key,row in comparison.items()},
            'relative_errors_scope':'nonzero reference entries only; no denominator floor',
            'zero_absolute_errors':{key:row['max_absolute_error_at_exact_zeros'] for key,row in comparison.items()},
            'comparison_policy':'nonzero relative <1e-10; exact-zero absolute <=32*eps_float64*max(1,maxabs(reference array))'}

def _valid_independent_metadata(ref):
    precision = ref.get('working_precision')
    return (ref.get('gram_identity') is True and type(precision) in (int,float) and
            np.isfinite(precision) and precision >= 100 and
            isinstance(ref.get('version'),str) and bool(ref['version']))

def _compare(M, ref):
    ours = {"nodes": grid.s_nodes(M), "weights": grid.quad_weights(M),
            "hilbert": F.hilbert_kernel(M),
            "kinematic": np.array([F.kinematic_factor(e, grid.s_nodes(M)) for e in (0, 1)]),
            "targets": np.array(list(C.printed_targets().values()))}
    errors = {}
    for key,value in ours.items():
        if key not in ref:
            raise ValueError(f'Mathematica reference array missing: {key}')
        errors[key] = _compare_values(value, ref[key])
    return errors


def reuse_audit(M, workdir, reference_report):
    """Recheck current arithmetic against an authenticated completed MMA control."""
    source,path = Path(reference_report).resolve(),Path(workdir).resolve()
    source_bytes = source.read_bytes()
    report = json.loads(source_bytes)
    script = Path(__file__).resolve().parents[3]/'scripts/mma/audit_2309.wls'
    digest = hashlib.sha256(script.read_bytes()).hexdigest()
    if (type(report.get('passed')) is not bool or type(report.get('returncode')) is not int or
            report['returncode'] != 0 or report.get('script_sha256') != digest or report.get('M') != M or
            not _valid_independent_metadata(report)):
        raise ValueError('Mathematica reuse requires authenticated successful independent execution of the current script')
    original = source.parent/'independent.json';reference_bytes=original.read_bytes()
    ref=json.loads(reference_bytes);reference_hash=hashlib.sha256(reference_bytes).hexdigest()
    if report.get('independent_sha256')!=reference_hash:
        raise ValueError('Mathematica reference output hash mismatch')
    errors = _compare(M,ref)
    if (not _valid_independent_metadata(ref) or
            any(ref.get(key)!=report.get(key) for key in ('working_precision','version','gram_identity')) or
            not all(row['passed'] for row in errors.values())):
        raise ValueError('Mathematica reference fails current formula comparison')
    if path.exists() and any(path.iterdir()):
        raise FileExistsError('Preserve existing Mathematica evidence; use a fresh directory')
    path.mkdir(parents=True,exist_ok=True)
    (path/'independent.json').write_bytes(reference_bytes)
    if source.read_bytes()!=source_bytes or original.read_bytes()!=reference_bytes:
        raise ValueError('Mathematica source evidence changed during recheck')
    result = {'passed':True,'returncode':0,'reused':True,'source_report':str(source),
              'source_report_sha256':hashlib.sha256(source_bytes).hexdigest(),
              'source_comparison_passed':report['passed'],'source_execution_returncode':report['returncode'],
              'independent_sha256':reference_hash,'script_sha256':digest,**_comparison_fields(errors),
              'working_precision':ref['working_precision'],'version':ref['version'],
              'gram_identity':True,'M':M,
              'scope':'same-script Mathematica output; current Python formulas rechecked, no solver certificate'}
    (path/'report.json').write_text(json.dumps(result,indent=2))
    return result


def audit(M, workdir):
    root = Path(__file__).resolve().parents[3]
    path = Path(workdir).resolve()
    path.mkdir(parents=True, exist_ok=True)
    script = root / "scripts/mma/audit_2309.wls"
    cmd = [str(root / "scripts/mma/wolfram.sh"), str(script),
           f"rw:{path}:/audit", "--", str(M), "/audit/independent.json"]
    cidfile = path / "mma.cid"
    with (path / "mma.log").open("w") as log:
        from .sdpb import stop_owned_process, wait_process
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
            env=dict(os.environ, MMA_CIDFILE=str(cidfile)), start_new_session=True)
        try:
            rc = wait_process(proc, MMA_TIMEOUT_SECONDS)
        finally:
            if proc.poll() is None:
                stop_owned_process(proc, cidfile, log)
    result = {"returncode": rc, "command": cmd,
              "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(), "passed": False}
    if rc == 0 and (path / "independent.json").exists():
        ref = json.loads((path / "independent.json").read_text())
        errors = _compare(M,ref)
        result.update(**_comparison_fields(errors), gram_identity=ref["gram_identity"],
                      version=ref["version"], working_precision=ref["working_precision"],
                      M=M,independent_sha256=hashlib.sha256((path/'independent.json').read_bytes()).hexdigest())
        result["passed"] = _valid_independent_metadata(ref) and all(row['passed'] for row in errors.values())
    (path / "report.json").write_text(json.dumps(result, indent=2))
    return result
