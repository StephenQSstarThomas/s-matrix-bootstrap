"""Saved-amplitude audits and accuracy budgets; no optimization or fitted inputs."""
from pathlib import Path
import hashlib
import json
import math
import re

import numpy as np


def rho_diagnostic(ph):
    """Apply the frozen single-point part of C6; a high-energy crossing is insufficient."""
    E, d, eta = (np.asarray(ph[k], float) for k in ("E_GeV", "delta_deg", "eta"))
    low = (E >= .28) & (E <= 1.2)
    bracket, crossing = None, None
    for i in range(len(E)-1):
        if low[i] and low[i+1] and d[i] < 90 <= d[i+1]:
            bracket = [float(E[i]), float(E[i+1])]
            crossing = float(E[i] + (90-d[i])*(E[i+1]-E[i])/(d[i+1]-d[i]))
            break
    minimum = float(np.min(eta[low])) if np.any(low) else None
    finite = bool(np.any(low) and np.all(np.isfinite(d[low])) and np.all(np.isfinite(eta[low])))
    passed = finite and crossing is not None and .795 <= crossing <= .845 and minimum >= .9
    return {"low_energy_crossing_GeV": crossing, "native_crossing_bracket_GeV": bracket,
            "min_eta_below_1p2GeV": minimum, "single_point_C6_pass": bool(passed),
            "registered_crossing_window_GeV": [.795, .845], "analysis_window_GeV": [.28, 1.2],
            "scope": "nearest-node phase lift and linear crossing readout; no pole or off-node certificate"}


def trace_report(path, precision):
    """Allow a live final JSON fragment, but never interpret it as terminal state."""
    path = Path(path)
    if not path.exists():
        return {"iterations_observed": 0}
    rows = [json.loads(s) for s in re.findall(r'\{[^\n]+\}', path.read_text())]
    if not rows:
        return {"iterations_observed": 0}
    condition = max(float(r["Q_cond_number"]) for r in rows)
    times = [float(r["iter_time"]) for r in rows]
    return {"iterations_observed": len(rows), "max_Q_condition_number": condition,
            "estimated_digits_remaining": (precision-math.log2(condition))*math.log10(2),
            "median_iteration_seconds": float(np.median(times)),
            "iteration_seconds_range": [min(times), max(times)],
            "last_gap": float(rows[-1]["gap"]),
            "scope": "conditioning estimate only; requires final convergence and same-problem refinement"}


def current_mechanism(checker, im, rho_hat):
    """Separate elastic saturation, two-pion saturation and Watson alignment.

    det(G_hat)=(1-|S|²)(rho_hat-|F|²)-|F-S F*|². No new constraint or objective.
    """
    from flint import arb, acb
    from .arbaudit import _enclosure, _kinematic_square
    rows,moments = [],[]
    for ell,wave in ((0,'S0'),(1,'P1')):
        twopion,spectral = [],[]
        for i,S in enumerate(checker.last_primary[wave]):
            F = acb(1+sum((checker.K[i][j]*im[ell][j] for j in range(checker.M)),arb(0)), im[ell][i])
            elastic = 1-(S*S.conjugate()).real
            excess = rho_hat[ell][i]-(F*F.conjugate()).real
            residual = F-S*F.conjugate()
            residual2 = (residual*residual.conjugate()).real
            k2 = _kinematic_square(ell,checker.x[i])
            twopion.append(k2*(F*F.conjugate()).real);spectral.append(k2*rho_hat[ell][i])
            rows.append({'wave':wave,'node':i,'below_s0':bool(checker.x[i]<=arb(3600)/49),
                         'elastic_slack':_enclosure(elastic),'two_pion_excess_scaled':_enclosure(excess),
                         'watson_residual_squared':_enclosure(residual2),
                         'gram_determinant_scaled':_enclosure(elastic*excess-residual2)})
        for n in ((0,1) if ell==0 else (-1,0)):
            weights=[arb.pi()*checker.w[i]*x**n if x<=arb(3600)/49 else arb(0) for i,x in enumerate(checker.x)]
            part=sum((w*v for w,v in zip(weights,twopion)),arb(0))
            total=sum((w*v for w,v in zip(weights,spectral)),arb(0))
            moments.append({'wave':wave,'n':n,'two_pion_moment':_enclosure(part),'spectral_moment':_enclosure(total),
                            'two_pion_fraction':_enclosure(part/total) if total>0 else None})
    return {'rows':rows,'moment_saturation':moments,'scope':'native-node identity in rho_hat=rho/k² units; moments in raw m_pi=1 units; diagnostic only',
            'phase_caveat':'F=0 does not define a form-factor phase; PSD alone need not saturate elasticity or spectral density'}


def audit_saved_currents(source_report, outdir, bits=384, error_target='1e-12', timeout=600, snapshot=None):
    """Replay just the two native current waves from canonical saved coordinates.

    The deadline is checked between wave evaluations; a single quadrature is
    allowed to finish. No all-wave, chiral, objective or dual replay is done.
    """
    import time
    from flint import arb, acb, arb_mat, ctx
    from .arbaudit import ArbAudit, _enclosure, _kinematic_square
    from .projector import Layout
    from .sdpb import write_json
    start = time.monotonic()
    source, outdir = Path(source_report).resolve(), Path(outdir).resolve()
    if outdir.exists() and (not outdir.is_dir() or any(outdir.iterdir())):
        raise FileExistsError('Preserve existing output; choose a fresh directory')
    if bits < 64 or not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('Require bits >= 64 and a positive finite timeout')
    ctx.prec = bits
    if not arb(error_target).is_finite() or not arb(error_target) > 0:
        raise ValueError('Require a positive finite angular error target')
    rec = json.loads(source.read_text())
    if (rec.get('solver') != 'SDPB' or rec.get('status') not in ('numerically_accepted', 'not_accepted')
            or 'verification' not in rec or 'sdpb' not in rec):
        raise ValueError('Select a completed leaf SDPB report with a saved solution')
    spec, pmp = rec['spec'], rec['pmp']
    precision = pmp.get('operator_precision', {})
    if (spec.get('scattering_prescription') != 'sine-cardinal' or
            precision.get('scattering_prescription') != 'sine-cardinal'):
        raise ValueError('currents requires the recorded sine-cardinal analytic family')
    if (not spec.get('uv') or spec.get('operator_dps', 17) <= 17 or
            precision.get('source_dps') != spec['operator_dps']):
        raise ValueError('currents requires a precise source run with saved UV variables')
    M, L, na = spec['M'], spec['L'], pmp['n_a']
    if any(type(v) is not int or v <= 0 for v in (M, L, na)):
        raise ValueError('Invalid saved amplitude layout')
    n = Layout(M).n
    if (pmp['n_vars'] != 1+na+4*M or pmp.get('reduce_basis') != spec['reduce_basis'] or
            (not spec['reduce_basis'] and na != n) or na > n):
        raise ValueError('Saved variable layout is inconsistent')
    paths = [source, source.parent/'pmp.json', source.parent/'solution.npz', source.parent/'out/y.txt']
    if snapshot:
        paths.append(Path(snapshot).resolve())
    if spec['reduce_basis']:
        paths.append(source.parent/'basis.npy')
    def digest(path):
        with path.open('rb') as fh:
            return hashlib.file_digest(fh, 'sha256').hexdigest()
    hashes = {str(path): digest(path) for path in paths}
    for path, expected in ((paths[1], pmp['sha256']), (paths[2], rec.get('solution_sha256'))):
        if hashes[str(path)] != expected:
            raise ValueError(f'Saved input hash mismatch: {path.name}')
    modules = ('accuracy.py', '__main__.py', 'crosscheck.py','arbaudit.py', 'sine.py', 'precision.py',
               'projector.py', 'grid.py', 'constraints.py')
    producers = {name: digest(Path(__file__).with_name(name)) for name in modules}
    from .crosscheck import _angle_identity
    identity = _angle_identity(rec,snapshot)
    for name in ('grid.py',):
        if rec.get('source_sha256', {}).get('src/smatrix_bootstrap/sdp/'+name) != producers[name]:
            raise ValueError(f'Recorded analytic family/layout source hash mismatch: {name}')
    lines = paths[3].read_text().splitlines()
    if not lines or lines[0].split() != [str(na+4*M), '1']:
        raise ValueError('Full y text header does not match the saved variable layout')
    values = [line.strip() for line in lines[1:] if line.strip()]
    if len(values) != na+4*M:
        raise ValueError('Full y text length does not match the saved variable layout')
    y = [arb(value) for value in values]
    if not all(value.is_finite() for value in y):
        raise ValueError('Nonfinite saved y value')
    rounded = np.array([float(value) for value in values])
    with np.load(paths[2], allow_pickle=False) as sol:
        shapes = {'y':(na+4*M,), 'a':(na,), 'c':(n,), 'ImF':(2,M), 'rho_hat':(2,M)}
        if any(key not in sol or sol[key].shape != shape for key, shape in shapes.items()):
            raise ValueError('Saved solution arrays do not match the variable layout')
        arrays = {'y':rounded, 'a':rounded[:na], 'ImF':rounded[na:na+2*M].reshape(2,M),
                  'rho_hat':rounded[na+2*M:].reshape(2,M)}
        if any(not np.array_equal(sol[key], value) for key, value in arrays.items()):
            raise ValueError('Full y text differs from the authenticated saved solution at binary64')
    if rec.get('y_sha256') and rec['y_sha256'] != hashes[str(paths[3])]:
        raise ValueError('Recorded full y text hash mismatch')
    basis = None
    if spec['reduce_basis']:
        if hashes[str(paths[-1])] != rec.get('basis', {}).get('sha256'):
            raise ValueError('Source coordinate basis hash mismatch')
        basis = np.load(paths[-1], allow_pickle=False)
        if basis.dtype != np.dtype('float64') or basis.shape != (n, na) or not np.all(np.isfinite(basis)):
            raise ValueError('Saved basis must be finite binary64 with the recorded dimensions')
    checker = ArbAudit(M, L, bits, prescription='sine-cardinal', source_dps=spec['operator_dps'])
    av = arb_mat([[value] for value in y[:na]])
    c = y[:na] if basis is None else [
        (arb_mat([[arb(float(value)) for value in row]])*av)[0,0] for row in basis]
    im = [y[na+e*M:na+(e+1)*M] for e in (0,1)]
    rh = [y[na+2*M+e*M:na+2*M+(e+1)*M] for e in (0,1)]
    nodes, S_store = [], {}
    checker.last_primary = {'S0':[], 'P1':[]}
    for ell, wave in ((0,'S0'), (1,'P1')):
        for i, x in enumerate(checker.x):
            if time.monotonic()-start > timeout:
                raise TimeoutError('Current diagnostic deadline reached between native wave evaluations')
            f, info = checker.family.wave(c, ell, ell, node=i, error_target=error_target)
            ctx.prec = bits
            S = 1+acb(0,1)*arb.pi()*(1-4/x).sqrt()*f
            F = acb(1+sum((checker.K[i][j]*im[ell][j] for j in range(M)),arb(0)),im[ell][i])
            k2 = _kinematic_square(ell, x)
            S_store[ell,i] = S
            checker.last_primary[wave].append(S)
            nodes.append({'wave':wave, 'node':i, 's':_enclosure(x),
                'E_GeV':_enclosure(arb('0.14')*x.sqrt()), 'below_s0':bool(x<=arb(3600)/49),
                'S':{'real':_enclosure(S.real), 'imag':_enclosure(S.imag)},
                'F':{'real':_enclosure(F.real), 'imag':_enclosure(F.imag)},
                'rho':_enclosure(k2*rh[ell][i]), 'rho_hat':_enclosure(rh[ell][i]),
                'k_squared':_enclosure(k2), 'quadrature':info})
    uv = checker._uv_audit(S_store, im, rh, spec['sr_caliber'], spec['eps_ff'], spec['m_q'],
                           frozen_at_s0=spec['ff_frozen_at_s0'])
    report = {'source_report':str(source), 'source_sha256':hashes[str(source)], 'input_sha256':hashes,
        'producer_sha256':producers, 'recorded_source_sha256':rec.get('source_sha256', {}),
        'source_identity':identity,
        'source_numerically_accepted':rec.get('accepted', False), 'spec':spec,
        'certified':False, 'optimization_performed':False, 'bits':bits, 'angular_error_target':error_target,
        'seconds':time.monotonic()-start, 'timeout_seconds':timeout,
        'deadline_scope':'checked between native wave evaluations; one quadrature may finish after deadline',
        'layout':{'n_c':n, 'n_a':na, 'n_y':len(y), 'basis_shape':list(basis.shape) if basis is not None else None,
                  'y_order':'a, ImF_S0, ImF_P1, rho_hat_S0, rho_hat_P1; normalized y0 omitted'},
        'coefficient_reconstruction':'full SDPB y decimal text in Arb; saved binary64 basis entries are exact dyadics',
        'y_integrity_scope':('full text matches the solve-time recorded hash' if rec.get('y_sha256') else
            'fresh full-text hash; authenticated NPZ y matches at binary64; no solve-time sub-binary64 integrity hash'),
        'nodes':nodes, 'current_mechanism':current_mechanism(checker, im, rh), 'uv_checks':uv,
        'uv_parts_imposed':spec['uv_parts'],
        'scope':'S0/P1 native-node current diagnostic in the same recorded sine-cardinal family; all UV families evaluated; '
                'no full scattering/chiral feasibility, support optimality or continuum certificate'}
    outdir.mkdir(parents=True, exist_ok=True)
    write_json(outdir/'report.json', report)
    return report


def compare_reports(first, second):
    """Empirical same-objective refinement, never a bound on all near-optimal points."""
    def contract(r):
        spec = {k: v for k, v in r["spec"].items() if k not in ("tag", "disk_mask")}
        return spec, r["direction"], r.get("fix_f00")
    if contract(first) != contract(second):
        raise ValueError("Accuracy comparison needs the same finite model, objective and fixed section")
    if not first.get("pmp", {}).get("sha256") or first["pmp"]["sha256"] != second.get("pmp", {}).get("sha256"):
        raise ValueError("Solver refinement requires identical emitted PMP coefficients")
    if first["spec"].get("reduce_basis") and (not first.get("basis", {}).get("sha256") or
            first["basis"]["sha256"] != second.get("basis", {}).get("sha256")):
        raise ValueError("Solver refinement requires the same saved amplitude basis")
    out = {"both_numerically_accepted": first.get("accepted", False) and second.get("accepted", False),
           "scope": "measured drift, not a uniqueness or near-optimal-face certificate", "waves": {}}
    for wave in ("S0", "S2", "P1"):
        a, b = first["observables"][wave], second["observables"][wave]
        if a["E_GeV"] != b["E_GeV"]:
            raise ValueError("Native grids differ")
        da, db = np.asarray(a["delta_deg"]), np.asarray(b["delta_deg"])
        sa = np.asarray(a["eta"])*np.exp(2j*np.deg2rad(da))
        sb = np.asarray(b["eta"])*np.exp(2j*np.deg2rad(db))
        mask = np.asarray(a["E_GeV"]) <= 1.2
        out["waves"][wave] = {"max_complex_S_drift": float(np.max(abs(sa[mask]-sb[mask]))),
                              "max_phase_drift_deg": float(np.max(abs(da[mask]-db[mask])))}
    out["rho"] = [rho_diagnostic(r["observables"]["P1"]) for r in (first, second)]
    return out


def audit_saved_run(source_report, outdir, bits=384, compare=None):
    """Independently rebuild native functions with Arb from the saved complete c.

    Exact-source interval verdicts are retained even when they disagree with
    floating numerical feasibility. Primal and optimality certificates differ.
    """
    from .arbaudit import ArbAudit
    from .sdpb import write_json
    source, outdir = Path(source_report).resolve(), Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if (outdir / "report.json").exists():
        raise FileExistsError("Preserve existing audit; choose a new directory")
    rec = json.loads(source.read_text())
    if "spec" not in rec or "verification" not in rec:
        raise ValueError("Select a completed leaf SDPB report with a saved solution")
    solpath = source.parent / "solution.npz"
    sol = np.load(solpath)
    s = rec["spec"]
    report = {"source_report": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
              "solution_sha256": hashlib.sha256(solpath.read_bytes()).hexdigest(),
              "audit_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "source_numerically_accepted": rec.get("accepted", False), "certified": False,
              "trace": trace_report(source.parent / "out/iterations.json", rec["settings"]["precision"]),
              "rho": rho_diagnostic(rec["observables"]["P1"])}
    files = [source, solpath, source.parent / "pmp.json", source.parent / "out/y.txt"]
    files += sorted((source.parent / "out").glob("x_*.txt"))
    files += [source.parent / "basis.npy"] if s["reduce_basis"] else []
    files += [Path(compare).resolve()] if compare else []
    report["input_sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    report["producer_sha256"] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                  for name in ("accuracy.py", "arbaudit.py", "dual.py", "projector.py")}
    if report["input_sha256"][str(source.parent / "pmp.json")] != rec["pmp"]["sha256"]:
        raise ValueError("Emitted PMP hash differs from the solved input recorded in the report")
    checker = ArbAudit(s["M"],s["L"],bits,prescription=s.get('scattering_prescription','mixed-pv'),
                       source_dps=s.get('operator_dps',40))
    from flint import arb, arb_mat
    text = (source.parent / "out/y.txt").read_text().splitlines()[1:]
    y = [arb(v.strip()) for v in text if v.strip()]
    na, M = rec["pmp"]["n_a"], s["M"]
    if s["reduce_basis"]:
        basis = np.load(source.parent / "basis.npy")
        av = arb_mat([[v] for v in y[:na]])
        c = [(arb_mat([[arb(float(v)) for v in row]])*av)[0, 0] for row in basis]
    else:
        c = y[:na]
    report["coefficient_reconstruction"] = "saved basis as exact binary64 constants times full SDPB y text in Arb"
    report["max_c_float_rounding_difference"] = max(float((v-arb(float(old))).abs_upper())
                                                   for v, old in zip(c, sol["c"]))
    im = [y[na+e*M:na+(e+1)*M] for e in (0, 1)] if s["uv"] else None
    rh = [y[na+2*M+e*M:na+2*M+(e+1)*M] for e in (0, 1)] if s["uv"] else None
    report["original_arithmetic"] = checker.audit(c, im, rh,
        chi_caliber=s["chi_caliber"] if s["chiral"] else None,
        eps_chi=s["eps_chi"], sr_caliber=s["sr_caliber"], eps_ff=s["eps_ff"],
        m_q=s["m_q"], ff_frozen_at_s0=s["ff_frozen_at_s0"])
    original = report["original_arithmetic"]
    if s['uv']:
        report['current_mechanism'] = current_mechanism(checker,im,rh)
    parts = {"gram": "gram_ok", "fesr": "fesr_ok", "ff": "form_factor_ok"}
    report["original_primal_certified"] = bool(s["B"] is None and original["unitarity_ok"] and
        (not s["chiral"] or original.get("chiral_verdict") == "certified_pass") and
        (not s["uv"] or (all(original[parts[p]] for p in s["uv_parts"]) and
                          all(v >= 0 for wave in rh for v in wave))))
    from .dual import replay_dual
    report["emitted_pmp_dual"] = replay_dual(source.parent, bits)
    if s["uv"]:
        report["uv_parts_imposed"] = s["uv_parts"]
        report["uv_note"] = "Independent audit evaluates all UV families; only imposed families belong to feasibility"
    if compare:
        other = json.loads(Path(compare).read_text())
        if s["reduce_basis"]:
            rec["basis"]["sha256"] = report["input_sha256"][str(source.parent / "basis.npy")]
            other["basis"]["sha256"] = hashlib.sha256((Path(compare).parent / "basis.npy").read_bytes()).hexdigest()
        other_pmp = Path(compare).parent / "pmp.json"
        digest = hashlib.sha256(other_pmp.read_bytes()).hexdigest()
        if digest != other["pmp"]["sha256"]:
            raise ValueError("Comparison PMP differs from the comparison solve record")
        report["input_sha256"][str(other_pmp.resolve())] = digest
        if s["reduce_basis"]:
            report["input_sha256"][str((Path(compare).parent / "basis.npy").resolve())] = other["basis"]["sha256"]
        report["refinement"] = compare_reports(rec, other)
    report["remaining_accuracy_requirements"] = ["accurate source-to-PMP coefficients", "independent original-model dual bound",
        "same-problem tighter-gap and higher-precision comparison", "near-optimal amplitude variation for C6/C7"]
    write_json(outdir / "report.json", report)
    return report
