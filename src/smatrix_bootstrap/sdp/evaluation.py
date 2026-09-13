"""Verified post-solve overlays for saved amplitudes; no optimization calls."""
from pathlib import Path
import hashlib
import json
import time

import numpy as np

from .observables import WAVES, arb_subthreshold_curves, subthreshold_grid, subthreshold_coverage

KERNELS = ('precision.py','sine.py','projector.py','legendreq.py','hilbert.py','grid.py')


def _digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def _data_digest(data):
    return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def _saved_amplitude(source):
    """Use precise reuse validation, then retain full y decimal text in Arb."""
    from flint import arb, arb_mat, ctx
    from .arbaudit import ArbAudit
    from .pmp import Pmp
    from .projector import Layout
    rec = json.loads(source.read_text())
    if (rec.get('solver') != 'SDPB' or rec.get('status') not in ('numerically_accepted','not_accepted')
            or 'verification' not in rec or 'sdpb' not in rec or rec['spec'].get('operator_dps',17) <= 17):
        raise ValueError('A completed precise SDPB leaf with a saved full solution is required')
    spec, pmp, root = rec['spec'], rec['pmp'], source.parent
    M, L, na = spec['M'], spec['L'], pmp['n_a']
    if any(type(v) is not int or v <= 0 for v in (M,L,na)):
        raise ValueError('Invalid source amplitude dimensions')
    n, ny = Layout(M).n, na+(4*M if spec['uv'] else 0)
    if (pmp['n_vars'] != ny+1 or pmp.get('reduce_basis') != spec['reduce_basis'] or
            na > n or (not spec['reduce_basis'] and na != n)):
        raise ValueError('Source variable layout mismatch')
    paths = [source,root/'pmp.json',root/'out/y.txt',root/'solution.npz',root/'disk_mask.npy']
    if spec['reduce_basis']:
        paths.append(root/'basis.npy')
        basis = np.load(paths[-1],allow_pickle=False)
        if basis.dtype != np.dtype('float64') or basis.shape != (n,na) or not np.all(np.isfinite(basis)):
            raise ValueError('Source basis must be finite binary64 with the recorded dimensions')
    hashes = {str(path):_digest(path) for path in paths}
    for name,expected in [('pmp.json',pmp['sha256']),('solution.npz',rec.get('solution_sha256')),
                           ('out/y.txt',rec.get('y_sha256'))]:
        if expected is None or hashes[str(root/name)] != expected:
            raise ValueError(f'Source input hash mismatch or missing solve-time hash: {name}')
    lines = (root/'out/y.txt').read_text().splitlines()
    if not lines or lines[0].split() != [str(ny),'1']:
        raise ValueError('Full y text header does not match the saved layout')
    # This restores only source kernels, saved basis and projection-plane rows.
    # It never writes a PMP, calls a solver, or replays a dual certificate.
    w = Pmp.from_saved(source,rec['direction'],rec.get('fix_f00'))
    sol = w.read_solution(root/'out/y.txt')
    if not np.all(np.isfinite(sol['y'])):
        raise ValueError('Nonfinite source y values')
    with np.load(root/'solution.npz',allow_pickle=False) as saved:
        shapes = {'y':(ny,), 'a':(na,), 'c':(n,)}
        if spec['uv']:
            shapes.update(ImF=(2,M),rho_hat=(2,M))
        if any(key not in saved or saved[key].shape != shape for key,shape in shapes.items()):
            raise ValueError('Saved solution array layout mismatch')
        if any(not np.array_equal(saved[key],sol[key]) for key in shapes if key != 'c'):
            raise ValueError('Full y text differs from the authenticated binary64 solution copy')
    checker = ArbAudit(M,L,max(256,w.precise.bits),prescription=spec['scattering_prescription'],
                       source_dps=spec['operator_dps'])
    ctx.prec = checker.bits
    av = arb_mat([[arb(value)] for value in sol['y_text'][:na]])
    cv = av if w.precise_basis is None else w.precise_basis*av
    c = [cv[i,0] for i in range(cv.nrows())]
    if not all(value.is_finite() for value in c):
        raise ValueError('Nonfinite reconstructed source amplitude')
    return rec, checker, c, hashes


def evaluate_subthreshold(source_report, outdir, timeout=600):
    """Re-evaluate the same saved amplitude on the complete Fig.5 display grid.

    Precision and analytic family are inherited from the source producer.
    Requires its recorded full-y hash, rather than authenticating only rounding.
    """
    from .sdpb import write_json
    start = time.monotonic()
    source, dest = Path(source_report).resolve(), Path(outdir).resolve()
    if dest.exists() and (not dest.is_dir() or any(dest.iterdir())):
        raise FileExistsError('Use a fresh subthreshold result directory')
    if not np.isfinite(timeout) or timeout <= 0:
        raise ValueError('A positive finite evaluation timeout is required')
    rec, checker, c, hashes = _saved_amplitude(source)
    old = rec.get('subthreshold',{})
    missing = subthreshold_coverage(old,old.get('s',[]))
    if missing:
        raise ValueError('Saved subthreshold display cannot be checked: '+missing)
    base = np.unique(np.r_[np.linspace(.05,3.95,40),np.asarray(old['s'],float)])
    data = arb_subthreshold_curves(checker,c,base,deadline=start+timeout)
    indices = np.searchsorted(data['s'],old['s'])
    if not np.array_equal(np.asarray(data['s'])[indices],old['s']):
        raise ValueError('Saved display samples were not retained exactly')
    comparison = {}
    for wave in WAVES:
        previous = np.asarray(old[wave],float)
        delta = np.abs(np.asarray(data[wave])[indices]-previous)
        tolerance = 64*np.finfo(float).eps*np.maximum(1.,np.abs(previous))
        comparison[wave] = {'samples':len(indices),'max_absolute_difference':float(delta.max()),
            'consistent':bool(np.all(delta <= tolerance)),
            'comparison_tolerance':'64 binary64 epsilon times max(1,abs(saved value)); display consistency only'}
    if not all(row['consistent'] for row in comparison.values()):
        raise ValueError('Re-evaluation disagrees with saved samples of the same amplitude')
    grid = data['grid_provenance']
    hashes[grid['reference_csv']] = grid['reference_csv_sha256']
    if any(_digest(path) != digest for path,digest in hashes.items()):
        raise ValueError('A saved input changed during subthreshold evaluation')
    report = {'kind':'subthreshold_evaluation', 'status':'complete', 'source_report':str(source),
        'source_sha256':hashes[str(source)], 'input_sha256':hashes, 'spec':rec['spec'],
        'source_numerically_accepted':rec.get('accepted',False), 'bits':checker.bits,
        'source_operator_dps':rec['spec']['operator_dps'], 'optimization_performed':False, 'certified':False,
        'producer_sha256':{name:_digest(Path(__file__).with_name(name)) for name in
            (*KERNELS,'evaluation.py','observables.py','arbaudit.py','precision_pmp.py','pmp.py')},
        'coefficient_reconstruction':'recorded full y decimal text in Arb times exact binary64 saved basis; rounded NPZ c unused',
        'old_sample_consistency':comparison, 'subthreshold':data, 'subthreshold_sha256':_data_digest(data),
        'seconds':time.monotonic()-start, 'timeout_seconds':timeout,
        'deadline_scope':'checked between individual evaluations; one quadrature may finish after the deadline',
        'scope':'same complete source amplitude, post-solve Fig.5 evaluation only; no new constraints, objective, '
                'selection, primal feasibility or dual/support certificate; C3 compares same-epsilon paper curves, not a linearity budget'}
    dest.mkdir(parents=True,exist_ok=True)
    if any(dest.iterdir()):
        raise FileExistsError('Subthreshold output directory became nonempty during evaluation')
    write_json(dest/'report.json',report)
    return report


def attach_subthreshold_overlays(records, reports):
    """Explicit verified in-memory overlays; the original solve reports stay intact."""
    output = [dict(record) for record in records]
    attached = set()
    for supplied in reports:
        path = Path(supplied).resolve()
        if path.is_dir():
            path = path/'report.json'
        replay = json.loads(path.read_text())
        if replay.get('kind') != 'subthreshold_evaluation' or replay.get('status') != 'complete':
            raise ValueError('A completed subthreshold evaluation report is required')
        source = str(Path(replay['source_report']).resolve())
        if source in attached:
            raise ValueError('Multiple subthreshold overlays for one source are ambiguous')
        attached.add(source)
        matches = [record for record in output if record.get('_source_report') == source]
        if len(matches) != 1:
            raise ValueError('Subthreshold overlay must match exactly one loaded source leaf')
        expected = replay.get('input_sha256',{})
        original = json.loads(Path(source).read_text())
        required = [Path(source),Path(source).parent/'pmp.json',Path(source).parent/'out/y.txt',
                    Path(source).parent/'solution.npz',Path(source).parent/'disk_mask.npy']
        if original['spec']['reduce_basis']:
            required.append(Path(source).parent/'basis.npy')
        if any(str(p) not in expected for p in required):
            raise ValueError('Subthreshold overlay is missing source input hashes')
        if replay.get('source_sha256') != expected[source] or any(_digest(p) != h for p,h in expected.items()):
            raise ValueError('Subthreshold overlay source input hash mismatch')
        for name in KERNELS:
            if replay.get('producer_sha256',{}).get(name) != _digest(Path(__file__).with_name(name)):
                raise ValueError('Subthreshold overlay source kernel changed: '+name)
        data = replay['subthreshold']
        grid, provenance = subthreshold_grid()
        if (_data_digest(data) != replay.get('subthreshold_sha256') or
                data.get('grid_provenance',{}).get('reference_csv_sha256') != provenance['reference_csv_sha256'] or
                data['grid_provenance'].get('evaluation_grid') != data['s'] or
                not np.all(np.isin(grid,data['s'])) or subthreshold_coverage(data,grid)):
            raise ValueError('Subthreshold overlay data/hash/reference grid mismatch')
        record = matches[0]
        if replay.get('spec') != record['spec'] or record['spec'] != original['spec']:
            raise ValueError('Subthreshold overlay belongs to a different source model')
        record['subthreshold'] = data
        record['subthreshold_evaluation'] = {'report':str(path),'sha256':_digest(path),
            'source_sha256':replay['source_sha256'],'scope':'verified in-memory evaluation overlay; original solve unchanged'}
    return output


def subthreshold_claim_tables(records):
    """C3 tables for existing objectives, with no choice among duplicate states."""
    from .claims import c3
    from .figures import accepted_support
    from .assembly import recorded_basis_key
    groups = {}
    for record in records:
        spec = record['spec']
        if not accepted_support(record) or not spec['chiral'] or spec['uv'] or spec['eps_chi'] not in (.002,.004,.006):
            continue
        identity = {'spec':{k:v for k,v in spec.items() if k not in ('eps_chi','tag','disk_mask')},
            'basis':recorded_basis_key(record), 'direction':record.get('result',{}).get('direction'),
            'fix_f00':record.get('fix_f00')}
        key = _data_digest(identity)[:16]
        groups.setdefault(key,[]).append(record)
    tables = {}
    for key, group in groups.items():
        by_eps = {record['spec']['eps_chi']:record for record in group}
        claim = ({'verdict':'not run','evidence':'multiple source amplitudes at one epsilon; explicit representative selection required'}
                 if len(by_eps) != len(group) else c3(by_eps))
        tables[key] = {'claim':claim, 'available_epsilon':sorted(by_eps),
            'source_reports':[record.get('_source_report') for record in group],
            'subthreshold_replays':[record['subthreshold_evaluation'] for record in group if 'subthreshold_evaluation' in record]}
    return tables
