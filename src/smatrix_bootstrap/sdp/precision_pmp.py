"""Preserve source precision across the PMP coordinate change and readback."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

import numpy as np
from flint import arb, arb_mat, ctx

from .precision import PrecisionRows

_CACHE = {}


def restore(cls, source_report, direction, fix_f00, face=None, functional=None, allow_face_source=False):
    """Reuse a completed precise PMP's constraints and coordinate basis.

    ``allow_face_source`` re-opens a face-diagnostic leaf for its own post-processing (readback and verification of
    the same PMP); chaining a new support from such a leaf stays refused because its PMP carries the extra slab block.
    """
    from .assembly import Operators
    from .spec import ModelSpec
    from .pmp import check_face, check_functional
    root = Path(source_report).resolve().parent
    source = json.loads(Path(source_report).read_text())
    if source.get('solver') != 'SDPB' or source['spec'].get('operator_dps',17)<=17:
        raise ValueError('Reuse requires a precise SDPB source report')
    # Constraint bytes and a newly computed objective must use the same source.
    producers = source.get('source_sha256',{})
    # The saved basis bypasses assembly; only the plane/source kernels are rebuilt.
    for name in ('precision.py','sine.py','projector.py','legendreq.py','hilbert.py','grid.py'):
        key = 'src/smatrix_bootstrap/sdp/'+name
        if producers.get(key) != hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest():
            raise ValueError(f'Source operator changed: {name}; rebuild constraints before optimizing')
    scheme = source['spec'].get('scattering_prescription','mixed-pv')
    if source['pmp'].get('operator_precision',{}).get('scattering_prescription') != scheme:
        raise ValueError('Source operator prescription metadata mismatch')
    if hashlib.sha256((root/'pmp.json').read_bytes()).hexdigest()!=source['pmp']['sha256']:
        raise ValueError('Source PMP hash mismatch')
    w = cls.__new__(cls)
    w.spec = ModelSpec(**source['spec'])
    w.direction,w.fix_f00,w.digits = direction,fix_f00,source['pmp']['digits']
    if source.get('face') is not None and not allow_face_source:
        raise ValueError('A face-diagnostic leaf is terminal (it carries an extra slab block); chain supports from its own source leaf')
    # a functional-only leaf (Watson step, node functional without a slab) has exactly the base blocks and may be chained
    w.face,w.functional = check_face(face),check_functional(functional,w.spec)
    w.basis = np.load(root/'basis.npy') if w.spec.reduce_basis else None
    if w.basis is not None and hashlib.sha256((root/'basis.npy').read_bytes()).hexdigest()!=source['basis']['sha256']:
        raise ValueError('Source coordinate basis hash mismatch')
    if w.spec.scattering_prescription=='sine-cardinal':
        from .sine import SineFamily
        w.precise=SineFamily(w.spec.M,w.spec.L,w.spec.operator_dps)
    else:w.precise=PrecisionRows(w.spec.M,w.spec.operator_dps)
    from types import SimpleNamespace
    from .projector import Layout,ells_for
    w.ops=SimpleNamespace(M=w.spec.M,L=w.spec.L,lay=Layout(w.spec.M),
        s=np.array([float(v.mid()) for v in w.precise.s]),
        index=[(I,e) for I in (0,1,2) for e in ells_for(I,w.spec.L)])
    w.n_a,w.n_vars,w.i_a = source['pmp']['n_a'],source['pmp']['n_vars'],1
    w.i_ImF = 1+w.n_a if w.spec.uv else None
    w.i_rho = w.i_ImF+2*w.spec.M if w.spec.uv else None
    w.n_rho = w.ops.lay.r2.stop-w.ops.lay.r1.start
    w.n_aux = {None:0,'linf':0,'l2':w.n_rho,'l4':2*w.n_rho}[w.spec.reg_norm]
    w.i_aux = w.n_vars-w.n_aux
    if w.n_vars != 1+w.n_a+(4*w.spec.M if w.spec.uv else 0)+w.n_aux or source['pmp'].get('n_aux',0)!=w.n_aux:
        raise ValueError('Saved variable layout is inconsistent')
    w.keep = np.load(root/'disk_mask.npy')
    if w.keep.shape != (3*w.spec.L*w.spec.M,) or w.keep.dtype != np.bool_:
        raise ValueError('Source disk mask layout mismatch')
    expected = (np.ones(len(w.ops.index)*w.spec.M,dtype=bool) if w.spec.disk_mask is None
                else np.asarray(w.spec.disk_mask,dtype=bool).copy())
    if expected.shape != w.keep.shape:raise ValueError('ModelSpec disk mask layout mismatch')
    if w.spec.uv and 'gram' in w.spec.uv_parts:
        for ell,I in ((0,0),(1,1)):
            a0 = w.ops.index.index((I,ell))
            expected[a0*w.spec.M:(a0+1)*w.spec.M] = False
    if not np.array_equal(w.keep,expected):
        raise ValueError('Source disk mask content differs from the ModelSpec effective mask')
    if w.basis is not None and w.basis.shape != (w.ops.lay.n,w.n_a):
        raise ValueError('Source basis dimensions mismatch')
    w.precise_basis = None if w.basis is None else arb_mat([[arb(float(x)) for x in row] for row in w.basis])
    xy = w.precise.projected([w.precise.rows(0,0,point=3)[0],w.precise.rows(1,1,point=3)[0]],w.precise_basis)
    w.f00,w.f11 = xy[0],xy[1]
    w.precision_info = source['pmp']['operator_precision']
    w.ops.basis_projection_error = source.get('basis',{}).get('normalized_projection_error')
    w.prepared = root,source
    return w


def write_prepared(w, path):
    """Preserve all source block coefficients; replace only the selected section/objective."""
    from .dual import _stream_pmp
    root,source = w.prepared
    blocks = _stream_pmp(root/'pmp.json'); header = next(blocks)
    count = source['pmp']['n_blocks']-(2 if source.get('fix_f00') is not None else 0)
    old_section = None if source.get('fix_f00') is None else w._row(y0=source['fix_f00'],a=-w.f00)
    sizes = {}
    with open(path,'w') as stream:
        stream.write('{"objective":'+json.dumps(w._num_vec(w.objective()))+',"normalization":')
        stream.write(json.dumps(header['normalization'])+',"PositiveMatrixWithPrefactorArray":[')
        n = 0
        actual = 0
        for j,block in enumerate(blocks):
            actual += 1
            if j>=count:
                if old_section is None or j-count>=2:
                    raise ValueError('Unexpected source PMP trailing blocks')
                expected = [[[[v] for v in w._num_vec((1 if j==count else -1)*old_section)]]]
                if block['polynomials'] != expected:
                    raise ValueError('Source trailing blocks are not the recorded fixed section')
                continue
            if n:stream.write(',')
            stream.write(json.dumps(block,separators=(',',':')))
            size=len(block['polynomials']);sizes[size]=sizes.get(size,0)+1;n+=1
        if actual != source['pmp']['n_blocks']:
            raise ValueError('Source PMP block count mismatch')
        if w.fix_f00 is not None:
            row = w._row(y0=w.fix_f00,a=-w.f00)
            for sign in (1,-1):
                block={'DampedRational':{'base':'1','constant':'1','poles':[]},
                       'polynomials':[[[[v] for v in w._num_vec(sign*row)]]]}
                stream.write(',' if n else '');stream.write(json.dumps(block,separators=(',',':')));n+=1
            sizes[1]=sizes.get(1,0)+2
        if w.face is not None:
            block={'DampedRational':{'base':'1','constant':'1','poles':[]},
                   'polynomials':[[[[v] for v in w._num_vec(w.face_row())]]]}
            stream.write(',' if n else '');stream.write(json.dumps(block,separators=(',',':')));n+=1
            sizes[1]=sizes.get(1,0)+1
        stream.write(']}')
    info = dict(source['pmp'])
    info.update(path=str(path),bytes=Path(path).stat().st_size,n_blocks=n,block_sizes=sizes,
                source_pmp_sha256=source['pmp']['sha256'],source_report=str(root/'report.json'))
    return info


def populate(w):
    """Prepare high-precision projected rows; cache one model for tip/ref/mid."""
    spec, ops, basis = w.spec, w.ops, w.basis
    if spec.sparsify:
        raise ValueError("Precise source assembly does not sparsify coefficients")
    p = getattr(ops,'precise_source',None) or PrecisionRows(spec.M,spec.operator_dps)
    ctx.threads = 8
    needed = w.keep.copy()
    if spec.uv and 'gram' in spec.uv_parts:
        for wave in ((0,0),(1,1)):
            pos = ops.index.index(wave)*spec.M
            needed[pos:pos+spec.M] = True
    key = (spec.M,spec.L,spec.operator_dps,spec.scattering_prescription,needed.tobytes(),
           None if basis is None else hashlib.sha256(basis.tobytes()).hexdigest())
    started = time.monotonic()
    cached = _CACHE.get(key)
    if cached is None:
        vm = None if basis is None else arb_mat([[arb(float(x)) for x in row] for row in basis])
        na = ops.lay.n if basis is None else basis.shape[1]
        H = np.zeros((len(needed),2,na),dtype=object)
        ids = np.flatnonzero(needed)
        for begin in range(0,len(ids),16):
            batch = ids[begin:begin+16]
            rows = []
            for index in batch:
                a,k = divmod(int(index),spec.M); I,ell = ops.index[a]
                rows.extend(ops.precise_h[index] if hasattr(ops,'precise_h') else p.rows(I,ell,node=k)*p.kap[k])
            result = p.projected(rows,vm).reshape(len(batch),2,na)
            H[batch] = result
            if begin%64 == 0 or begin+16 >= len(ids):
                print(f"precise rows {min(begin+16,len(ids))}/{len(ids)} "
                      f"in {time.monotonic()-started:.1f}s",flush=True)
        if hasattr(ops,'precise_sub'):sub=ops.precise_sub
        else:
            sub = [p.rows(0,0,point=3)[0],p.rows(1,1,point=3)[0]]
            for sj in (arb(1)/2,arb(1),arb(3)/2,arb(2)):
                f0,f1,f2 = (p.rows(I,ell,point=sj)[0] for I,ell in ((0,0),(1,1),(2,0)))
                sub.extend([f0-3*(2*sj-1)/(sj-4)*f1, f2-3*(2-sj)/(sj-4)*f1])
        sub = p.projected(sub,vm)
        cached = H,sub,vm
        _CACHE.clear(); _CACHE[key] = cached
        reused = False
    else:
        reused = True
    H,sub,vm = cached
    lam2 = np.array([arb(float(v)) for v in ops.row_scale],dtype=object)
    lam = np.array([v.sqrt() for v in lam2],dtype=object)
    w.P_re,w.P_im,w.R_im = H[:,0]/lam[:,None],H[:,1]/lam[:,None],H[:,1]/lam2[:,None]
    w.f00,w.f11,w.chi = sub[0],sub[1],sub[2:]
    w.gram = {}
    for e,I in ((0,0),(1,1)):
        start = ops.index.index((I,e))*spec.M
        w.gram[e] = (H[start:start+spec.M,0],H[start:start+spec.M,1])
    w.precise,w.precise_basis = p,vm
    w.precision_info = {"source_dps":spec.operator_dps,"working_bits_with_guard":p.bits,
        'scattering_prescription':spec.scattering_prescription,
        "projection_seconds":time.monotonic()-started,"cache_reused":reused,"assembly_threads":8,
        "basis_scope":"chosen binary64 coordinate matrix; source rows and projection rebuilt in Arb"}


def verify(w, sol, tolerance=1e-8):
    """Apply the existing numerical tolerances to independently evaluated source functions."""
    from .arbaudit import ArbAudit
    from .formfactor import gram_block
    from .observables import phase_from_S,crossing_energy,modulus_peak,WAVES
    spec,M = w.spec,w.spec.M
    checker = ArbAudit(M,spec.L,max(256,w.precise.bits),prescription=spec.scattering_prescription,
                       source_dps=spec.operator_dps)
    y = [arb(v) for v in sol['y_text']]; av = arb_mat([[v] for v in y[:w.n_a]])
    cv = av if w.precise_basis is None else w.precise_basis*av
    c = [cv[i,0] for i in range(cv.nrows())]
    im = [y[w.n_a+e*M:w.n_a+(e+1)*M] for e in (0,1)] if spec.uv else None
    rh = [y[w.n_a+2*M+e*M:w.n_a+2*M+(e+1)*M] for e in (0,1)] if spec.uv else None
    audit = checker.audit(c,im,rh,chi_caliber=spec.chi_caliber if spec.chiral else None,
        eps_chi=spec.eps_chi,sr_caliber=spec.sr_caliber,eps_ff=spec.eps_ff,
        m_q=spec.m_q,ff_frozen_at_s0=spec.ff_frozen_at_s0,sr_free=spec.sr_free,eps_sr=spec.eps_sr,
        eps_ff_s0=spec.eps_ff_s0,eps_ff_p1=spec.eps_ff_p1)
    lower = lambda value: arb(value['lower'])
    midpoint = lambda value: float((arb(value['lower'])+arb(value['upper']))/2)
    violations = [-lower(r['slack']) for r in audit['unitarity']]
    active = [v/lower(r['h_squared']) for r,v in zip(audit['unitarity'],violations)
              if lower(r['h_squared']) > arb('1e-12')]
    worst, relative = max(violations),max(active,default=arb(-1))
    unitary = bool(worst <= arb('1e-8') and relative <= arb('1e-6'))
    checks = {'unitarity':unitary}
    out = {'unitarity':{'feasible':unitary,'max_absolute_violation_all':float(worst.upper()),
        'max_relative_violation_active':float(relative.upper()),'n_rows':len(violations)},
        'certified':False,'primal_feasible':False,'constraint_checks':checks,'source_audit':audit,
        'verification_tolerance':tolerance,'scope':'independent Arb source evaluation with unchanged numerical tolerances',
        'canonical_solution':'full SDPB y text and saved coordinate basis; float c is a convenience copy'}
    if spec.B is not None:
        raise ValueError('Precise mainline verification has no legacy density ball B')
    lay = w.ops.lay
    rho = [c[i] for i in list(range(lay.r1.start,lay.r1.stop))+list(range(lay.r2.start,lay.r2.stop))]
    linf = max(float(v.abs_upper()) for v in rho)
    out['density_norms'] = {'rho_linf':linf,
        'rho_l2':float(sum((v*v for v in rho),arb(0)).sqrt().upper()),
        'rho_l4':float(sum((v**4 for v in rho),arb(0)).root(4).upper()),
        'scope':'upper bounds of Arb enclosures of the double-density node values'}
    if spec.reg_norm is not None:
        used = out['density_norms']['rho_'+spec.reg_norm]
        checks['regulariser'] = used <= spec.reg_bound*(1+tolerance)
        out['regulariser'] = {'norm':spec.reg_norm,'bound':spec.reg_bound,'norm_used':used,'rho_linf':linf,
            'active_fraction':float(np.mean([float(v.abs_upper()) >= .99*spec.reg_bound for v in rho])) if spec.reg_norm=='linf' else None,
            'active':bool(used >= .99*spec.reg_bound)}
    if spec.chiral:
        scale = arb(str(spec.eps_chi))
        budget = scale*arb(str(tolerance)) if spec.chi_caliber=='chi-a' else scale**2*arb(str(2*tolerance+tolerance**2))
        checks['chiral'] = all(lower(r['slack']) >= -budget for r in audit['chiral']['constraints'])
        out['chiral'] = audit['chiral']
    if spec.uv:
        checks['rho_nonnegative'] = all(v >= -arb(str(tolerance)) for row in rh for v in row)
        if 'gram' in spec.uv_parts:
            minimum = float('inf')
            for e,wave in ((0,'S0'),(1,'P1')):
                for i in range(M):
                    fr = 1+sum((checker.K[i][j]*im[e][j] for j in range(M)),arb(0))
                    S = checker.last_primary[wave][i]
                    G = gram_block(complex(float(S.real.mid()),float(S.imag.mid())),
                        complex(float(fr.mid()),float(im[e][i].mid())),float(rh[e][i].mid()))
                    scale = np.sqrt(np.maximum(abs(np.diag(G)),1.))
                    minimum = min(minimum,float(np.linalg.eigvalsh(G/scale[:,None]/scale[None,:]).min()))
            checks['gram'] = minimum >= -tolerance
            out['gram'] = {'min_equilibrated_eigenvalue':minimum,'strict_counts':audit['gram_counts']}
        if 'fesr' in spec.uv_parts:
            checks['fesr'] = all(lower(r['slack']) >= -lower(r['tolerance'])*arb(str(tolerance))
                                 for r in audit['fesr'] if not r.get('free'))
            out['fesr'] = {'rows':audit['fesr']}
        if 'ff' in spec.uv_parts:
            checks['ff'] = all(lower(r['slack']) >= -lower(r['cap'])*arb(str(2*tolerance+tolerance**2)) for r in audit['form_factor'])
            out['form_factor'] = {'rows':audit['form_factor']}
    for name in ('f00','f11'):
        out[name+'_3'] = midpoint(audit['projection_at_3'][name])
    if w.functional is None:
        obj = sum((arb(str(d))*(arb(audit['projection_at_3'][name]['lower'])+arb(audit['projection_at_3'][name]['upper']))/2
                   for d,name in zip(w.direction,('f00','f11'))),arb(0))
    else:
        f = w.functional
        if f['kind']=='watson':
            from .watson import weight_of
            raw = arb(0)
            for wave in f['waves']:
                for i,((tr,ti),k) in enumerate(zip(f['targets'][wave],f['nodes'])):
                    S = checker.last_primary[wave][k]
                    raw = raw + arb(str(weight_of(f,wave,i)))*(arb(str(tr))*S.imag + arb(str(ti-1.0))*(1-S.real))
        elif f['kind'] in ('ImKH','ImS'):
            S = checker.last_primary[f['wave']][f['node']]
            raw = 1-S.real if f['kind']=='ImKH' else S.imag
        elif f['kind']=='SRmom':
            row = next(r for r in audit['fesr'] if r['wave']==f['wave'] and r['n']==f['node'])
            raw = (arb(row['moment_interval']['lower'])+arb(row['moment_interval']['upper']))/2
        else:
            raw = y[w.n_a+(2*M if f['kind']=='rho' else 0)+f['ell']*M+f['node']]
        out['functional'] = dict(f,value=float(raw.mid()),
            scope='Arb re-evaluation of the registered node functional from the same y')
        if f['kind']=='watson':
            from .watson import objective_value
            Sf = {wave:[complex(float(v.real.mid()),float(v.imag.mid())) for v in checker.last_primary[wave]] for wave in f['waves']}
            out['functional']['value_float_check'] = objective_value(f,Sf)
        obj = raw if f['sense']=='max' else -raw
    out['objective_recomputed'] = float(obj.mid())
    out['c_norm_inf'] = max(float(v.abs_upper()) for v in c)
    if w.fix_f00 is not None:
        out['section_residual'] = abs(out['f00_3']-w.fix_f00)
        checks['section'] = out['section_residual'] <= tolerance*max(1,abs(w.fix_f00))
    if w.face is not None:
        d0,d1 = w.face['direction']; floor = w.face['value']-w.face['margin']
        here = d0*out['f00_3']+d1*out['f11_3']
        out['face'] = dict(w.face,floor=floor,value_here=here,slack=here-floor)
        checks['face'] = here >= floor-tolerance*max(1.,abs(floor))
    out['primal_feasible'] = bool(all(checks.values()))
    out['observables'] = {}
    for wave in WAVES:
        S = [complex(float(v.real.mid()),float(v.imag.mid())) for v in checker.last_primary[wave]]
        ph = phase_from_S(w.ops.s,S)
        out['observables'][wave] = {k:ph[k].tolist() for k in ('E_GeV','eta','delta_deg')}
        out['observables'][wave].update(first_zero_node=ph['first_zero_node'],phase_convention=ph['phase_convention'],
            crossing_90_GeV=crossing_energy(ph),modulus_peak_GeV=modulus_peak(ph),
            min_eta_below_1p2GeV=float(np.min(ph['eta'][ph['E_GeV']<=1.2])))
    from .observables import arb_subthreshold_curves
    out['subthreshold'] = arb_subthreshold_curves(checker,c)
    return out
