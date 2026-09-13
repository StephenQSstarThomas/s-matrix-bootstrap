"""Freeze nearby UV sections and authenticate their completed saved amplitudes."""
from __future__ import annotations

from collections import deque
from dataclasses import fields
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path

from .amplitude_export import _source
from .assembly import basis_identity
from .constraints import chiral_reference_point
from .dual import _stream_pmp
from .ir_selection import _finite, _sha
from .spec import ModelSpec

PLAN='uv-nearby-plan-v1'
RECEIPT='uv-nearby-selection-v1'
RULE='x_near=x_ref+min(0.01*x_ref,(x_tip-x_ref)/2); maximize f11 on ref/near sections'
ROLE_MAP={'tip':'tip','ref':'ref','mid':'near'}
KERNELS=('precision.py','sine.py','projector.py','legendreq.py','hilbert.py','grid.py')
PREF={'base':'1','constant':'1','poles':[]}


def _fresh(outdir):
    dest=Path(outdir).resolve()
    if dest.exists() and any(dest.iterdir()):raise FileExistsError('Preserve prior evidence; use a fresh directory')
    return dest


def _write(dest,data):
    dest.mkdir(parents=True,exist_ok=True)
    with (dest/'report.json').open('x') as stream:json.dump(data,stream,indent=2,allow_nan=False)
    return data


def _entry(path,dest):
    return {'path':os.path.relpath(Path(path).resolve(),dest),'sha256':_sha(path)}


def _resolve(entry,root):
    if not isinstance(entry,dict) or set(entry)!={'path','sha256'}:
        raise ValueError('A relative source path and complete hash are required')
    value=entry['path']
    if not isinstance(value,str) or Path(value).is_absolute():raise ValueError('Sources must use relative paths')
    path=(root/value).resolve()
    if _sha(path)!=entry['sha256']:raise ValueError('Selection source hash changed')
    return path


def _vector(values,n):
    if not isinstance(values,list) or len(values)!=n:raise ValueError('PMP vector layout mismatch')
    return [Fraction(str(v)) for v in values]


def _program(path,r,plane,role,fixed):
    """Verify objective and both section blocks before excluding them from identity."""
    pmp=path.parent/'pmp.json';digest=_sha(pmp);n=r['pmp']['n_vars']
    if digest!=r['pmp']['sha256']:raise ValueError('PMP hash mismatch')
    blocks=_stream_pmp(pmp);header=next(blocks)
    norm=_vector(header['normalization'],n)
    if norm!=[1]+[0]*(n-1):raise ValueError('PMP normalization mismatch')
    obj=_vector(header['objective'],n)
    if obj!=_vector(plane['f00' if role=='tip' else 'f11'],n):
        raise ValueError('Actual PMP objective is not the required source partial wave')
    tail=deque();h=hashlib.sha256();count=0;base_count=0
    def accept(block):
        data=json.dumps(block,sort_keys=True,separators=(',',':')).encode()
        h.update(len(data).to_bytes(8,'big'));h.update(data)
    for block in blocks:
        count+=1;tail.append(block)
        if len(tail)>(2 if fixed is not None else 0):accept(tail.popleft());base_count+=1
    if count!=r['pmp']['n_blocks'] or base_count<1:raise ValueError('PMP block count mismatch')
    if fixed is not None:
        expected=[Fraction(str(fixed))*a-b for a,b in zip(norm,_vector(plane['f00'],n))]
        if len(tail)!=2:raise ValueError('Both fixed-section blocks are required')
        for sign,block in zip((1,-1),tail):
            pol=block.get('polynomials',[])
            if (block.get('DampedRational')!=PREF or len(pol)!=1 or len(pol[0])!=1
                    or len(pol[0][0])!=n or any(not isinstance(v,list) or len(v)!=1 for v in pol[0][0])):
                raise ValueError('Malformed fixed-section block or prefactor')
            if _vector([v[0] for v in pol[0][0]],n)!=[sign*v for v in expected]:
                raise ValueError('PMP fixed section is not the frozen source f00 row')
    if _sha(pmp)!=digest:raise ValueError('PMP changed during validation')
    return {'sha256':digest,'physical_blocks_sha256':h.hexdigest(),'physical_block_count':base_count,
            'normalization':[str(v) for v in norm]},obj


def _planes(path):
    """Rebuild only two plane rows in the saved basis; no constraints or SVD."""
    from .pmp import Pmp
    w=Pmp.from_saved(path,(1.,0.),None)
    return {key:[str(Fraction(str(v))) for v in w._num_vec(w._row(a=row))]
            for key,row in (('f00',w.f00),('f11',w.f11))}


def _leaf(path,dest,role,fixed=None,plane=None):
    path=Path(path).resolve()
    try:r,layout,basis,tokens,hashes=_source(path)
    except (KeyError,TypeError,IndexError,AttributeError,OSError) as exc:
        raise ValueError('Incomplete saved amplitude source') from exc
    spec=r['spec'];v=r['verification'];cv=r.get('convergence',{})
    if set(spec)!={f.name for f in fields(ModelSpec)}:raise ValueError('Complete ModelSpec required')
    if (spec['uv'] is not True or spec['chiral'] is not True or spec['B'] is not None
            or set(spec['uv_parts'])!={'gram','fesr','ff'} or spec['scattering_prescription']!='sine-cardinal'):
        raise ValueError('Selection requires a full UV/no-B/chiral sine source')
    required=('unitarity','chiral','rho_nonnegative','gram','fesr','ff')+(('section',) if fixed is not None else ())
    if (r['status']!='numerically_accepted' or r.get('accepted') is not True
            or v.get('primal_feasible') is not True or cv.get('solver_optimal') is not True
            or v.get('source_audit',{}).get('scattering_prescription')!='sine-cardinal'
            or any(v.get('constraint_checks',{}).get(k) is not True for k in required)):
        raise ValueError('Complete numerical and source acceptance is required')
    kernels={name:r.get('source_sha256',{}).get('src/smatrix_bootstrap/sdp/'+name) for name in KERNELS}
    if any(value!=_sha(Path(__file__).with_name(name)) for name,value in kernels.items()):
        raise ValueError('Source partial-wave kernel changed')
    direction=[1.,0.] if role=='tip' else [0.,1.]
    if r.get('direction')!=direction or ((r.get('fix_f00') is None)!=(fixed is None)):
        raise ValueError('Source role direction or section differs from the plan')
    if fixed is not None and _finite(r['fix_f00'],'section')!=fixed:
        raise ValueError('Source fixed section differs from the frozen value')
    identity=basis_identity(r,path)
    if not identity['verified']:raise ValueError('Actual saved basis is not authenticated')
    if spec['reduce_basis'] and identity['path']!=str(path.parent/'basis.npy'):
        raise ValueError('A local basis is required for portable selection')
    plane=_planes(path) if plane is None else plane
    program,objective=_program(path,r,plane,role,fixed)
    x=_finite(v.get('f00_3'),'f00');y=_finite(v.get('f11_3'),'f11')
    values=[Fraction(1)]+[Fraction(t) for t in tokens]
    if any(abs(float(sum(a*b for a,b in zip(_vector(plane[k],len(values)),values)))-q)>1e-8*max(1.,abs(q))
           for k,q in (('f00',x),('f11',y))):raise ValueError('Saved projection does not match complete y and source plane')
    source_obj=_finite(v.get('objective_recomputed'),'source objective');target=x if role=='tip' else y
    lo=_finite(cv.get('candidate_lower'),'support lower');hi=_finite(cv.get('candidate_upper'),'support upper')
    gap=_finite(cv.get('absolute_gap'),'gap');readback=_finite(cv.get('objective_readback_error'),'readback error')
    actual=float(sum(a*b for a,b in zip(objective,values)))
    if (hi<lo or gap<0 or readback<0 or abs(source_obj-target)>1e-8*max(1.,abs(target))
            or abs(lo-source_obj)>1e-8*max(1.,abs(lo)) or abs(actual-source_obj)>1e-8*max(1.,abs(actual))
            or abs(hi-lo-gap)>1e-12*max(1.,abs(lo),abs(hi))):
        raise ValueError('Inconsistent saved objective/support/projection values')
    if fixed is not None and abs(x-fixed)>1e-8:raise ValueError('Source does not satisfy the frozen section')
    files={'report':path,'pmp':path.parent/'pmp.json','y':path.parent/'out/y.txt',
           'solution':path.parent/'solution.npz','process':path.parent/'sdpb_process.json',
           'mask':path.parent/'disk_mask.npy'}
    if spec['reduce_basis']:files['basis']=path.parent/'basis.npy'
    inputs={key:_entry(value,dest) for key,value in files.items()}
    if any(_sha(p)!=digest for p,digest in hashes.items()):raise ValueError('Source changed during validation')
    xr,yr=chiral_reference_point();distance=math.hypot(x-xr,y-yr)
    error=gap+readback+abs(actual-source_obj)+(abs(x-fixed) if fixed is not None else 0)+32*math.ulp(max(1.,abs(x),abs(y)))
    meta={'inputs':inputs,'role':role,'direction':direction,'fix_f00':fixed,'projection':[x,y],
          'distance':distance,'distance_over_xref':distance/xr,'distance_error_budget':error,
          'distance_interval':[max(0.,distance-error),distance+error],
          'basis_key':identity['key'],'pmp_identity':program,'source_kernel_sha256':kernels}
    return meta,r,plane


def _plan_data(tip_report,dest):
    tip,r,plane=_leaf(tip_report,dest,'tip');xr,yr=chiral_reference_point();xt=tip['projection'][0]
    if xt<=xr:raise ValueError('Accepted UV tip must lie strictly to the right of x_ref')
    exact=Fraction(str(xr))+min(Fraction(str(xr))/100,(Fraction(str(xt))-Fraction(str(xr)))/2)
    near=float(exact)
    if not xr<near<xt:raise ValueError('No representable distinct interior nearby section')
    return {'rule':RULE,'policy':'uv-near-ref-v1','target':[xr,yr],'x_tip':xt,'x_ref':xr,'x_near':near,
            'x_near_exact_rational':str(exact),'actual_section_tokens':{'ref':repr(xr),'near':repr(near)},
            'spec':r['spec'],'tip':tip,'source_plane':plane,'phase_data_used':False,
            'two_dimensional_radius':None,'legacy_registered_chain_complete':False,
            'scope':'Declared nearby x sections; full geometric distances are diagnostic, not a near-black qualification or author point recovery.'},r


def plan_uv(tip_report,outdir):
    """Freeze nearby x values from one completed tip, before new section supports."""
    dest=_fresh(outdir);data,_=_plan_data(tip_report,dest)
    return _write(dest,{'schema':PLAN,'created_utc':datetime.now(timezone.utc).isoformat(),
                        'plan':data,'producer_sha256':_sha(__file__)})


def _verify_plan(path):
    path=Path(path).resolve();doc=json.loads(path.read_text())
    if set(doc)!={'schema','created_utc','plan','producer_sha256'} or doc['schema']!=PLAN:
        raise ValueError('Unknown or malformed UV plan')
    if doc['producer_sha256']!=_sha(__file__):raise ValueError('UV plan producer identity changed')
    datetime.fromisoformat(doc['created_utc'])
    tip_path=_resolve(doc['plan']['tip']['inputs']['report'],path.parent)
    data,r=_plan_data(tip_path,path.parent)
    if data!=doc['plan']:raise ValueError('UV plan or its complete tip identity changed')
    return data,r,tip_path


def _evaluate(plan_report,ref_report,near_report,dest):
    plan_path=Path(plan_report).resolve();plan,r0,tip_path=_verify_plan(plan_path)
    paths={'tip':tip_path,'ref':Path(ref_report).resolve(),'near':Path(near_report).resolve()}
    if len(set(paths.values()))!=3:raise ValueError('Three distinct role sources are required')
    tip=dict(plan['tip']);tip['inputs']={key:_entry(_resolve(value,plan_path.parent),dest) for key,value in tip['inputs'].items()}
    points={'tip':tip};reports={'tip':r0}
    for role in ('ref','near'):
        meta,r,_=_leaf(paths[role],dest,role,plan['x_'+role],plan['source_plane'])
        if r['spec']!=plan['spec'] or meta['basis_key']!=tip['basis_key']:
            raise ValueError('UV representatives differ in full ModelSpec or actual basis')
        for key in ('physical_blocks_sha256','physical_block_count','normalization'):
            if meta['pmp_identity'][key]!=tip['pmp_identity'][key]:raise ValueError('UV physical PMP blocks/normalization differ')
        if meta['inputs']['mask']['sha256']!=tip['inputs']['mask']['sha256']:
            raise ValueError('UV effective disk masks differ')
        points[role]=meta;reports['mid' if role=='near' else role]=r
    data={'plan':_entry(plan_path,dest),'policy':plan['policy'],'rule':RULE,'target':plan['target'],
          'points':points,'role_mapping':ROLE_MAP,'spec':plan['spec'],'phase_data_used':False,
          'two_dimensional_radius':None,'legacy_registered_chain_complete':False,
          'scope':plan['scope'],'distance_uncertainty_scope':'Numerical gap/readback envelope, not a rigorous dual enclosure.',
          'chronology_scope':'Plan has an explicit creation time; source hashes do not by themselves prove blind preregistration.'}
    return data,reports,paths


def select_uv(plan_report,ref_report,near_report,outdir):
    """Authenticate all three roles without inspecting their phase or linearity data."""
    dest=_fresh(outdir);data,_,_=_evaluate(plan_report,ref_report,near_report,dest)
    return _write(dest,{'schema':RECEIPT,'selection':data,'producer_sha256':_sha(__file__)})


def verify_uv_selection(receipt):
    """Revalidate relative sources and return report data plus explicit role metadata."""
    path=Path(receipt).resolve();doc=json.loads(path.read_text())
    if set(doc)!={'schema','selection','producer_sha256'} or doc['schema']!=RECEIPT:
        raise ValueError('Unknown or malformed UV selection receipt')
    if doc['producer_sha256']!=_sha(__file__):raise ValueError('UV receipt producer identity changed')
    saved=doc['selection'];plan=_resolve(saved['plan'],path.parent)
    ref=_resolve(saved['points']['ref']['inputs']['report'],path.parent)
    near=_resolve(saved['points']['near']['inputs']['report'],path.parent)
    data,reports,paths=_evaluate(plan,ref,near,path.parent)
    if data!=saved:raise ValueError('UV receipt or any selected source changed')
    return {'reports':reports,'role_mapping':dict(ROLE_MAP),
            'source_paths':{('mid' if k=='near' else k):str(v) for k,v in paths.items()},
            'metadata':dict(data,revalidated=True,receipt_sha256=_sha(path),verification_producer_sha256=_sha(__file__))}
