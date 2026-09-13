"""Geometric selection of two completed IR section endpoints, without phase inputs."""
from __future__ import annotations

from collections import deque
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path

from .assembly import basis_identity
from .constraints import chiral_reference_point
from .dual import _stream_pmp

SCHEMA='ir-section-endpoint-selection-v1'
RULE='closest of these two endpoints by unweighted Euclidean distance; not global closest'
PRIOR='Upper-proxy phase diagnostics were already known before this supplement; selection uses no phase/linearity data and does not establish blind preregistration.'


def _sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def _finite(value,name):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
        raise ValueError('Finite numeric '+name+' required')
    return float(value)


def distance_order(distances,errors):
    """Conservative distance envelopes; a tie/overlap never chooses an endpoint."""
    if len(distances)!=2 or len(errors)!=2:raise ValueError('Two distances and error budgets required')
    if any(not math.isfinite(v) or v<0 for v in (*distances,*errors)):
        raise ValueError('Distances and budgets must be finite and nonnegative')
    intervals=[[max(0.,d-e),d+e] for d,e in zip(distances,errors)]
    selected=('upper' if intervals[0][1]<intervals[1][0] else
              'lower' if intervals[1][1]<intervals[0][0] else None)
    return {'selected':selected,'distance_intervals':intervals}


def _pmp(path,record):
    """Bind complete normalized constraint blocks and the actual opposite objectives."""
    pmp=path.parent/'pmp.json';meta=record['pmp'];digest=_sha(pmp)
    if digest!=meta.get('sha256'):raise ValueError('Saved PMP hash mismatch')
    blocks=_stream_pmp(pmp);header=next(blocks);n=meta.get('n_vars')
    if type(n) is not int or n<2:raise ValueError('Invalid PMP variable count')
    objective=[Fraction(str(v)) for v in header['objective']]
    norm=[Fraction(str(v)) for v in header['normalization']]
    if len(objective)!=n or norm!=[Fraction(1)]+[Fraction(0)]*(n-1):
        raise ValueError('PMP variable layout/normalization mismatch')
    h=hashlib.sha256();tail=deque(maxlen=2);count=0
    for block in blocks:
        data=json.dumps(block,sort_keys=True,separators=(',',':')).encode()
        h.update(len(data).to_bytes(8,'big'));h.update(data);tail.append(block);count+=1
    if _sha(pmp)!=digest:raise ValueError('PMP changed while validating its blocks')
    if count!=meta.get('n_blocks') or count<2:raise ValueError('PMP block count mismatch')
    def scalar(block):
        pol=block['polynomials']
        if len(pol)!=1 or len(pol[0])!=1 or len(pol[0][0])!=n or any(len(v)!=1 for v in pol[0][0]):
            raise ValueError('Expected two scalar fixed-section blocks at the end')
        return [Fraction(str(v[0])) for v in pol[0][0]]
    a,b=(scalar(block) for block in tail)
    if a[0]!=Fraction(str(record['fix_f00'])) or any(x+y for x,y in zip(a,b)):
        raise ValueError('Saved fixed-section blocks do not match the reported section')
    return {'sha256':digest,'constraint_blocks_sha256':h.hexdigest(),'n_blocks':count,
            'normalization':[str(v) for v in norm]},objective


def _endpoint(path,sign,receipt_dir):
    """Only source identity, projection, acceptance and numerical support data are read."""
    path=Path(path).resolve();raw=path.read_bytes();r=json.loads(raw)
    if not isinstance(r,dict):raise ValueError('Endpoint leaf must be a JSON object')
    spec=r.get('spec',{});v=r.get('verification',{});cv=r.get('convergence',{})
    if not all(isinstance(obj,dict) for obj in (spec,v,cv,r.get('pmp'))):
        raise ValueError('Complete endpoint model, acceptance and PMP metadata required')
    if not {'M','L','uv','chiral','B','reduce_basis','scattering_prescription'}.issubset(spec):
        raise ValueError('Incomplete endpoint ModelSpec')
    if (r.get('solver')!='SDPB' or r.get('status')!='numerically_accepted' or r.get('accepted') is not True
            or v.get('primal_feasible') is not True or cv.get('solver_optimal') is not True):
        raise ValueError('Completed numerically accepted endpoint leaves are required')
    if spec.get('uv') is not False or spec.get('chiral') is not True or spec.get('B') is not None:
        raise ValueError('Selection requires IR/no-B/chiral endpoints')
    if spec.get('scattering_prescription')!='sine-cardinal' or v.get('source_audit',{}).get('scattering_prescription')!='sine-cardinal':
        raise ValueError('A common sine-cardinal source acceptance is required')
    if any(v.get('constraint_checks',{}).get(k) is not True for k in ('unitarity','chiral','section')):
        raise ValueError('Endpoint lacks full numerical/source constraint acceptance')
    if r.get('direction')!=[0.,float(sign)]:raise ValueError('Endpoint direction must be [0,+1] or [0,-1] in its declared role')
    xr,yr=chiral_reference_point();x=_finite(v.get('f00_3'),'f00');y=_finite(v.get('f11_3'),'f11')
    fixed=_finite(r.get('fix_f00'),'fixed section')
    if abs(fixed-xr)>1e-12 or abs(x-xr)>1e-8:raise ValueError('Both endpoints must impose and satisfy x_ref')
    objective=_finite(v.get('objective_recomputed'),'source objective')
    lo=_finite(cv.get('candidate_lower'),'support lower');hi=_finite(cv.get('candidate_upper'),'support upper')
    gap=_finite(cv.get('absolute_gap'),'support gap');readback=_finite(cv.get('objective_readback_error'),'objective readback error')
    if hi<lo or gap<0 or readback<0 or abs(objective-sign*y)>1e-8:
        raise ValueError('Inconsistent numerical support/projection data')
    if abs(lo-objective)>1e-8*max(1.,abs(lo)) or abs((hi-lo)-gap)>1e-12*max(1.,abs(lo),abs(hi)):
        raise ValueError('Support bounds/gap do not match the accepted projection')
    identity=basis_identity(r,path)
    if not identity['verified']:raise ValueError('Actual endpoint basis is not authenticated')
    program,vector=_pmp(path,r)
    relative=lambda p:os.path.relpath(Path(p).resolve(),receipt_dir)
    inputs={'report':{'path':relative(path),'sha256':hashlib.sha256(raw).hexdigest()},
            'pmp':{'path':relative(path.parent/'pmp.json'),'sha256':program['sha256']}}
    if spec.get('reduce_basis'):
        local_basis=path.parent/'basis.npy'
        if not local_basis.is_file() or _sha(local_basis)!=identity['actual_sha256']:
            raise ValueError('A matching local basis.npy is required for a portable receipt')
        inputs['basis']={'path':relative(local_basis),'sha256':identity['actual_sha256']}
    roundoff=32*math.ulp(max(1.,abs(x),abs(y),abs(lo),abs(hi)))
    error=max(gap,hi-lo)+abs(objective-lo)+readback+abs(x-xr)+roundoff
    result={'inputs':inputs,'projection':[x,y],'direction':[0.,float(sign)],'fix_f00':fixed,
            'distance':math.hypot(x-xr,y-yr),'distance_error_budget':error,
            'numerical_support':{'lower':lo,'upper':hi,'gap':gap,'readback_error':readback},
            'basis_key':identity['key'],'pmp_identity':program}
    return result,spec,vector


def evaluate_pair(upper_report,lower_report,receipt_dir):
    """Compare authenticated input programs before making the geometric decision."""
    receipt_dir=Path(receipt_dir).resolve()
    upper,spec,uvec=_endpoint(upper_report,1,receipt_dir)
    lower,other,lvec=_endpoint(lower_report,-1,receipt_dir)
    if spec!=other:raise ValueError('Endpoint ModelSpecs differ')
    if upper['basis_key']!=lower['basis_key']:raise ValueError('Endpoint retained bases differ')
    for key in ('constraint_blocks_sha256','normalization','n_blocks'):
        if upper['pmp_identity'][key]!=lower['pmp_identity'][key]:raise ValueError('Endpoint PMP constraints/normalization differ')
    if len(uvec)!=len(lvec) or not any(uvec) or any(a+b for a,b in zip(uvec,lvec)):
        raise ValueError('Endpoint PMP objectives must be nonzero exact opposites')
    if lower['projection'][1]>upper['projection'][1]+upper['distance_error_budget']+lower['distance_error_budget']:
        raise ValueError('Endpoint order conflicts with the numerical support error budgets')
    decision=distance_order([upper['distance'],lower['distance']],
                            [upper['distance_error_budget'],lower['distance_error_budget']])
    for row,interval in zip((upper,lower),decision['distance_intervals']):row['distance_interval']=interval
    return {'schema':SCHEMA,'rule':RULE,'target':list(chiral_reference_point()),'metric':'unweighted Euclidean',
            'status':'selected' if decision['selected'] else 'ambiguous','selected':decision['selected'],
            'endpoints':{'upper':upper,'lower':lower},'spec':spec,
            'uncertainty_scope':'conservative numerical support-gap/readback/section-residual envelope; not a rigorous dual enclosure',
            'phase_data_used':False,'prior_information':PRIOR,'global_closest_claim':False}


def select_ir(upper_report,lower_report,outdir):
    dest=Path(outdir).resolve()
    if dest.exists() and any(dest.iterdir()):raise FileExistsError('Preserve prior selection; use a fresh directory')
    data=evaluate_pair(upper_report,lower_report,dest)
    receipt={'selection':data,'producer_sha256':_sha(__file__)}
    dest.mkdir(parents=True,exist_ok=True)
    from .sdpb import write_json
    write_json(dest/'report.json',receipt)
    return receipt


def verify_ir_selection(receipt_path,selected_report):
    """Recompute the entire selection from relative sources, not just its chosen name."""
    receipt_path=Path(receipt_path).resolve();receipt=json.loads(receipt_path.read_text());saved=receipt['selection']
    if saved.get('schema')!=SCHEMA:raise ValueError('Unknown IR selection receipt schema')
    sources={}
    for name in ('upper','lower'):
        value=saved['endpoints'][name]['inputs']['report']['path']
        if not isinstance(value,str) or Path(value).is_absolute():raise ValueError('Selection sources must use relative paths')
        sources[name]=(receipt_path.parent/value).resolve()
    current=evaluate_pair(sources['upper'],sources['lower'],receipt_path.parent)
    if current!=saved:raise ValueError('IR selection receipt or its source inputs changed')
    if current['status']!='selected':raise ValueError('Ambiguous receipt does not select an IR amplitude')
    if Path(selected_report).resolve()!=sources[current['selected']]:raise ValueError('IR report is not the receipt-selected endpoint')
    return {'receipt_sha256':_sha(receipt_path),'selection':current,'revalidated':True,
            'verification_producer_sha256':_sha(__file__)}
