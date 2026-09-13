"""Deliver registered IR/UV comparisons and explicitly partial ref diagnostics."""
from pathlib import Path
import hashlib
import json

import numpy as np

from . import claims, constraints as C
from .accuracy import rho_diagnostic
from .figures import ENERGY_AXIS_NOTE, REF, load_csv, phase_comparison
from .assembly import basis_identity

POINTS = ('tip','ref','mid')
IR_REFERENCES={wave:f'figure7_{wave.lower()}_phases.csv' for wave in ('S0','S2','P1')}


def load_chain(path):
    path = Path(path).resolve()
    doc = json.loads(path.read_text())
    chain = doc.get('reports',{})
    if set(chain) != set(POINTS):
        raise ValueError('The physical comparison requires complete tip/ref/mid results')
    _accepted_sources(chain)
    xr = C.chiral_reference_point()[0]
    xt = chain['tip']['verification']['f00_3']
    for name,r in chain.items():
        direction = [1.,0.] if name=='tip' else [0.,1.]
        target = None if name=='tip' else xr if name=='ref' else (xt+xr)/2
        if list(r['direction']) != direction or ((r.get('fix_f00') is None) != (target is None)):
            raise ValueError('Representative objective differs from the frozen rule')
        if target is not None and abs(r['fix_f00']-target)>1e-12:
            raise ValueError('Representative section differs from the frozen rule')
    return chain


def _accepted_sources(chain):
    """Shared numerical/source acceptance and analytic-family checks."""
    for name,r in chain.items():
        v = r.get('verification',{})
        if not (r.get('solver')=='SDPB' and r.get('status')=='numerically_accepted' and
                r.get('accepted') and v.get('primal_feasible') and
                r.get('convergence',{}).get('solver_optimal') and 'source_audit' in v):
            raise ValueError(f'{name} lacks full source acceptance and numerical support convergence')
        scheme = r['spec'].get('scattering_prescription','mixed-pv')
        if scheme != 'sine-cardinal' or v['source_audit'].get('scattering_prescription') != scheme:
            raise ValueError('This delivery route requires the declared sine-cardinal family with matching source audit')


def load_reference(path,direction=(0.,1.)):
    """Read exactly the registered ref from a completed leaf or reports.ref."""
    doc=json.loads(Path(path).resolve().read_text())
    if not isinstance(doc,dict):raise ValueError('A completed ref leaf or reports.ref is required')
    r=doc.get('reports',{}).get('ref') if isinstance(doc.get('reports'),dict) else doc
    if (not isinstance(r,dict) or r.get('solver')!='SDPB' or r.get('status')!='numerically_accepted'
            or r.get('accepted') is not True or not isinstance(r.get('spec'),dict)
            or not isinstance(r.get('verification'),dict)):
        raise ValueError('Reference-only comparison requires a completed accepted ref leaf')
    if (not isinstance(r.get('convergence'),dict) or r['convergence'].get('solver_optimal') is not True
            or r['verification'].get('primal_feasible') is not True
            or not isinstance(r['verification'].get('source_audit'),dict)):
        raise ValueError('Reference lacks full source acceptance and numerical support convergence')
    if not {'M','L','uv','chiral','B','uv_parts'}.issubset(r['spec']):
        raise ValueError('Reference lacks the complete source model specification')
    if not isinstance(r['spec']['uv_parts'],(list,tuple)):
        raise ValueError('Reference UV sector specification is malformed')
    _accepted_sources({'ref':r})
    v=r['verification'];checks=v.get('constraint_checks',{})
    required={'unitarity','chiral','section'}
    if r['spec'].get('uv'):required|={'rho_nonnegative','gram','fesr','ff'}
    if not isinstance(checks,dict) or any(checks.get(k) is not True for k in required):
        raise ValueError('Reference lacks full numerical/source constraint acceptance')
    xr=C.chiral_reference_point()[0]
    for value,tolerance in ((r.get('fix_f00'),1e-12),(v.get('f00_3'),1e-8)):
        if not isinstance(value,(int,float)) or not np.isfinite(value) or abs(value-xr)>tolerance:
            raise ValueError('Reference must impose and satisfy f00(3)=x_ref')
    if r.get('direction')!=list(direction):raise ValueError('Reference direction differs from its verified selection rule')
    obs=r.get('observables',{})
    for wave in ('S0','S2','P1'):
        try:energy,phase,eta=(np.asarray(obs[wave][k],float) for k in ('E_GeV','delta_deg','eta'))
        except (KeyError,TypeError,ValueError) as exc:raise ValueError('Complete native phase/eta observations required') from exc
        if (energy.ndim!=1 or len(energy)<2 or phase.shape!=energy.shape or eta.shape!=energy.shape
                or not all(np.all(np.isfinite(a)) for a in (energy,phase,eta))
                or not np.all(np.diff(energy)>0)):
            raise ValueError('Malformed native phase/eta observations')
    return {'ref':r}


def c4_source(observables):
    """Fig.7 comparison budgets, frozen independently of the new lower-point phases."""
    rows=[]
    for wave,name in IR_REFERENCES.items():
        path=Path(REF)/name
        if not path.exists():return {'verdict':'not run','evidence':'Fig.7 source references unavailable'}
        ref=load_csv(name);m=ref.get('group',np.array([]))=='ir_magenta'
        if not np.any(m):return {'verdict':'not run','evidence':'Fig.7 ir_magenta group unavailable'}
        x,y=ref['energy_gev'][m],ref['phase_deg'][m];order=np.argsort(x);x,y=x[order],y[order]
        obs=observables.get(wave,{})
        E=np.asarray(obs.get('E_GeV',[]),float);d=np.asarray(obs.get('delta_deg',[]),float)
        if (E.ndim!=1 or len(E)<2 or d.shape!=E.shape or not np.all(np.isfinite(E))
                or not np.all(np.isfinite(d)) or not np.all(np.diff(E)>0)):
            return {'verdict':'not run','evidence':'Complete finite native phase arrays required'}
        targets=x[x<=1.2]
        if not len(targets) or E[-1]<targets[-1]:
            return {'verdict':'not run','evidence':'Candidate does not cover the final published IR comparison energy'}
        lo,hi=max(E[0],x[0]),targets[-1];used=(E>=lo)&(E<=hi)
        public=(x>=lo)&(x<=hi)
        if not np.any(used) or not np.any(public):return {'verdict':'not run','evidence':'No common published energy interval'}
        rms=float(np.sqrt(np.mean((d[used]-np.interp(E[used],x,y))**2)))
        end=np.flatnonzero(public)[-1];delta=float(np.interp(x[end],E,d)-y[end])
        rows.append({'wave':wave,'n':int(used.sum()),'common_energy_GeV':[float(lo),float(hi)],
                     'rms_deg':rms,'endpoint_energy_GeV':float(x[end]),'endpoint_difference_deg':delta,
                     'source_csv':name,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                     'pass':rms<=10. and abs(delta)<=10.})
    return {'verdict':'pass' if all(r['pass'] for r in rows) else 'FAIL','rows':rows,
            'rms_budget_deg':10.,'endpoint_budget_deg':10.,'reference_group':'ir_magenta',
            'scope':'Source-based Fig.7 phase comparison; operational 10-degree budgets, not paper statistical uncertainty; no reference extrapolation to 1.2 GeV.',
            'amendment_scope':'Corrects the legacy C4 gate that rejects the published IR curve; frozen before inspection of new lower-point phases.'}


def c7_source(points):
    """Original C7 bands/RMS, evaluated at each published wave/color endpoint."""
    if not points or not set(points).issubset(POINTS):return {'verdict':'not run','evidence':'Known UV representative names required'}
    rows=[]
    for name in (p for p in POINTS if p in points):
        for wave,filename,band in (('S0','figure10_s0_phases.csv',(85.,110.)),('S2','figure10_s2_phases.csv',(-40.,-15.))):
            ref=load_csv(filename);m=ref['group']==claims.PAIRING[name]
            x,y=ref['energy_gev'][m],ref['phase_deg'][m];order=np.argsort(x);x,y=x[order],y[order]
            obs=points[name].get(wave,{})
            E=np.asarray(obs.get('E_GeV',[]),float);d=np.asarray(obs.get('delta_deg',[]),float)
            if (not len(x) or E.ndim!=1 or len(E)<2 or d.shape!=E.shape or not np.all(np.isfinite(E))
                    or not np.all(np.isfinite(d)) or not np.all(np.diff(E)>0)):
                return {'verdict':'not run','evidence':'Complete finite phase curves and published groups required'}
            targets=x[x<=1.2]
            if not len(targets) or E[-1]<targets[-1]:
                return {'verdict':'not run','evidence':'Candidate does not cover the final published UV comparison energy'}
            lo,hi=max(x[0],E[0]),targets[-1];used=(E>=lo)&(E<=hi);public=(x>=lo)&(x<=hi)
            if not np.any(used) or not np.any(public):return {'verdict':'not run','evidence':'No common published energy interval'}
            rms=float(np.sqrt(np.mean((d[used]-np.interp(E[used],x,y))**2)))
            end=np.flatnonzero(public)[-1];value=float(np.interp(x[end],E,d))
            rows.append({'point':name,'wave':wave,'n':int(used.sum()),'rms_deg':rms,
                         'endpoint_energy_GeV':float(x[end]),'endpoint_phase_deg':value,'endpoint_band_deg':list(band),
                         'source_group':claims.PAIRING[name],'source_csv':filename,
                         'source_sha256':hashlib.sha256((Path(REF)/filename).read_bytes()).hexdigest(),
                         'pass':band[0]<=value<=band[1] and rms<=10.})
    return {'verdict':('pass' if all(r['pass'] for r in rows) else 'FAIL') if set(points)==set(POINTS) else 'not run',
            'rows':rows,'rms_budget_deg':10.,'requires_all_three_points':True,
            'scope':'Original C7 endpoint bands and RMS budget; endpoint is each wave/color last published in-range energy, with no extrapolation to 1.196 GeV. No eta inferred from phase references.'}


def ir_linearity(record):
    """Describe deviation from a Weinberg line anchored to this saved f00(3)."""
    out={'complete_C3':False,'applied_threshold':None,
         'scope':'Descriptive RMS over saved subthreshold samples. C3 8% compares matching-epsilon Fig.5 curves, not a Weinberg-line tolerance; no near-Weinberg qualification is inferred.'}
    sub=record.get('subthreshold',{});s=np.asarray(sub.get('s',[]),float)
    x=record.get('verification',{}).get('f00_3')
    if s.ndim!=1 or not len(s) or not isinstance(x,(float,int)) or not np.isfinite(x) or x==0:
        return dict(out,status='not available',reason='Saved subthreshold curves and f00(3) required')
    used=np.isfinite(s)&(s>0)&(s<4)
    baseline=np.linspace(.05,3.95,40)
    common=np.all(np.isin(baseline,s))
    if common:used=np.isin(s,baseline)
    if not np.any(used):return dict(out,status='not available',reason='No saved samples in 0<s<4')
    lines={'S0':(2*s-1)*x/5,'P1':(s-4)*x/15,'S2':(2-s)*x/5};waves={}
    for wave,line in lines.items():
        values=np.asarray(sub.get(wave,[]),float)
        if values.shape!=s.shape or not np.all(np.isfinite(values[used])):
            return dict(out,status='not available',reason='Incomplete finite subthreshold samples')
        rms=float(np.sqrt(np.mean((values[used]-line[used])**2)))
        waves[wave]={'rms':rms,'rms_over_f00_3':rms/abs(x)}
    return dict(out,status='evaluated',f00_3=x,n=int(used.sum()),waves=waves,
                metric_grid='common legacy40 evaluation samples' if common else 'available saved grid; compare only matching grids',
                s_range=[float(s[used].min()),float(s[used].max())],used_for_selection=False)


def deliver_ir(source_report,selection_report,outdir):
    """Deliver one geometrically selected IR state without borrowing a UV amplitude."""
    from .ir_selection import verify_ir_selection
    from .figures import ir_profile_artifacts
    receipt=verify_ir_selection(selection_report,source_report)
    chosen=receipt['selection']['selected'];endpoint=receipt['selection']['endpoints'][chosen]
    state=load_reference(source_report,endpoint['direction'])['ref']
    if state['spec']['uv']:raise ValueError('Selected IR figures require an IR amplitude')
    dest=Path(outdir)
    if dest.exists() and any(dest.iterdir()):raise FileExistsError('Use a fresh IR figure directory')
    result={'kind':'selected_IR_delivery','source_report':str(Path(source_report).resolve()),
            'source_sha256':hashlib.sha256(Path(source_report).read_bytes()).hexdigest(),
            'selection':receipt,'selected':chosen,'spec':state['spec'],
            'C4_source':c4_source(state['observables']),'linearity_diagnostic':ir_linearity(state),
            'complete_C3':False,'complete_UV_comparison':False,
            'scope':'One selected complete IR amplitude. Near-global geometry, all Fig5 epsilon curves and UV claims are separate.',
            'energy_axis_note':ENERGY_AXIS_NOTE}
    dest.mkdir(parents=True,exist_ok=True)
    _plots({'ref':state['observables']},{},dest,paper_style=True,chosen=chosen)
    ir_profile_artifacts(dest,state)
    (dest/'report.json').write_text(json.dumps(result,indent=2))
    return result


def deliver(ir_path, uv_path, outdir,reference_only=False,ir_selection_report=None):
    """Require the complete common amplitude contract before comparing its phases."""
    loader=load_reference if reference_only else load_chain
    selection=None
    if ir_selection_report:
        from .ir_selection import verify_ir_selection
        selection=verify_ir_selection(ir_selection_report,ir_path)
        chosen=selection['selection']['selected'];endpoint=selection['selection']['endpoints'][chosen]
        ir=load_reference(ir_path,endpoint['direction'])
    else:ir=loader(ir_path)
    uv_selection=None
    if json.loads(Path(uv_path).read_text()).get('schema')=='uv-nearby-selection-v1':
        if reference_only:raise ValueError('A complete UV nearby receipt is not a reference-only comparison')
        from .uv_selection import verify_uv_selection
        uv_selection=verify_uv_selection(uv_path);uv=uv_selection['reports']
    else:uv=loader(uv_path)
    contract = None
    for is_uv,chain in ((False,ir),(True,uv)):
        sector_contract = None
        for r in chain.values():
            spec = r['spec']
            if bool(spec['uv']) != is_uv or not spec['chiral'] or spec['B'] is not None:
                raise ValueError('Expected the same no-B chiral problem, with and without full UV')
            if is_uv and set(spec['uv_parts']) != {'gram','fesr','ff'}:
                raise ValueError('The gauge claim requires all UV constraints')
            common = {k:v for k,v in spec.items() if k not in ('uv','tag','disk_mask')}
            if sector_contract is not None and common != sector_contract:
                raise ValueError('Representatives mix different contracts within one sector')
            sector_contract = common
            common = {k:v for k,v in common.items() if k not in
                      ('uv_parts','sr_caliber','eps_ff','m_q','ff_frozen_at_s0')}
            if contract is None:
                contract = common
            if common != contract:
                raise ValueError('IR/UV results mix different physical or numerical contracts')
    ip = {name:r['observables'] for name,r in ir.items()}
    up = {name:r['observables'] for name,r in uv.items()}
    bases={sector:{name:(basis_identity(r,uv_selection['source_paths'][name]) if sector=='UV' and uv_selection else
                        basis_identity(r,path,name)) for name,r in chain.items()}
           for sector,chain,path in (('IR',ir,ir_path),('UV',uv,uv_path))}
    identities=[v for sector in bases.values() for v in sector.values()]
    same_basis=all(v['verified'] for v in identities) and len({v['key'] for v in identities})==1
    pairs={name:(ir['ref'],row) for name,row in uv.items()} if selection else {name:(ir[name],uv[name]) for name in ir.keys() & uv.keys()}
    matches={name:a['direction']==b['direction'] and a.get('fix_f00')==b.get('fix_f00') for name,(a,b) in pairs.items()}
    matched=bool(matches) and all(matches.values());causal=same_basis and matched
    legacy=claims.c4(ip);source_c4=c4_source(ir['ref']['observables']) if selection else None
    legacy_c7=claims.c7(up);source_c7=c7_source(up) if selection else None
    out = {'contract':contract,'uv_contract':sector_contract,
           'basis_identity':bases,'same_basis':same_basis,'causal_comparison_supported':causal,
           'matched_objectives':matches,'matched_objective_causal_comparison_supported':causal,
           'basis_comparison_scope':('same retained coordinate space; selection objective/section identity is reported separately' if same_basis else
               'retained coordinate spaces differ or are unverified; comparison is diagnostic, UV effects are not isolated'),
           'comparison_mode':('paper-style-reference-only' if reference_only else 'paper-style-selected-ir') if selection else
                             'reference-only' if reference_only else 'registered-three-representative',
           'incomplete_registered_chain':bool(reference_only),
           'available_points':{'IR':list(ir),'UV':list(uv)},
           'required_points':{'IR':['ref'] if selection else list(POINTS),'UV':list(POINTS)},
           'missing_registered_points':[name for name in POINTS if (not selection and name not in ir) or name not in uv],
           'C4':source_c4 if selection else legacy,'C4_legacy':legacy,'C6':claims.c6(up),
           'C7':source_c7 if selection else legacy_c7,'C7_legacy':legacy_c7,
           'rho':{sector:{name:rho_diagnostic(r['observables']['P1']) for name,r in chain.items()}
                  for sector,chain in (('IR',ir),('UV',uv))},
           'scope':('Reference-only finite-node IR/UV diagnostic; registered tip/ref/mid chain incomplete; C6/C7 not run'
                    if reference_only else 'Registered finite-node IR/UV contrast, not completion of all Fig.3–11 claims'),
           'selection_scope':selection['selection']['rule'] if selection else 'registered section proxy, not a proved closest point to physical f_pi',
           'energy_axis_note':ENERGY_AXIS_NOTE,
           'nominal_resolution':contract['M']==50 and contract['L']==10,
           'input_sha256':{str(Path(p).resolve()):hashlib.sha256(Path(p).read_bytes()).hexdigest()
                           for p in (ir_path,uv_path)}}
    out['main_contrast_pass'] = same_basis and not reference_only and out['nominal_resolution'] and all(out[k]['verdict']=='pass' for k in ('C4','C6','C7'))
    if selection:
        out.update(ir_selection=selection,C4_source=source_c4,C7_source=source_c7,selected_ir_diagnostics=ir_linearity(ir['ref']),
                   paper_state_comparison_complete=not reference_only,paper_state_comparison_pass=out['main_contrast_pass'],
                   scope='Paper-style selected IR versus UV states; matched-objective causal attribution, sampled linearity and full C3 completion remain separate.')
        out['input_sha256'][str(Path(ir_selection_report).resolve())]=selection['receipt_sha256']
        out['C6']['eta_scope']='The eta>=0.9 gate remains an existing native-solution diagnostic; published phase curves do not supply numerical eta evidence.'
    if uv_selection:
        out.update(uv_selection=uv_selection['metadata'],uv_role_mapping=uv_selection['role_mapping'],
                   legacy_registered_chain_complete=False,incomplete_registered_chain=True,
                   comparison_mode='paper-nearby-uv-selected-ir' if selection else 'paper-nearby-uv',
                   scope='UV states use validated frozen nearby sections; the legacy midpoint rule is not claimed. Physical comparison, proximity and full Fig.3-11 completion remain separate.',
                   uv_selection_scope='The internal mid slot denotes the frozen near section; it is not the legacy arithmetic midpoint or an authenticated author coordinate.')
        out['C6']['role_mapping']=uv_selection['role_mapping'];out['C7']['role_mapping']=uv_selection['role_mapping']
    dest = Path(outdir)
    if dest.exists() and any(dest.iterdir()):
        raise FileExistsError('Use a fresh contrast directory to preserve prior evidence')
    dest.mkdir(parents=True,exist_ok=True)
    (dest/'report.json').write_text(json.dumps(out,indent=2))
    _plots(ip,up,dest,same_basis,paper_style=bool(selection),chosen=None if not selection else chosen,matched=matched,
           uv_labels=None if uv_selection is None else uv_selection['role_mapping'])
    return out


def _plots(ir,uv,dest,same_basis=True,paper_style=False,chosen=None,matched=True,uv_labels=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .figures import load_csv
    colors = {'tip':'#b2182b','mid':'#ef8a62','ref':'#b35886'}
    names = {'S0':'figure10_s0_phases.csv','S2':'figure10_s2_phases.csv','P1':'figure9_p1_phases.csv'}
    for field,suffix,ylabel in (('delta_deg','phases','phase shift (degrees)'),('eta','inelasticity',r'$|S|$')):
        fig,axes = plt.subplots(1,3,figsize=(12,3.6),sharex=True)
        for ax,wave in zip(axes,('S0','S2','P1')):
            for name in (point for point in POINTS if point in uv):
                obs = uv[name][wave]; E=np.array(obs['E_GeV']); mask=E<=1.2
                ax.plot(E[mask],np.array(obs[field])[mask],color=colors[name],label='UV '+(uv_labels or {}).get(name,name))
            obs=ir['ref'][wave];E=np.array(obs['E_GeV']);mask=E<=1.2
            ax.plot(E[mask],np.array(obs[field])[mask],'k--',label='IR '+('selected '+chosen if chosen else 'ref'))
            if field=='delta_deg':
                if uv:ref=load_csv(names[wave])
                for name in (point for point in POINTS if point in uv):
                    mask=ref['group']==claims.PAIRING[name]
                    ax.plot(ref['energy_gev'][mask],ref['phase_deg'][mask],':',color=colors[name],lw=.9)
                if paper_style and (Path(REF)/IR_REFERENCES[wave]).exists():
                    ref=load_csv(IR_REFERENCES[wave]);mask=ref['group']=='ir_magenta'
                    ax.plot(ref['energy_gev'][mask],ref['phase_deg'][mask],'-.',color='.35',lw=.9,label='Fig.7 IR')
            else:
                ax.axhline(1,color='.6',lw=.7)
            ax.set(title=wave,xlabel='energy (GeV)',xlim=(.28,1.2),ylabel=ylabel)
            ax.grid(alpha=.15)
        axes[0].legend(frameon=False,fontsize=8)
        title=('Reference-only IR/UV diagnostic; registered chain incomplete' if set(uv)!=set(POINTS)
               else 'Same-input IR/UV comparison')
        if not same_basis:title='IR/UV diagnostic: retained bases differ or are unverified'
        elif paper_style:title='Paper-style IR/UV state comparison'+(' (UV ref only)' if set(uv)!=set(POINTS) else '')
        if paper_style and not matched:title+='; different selection objectives'
        if uv_labels:title+='; frozen nearby UV sections'
        if not uv:title='Selected IR state; no UV comparison'
        note=('; dotted UV / dash-dot IR: published phases' if uv else '; dash-dot: published IR phase') if field=='delta_deg' else '; horizontal line: unitarity limit'
        fig.suptitle(title+note,fontsize=10)
        fig.tight_layout()
        fig.savefig(dest/f'{suffix}.pdf');fig.savefig(dest/f'{suffix}.png',dpi=180)
        plt.close(fig)
