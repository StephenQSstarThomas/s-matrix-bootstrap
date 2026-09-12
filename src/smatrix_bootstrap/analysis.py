"""Predeclared representatives and scientific curve comparisons."""
from flint import arb
from pathlib import Path
from scipy.spatial import ConvexHull, QhullError
import csv
import matplotlib.pyplot as plt
import numpy as np
import shutil
from . import read_json, write_json
XREF = 5/(16*np.pi**2*(92/140)**2)
COLORS = dict(tip='red',mid='#ff5555',ref='#ff9999')
WAVES = ['S0','S2','P1']

def published_inner_vertices(points):
    positive=[(i,p) for i,p in enumerate(points) if p[0]>=0]
    negative=[(i,p) for i,p in enumerate(points) if p[0]<0]
    vertices=[[arb(0),arb(0)]]+[p for i,p in positive]
    mixtures=[dict(zero_amplitude=True)]+[dict(source_point=i) for i,p in positive]
    for i,p in positive:
        if not p[0]>0:continue
        for j,n in negative:
            left,right=p[0].lower(),n[0].lower()
            weight=float(np.nextafter(float((-right/(left-right)).upper()),np.inf))
            if not 0<weight<1:continue
            t=arb(weight);v=[t*a+(1-t)*b for a,b in zip(p,n)]
            if v[0]>=0:
                vertices.append(v);mixtures.append(dict(source_points=[i,j],positive_weight=weight))
    return vertices,mixtures

def energy_rows(job):
    from .kernels import source_rows
    from . import decode_real_ball
    if len(job) not in (7,8,9,10):raise ValueError('An energy job has 7 to 10 fields')
    M,L,bits,order,subtracted,node,encoded=job[:7];selected=job[7] if len(job)>7 else None
    source=source_rows(M,L,bits,order,subtracted,prescription=job[8] if len(job)>=9 else 'pv-midpoint',registry=job[9] if len(job)==10 else None)
    s=source.x[node] if node is not None else decode_real_ball(encoded)
    real,imag=[],[]
    waves=source.waves if selected is None else selected;spins={ell for I,ell in waves}
    for I,ell in waves:
        row=source.row(s,ell,I,node=node,spins=spins)
        real.append([float(v.real.mid()) for v in row]);imag.append([float(v.imag.mid()) for v in row])
    return real,imag

from .gauge import initialize

def joint_support_points(requests,isig=None,require_center=False):
    from .run import _path, _signature
    points=[];rejected=[];reference_signature=None;seen=set()
    for request in requests:
        path=_path(request);path=path/'report.json' if path.is_dir() else path
        try:
            if path in seen:continue
            seen.add(path);r=read_json(path);b=r['outer'];d=np.asarray(r.get('support_direction',r['parameters']['direction']),float);xy=np.asarray(r['targets'],float)
            if r.get('objective_kind','projection')!='projection':raise ValueError('A Watson functional is not a projection support')
            if r.get('status') in ('inconclusive','running') or r.get('timeout') or r.get('exit_code',0)!=0:raise ValueError('Support run did not finish successfully')
            if r['mode']!='gauge' or not r.get('joint_feasible') or not b.get('primal_feasible') or not b.get('joint_primal_feasible') or not r.get('all_original_constraints_checked'):
                raise ValueError('No verified joint point')
            if not r['solver'].get('all_amplitude_directions_retained'):raise ValueError('Full amplitude support search required')
            if require_center and r['solver'].get('representative_center_converged') is not True:raise ValueError('Converged representative center required')
            if d.shape!=(2,) or xy.shape!=(2,) or not np.any(d) or not np.isfinite(np.r_[d,xy,b['lower'],b['upper']]).all() or b['lower']>b['upper'] or not b.get('enclosure'):
                raise ValueError('Finite nonzero support direction and valid bounds required')
            sig=_signature(r)
            if isig is not None:
                if sig['chiral_norm']!=isig.get('chiral_norm','separate-l2'):raise ValueError('IR/UV chiral norm mismatch')
                for key in ('M','L','density_limit','infinity','unitarity_scope','prescription','preparation','scattering_samples'):
                    if sig[key]!=isig[key]:raise ValueError('IR/UV mismatch: '+key)
            if sig['chiral_tolerance']!=.002 or reference_signature is not None and sig!=reference_signature:raise ValueError('Joint input mismatch')
            for name in ('coefficients.json','current.json','joint.npz'):
                if not (path.parent/name).is_file():raise ValueError('Missing '+name)
            reference_signature=sig;points.append(dict(target=xy.tolist(),direction=d.tolist(),lower=b['lower'],upper=b['upper'],
                gap=b['upper']-b['lower'],upper_enclosure=b['enclosure'],report=str(path),coefficients=str(path.parent/'coefficients.json'),
                current=str(path.parent/'current.json'),joint=str(path.parent/'joint.npz'),joint_feasible=True,
                support_optimality_certified=r.get('support_optimality_certified',False),
                representative_center_converged=r['solver'].get('representative_center_converged') is True))
        except (OSError,ValueError,KeyError,TypeError) as error:rejected.append(dict(report=str(path),reason=str(error)))
    if not points:raise ValueError('No qualified full-space joint supports: '+str(rejected))
    return reference_signature,points,rejected


def gauge_regions(args):
    from . import gauge_region_report
    from .run import _path
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    baseline=read_json(args.region_summary);ir=next(r for r in baseline['regions'] if r['mode']=='chiral' and r['epsilon']==.002)
    reference_signature,points,rejected=joint_support_points(args.support_runs,baseline['model_signature'])
    iin=_hull(ir['geometry']['inner_vertices']);iout=_hull(ir['geometry']['outer_vertices'])
    if not ir['geometry']['outer_closed']:raise ValueError('Closed valid IR outer region required')
    inner=_clip(_hull([p['target'] for p in points]),[-1,0],0);outer=iout.copy()
    for p in points:outer=_clip(outer,p['direction'],p['upper'])
    if not len(outer):raise ArithmeticError('Joint support lines contradict the retained IR outer region')
    xs=[]
    if len(inner)>=3:
        lo=max(0.,inner[:,0].min(),iin[:,0].min());hi=min(inner[:,0].max(),iin[:,0].max())
        if lo<hi:xs=sorted(set(np.linspace(lo,hi,33).tolist()+([XREF] if lo<=XREF<=hi else [])))
    cuts=[]
    for x in xs:
        a,b,c,d=(_section(p,x) for p in (iin,iout,inner,outer))
        if any(v is None for v in (a,b,c,d)):continue
        upper=[a[1]-d[1],b[1]-c[1]];lower=[d[0]-a[0],c[0]-b[0]]
        cuts.append(dict(x=x,IR_inner=a,IR_outer=b,UV_inner=c,UV_outer=d,upper_shrink_interval=upper,
            lower_rise_interval=lower,upper_minus_lower_interval=[upper[0]-lower[1],upper[1]-lower[0]]))
    for cut in cuts:cut['upper_change_larger']=_interval_verdict(cut['upper_minus_lower_interval'])
    ref_cut=next((c for c in cuts if c['x']==XREF),None)
    from .kernels import certify_reference_section
    section_proof=certify_reference_section(ir['points'],points,XREF) if ref_cut is not None else None
    halfplanes=[dict(direction=p['direction'],upper=p['upper'],source=p['report']) for p in ir['points']]
    halfplanes += [dict(direction=[-1.,0.],upper=0.,source='Published x>=0 display window')]
    halfplanes += [dict(direction=p['direction'],upper=p['upper'],source=p['report']) for p in points]
    result=gauge_region_report(reference_signature,points,rejected,str(_path(args.region_summary)),(iin,iout,inner,outer),cuts,ref_cut,section_proof,halfplanes)
    from .io import plot_gauge_region
    plot_gauge_region(out,(iin,iout,inner,outer),points,XREF)
    write_json(out/'regions.json',result);return result

def gauge_phases(args):
    from . import save_figure
    from . import workspace; from .run import _path
    from . import _phase_data
    from .model import watson_diagnostics, native_peaks
    from .certificates import native_rho_certificate,ir_projection_exclusion
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True);profiles=[];reference_signature=None
    if len(args.profile_runs)!=3:raise ValueError('Exactly three preselected evaluations required')
    for path in args.profile_runs:
        d,x,delta,eta=_phase_data(path);sel=d['selection'];role=sel['role'];sig=dict(sel['model_signature'],chiral_norm=sel['model_signature'].get('chiral_norm','separate-l2'));sel=dict(sel,model_signature=sig)
        if role not in COLORS or reference_signature is not None and sig!=reference_signature or sel.get('coefficients_mixed') is not False:
            raise ValueError('Require three unchanged same-model selections')
        reference_signature=sig;hits=np.flatnonzero((delta[:-1,2]<90)&(delta[1:,2]>=90));crossing=None
        if len(hits):
            j=int(hits[0]);mass=x[j]+(90-delta[j,2])/(delta[j+1,2]-delta[j,2])*(x[j+1]-x[j])
            crossing=dict(native_energy_bracket_gev=x[j:j+2],native_phase_bracket_degrees=delta[j:j+2,2],
                linear_in_energy_gev=float(mass),relative_to_770mev=float(mass/.770-1),
                interpretation='First native upward 90° bracket; interpolation is descriptive, not a pole.')
        intensity=abs(eta[:,2]*np.exp(2j*np.radians(delta[:,2]))-1)**2/4
        profiles.append(dict(role=role,P1_intensity=intensity,P1_intensity_peaks=native_peaks(x,intensity),evaluation=str(_path(path)/'evaluation.json'),coefficients=d['coefficients'],selection=sel,
            energy_gev=x,phase_degrees=delta,eta=eta,unitarity=d['unitarity']['status'],maximum_eta_upper=d['maximum_eta'],
            P1_first_upward_90=crossing,P1_crossing_observed=crossing is not None,watson=watson_diagnostics(d)))
    if {p['role'] for p in profiles}!=set(COLORS):raise ValueError('tip, mid and ref are each required')
    profiles.sort(key=lambda p:list(COLORS).index(p['role']));x=profiles[0]['energy_gev']
    if any(not np.array_equal(x,p['energy_gev']) for p in profiles):raise ValueError('Require common native energy grid')
    for p in profiles:p['native_rho_certificate']=native_rho_certificate(read_json(p['evaluation']),args.bits)
    angles=np.stack([p['phase_degrees'] for p in profiles]);spread=np.ptp(angles,axis=0);separation={}
    circular=np.max(abs((angles[:,None]-angles[None,:]+90)%180-90),axis=(0,1))
    S=np.stack([p['eta'] for p in profiles])*np.exp(2j*np.radians(angles))
    distance=np.max(abs(S[:,None]-S[None,:]),axis=(0,1))
    for j,wave in enumerate(WAVES):
        k=int(np.argmax(spread[:,j]));separation[wave]=dict(pointwise_range_degrees=spread[:,j],maximum_degrees=float(spread[k,j]),
            energy_at_maximum_gev=float(x[k]),range_at_first_physical_node=float(spread[1,j]),range_at_last_node=float(spread[-1,j]),
            mod_pi_pairwise_max_degrees=circular[:,j],complex_S_pairwise_max_distance=distance[:,j],
            mod_pi_maximum_degrees=float(max(circular[:,j])),mod_pi_at_last_node=float(circular[-1,j]),
            complex_S_maximum=float(max(distance[:,j])),complex_S_at_last_node=float(distance[-1,j]))
    refs=[]
    if args.phase_reference:
        with Path(args.phase_reference).open() as stream:refs=[r for r in csv.DictReader(stream) if r['role'] in ('experiment','phenomenology')]
    ir_path=getattr(args,'baseline_profile',None);baseline=None
    if ir_path is not None:
        d,xb,db,eb=_phase_data(ir_path);bs=d['selection']['model_signature']
        keys=('M','L','preparation','prescription','density_limit','chiral_norm','infinity','unitarity_scope','tail_conditions_applied')
        if any(bs.get(k,'separate-l2' if k=='chiral_norm' else False if k=='tail_conditions_applied' else None)!=sig.get(k,False if k=='tail_conditions_applied' else None) for k in keys) or d['selection']['epsilon']!=sig['chiral_tolerance']:raise ValueError('IR baseline scattering/chiral inputs differ')
        if (Path(d['coefficients']).parent/'current.json').exists():raise ValueError('IR control must not contain joint current inputs')
        intensity=abs(eb[:,2]*np.exp(2j*np.radians(db[:,2]))-1)**2/4
        baseline=dict(evaluation=str(ir_path/'evaluation.json'),native_rho_certificate=native_rho_certificate(d,args.bits),coefficients=d['coefficients'],energy_gev=xb,phase_degrees=db,eta=eb,P1_intensity=intensity,P1_intensity_peaks=native_peaks(xb,intensity))
        baseline['UV_projection_exclusion']=ir_projection_exclusion(d['selection'],profiles,args.bits)
    paper={}
    for wave,filename in zip(WAVES,('figure10_s0_phases.csv','figure10_s2_phases.csv','figure9_p1_phases.csv')):
        source=workspace()/'references'/filename
        with source.open() as stream:records=list(csv.DictReader(stream))
        curves={group:np.asarray(sorted((float(r['energy_gev']),float(r['phase_deg'])) for r in records if r['group']==group))
            for group in ('red','pink','light_pink')}
        if any(v.ndim!=2 or v.shape[1]!=2 or not np.isfinite(v).all() or np.any(np.diff(v[:,0])<=0) for v in curves.values()):
            raise ValueError('Require three ordered v3 marker curves')
        paper[wave]=dict(csv=str(source),metadata=str(source.with_suffix('.metadata.json')),curves=curves,
            scope='Three paper amplitudes; no error band or extrapolation')
    for name,channels in [('fig9',[2]),('fig10',[0,1])]:
        fig,axes=plt.subplots(1,len(channels),figsize=(7*len(channels),4.8),layout='constrained');axes=np.atleast_1d(axes)
        for ax,j in zip(axes,channels):
            for k,curve in enumerate(paper[WAVES[j]]['curves'].values()):
                ax.plot(*curve.T,'--',color='#9c8066',alpha=.65,lw=.7,label='v3: three selected amplitudes' if k==0 else None,zorder=1)
            for role,style,label in [('phenomenology','--','PDF phenomenology'),('experiment','.','PDF experiment')]:
                pts=sorted((float(r['energy_gev']),float(r['delta_deg'])) for r in refs if r['wave']==WAVES[j] and r['role']==role)
                if pts:ax.plot(*np.asarray(pts).T,style,color='.55',ms=4,lw=1,label=label)
            if baseline is not None:ax.plot(baseline['energy_gev'],baseline['phase_degrees'][:,j],':',color='purple',lw=1,label='C: chiral only')
            for p in profiles:ax.plot(x,p['phase_degrees'][:,j],'.-',color=COLORS[p['role']],ms=4,lw=.8,label=p['role'])
            if j==2:ax.axhline(90,color='.5',ls='--',lw=.6)
            ax.set(xlim=(.28,1.2),xlabel='E [GeV]',ylabel='δ [degrees]',title=WAVES[j]);ax.grid(alpha=.15);ax.legend(fontsize=8)
        fig.suptitle('Fig.9: QCD-coupled P1' if name=='fig9' else 'Fig.10: the same three amplitudes, S0 and S2')
        save_figure(fig,out/f'{name}',180)
    fig,axes=plt.subplots(1,3,figsize=(14,4.2),layout='constrained')
    for j,ax in enumerate(axes):
        for p in profiles:ax.plot(x,np.asarray(p['eta'])[:,j],'.-',color=COLORS[p['role']],ms=4,lw=.8,label=p['role'])
        ax.axhline(1,color='.5',ls='--',lw=.6);ax.set(xlim=(.28,1.2),ylim=(0,1.03),xlabel='E [GeV]',ylabel='η = |S|',title=WAVES[j]);ax.grid(alpha=.15);ax.legend(fontsize=8)
    fig.suptitle('The same three amplitudes: native inelasticities')
    save_figure(fig,out/f'inelasticities',180)
    from .spectra import plot_spectral_profiles
    plot_spectral_profiles(out,profiles,baseline)
    result=dict(status='gauge_phases_completed',model_signature=sig,profiles=profiles,three_curve_separation=separation,
        chiral_only_baseline=baseline,reference=str(args.phase_reference) if args.phase_reference else None,paper_profiles=paper,
        paper_comparison_scope='Dashed: three v3 selections; no error band or selection input.',
        reference_roles=['experiment','phenomenology'],phase_convention='delta=unwrap(arg S)/2 from threshold; eta=|S|',
        selection_uses_phase_data=any(p['selection'].get('watson_iteration',0)>0 for p in profiles),continuum_certified=False,off_node_physical_continuation=False,
        watson_iterations={p['role']:p['selection'].get('watson_iteration',0) for p in profiles},experimental_phase_fitting=False,
        winding_scope='Native delta mod pi only; missing crossing does not exclude poles. Circular distances remove pi shifts.',
        comparison_scope='Native samples; no fitting, extrapolation or exact 6% target.',
        reference_scope='PDF markers/phenomenology; experimental uncertainty/covariance excluded.',
        baseline_scope='Explicit IR control; same model, possibly different objective.')
    write_json(out/'phases.json',result);return result

def select_gauge(args):
    from . import gauge_selection_report
    from .run import _path
    if getattr(args,'selection_reference',None) is not None:
        raise ValueError('Paper marker coordinates are comparison outputs, not scientific selection inputs')
    direct=bool(getattr(args,'support_runs',[]));required=getattr(args,'require_center',False)
    data={} if direct else read_json(args.region_summary)
    requests=args.support_runs if direct else [p['report'] for p in data['points']]
    sig,points,rejected=joint_support_points(requests,require_center=required)
    if direct and (len(points)!=3 or not all(p['support_optimality_certified'] for p in points)):
        raise ValueError('Three supports meeting their gap budgets required: '+str(rejected))
    if not direct and sig!=data['model_signature']:raise ValueError('Region/support model mismatch')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    tips=[p for p in points if p['direction'][0]>0 and p['direction'][1]==0]
    if not tips:raise ValueError('A +x joint support is required for the tip')
    frozen_dir=args.coefficients.parent if args.coefficients is not None else out/'tip'
    frozen_path=frozen_dir/'selection.json';frozen=read_json(frozen_path) if frozen_path.is_file() else None
    if frozen is not None:frozen=dict(frozen,model_signature=dict(frozen['model_signature'],chiral_norm=frozen['model_signature'].get('chiral_norm','separate-l2')))
    if frozen is None:tip=max(tips,key=lambda p:p['target'][0])
    else:
        matches=[p for p in tips if [_path(v) for v in frozen.get('source_reports',[])]==[_path(p['report'])]]
        if len(matches)!=1 or frozen.get('model_signature')!=sig:raise ValueError('Frozen tip needs the same support set')
        tip=matches[0]
        if [_path(v) for v in frozen.get('source_coefficients',[])]!=[_path(tip['coefficients'])] or not np.array_equal(frozen['target'],tip['target']):
            raise ValueError('Frozen tip changed; phase-based reselection forbidden')
    refx=XREF;mid=(XREF+tip['target'][0])/2;template=None
    rule='Tip +x; nearest audited upper points to xref and (xref+x_tip)/2; no phase ranking/mixing'
    def different(p,chosen):return all(np.max(abs(np.asarray(p['target'])-q['target']))>1e-12 for q in chosen)
    upper=[p for p in points if p['direction'][1]>0 and different(p,[tip])]
    if len(upper)<2:raise ValueError('Two distinct actual upper-branch supports are required')
    ref=min(upper,key=lambda p:(abs(p['target'][0]-refx),p['gap'],p['report']))
    remaining=[p for p in upper if different(p,[ref])]
    if not remaining:raise ValueError('Three distinct projected representatives required')
    middle=min(remaining,key=lambda p:(abs(p['target'][0]-mid),p['gap'],p['report']))
    selected=[];pcount=1+2*sig['M']+sig['M']**2+sig['M']*(sig['M']+1)//2
    for role,point,targetx in [('tip',tip,tip['target'][0]),('mid',middle,mid),('ref',ref,refx)]:
        reuse=role=='tip' and frozen is not None;dest=out/role;dest.mkdir(exist_ok=reuse);source=Path(point['report']).parent
        c=read_json(source/'coefficients.json')['coefficients']
        with np.load(source/'joint.npz') as saved:
            full=saved['point'].copy()
            if full.shape!=(pcount+4*sig['M'],) or not np.array_equal(full[:pcount],c):raise ValueError('Complete joint/amplitude mismatch')
        if reuse:
            if any(read_json(frozen_dir/name)!=read_json(source/name) for name in ('coefficients.json','current.json')):
                raise ValueError('Frozen tip amplitude/current differs')
            with np.load(frozen_dir/'joint.npz') as saved:
                if not np.array_equal(saved['point'],full):raise ValueError('Frozen joint coefficients changed')
        for name in ('coefficients.json','current.json','joint.npz'):shutil.copyfile(source/name,dest/name)
        sel=dict(role=role,color=COLORS[role],rule=rule,selection_template=template,
            reference=[XREF,-XREF/15],target=point['target'],target_x=targetx,x_deviation=point['target'][0]-targetx,epsilon=.002,
            source_reports=[point['report']],source_coefficients=[point['coefficients']],source_joint=point['joint'],
            model_signature=sig,direction=point['direction'],support_lower=point['lower'],support_upper=point['upper'],support_gap=point['gap'],
            support_optimality_certified=point['support_optimality_certified'],representative_center_converged=point['representative_center_converged'],
            boundary_y_distance_upper=point['gap']/point['direction'][1] if point['direction'][1]>0 else None,exact_extremizer_claimed=False,
            joint_feasible=True,selection_uses_phase_data=False,selection_uses_paper_marker_coordinates=False,coefficients_mixed=False,finite_only=True)
        if reuse:sel=frozen
        write_json(dest/'selection.json',sel)
        selected.append(dict(directory=str(dest.resolve()),reused_frozen_selection=reuse,**sel))
    inner=np.empty((0,2)) if direct else np.asarray(data['inner_vertices']);outer=np.empty((0,2)) if direct else np.asarray(data['outer_vertices'])
    result=gauge_selection_report(sig,selected,None if direct else str(_path(args.region_summary)),inner,outer,frozen,rule,template)
    result.update(source_support_reports=[p['report'] for p in points],rejected=rejected,representative_center_required=required,
        direct_support_selection=direct,complete_region_claimed=False,plots_created=not direct)
    if direct:result.update(reference_x_in_inner_projection=None,reference_x_in_outer_projection=None)
    else:
        from .io import plot_gauge_selections
        plot_gauge_selections(out,data,selected,XREF)
    write_json(out/'selection.json',result);return result

def _hull(points):
    p=np.unique(np.asarray(points,float).reshape(-1,2),axis=0)
    if len(p)<3:return p
    q=np.column_stack((p[:,0]/XREF,(p[:,1]+p[:,0]/15)/.002))
    try:return p[ConvexHull(q).vertices]
    except QhullError:return p[np.argsort(p[:,0])][[0,-1]]

def _clip(p,d,h):
    p=np.asarray(p,float).reshape(-1,2);d=np.asarray(d,float);out=[]
    if not len(p):return p
    for a,b in zip(p,np.roll(p,-1,axis=0)):
        va,vb=d@a-h,d@b-h
        if va<=0:out.append(a)
        if (va<0<vb) or (vb<0<va):out.append(a+va/(va-vb)*(b-a))
    return np.asarray(out,float).reshape(-1,2)

def _section(p,x):
    p=np.asarray(p,float).reshape(-1,2);ys=[]
    for a,b in zip(p,np.roll(p,-1,axis=0)):
        if a[0]==b[0]:
            if x==a[0]:ys.extend([a[1],b[1]])
        elif min(a[0],b[0])<=x<=max(a[0],b[0]):ys.append(a[1]+(x-a[0])*(b[1]-a[1])/(b[0]-a[0]))
    return [min(ys),max(ys)] if ys else None

def _interval_verdict(interval):
    return True if interval[0]>0 else False if interval[1]<0 else None

def compare_phase(x,delta,points):
    xp,yp=np.asarray(points).T;mask=(x>=xp.min())&(x<=xp.max());error=delta[mask]-np.interp(x[mask],xp,yp)
    return dict(samples=int(mask.sum()),energy_range=[float(x[mask].min()),float(x[mask].max())],mean_absolute_degrees=float(np.mean(abs(error))),rms_degrees=float(np.sqrt(np.mean(error**2))),maximum_absolute_degrees=float(max(abs(error))))

def phase_shifts(f,energies):
    s=np.asarray(energies,float);f=np.asarray(f,complex)
    if s.ndim!=1 or f.ndim!=2 or len(s)!=len(f) or s[0]!=4 or np.any(np.diff(s)<=0):raise ValueError('Ordered physical grid starting at threshold required')
    S=1+1j*np.pi*np.sqrt(1-4/s)[:,None]*f
    if not np.isfinite(S).all() or np.any(abs(S)==0):raise ValueError('A finite nonzero S is required to define its phase')
    return np.unwrap(np.angle(S),axis=0)/2,np.abs(S)


from .spectra import ff_guided_phase
