"""Typed records and saved figure rendering."""
from flint import arb, ctx
from pathlib import Path
from scipy.spatial import ConvexHull
import csv
import json
import matplotlib.pyplot as plt
import numpy as np
from . import read_cflat, read_json, save_figure, write_json, AMPLITUDE_COORDINATES
CONFIGURATIONS = ((50,8),(50,10),(50,12),(45,10),(60,10))

def calculation_inputs(args):
    """Snapshot directly consumed operators, points and phase inputs."""
    from . import digest
    paths=[];names=('report.json','coefficients.json','joint.npz','objective.npz','selection.json','current.json','watson_update.json','barrier_state.npz','phase_I_state.npz')
    for folder,files in ((args.preparation,('report.json','amplitude.npz')),
            (args.current_preparation,('report.json','current_data.npz','gram_linear.npz','moment_linear.npz','ff_cap_linear.npz'))):
        if folder is not None:paths.extend(folder/name for name in files)
    sources=[p for p in (args.coefficients,args.interior_coefficients) if p is not None]
    profiles=list(getattr(args,'profile_runs',[]))
    if getattr(args,'baseline_profile',None):profiles.append(args.baseline_profile)
    for folder in profiles:
        path=Path(folder)/'evaluation.json';paths.append(path)
        if path.is_file():
            data=read_json(path);sources.append(Path(data['coefficients']))
            current=data.get('selection',{}).get('model_signature',{}).get('current_preparation')
            if current:paths.append(Path(current)/'current_data.npz')
    for source in sources:paths.append(source);paths.extend(source.parent/name for name in names)
    anchor=getattr(args,'probe_anchor',None)
    if anchor is not None:
        anchor=anchor/'report.json' if anchor.is_dir() else anchor;paths.append(anchor);paths.extend(anchor.parent/name for name in names)
    if getattr(args,'command',None)=='support-probe' and args.coefficients is not None:paths.append(args.coefficients.parent/'probe_goal.json')
    if getattr(args,'command',None) in ('support-probe','select-gauge','gauge-regions'):
        for source in args.support_runs:
            folder=Path(source);folder=folder.parent if folder.is_file() else folder
            paths.extend(folder/name for name in names)
    if getattr(args,'command',None)=='gauge-regions':paths.append(args.region_summary)
    for path in list(paths):
        if path.name=='objective.npz' and path.is_file():
            with np.load(path) as data:source=Path(json.loads(str(data['metadata']))['source_coefficients'])
            paths.append(source);paths.extend(source.parent/name for name in names)
    if getattr(args,'phase_reference',None):paths.append(args.phase_reference)
    if getattr(args,'command',None)=='gauge-phases':paths.extend(Path('references')/f for f in ('figure9_p1_phases.csv','figure10_s0_phases.csv','figure10_s2_phases.csv'))
    return {str(p.resolve()):dict(sha256=digest(p)) for p in paths if p.is_file()}


def select_chiral_amplitude(args):
    if args.support_runs:
        from .scattering import select_chiral_support
        return select_chiral_support(args)
    from .imaginary import support_outer
    data=read_json(args.region_summary);sig=data['model_signature']
    if sig['prescription'] not in AMPLITUDE_COORDINATES or sig['unitarity_scope']!='sampled' or sig['infinity']!='free':raise ValueError('Declared sampled/free model required')
    region=next(g for g in data['regions'] if g['epsilon']==args.chiral_tolerance)
    points=region['points'];P=np.array([p['target'] for p in points]);hull=ConvexHull(P).vertices
    xref=5/(16*np.pi**2*(92/140)**2);reference=np.array([xref,-xref/15]);choices=[]
    for i,j in zip(hull,np.roll(hull,-1)):
        dx=P[j,0]-P[i,0]
        if dx==0:continue
        t=(xref-P[i,0])/dx
        if 0<=t<=1:choices.append(((1-t)*P[i,1]+t*P[j,1],int(i),int(j),float(t)))
    if not choices:raise ValueError('Reference coupling outside the saved boundary hull')
    _,i,j,t=max(choices);paths=[points[k]['coefficients'] for k in (i,j)]
    vectors=[read_cflat(p,sig['M'],prescription=sig['prescription'])[0] for p in paths]
    c=np.asarray((1-np.longdouble(t))*vectors[0].astype(np.longdouble)+np.longdouble(t)*vectors[1],float)
    with np.load(Path(sig['preparation'])/'amplitude.npz') as saved:H=saved['rows'];energies=saved['energies']
    n=(H.shape[0]-11)//2;kap=np.repeat(np.pi*np.sqrt(1-4/energies),3*sig['L'])
    audit=support_outer(H,kap,np.zeros(n),np.zeros(n),np.zeros(8),[0.,0.],sig['density_limit'],args.chiral_tolerance,bits=args.bits,coefficients=c,waves_per_isospin=sig['L'],chiral_norm=sig.get('chiral_norm','separate-l2'))
    if not audit['primal_feasible']:raise ArithmeticError('Rounded boundary interpolation failed the original finite constraints')
    target=np.asarray(H[-2:],np.longdouble)@c.astype(np.longdouble)
    with ctx.workprec(args.bits):
        xx,yy=map(arb,audit['targets']);bounds=[]
        for p in points:
            a,b=p['direction']
            if b>0:bounds.append(((arb(p['upper'])-arb(a)*xx)/arb(b)-yy,p['report']))
        distance,plane=min(bounds,key=lambda v:float(v[0].upper()))
        if distance<0:raise ArithmeticError('A feasible selected point contradicts a global upper support')
        raw=float(np.nextafter(float(distance.upper()),np.inf));scaled=float(np.nextafter(float((distance/arb(args.chiral_tolerance)).upper()),np.inf))
    budget=region['geometry']['budgets']['distance'];local_ready=scaled<=budget
    if not region['geometry']['geometry_ready'] and not local_ready:raise ValueError(f'Upper section at xref needs refinement: scaled gap {scaled:g} > {budget:g}')
    from . import chiral_selection_record
    selection=chiral_selection_record(reference,target,args.chiral_tolerance,[points[i],points[j]],t,P[hull],region['geometry'],audit,sig,(raw,scaled,budget,plane,local_ready))
    write_json(args.output/'coefficients.json',dict(coordinates=AMPLITUDE_COORDINATES[sig['prescription']],prescription=sig['prescription'],coefficients=c.tolist()))
    write_json(args.output/'selection.json',selection)
    return dict(status='selected',model_signature=sig,selection=selection,sampled_feasible=True,optimization_performed=False)

def plot_profiles(args):
    from .spectra import plot_profiles as render
    return render(args)


def plot_regions(directory, reference_directory):
    from scipy.spatial import QhullError
    directory, reference_directory = Path(directory), Path(reference_directory)
    data = json.loads((directory/'regions.json').read_text())
    if not any(p.get('feasible') is True for r in data['regions'] for p in r['points']):
        raise ValueError('Plot requires computed supports')
    reference = json.loads((reference_directory/'report.json').read_text())
    refs = {name: list(csv.DictReader((reference_directory/f'{name}_points.csv').open()))
            for name in ('purSplot', 'chiplot')}
    palette = {v['epsilon']: v['color'] for v in reference['curves']['chiplot']['series']}
    palette.update({.001: 'orangered', .0006: 'mediumpurple', .0002: 'saddlebrown'})
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    def curve(ax, points, color, label, shear=False, outer=False):
        points = np.asarray(points, dtype=float).reshape(-1, 2).copy()
        if not len(points): return
        if shear: points[:, 1] += points[:, 0]/15
        if len(points) >= 3:
            try:
                hull = points[ConvexHull(points).vertices]
                ax.plot(*np.vstack([hull, hull[0]]).T, color=color,
                        linestyle='--' if outer else '-', linewidth=.8, label=label)
                if not outer: ax.fill(*hull.T, color=color, alpha=.06)
            except QhullError: ax.scatter(*points.T, s=12, color=color, label=label)
        else: ax.scatter(*points.T, s=12, color=color, label=label)
    comparisons = []
    for region in data['regions']:
        pure, epsilon = region['mode'] == 'pure', region.get('epsilon')
        name, color = ('purSplot', 'royalblue') if pure else ('chiplot', palette.get(epsilon,'tab:cyan'))
        points = [p for p in region['points'] if p.get('feasible') is True
                  and np.isfinite(p['target']).all()]
        inner = np.array([p['target'] for p in points], dtype=float).reshape(-1, 2)
        geometry=region['geometry'];inner=np.asarray(geometry['inner_vertices']).reshape(-1,2);outer=geometry['outer_vertices'] if geometry['outer_closed'] else []
        published = np.array([(float(p['x']), float(p['y'])) for p in refs[name]
            if p['role'] == 'boundary_sample' and (pure or float(p['epsilon']) == epsilon)]).reshape(-1,2)
        for ax, shear in [(axes[0], False)] if pure else [(axes[1], False), (axes[2], True)]:
            xy = published.copy()
            if shear: xy[:, 1] += xy[:, 0]/15
            ax.scatter(*xy.T, s=3, color=color, alpha=.35, label='Paper' if pure else None)
            label = 'Computed inner' if pure else f'ε={epsilon:g}'
            curve(ax, inner, color, label, shear); curve(ax, outer, color, 'Outer' if pure else None, shear, True)
        gaps = [p['upper']-p['lower'] for p in points]
        comparisons.append(dict(mode=region['mode'], epsilon=epsilon, feasible_points=len(inner),
            outer_closed=geometry['outer_closed'],max_support_gap=max(gaps, default=None),
            direction_coverage=geometry['direction_coverage'],
            computed_min=inner.min(axis=0).tolist() if len(inner) else None,
            computed_max=inner.max(axis=0).tolist() if len(inner) else None,
            paper_min=published.min(axis=0).tolist() if len(published) else None,
            paper_max=published.max(axis=0).tolist() if len(published) else None))
    chiral = reference['curves']['chiplot']
    black = next(p for p in chiral['highlights'] if p['color'] == 'black')
    for ax, shear in [(axes[1], False), (axes[2], True)]:
        if not any(r['mode']=='chiral' and r['points'] for r in data['regions']):
            ax.text(.5,.5,'No chiral supports',ha='center',transform=ax.transAxes);continue
        ax.axline((0, 0), slope=0 if shear else -1/15, color='.4', linewidth=.6)
        for p, marker, label in [(black, 'o', 'PDF black dot'), (chiral['formula_reference_mpi140_fpi92'], '+', '140/92 formula')]:
            ax.scatter(p['x'], p['y']+(p['x']/15 if shear else 0), color='black', marker=marker, label=label)
    if not any(r['mode']=='pure' and r['points'] for r in data['regions']):
        axes[0].text(.5,.5,'No pure-S supports',ha='center',transform=axes[0].transAxes)
    for ax, title in zip(axes, ['Fig.3: pure S', 'Fig.4: chiral supports', 'Chiral widths: y+x/15']):
        ax.set(title=title, xlabel=r'$f^0_0(3)$'); ax.grid(alpha=.15)
        if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=6)
    axes[0].set_ylabel(r'$f^1_1(3)$'); axes[1].set_ylabel(r'$f^1_1(3)$'); axes[2].set_ylabel(r'$f^1_1(3)+f^0_0(3)/15$')
    fig.suptitle('Conditional model; finite sampling — solid: feasible hull, dashed: support outer bound')
    save_figure(fig,directory/'regions')
    summary = dict(regions=comparisons,nesting=data['nesting'],status=data['status'],reference_semantics='PDF marker centers; finite plotted samples',
                   reference_points=dict(pdf=black, formula=chiral['formula_reference_mpi140_fpi92']))
    (directory/'comparisons.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary

def summarize_regions(directory,support_runs,reference_directory,nodes=50,chiral_norm='separate-l2'):
    from . import workspace
    from . import decode_real_ball
    from .run import encode_real_ball; from .kernels import region_geometry
    from flint import arb_mat
    directory=Path(directory);groups={};rejected=[];signature=None;seen=set();target_rows={}
    if type(nodes) is not int or nodes<1:raise ValueError('Positive nodes required')
    identity=lambda value:None if value is None else str((workspace()/Path(value)).resolve())
    for requested in support_runs:
        path=Path(identity(requested));path=path/'report.json' if path.is_dir() else path
        reasons=[];entry=dict(report=str(path));coefficient_file=path.parent/'coefficients.json'
        if path in seen:
            rejected.append(dict(entry,reasons=['duplicate_report']));continue
        seen.add(path)
        try:
            report=read_json(path);parameters=report.get('parameters',{});outer=report.get('outer',{});diagnostic=report.get('high_spin_necessary',{})
            scope=report.get('unitarity_scope');fg=report.get('global_high_spin_enforced',False);tail=report.get('tail_conditions_applied',False)
            if scope!='sampled' or fg or tail or report.get('infinity')!='free' or report.get('prescription') not in AMPLITUDE_COORDINATES:raise ValueError('Require PV sampled/free or analytic-cardinal sampled/free')
            if report.get('command') not in ('boundary','dual'):reasons.append('not_a_support_report')
            if report.get('M')!=nodes:reasons.append(f'wrong_M: expected {nodes}, found {report.get("M")}')
            if outer.get('primal_feasible') is not True or not outer.get('enclosure'):reasons.append('unverified_outer_bound')
            if report.get('sampled_feasible') is not True:reasons.append('sampled_feasible is not true')
            if report.get('status') in ('running','inconclusive') or report.get('timeout') or report.get('exit_code',0)!=0:
                reasons.append('run_did_not_complete_successfully')
            if report.get('source_changes') or report.get('input_changes'):reasons.append('recorded_source_or_input_changes')
            if not coefficient_file.is_file():reasons.append('missing_complete_coefficients')
            if reasons:
                rejected.append(dict(entry,reasons=reasons));continue
            mode=report['mode'];epsilon=None if mode=='pure' else float(report['chiral_tolerance'])
            if mode not in ('pure','chiral') or epsilon is not None and (not np.isfinite(epsilon) or epsilon<=0):
                raise ValueError('Invalid mode or epsilon')
            if mode=='chiral' and report.get('chiral_norm','separate-l2')!=chiral_norm:raise ValueError('Different chiral norm')
            if not parameters.get('preparation'):raise ValueError('Missing preparation')
            current=dict(chiral_norm=chiral_norm,M=nodes,L=int(report['L']),density_limit=float(report['density_limit']),infinity=report['infinity'],
                unitarity_scope=scope,global_high_spin_enforced=fg,tail_conditions_applied=tail,prescription=report['prescription'],
                preparation=identity(parameters['preparation']),additional_constraints=identity(parameters.get('additional_constraints')),
                sampling_nodes=report.get('sampling_nodes'),scattering_samples=report.get('scattering_samples'))
            if current['L']!=report['L'] or current['L']<1 or not np.isfinite(current['density_limit']) or current['density_limit']<=0 or current['infinity']!='free':
                raise ValueError('Invalid model identity')
            target=np.asarray(report['targets'],float);direction=np.asarray(report.get('support_direction',parameters['direction']),float)
            bound=arb(outer['enclosure'])
            if not bound.is_finite():raise ValueError('Finite Arb upper enclosure required')
            lower,upper=np.nextafter(float(outer['lower']),-np.inf),np.nextafter(max(float(outer['upper']),float(bound.upper())),np.inf)
            if target.shape!=(2,) or direction.shape!=(2,) or not np.isfinite(np.r_[target,direction,lower,upper]).all() or np.hypot(*direction)==0 or lower>upper:
                raise ValueError('Invalid target, direction or bounds')
            vector,declared=read_cflat(coefficient_file,nodes)
            if declared['prescription']!=current['prescription']:raise ValueError('Report/coefficient prescription mismatch')
            if signature is not None and current!=signature:
                different=[key for key in signature if current[key]!=signature[key]]
                rejected.append(dict(entry,reasons=['incomparable_model: '+', '.join(different)],model=current));continue
            with ctx.workprec(384):
                preparation=current['preparation']
                if preparation not in target_rows:
                    with np.load(Path(preparation)/'amplitude.npz',allow_pickle=False) as saved:rows=saved['rows'][-2:].copy()
                    if rows.shape!=(2,len(vector)):raise ValueError('Invalid target rows')
                    target_rows[preparation]=arb_mat([[arb(float(x)) for x in row] for row in rows])
                balls=(target_rows[preparation]*arb_mat(len(vector),1,list(map(float,vector)))).entries()
                projected=sum((arb(float(x))*y for x,y in zip(direction,balls)),arb(0))
                if projected<arb(float(lower)) or projected>arb(float(upper)):raise ValueError('Coefficients contradict support bounds')
                target=[float(v.mid()) for v in balls];enclosures=list(map(encode_real_ball,balls))
            if signature is None:signature=current
            key=(mode,epsilon)
            group=groups.setdefault(key,dict(mode=mode,epsilon=epsilon,chiral_norm=chiral_norm,points=[]))
            group['points'].append(dict(target=target,target_enclosures=enclosures,direction=direction.tolist(),lower=float(lower),upper=float(upper),
                gap=float(np.nextafter(upper-lower,np.inf)),pointwise_tolerance_met=report.get('support_optimality_certified') is True and upper-lower<=parameters.get('gap',1e-4),
                diagnostics=dict(high_spin={k:diagnostic.get(k) for k in ('status','violations')},tail_passed=report.get('analytic_tail_conditions_passed'),physical_unitarity_accepted=report.get('physical_unitarity_accepted',False)),
                feasible=True,report=str(path),coefficients=str(coefficient_file)))
        except (OSError,ValueError,TypeError,KeyError) as error:
            rejected.append(dict(entry,reasons=['invalid_record: '+str(error)]))
    contradictions=[];ordered=sorted(groups.values(),key=lambda g:float('inf') if g['epsilon'] is None else g['epsilon'])
    with ctx.workprec(384):
        for j,small in enumerate(ordered):
            for large in ordered[j:]:
                for p in small['points']:
                    for q in large['points']:
                        d=list(map(arb,q['direction']));s=list(map(arb,p['direction']));claims=[sum((x*decode_real_ball(y) for x,y in zip(d,p['target_enclosures'])),arb(0))-q['upper']]
                        if (s[0]*d[1]-s[1]*d[0]).is_zero() and sum((x*y for x,y in zip(d,s)),arb(0))>0:
                            k=0 if not d[0].is_zero() else 1;claims.append(arb(p['lower'])*d[k]/s[k]-q['upper'])
                        if any(v>0 for v in claims):contradictions.append(dict(smaller=p['report'],larger=q['report'],positive_margin=[v.str(20) for v in claims if v>0]))
    for group in ordered:
        group['published_positive_x']=group['epsilon'] is not None;group['geometry']=region_geometry(group)
    six={.0002,.0006,.001,.002,.004,.006}<=set(g['epsilon'] for g in ordered)
    from . import region_summary_record
    data=region_summary_record(nodes,signature,groups,rejected,ordered,contradictions,six)
    if data['B2_geometry_ready'] and data['B3_geometry_ready']:data['status']='geometry_ready'
    write_json(directory/'regions.json',data)
    if groups:
        comparison=plot_regions(directory,reference_directory)
        for group,item in zip(data['regions'],comparison['regions']):
            group.update(direction_coverage=item['direction_coverage'],outer_closed=item['outer_closed'])
        data['plots_created']=True;write_json(directory/'regions.json',data)
    return data

def plot_gauge_region(out,polygons,points,xref):
    iin,iout,inner,outer=polygons
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for ax,shear in zip(axes,(False,True)):
        for p,color,label,style in [(iin,'forestgreen','IR feasible hull','-'),(iout,'forestgreen','IR outer','--'),(inner,'darkturquoise','UV feasible hull','-'),(outer,'darkturquoise','UV outer','--')]:
            if not len(p):continue
            q=np.vstack([p,p[0]]).copy()
            if shear:q[:,1]+=q[:,0]/15
            ax.plot(*q.T,color=color,lw=1,ls=style,label=label)
            if style=='-' and len(p)>=3:ax.fill(*q.T,color=color,alpha=.07)
        q=np.asarray([p['target'] for p in points]);ax.scatter(q[:,0],q[:,1]+(q[:,0]/15 if shear else 0),s=15,color='darkturquoise')
        ax.scatter(xref,0 if shear else -xref/15,color='black',s=25,label='140/92 reference')
        ax.axline((0,0),slope=0 if shear else -1/15,color='.4',lw=.7);ax.set_xlim(left=0)
        ax.set(xlabel=r'$f^0_0(3)$',ylabel=r'$f^1_1(3)+f^0_0(3)/15$' if shear else r'$f^1_1(3)$');ax.grid(alpha=.15);ax.legend(fontsize=7)
    fig.suptitle('Fig.8: fixed chiral setting, with current/QCD constraints')
    save_figure(fig,out/f'regions',180)

def plot_gauge_selections(out,data,selections,xref):
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    for key,color,label,style in [('IR_inner_vertices','forestgreen','IR feasible hull','-'),('inner_vertices','darkturquoise','UV feasible hull','-'),('outer_vertices','darkturquoise','UV outer','--')]:
        p=np.asarray(data[key]);p=np.vstack([p,p[0]]);ax.plot(*p.T,color=color,label=label,ls=style,lw=1)
    for sel in selections:ax.scatter(*sel['target'],color=sel['color'],s=45,label=sel['role'],zorder=5)
    ax.scatter(xref,-xref/15,color='black',s=30,label='140/92 reference');ax.axline((0,0),slope=-1/15,color='.4',lw=.7)
    ax.set(xlim=(0,None),xlabel=r'$f^0_0(3)$',ylabel=r'$f^1_1(3)$',title='Fig.8: three representatives fixed before phases');ax.grid(alpha=.15);ax.legend(fontsize=8)
    save_figure(fig,out/f'fig8_selected',180)

def plot_resolution_profiles(*args,**kwargs):
    from .spectra import plot_resolution_profiles as plot
    return plot(*args,**kwargs)


def watson_lineage(args,target,functional,write_json,read_json,result):
    from . import digest
    source=Path(functional['source_coefficients']).parent/'selection.json'
    if not source.is_file():return
    old=read_json(source);record=dict(old)
    sig=old['model_signature'];norm=getattr(args,'chiral_norm','separate-l2');M=sig['M'];bound=args.density_limit or 100*(M*M+M*(M+1)//2)
    if any(Path(sig[k]).resolve()!=getattr(args,k).resolve() for k in ('preparation','current_preparation')) or sig['density_limit']!=bound or sig['chiral_tolerance']!=args.chiral_tolerance or sig.get('chiral_norm','separate-l2')!=norm:raise ValueError('Geometric role requires the same model')
    record['model_signature']=dict(sig,chiral_norm=norm)
    for key in ('direction','support_lower','support_upper','support_gap','boundary_y_distance_upper','support_optimality_certified','target_x','x_deviation','resolution_comparison'):record.pop(key,None)
    record.update(target=target.tolist(),initial_target=old.get('initial_target',old['target']),initial_selection=old.get('initial_selection',str(source)),
        watson_iteration=old.get('watson_iteration',0)+1,watson_reference=functional['source_coefficients'],
        source_reports=[str(args.output/'report.json')],source_coefficients=[str(args.output/'coefficients.json')],source_joint=str(args.output/'joint.npz'),
        rule='2403 Watsonian from a fixed initial role; projection coordinates released',
        uses_previous_bootstrap_form_factor_phase=True,initial_selection_uses_phase_data=False,selection_uses_phase_data=True,
        selection_uses_experimental_data=False,coefficients_mixed=False)
    record.update(objective_kind='watson',joint_feasible=bool(result.get('joint_feasible')),
        representative_center_converged=bool(result.get('joint_feasible') and result.get('solver',{}).get('representative_center_converged') is True),
        support_optimality_certified=bool(result.get('support_optimality_certified')))
    for name,suffix in (('coefficient','coefficients.json'),('current','current.json'),('joint','joint.npz')):record[name+'_sha256']=digest(args.output/suffix)
    write_json(args.output/'selection.json',record)
