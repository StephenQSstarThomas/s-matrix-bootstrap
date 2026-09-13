"""Explicit analytic scattering samples, retaining native current identities."""
from pathlib import Path
import numpy as np
from flint import arb,acb,arb_mat,ctx
from . import read_json,write_json,digest


def prepared_sampling(data,M,L):
    """Return every physical row label; native constraints must remain first."""
    from .operators import midpoint_grid
    from .kernels import canonical_waves
    with ctx.workprec(max(256,ctx.prec)):native=np.repeat([float(s) for s in midpoint_grid(M)[0]],3*L)
    waves=np.tile(canonical_waves(L),(M,1))
    ss=np.asarray(data.get('sample_energies',native),float);ww=np.asarray(data.get('sample_waves',waves),int)
    n=(len(data['rows'])-11)//2
    if ss.shape!=(n,) or ww.shape!=(n,2) or n<3*M*L or not np.isfinite(ss).all() or np.any(ss<=4):raise ValueError('Finite physical sample metadata required')
    if not np.array_equal(ss[:len(native)],native) or not np.array_equal(ww[:len(native)],waves):raise ValueError('Exact native sampling prefix required')
    if any(I not in (0,1,2) or ell<0 or ell%2!=int(I==1) for I,ell in ww):raise ValueError('Allowed physical partial waves required')
    if len({(float(s),int(I),int(l)) for s,(I,l) in zip(ss,ww)})!=n:raise ValueError('Duplicate scattering sample')
    return ss,ww


def prepare_additional(args,progress):
    from .basis import CardinalSourceRows
    from .quadrature import polynomial_rows
    from .certificates import verify_scattering_layout
    parent=Path(args.preparation);r=read_json(parent/'report.json');meta=r['amplitude_model'];M,L=args.nodes,args.waves
    if (meta['prescription'],meta['M'],meta['L'])!=('analytic-cardinal',M,L):raise ValueError('Additional rows require the same analytic M/L preparation')
    source=CardinalSourceRows(M,L,args.bits,args.angular_order,registry=args.source_registry)
    if source.metadata['source_registry_sha256']!=meta.get('source_registry_sha256'):raise ValueError('Same analytic source registry required')
    with np.load(parent/'amplitude.npz') as d:
        ss,ww=prepared_sampling(d,M,L);H=d['rows'];E=d['radius_upper'];labels=d['coefficient_labels'];nodes=d['energies'];waves=d['waves']
    n=len(ss);p=H.shape[1]
    if E.shape!=H.shape or not np.isfinite(H).all() or not np.isfinite(E).all() or np.any(E<0):raise ValueError('Complete analytic coefficient enclosures required')
    targets=sorted(set(map(float,args.unitarity_energies)))
    if any(not np.isfinite(s) or s<=4 for s in targets):raise ValueError('Additional energies must be finite and strictly above threshold')
    requested=[(0,0),(1,1),(2,0)] if args.primary_waves else source.waves
    existing={(float(s),int(I),int(l)) for s,(I,l) in zip(ss,ww)};real=[];imag=[];er=[];ei=[];added=[];accuracy={}
    with ctx.workprec(args.bits):
        for s in targets:
            missing=[w for w in requested if (s,*w) not in existing]
            if not missing:continue
            rows,bounds=polynomial_rows(arb(s),M,missing,order=args.angular_order);accuracy[str(s)]=bounds
            for I,l in missing:
                row=rows[I,l];parts=[];radii=[]
                for part in ('real','imag'):
                    vals=[getattr(v,part) for v in row];mid=np.array([float(v.mid()) for v in vals])
                    rad=[float(np.nextafter(float(abs(v-arb(float(x))).upper()),np.inf)) if not v.is_zero() else 0. for v,x in zip(vals,mid)]
                    parts.append(mid);radii.append(rad)
                real.append(parts[0]);imag.append(parts[1]);er.append(radii[0]);ei.append(radii[1]);added.append((s,I,l))
            progress(stage='additional_analytic_rows',energy=s,waves=len(missing))
    if not added:raise ValueError('No new physical samples requested')
    H=np.vstack((H[:n],real,H[n:2*n],imag,H[2*n:]));E=np.vstack((E[:n],er,E[n:2*n],ei,E[2*n:]))
    ss=np.r_[ss,[s for s,_,_ in added]];ww=np.vstack((ww,[(I,l) for _,I,l in added]));N=len(ss)
    verify_scattering_layout(H,[arb.pi()*(1-4/arb(float(s))).sqrt() for s in ss],M,L)
    np.savez_compressed(args.output/'amplitude.npz',rows=H,radius_upper=E,energies=nodes,waves=waves,
        coefficient_labels=labels,sample_energies=ss,sample_waves=ww)
    model=dict(meta,physical_samples=N,analytic_constraint_grid='Native prefix plus explicit additional samples',
        additional_scattering_samples=[dict(s=float(s),isospin=int(I),ell=int(l)) for s,(I,l) in zip(ss[3*M*L:],ww[3*M*L:])])
    inputs={str(q.resolve()):dict(sha256=digest(q)) for q in (parent/'report.json',parent/'amplitude.npz')}
    inputs.update(source.inputs)
    return dict(status='prepared',amplitude_model=model,inputs=inputs,restrictions=r['restrictions'],
        operator=dict(path=str(args.output/'amplitude.npz'),shape=list(H.shape),storage='Float64 plus outward analytic row radii'),
        extra_rows_added=len(added),native_rows_unchanged=True,quadrature_accuracy=accuracy,
        optimization_performed=False,old_feasibility_transferred=False,maximum_operator_radius_upper=float(E.max()))


def replay_support_data(args,prior,M,n,high):
    """Read a fixed objective and the latest duals; never generate a new Q."""
    from . import zero_joint_duals
    import json
    kind=prior.get('objective_kind','projection');objective=functional=None
    if kind not in ('projection','watson'):raise ValueError('Unknown saved objective')
    if kind=='watson':
        with np.load(args.coefficients.parent/'objective.npz') as saved:
            objective=saved['coefficients'].copy();functional=json.loads(str(saved['metadata']))
        if objective.ndim!=1 or not np.isfinite(objective).all() or functional!=prior.get('functional'):raise ValueError('Saved Watson objective identity differs')
    same=all(prior.get('parameters',{}).get(k) is not None and Path(prior['parameters'][k]).resolve()==getattr(args,k).resolve() for k in ('preparation','current_preparation'))
    clean=not (prior.get('source_changes') or prior.get('input_changes') or prior.get('timeout') or prior.get('exit_code',0) or prior.get('status') in ('running','inconclusive'))
    retained=bool(prior.get('support_optimality_certified') and same and clean)
    duals=zero_joint_duals(M,n,high);normal=prior['support_direction'] if retained and kind=='projection' else [0,0]
    if same and clean:
        with np.load(args.coefficients.parent/'joint.npz') as saved:duals={k:saved[k] for k in list(duals)+['asymptotic_kR','asymptotic_kI'] if k in saved}
    return retained,normal,duals,objective,functional


def analytic_joint_audit(args,progress):
    """Primal certificate including the imported analytic row enclosures."""
    from . import read_cflat,current_data,UVConfig
    from dataclasses import asdict
    from .operators import midpoint_grid
    from .kernels import pv_matrix
    from .model import chiral_slices,unitarity_margin,current_kinematic_squares
    from .imaginary import density_fourth_power
    ctx.prec=args.bits;meta=read_json(args.preparation/'report.json')['amplitude_model']
    M,L=meta['M'],meta['L']
    if meta['prescription']!='analytic-cardinal':raise ValueError('Analytic preparation required')
    c,_=read_cflat(args.coefficients,M,prescription='analytic-cardinal');p=len(c)
    current=current_data(args.current_preparation,args.preparation)
    from .endpoints import complete_imf,validate_endpoint_model
    endpoint=validate_endpoint_model(args,current);order=current['metadata'].get('ff_endpoint_order',0)
    if current['metadata']['uv']['config']!=asdict(UVConfig()):
        raise ValueError('This analytic audit implements the frozen printed/raw/hard original inputs')
    prior=read_json(args.coefficients.parent/'report.json')
    if prior['chiral_norm']!=args.chiral_norm or prior['density_limit']!=args.density_limit or prior['chiral_tolerance']!=args.chiral_tolerance:
        raise ValueError('Analytic audit constraints differ from candidate')
    with np.load(args.coefficients.parent/'joint.npz') as data:point=data['point'].copy()
    original_point=point.copy()
    if not np.array_equal(point[:p],c):raise ValueError('Joint NPZ and coefficients disagree')
    with np.load(args.preparation/'amplitude.npz') as data:
        H=data['rows'];E=data['radius_upper'];sample_energies,sample_waves=prepared_sampling(data,M,L);n=len(sample_energies)
    if H.shape!=E.shape or H.shape!=(2*n+11,p) or not np.isfinite(E).all() or np.any(E<0):
        raise ValueError('Complete outward analytic coefficient bounds required')
    retained,normal,duals,objective,functional=replay_support_data(args,prior,M,n,len(current['arrays']['high_energy_indices']))
    if objective is not None:
        from shutil import copyfile
        copyfile(args.coefficients.parent/'objective.npz',args.output/'objective.npz')
    cc=arb_mat(p,1,list(map(arb,c)));values=[]
    for lo in range(0,len(H),64):
        B=arb_mat([[arb(float(v),float(e)) for v,e in zip(h,r)] for h,r in zip(H[lo:lo+64],E[lo:lo+64])])
        values.extend((B*cc).entries());progress(stage='analytic_primal_rows',completed=min(lo+64,len(H)),total=len(H))
    s,w=midpoint_grid(M);K=pv_matrix(M);s0=(arb('1.2')/arb('.14'))**2
    f=[acb(values[j],values[n+j]) for j in range(n)]
    exact_energies=[s[j//(3*L)] if j<3*M*L else arb(float(x)) for j,x in enumerate(sample_energies)]
    margins=[unitarity_margin(v,x) for v,x in zip(f,exact_energies)]
    chi=[arb(str(args.chiral_tolerance))**2-sum((values[2*n+i]**2 for i in range(g.start,g.stop)),arb(0)) for g in chiral_slices(args.chiral_norm)]
    density=arb(str(args.density_limit))**4-density_fourth_power(cc.entries(),M)
    im=[complete_imf(list(map(arb,point[p+ch*M:p+(ch+1)*M])),endpoint,order or 1) for ch in range(2)]
    rho=[list(map(arb,point[p+(2+ch)*M:p+(3+ch)*M])) for ch in range(2)]
    FF=[];gm=[];fm=[];gram_pass=[];gram_methods=[];lifts=[];mq=(arb('.004')+arb('.0073'))/(2*arb('.14'))
    caps=[2*mq**2*arb('.00006'),arb('.00003')];high=[j for j,x in enumerate(s) if x>s0]
    for ch in range(2):
        vals=[];channel=[]
        for j,x in enumerate(s):
            v=acb(1+sum((K[j,b]*im[ch][b] for b in range(M)),arb(0)),im[ch][j])
            k2=current_kinematic_squares(x)[ch];v*=k2.sqrt();vals.append(v)
            S=1+acb(0,arb.pi()*(1-4/x).sqrt())*f[j*3*L+ch*L]
            delta=unitarity_margin(f[j*3*L+ch*L],x);floor=None
            if delta>0:
                floor=abs(v)**2+abs(v-S*v.conjugate())**2/delta
                if j in high and not rho[ch][j]>=floor:
                    before=float(rho[ch][j]);after=max(before,float(np.nextafter(float(floor.upper()),np.inf)))
                    rho[ch][j]=arb(after);point[p+(2+ch)*M+j]=after
                    lifts.append(dict(channel=ch,node=j,before=before,after=after,Schur_floor=floor.str(25)))
            G=arb_mat([[1+S.real,S.imag,arb(2).sqrt()*v.real],
                [S.imag,1-S.real,arb(2).sqrt()*v.imag],
                [arb(2).sqrt()*v.real,arb(2).sqrt()*v.imag,rho[ch][j]]])
            q=[G[a,a] for a in range(3)]+[G[a,a]*G[b,b]-G[a,b]**2 for a in range(3) for b in range(a)]+[G.det()]
            channel.append(q)
            principal=all(t>=0 for t in q);schur=floor is not None and delta>0 and rho[ch][j]>=floor
            gram_pass.append(bool(principal or schur));gram_methods.append('principal_minors' if principal else 'strict_Schur' if schur else 'unresolved')
            if j in high:fm.append(caps[ch]-abs(v)**2)
        FF.append(vals);gm.append(channel)
    powers=(0,1,-1,0);targets=[arb('3.09e-8')*(arb('27.38')/2+arb('.61')),
        arb('3.09e-8')*arb('27.38')/3,arb('4.34e-6')*arb('13.26'),
        arb('4.34e-6')*(arb('13.26')/2-arb('.41'))]
    moments=[sum((arb.pi()*w[j]*s[j]**a*rho[t//2][j] for j in range(M) if s[j]<=s0),arb(0)) for t,a in enumerate(powers)]
    mm=[arb('.002')-abs(v-t*s0**(a+2)) for v,t,a in zip(moments,targets,powers)]
    checks=margins+chi+[density]+fm+mm
    if current['metadata'].get('asymptotic_unitarity'):
        from .endpoints import asymptotic_margins
        checks+=asymptotic_margins(c,M)
    passed=all(v>=0 for v in checks) and all(gram_pass) and (not endpoint or cc[0,0].is_zero())
    violations=[dict(index=j,s=float(sample_energies[j]),isospin=int(sample_waves[j,0]),ell=int(sample_waves[j,1]),margin=v.str(25)) for j,v in enumerate(margins) if v<0]
    result=dict(status='analytic_joint_primal_certified' if passed else 'analytic_joint_primal_violated' if violations else 'analytic_joint_primal_unresolved',
        analytic_primal_feasible=passed,M=M,L=L,amplitude_model=meta,
        strict_scattering_violations=violations,
        counts=dict(scattering=n,current_grams=2*M,chiral=len(chi),FESR=4,FF_caps=len(fm),density=1,endpoint_equalities=1+2*order if endpoint else 0,asymptotic_inequalities=5 if current['metadata'].get('asymptotic_unitarity') else 0),
        endpoint_completion='Exact T0=0 and algebraically derived FF endpoint entries; printed dependent floats are display values' if endpoint else None,
        minimum_scattering_margin=min(margins).str(25),chiral_margins=[v.str(25) for v in chi],
        chiral_margin_lower=[v.lower().str(25) for v in chi],density_margin=density.str(25),FESR_margins=[v.str(25) for v in mm],FF_margins=[v.str(25) for v in fm],
        Gram_minors=[[[v.str(25) for v in q] for q in channel] for channel in gm],Gram_certified=gram_pass,Gram_certificate_methods=gram_methods,
        high_spectrum_lifts=lifts,Gram_criterion='Principal minors or strict elastic block and rigorous Schur floor',
        scope='All retained nodes: analytic row enclosures, exact native kinematics/K, printed FESR decimals and hard-midpoint weights; no continuum certification',
        amplitude_coefficients_modified=False,current_spectra_modified=bool(lifts),optimization_performed=False,
        inputs={str(q.resolve()):dict(sha256=digest(q)) for q in
            (args.preparation/'report.json',args.preparation/'amplitude.npz',args.coefficients,args.coefficients.parent/'joint.npz')})
    if endpoint:
        from .endpoints import asymptotic_margins
        tail=asymptotic_margins(c,M);result['asymptotic_unitarity']=dict(necessary_margins=[v.str(25) for v in tail],strict_eventual_conditions_passed=all(v>0 for v in tail),strict_obstruction=any(v<0 for v in tail),scope='Leading large-s coefficients for each fixed partial wave; T0=0 required, no finite crossover bound')
    result.update(objective_kind=prior.get('objective_kind','projection'),functional=functional,source_coefficients=str(args.coefficients))
    write_json(args.output/'source_audit.json',result)
    if passed:
        from .certificates import joint_audit
        from . import joint_result,watson_lineage
        from .gauge import center_acceptance
        kap=np.pi*np.sqrt(1-4/sample_energies)
        args.resolved_prescription='analytic-cardinal'
        point,finite=joint_audit(H,kap,current,point,M,L,args.density_limit,args.chiral_tolerance,
            normal,duals,args.bits,chiral_norm=args.chiral_norm,objective=objective)
        if not finite['primal_feasible']:raise ArithmeticError('Analytic lift failed the original finite constraints')
        identity=center_acceptance(original_point,point,M,current)
        method=dict(solver='Analytic Schur extension of unbounded high spectra',support_search=False,
            support_direction=normal,all_amplitude_directions_retained=True,original_variables=p+4*M,
            source_support_report=str(args.coefficients.parent/'report.json') if retained else None,
            center_acceptance=identity,representative_center_converged=bool(retained and identity['eligible'] and prior['solver'].get('representative_center_converged') is True))
        saved=joint_result(args,H,kap,current,point,duals,finite,M,L,args.density_limit,method,objective,functional)
        result=dict(saved,**{k:v for k,v in result.items() if k not in saved})
        result.update(status='analytic_joint_primal_certified',support_optimality_certified=bool(retained and saved['support_optimality_certified']),
            support_scope='Replayed float64 operator bound; analytic primal certified separately')
        if objective is not None:
            from .model import joint_observable_change
            with np.load(Path(functional['source_coefficients']).parent/'joint.npz') as old:previous=old['point']
            write_json(args.output/'watson_update.json',joint_observable_change(H,kap,current,previous,point,M,L))
            watson_lineage(args,np.asarray(result['targets']),functional,write_json,read_json,result)
    return result
