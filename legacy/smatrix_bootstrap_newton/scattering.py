"""Feasible Newton support paths for the scattering-only constraints."""
from concurrent.futures import ThreadPoolExecutor
from scipy import linalg
import numpy as np
import time
from . import extended_product, extended_sparse_blocks, extended_triangular

def barrier_support(A, kappa, M, L, bound, direction, epsilon, z, progress, seconds=600, start_mu=1e-2, gap=1e-4, qr_reuse=6, raw=False, tail=None, keep_extended=False, capture=None, chiral_norm='separate-l2', fixed_x=None, require_center=False, chiral_weight=1.,asymptotic=None):
    from . import support_gradient_width,radial_coordinates
    from scipy.linalg import qr
    from .quotient import barrier_terms
    from .model import chiral_slices
    groups=chiral_slices(chiral_norm)
    from .endpoints import augment_asymptotic,asymptotic_slope,asymptotic_duals
    n, free = (A.shape[0]-11)//2, A.shape[1]-(M*M+M*(M+1)//2)
    kap = np.asarray(kappa,float)
    if len(kap)!=n:kap=np.repeat(kap,3*L)
    if len(kap) != n or free not in (2*M, 2*M+1): raise ValueError('Physical rows or free coordinates disagree')
    if not raw and free!=2*M and not (free==2*M+1 and getattr(asymptotic,'bounded_constant',False)):raise ValueError('Require proved bounded absorptive coordinates')
    re, im = kap[:,None]*A[:n], kap[:,None]*A[n:2*n]
    t = np.minimum(1., bound*np.linalg.norm(im[:,free:],ord=4/3,axis=1))
    t[np.arange(M)*3*L]=1; t[np.arange(M)*3*L+2*L]=1
    t[t<=0]=1.  # Zero-density rows still define disks.
    inverse_kappa, dual_scale = (1/np.asarray(kap,np.longdouble), np.asarray(kap,np.longdouble)) if raw else (np.ones(n,np.longdouble), 1.)
    R, I, C = (np.asarray(A[:n] if raw else re,np.longdouble)/np.sqrt(np.asarray(t[:,None],np.longdouble)),
        np.asarray(A[n:2*n] if raw else im,np.longdouble)/np.asarray(t[:,None],np.longdouble), np.asarray(A[2*n:2*n+8],np.longdouble))
    if tail is not None:raise ValueError('Analytic sine tails are retired; use the complete PV sampled problem')
    t=np.asarray(t,np.longdouble)
    if not np.isfinite(chiral_weight) or chiral_weight<=0:raise ValueError('Positive chiral barrier weight required')
    def weighted_terms(*a,**kw):
        result=barrier_terms(*a,**kw)
        if result is None:return None
        result['value']-=(chiral_weight-1)*np.log(result['chiral_slack']).sum()
        if 'features' in result and epsilon is not None:
            at=n+int(np.sum(np.any(R!=0,axis=1)))+int(np.sum(t>0))
            count=sum(1+g.stop-g.start for g in groups);scale=np.sqrt(chiral_weight)
            result['features'][at:at+count]*=scale;result['rhs'][at:at+count]*=scale
        return result if asymptotic is None else augment_asymptotic(asymptotic,result,a[0],'features' in result)
    c = np.asarray(direction,np.longdouble)@np.asarray(A[-2:],np.longdouble)
    z = np.asarray(z, dtype=np.longdouble)
    section=fixed_x is not None;X=np.asarray(A[-2],np.longdouble);v=-X[1:]/X[0] if section else None
    base_direction=np.asarray(direction,float)
    if section and abs(X@z-fixed_x)>1e-10:raise ValueError('A feasible seed on the requested section is required')
    pool=ThreadPoolExecutor(max_workers=8) if len(z)>1000 else None
    started, iterations, history = time.monotonic(), 0, []
    terms = weighted_terms(z,R,I,t,C,epsilon,bound,free,False,inverse_kappa=inverse_kappa,chiral_norm=chiral_norm)
    if terms is None: raise ArithmeticError('The coherent analytic seed must be strictly feasible')
    state,best=terms,None;centered=None
    T, age, qr_count, cg_count, dec = None, qr_reuse, 0, 0, np.inf
    def precondition(v, transpose=False):
        return (extended_triangular(T,scales*v,True) if transpose else
                scales*extended_triangular(T,v))
    relative_width=8*np.finfo(np.longdouble).eps
    def radial_center(z,state,mu):
        x,y,chi=state['x'],state['y'],state['chi']
        qa,qb=2*inverse_kappa*y,x*x+t*y*y
        cn=np.asarray([np.sum(chi[g]**2) for g in groups]) if epsilon is not None else np.empty(0)
        def at(a):
            shifted=dict(x=a*x,y=a*y,slack=a*(qa-a*qb),chi=a*chi,
                chiral_slack=epsilon**2-a*a*cn if epsilon is not None else np.empty(0))
            future=weighted_terms(a*z,R,I,t,C,epsilon,bound,free,False,shifted,inverse_kappa=inverse_kappa,chiral_norm=chiral_norm)
            if future is None:return None,np.inf
            slope=-np.sum((qa-2*a*qb)/future['slack'])+chiral_weight*np.sum(2*a*cn/future['chiral_slack'])
            slope+=future['density_gradient']@(z[free:]/bound)-c@z/mu
            if asymptotic is not None:slope+=asymptotic_slope(asymptotic,future,z)
            return future,slope
        lo,hi,best=np.longdouble(0),np.longdouble(1),None;future,slope=at(hi)
        while future is not None and slope<0:
            lo,best=hi,future;hi*=2
            if not np.isfinite(hi):raise ArithmeticError('Radial minimum was not finitely bracketed')
            future,slope=at(hi)
        for _ in range(70):
            mid=(lo+hi)/2;future,slope=at(mid)
            if future is None or slope>=0:hi=mid
            else:lo,best=mid,future
            if hi-lo<=relative_width*hi:break
        if best is None:raise ArithmeticError('A positive interior radial step was not resolved')
        return lo,lo*z,best
    mu_factor=.5 if not raw else .1
    for mu in [start_mu*mu_factor**j for j in range(100) if start_mu*mu_factor**j >= 1e-11]:
        age=qr_reuse;central_gradient=central_step=None
        while time.monotonic()-started<=seconds:
            terms = weighted_terms(z,R,I,t,C,epsilon,bound,free,state=state,inverse_kappa=inverse_kappa,chiral_norm=chiral_norm)
            F, b, method = terms['features'], terms['rhs'], 'pcg'
            full_features=F
            cost=c
            if section:F=F[:,1:]+F[:,0,None]*v;cost=c[1:]+c[0]*v
            at=n+int(np.sum(np.any(R!=0,axis=1)))+int(np.sum(t>0))
            radial=at+np.r_[0,np.cumsum([1+g.stop-g.start for g in groups])[:-1]] if epsilon is not None else []
            F,cost,restore=radial_coordinates(F,cost,radial,free-int(section))
            products=extended_sparse_blocks(F)
            if T is None or age >= qr_reuse or dec > 1:
                scales=1/np.sqrt(np.einsum('ij,ij->j',F,F))
                T=qr(np.asarray(F,float)*np.asarray(scales,float),mode='r',overwrite_a=True,check_finite=False)[0][:F.shape[1]].copy()
                age, qr_count, method = 0, qr_count+1, 'qr'
            adjoint=lambda v:extended_product(products,v,True,pool)
            gradient=np.asarray(cost,np.longdouble)/mu+adjoint(b)
            right=precondition(gradient,True)
            step=np.zeros_like(right); residual=right.copy(); search=residual.copy()
            rr=residual@residual; tolerance=rr*1e-14
            for krylov in range(64 if method=='qr' else 12):
                action=precondition(adjoint(extended_product(products,precondition(search),pool=pool)),True)
                denominator=search@action
                if not denominator>0:break
                length=rr/denominator; step+=length*search; residual-=length*action
                updated=residual@residual
                if updated<=tolerance:break
                search=residual+(updated/rr)*search; rr=updated
            dz=precondition(step); dec=float(gradient@dz)
            applied=extended_product(products,dz,pool=pool); curvature=float(applied@applied); cg_count+=krylov+1
            linear_residual=float(np.linalg.norm(precondition(adjoint(applied)-gradient,True))/max(np.linalg.norm(right),np.finfo(float).tiny))
            if linear_residual>1e-6:
                if method=='qr':
                    progress(stage='barrier_inexact_direction',mu=mu,iteration=iterations,linear_residual=linear_residual,decrement=dec,curvature=curvature)
                else:
                    progress(stage='barrier_refactor',mu=mu,linear_residual=linear_residual)
                    T=None;age=qr_reuse;continue
            dz=restore(dz)
            if section:dz=np.r_[v@dz,dz]
            near=0<=dec<1e-6 and 0<=curvature<1e-6 and linear_residual<=1e-6
            if near:
                full_gradient=c/mu+np.einsum('ij,i->j',full_features,b)
                if section:
                    full_gradient-=full_gradient[0]*X/X[0]
                residual=mu*full_gradient
                width=(support_gradient_width(A,kap,residual,M,L,free,bound) if raw else
                    2*np.sum(abs(residual[:free]))+2*bound*np.linalg.norm(residual[free:],ord=4/3))
                if width<=gap/4:
                    central_gradient=mu*full_gradient.copy();central_step=dz.copy()
                    progress(stage='central_converged',mu=mu,decrement=dec,curvature=curvature,linear_residual=linear_residual,gradient_width=width)
                    break
            age += 1
            dx,dy,dchi=R@dz,I@dz,C@dz
            x,y=terms['x'],terms['y']
            q1,q2=2*(inverse_kappa-t*y)*dy-2*x*dx,dx*dx+t*dy*dy
            chi1=np.asarray([-2*np.sum(terms['chi'][g]*dchi[g]) for g in groups]) if epsilon is not None else np.empty(0)
            chi2=np.asarray([np.sum(dchi[g]**2) for g in groups]) if epsilon is not None else np.empty(0)
            def on_ray(a):
                shifted=dict(x=x+a*dx,y=y+a*dy,slack=terms['slack']+a*q1-a*a*q2,
                    chi=terms['chi']+a*dchi if epsilon is not None else np.empty(0),
                    chiral_slack=terms['chiral_slack']+a*chi1-a*a*chi2)
                future=weighted_terms(z+a*dz,R,I,t,C,epsilon,bound,free,False,shifted,inverse_kappa=inverse_kappa,chiral_norm=chiral_norm)
                if future is None:return None,np.inf
                slope=-np.sum((q1-2*a*q2)/future['slack'])-chiral_weight*np.sum((chi1-2*a*chi2)/future['chiral_slack'])
                slope+=future['density_gradient']@(dz[free:]/bound)-(c@dz)/mu
                if asymptotic is not None:slope+=asymptotic_slope(asymptotic,future,dz)
                return future,slope
            if near:
                future,_=on_ray(np.longdouble(1))
                if future is None:raise ArithmeticError('Near-center Newton step left the strict interior')
                z,state=z+dz,future;iterations+=1
                progress(stage='central_newton_step',mu=mu,iteration=iterations,gradient_width=width,decrement=dec)
                continue
            _,g0=on_ray(0.)
            if not np.isfinite([g0,dec,curvature,linear_residual]).all() or not g0<0 or linear_residual>.1 or abs(curvature+g0)>.1*max(curvature,abs(g0)):
                progress(stage='non_descent' if not g0<0 else 'linear_precision_failure',mu=mu,directional_gradient=float(g0),linear_residual=linear_residual,directional_curvature=curvature,linear_solver=method)
                age=qr_reuse
                if method=='pcg':continue
                if capture is not None and hasattr(capture,'failure'):capture.failure(z,mu,state,linear_residual=linear_residual,decrement=dec,curvature=curvature,iteration=iterations)
                dec=np.inf;break
            lo,hi,ray_best=np.longdouble(0),np.longdouble(1),terms
            future,slope=on_ray(hi)
            while future is not None and slope<0 and hi<1e12:
                lo,ray_best=hi,future
                hi*=2;future,slope=on_ray(hi)
            for _ in range(100):
                alpha=(lo+hi)/2;future,slope=on_ray(alpha)
                if future is None or slope>=0:hi=alpha
                else:lo,ray_best=alpha,future
                if hi-lo<=relative_width*abs(hi):break
            alpha=lo
            line_decrease=terms['value']-ray_best['value']+alpha*(c@dz)/mu
            if linear_residual>1e-6 and not line_decrease>0:
                if capture is not None and hasattr(capture,'failure'):capture.failure(z,mu,state,linear_residual=linear_residual,decrement=dec,curvature=curvature,iteration=iterations)
                raise ArithmeticError('Inexact scattering direction failed actual barrier descent')
            z,state=z+alpha*dz,ray_best
            radial_scale=np.longdouble(1)
            if not section:radial_scale,z,state=radial_center(z,state,mu)
            iterations += 1
            progress(stage='barrier', mu=mu, iteration=iterations, objective=float(c@z),
                decrement=dec, step=float(alpha), elapsed=time.monotonic()-started, linear_solver=method,
                directional_gradient=float(g0),directional_curvature=curvature,linear_residual=linear_residual,krylov_iterations=krylov+1,radial_scale=float(radial_scale),line_merit_decrease=float(line_decrease),inexact_direction=linear_residual>1e-6)
            if method == 'pcg' and alpha < .01: age = qr_reuse
            if method == 'pcg' and alpha <= 1e-12 and time.monotonic()-started <= seconds: continue
            if dec < 1e-6 or alpha <= 1e-12 or time.monotonic()-started > seconds: break
        terms = state
        kr = 2*mu*terms['x'][:n]/(dual_scale*np.sqrt(t[:n])*terms['slack'][:n])
        ki = 2*mu*(t[:n]*terms['y'][:n]-inverse_kappa[:n])/(dual_scale*t[:n]*terms['slack'][:n])
        effective_c=np.asarray(c,np.longdouble).copy();extra_support=np.longdouble(0);extra_duals={}
        if 'asymptotic' in terms:
            q=terms['asymptotic'];u=2*mu*q['x']/q['slack'];v_tail=-2*mu/q['slack']
            effective_c-=u@asymptotic.asymptotic_R+v_tail@asymptotic.asymptotic_I
            extra_support=-np.sum(u*u/(2*v_tail));extra_duals=asymptotic_duals(asymptotic,terms,mu)
        yy = (np.concatenate([2*mu*chiral_weight*terms['chi'][g]/slack for g,slack in zip(groups,terms['chiral_slack'])])
              if epsilon is not None else np.zeros(8))
        normal=base_direction.copy();lam=np.longdouble(0)
        if section:
            r0=effective_c[0]-kr@re[:,0]-ki@im[:,0]-yy@C[:,0]
            lam=-r0/X[0];effective_c=effective_c+lam*X;normal[0]+=float(lam)
        if capture:kr0,ki0=kr.copy(),ki.copy()
        if raw:
            from .quotient import raw_support_residual
            kr,ki,density_residual = raw_support_residual(A,kap,kr,ki,yy,effective_c,M,L,free)
            free_bound = 0.
        else:
            residual = effective_c-kr@re-ki@im-yy@C
            density_residual = residual[free:]
            free_bound = 2*np.maximum(residual[:2*M],0).sum()
            if free == 2*M+1: free_bound += abs(residual[2*M])
        support = np.hypot(kr,ki)
        disk = np.where(ki<0,kr*kr/(support-ki),support+ki).sum()
        chi_support = (epsilon*sum(np.linalg.norm(yy[g]) for g in groups) if epsilon is not None else 0)
        outer = disk+chi_support+free_bound+bound*np.linalg.norm(density_residual,ord=4/3)+extra_support
        conditional=outer-lam*fixed_x if section else outer
        if capture:capture(mu=mu,z=z,t=t,c=c,support_direction=normal,fixed_x=np.nan if fixed_x is None else float(fixed_x),inverse_kappa=inverse_kappa,state=terms,kR_raw=kr0,kI_raw=ki0,kR=kr,kI=ki,dual_y=yy,rho_residual=density_residual,mu_gradient=central_gradient,newton_step=central_step,outer_estimate=float(conditional),
            iteration=iterations,decrement=dec,curvature=curvature if central_gradient is not None else np.nan,linear_residual=linear_residual if central_gradient is not None else np.nan,**extra_duals)
        point = dict(mu=mu,objective=float(c@z),outer_estimate=float(conditional),gap_estimate=float(conditional-c@z))
        history.append(point); progress(stage='central_point',**point)
        if central_gradient is not None:centered=z.copy()
        if best is None or conditional<best[0]: best=(conditional,kr.copy(),ki.copy(),yy.copy(),normal.copy(),extra_duals)
        chosen=centered if require_center else z
        if chosen is not None and float(best[0]-c@chosen)<gap or time.monotonic()-started>seconds or dec>1e-3:break
    if pool is not None:pool.shutdown()
    if require_center and centered is not None:z=centered
    return np.asarray(z,np.longdouble if keep_extended else float),np.asarray(best[1],float),np.asarray(best[2],float),np.asarray(best[3],float),dict(
        solver='analytic density barrier with feasible Newton steps',iterations=iterations,seconds=time.monotonic()-started,
        mu_reduction=mu_factor,chiral_barrier_weight=float(chiral_weight),radial_chiral_preconditioned=True,support_direction=best[4].tolist(),fixed_x=fixed_x,require_center=require_center,representative_center_converged=centered is not None if require_center else None,fixed_x_deviation=None if not section else float(X@z-fixed_x),
        qr_factorizations=qr_count,pcg_steps=cg_count,qr_reuse=qr_reuse,raw_coordinates=raw,matrix_threads=8 if pool else 1,
        krylov_storage='longdouble CSR exact nonzeros; dense QR',
        tail_weights={'kR':[],'kI':[]},tail_constraints_applied=False,asymptotic_duals=best[5],
        history=history,feasible_path=True,all_density_directions_retained=True)


def select_chiral_support(args):
    """Lock a converged fixed-section IR support before inspecting phases."""
    from pathlib import Path
    from . import read_json,write_json,read_cflat,digest
    from .ir import audit_ir,center_identity
    from .sampling import prepared_sampling
    if len(args.support_runs)!=1:raise ValueError('One fixed-section IR support required')
    path=Path(args.support_runs[0]);r=read_json(path/'report.json');p=r['parameters']
    if r['mode']!='chiral' or not r.get('sampled_feasible') or not r.get('support_optimality_certified') or not r['solver'].get('representative_center_converged'):
        raise ValueError('Feasible supported and converged IR representative required')
    if r['chiral_norm']!=args.chiral_norm or r['chiral_tolerance']!=args.chiral_tolerance:
        raise ValueError('Different chiral selection inputs')
    xref=5/(16*np.pi**2*(92/140)**2);direction=r['support_direction'];M,L=r['M'],r['L']
    if p.get('fixed_x')!=xref or direction[1]<=0:raise ValueError('Preregistered Weinberg upper section required')
    preparation=Path(p['preparation']).resolve();c,_=read_cflat(path/'coefficients.json',M,prescription=r['prescription'])
    for name in ('report.json','amplitude.npz'):
        source=preparation/name;expected=r.get('inputs',{}).get(str(source))
        if expected and expected['sha256']!=digest(source) or r['infinity']=='zero' and not expected:raise ValueError('IR preparation provenance missing or changed')
    with np.load(preparation/'amplitude.npz') as data:H=data['rows'];ss,_=prepared_sampling(data,M,L)
    n=len(ss);kap=np.pi*np.sqrt(1-4/ss);tail=bool(r.get('tail_conditions_applied',False))
    from .kernels import density_labels
    ds=np.array([2. if f=='rho2' and i!=j else 1. for f,i,j in density_labels(M)[1+2*M:]])
    center=center_identity(path,c,M,ds)
    with np.load(path/'candidate.npz') as data:
        if not np.array_equal(c,data['coefficients']):raise ValueError('IR candidate/JSON mismatch')
        duals={k:data[k].copy() for k in ('kR','kI','y','asymptotic_kR','asymptotic_kI') if k in data}
    audit=audit_ir(H,kap,c,M,L,r['density_limit'],args.chiral_tolerance,direction,duals,
        zero=r['infinity']=='zero',asymptotic=tail,bits=args.bits,chiral_norm=args.chiral_norm)
    if not audit['primal_feasible']:raise ArithmeticError('Saved IR amplitude fails original constraints')
    target=np.asarray(H[-2:],np.longdouble)@c
    if abs(target[0]-xref)>1e-10:raise ValueError('IR C is off the preregistered section')
    distance=(audit['upper']-audit['lower'])/direction[1]
    if distance/args.chiral_tolerance>.01:raise ValueError('Local boundary needs refinement')
    sig={k:r[k] for k in ('M','L','density_limit','infinity','unitarity_scope','prescription','chiral_norm','scattering_samples')}
    sig.update(preparation=str(preparation),global_high_spin_enforced=False,tail_conditions_applied=tail,additional_constraints=None,sampling_nodes=len(np.unique(ss)))
    if n>3*M*L:sig['additional_constraints']=dict(kind='explicit-analytic-samples',preparation=str(preparation))
    selection=dict(rule='Converged upper IR support at preregistered Weinberg x before phase evaluation',reference=[xref,-xref/15],
        target=list(map(float,target)),center_identity=center,coefficient_sha256=digest(path/'coefficients.json'),epsilon=args.chiral_tolerance,model_signature=sig,source_reports=[str((path/'report.json').resolve())],
        coefficients_mixed=False,selection_uses_phase_data=False,representative_center_converged=True,
        exact_extremizer_claimed=False,regional_geometry_ready=False,local_upper_boundary_ready=True,
        local_upper_boundary=dict(distance_upper=distance,scaled_distance_upper=distance/args.chiral_tolerance,budget=.01),original_H_feasibility=audit)
    write_json(args.output/'coefficients.json',read_json(path/'coefficients.json'));write_json(args.output/'selection.json',selection)
    return dict(status='selected',selection=selection,model_signature=sig,sampled_feasible=True,optimization_performed=False,
        inputs={str((path/f).resolve()):dict(sha256=digest(path/f)) for f in ('report.json','coefficients.json')})


def center_capture(output,M,coordinates,prescription,changed,density_scales,raw_map=None,numerical_coordinates=None):
    """Save Newton centers without changing their numerical coordinates."""
    from . import write_json
    from .quotient import subtract_amplitude_coordinates
    free=1+2*M
    numerical_coordinates=numerical_coordinates or ('absorptive-zero-T0' if raw_map is not None else 'subtracted' if changed else 'unsubtracted')
    center_best=[float('inf')]
    def raw(z):
        if raw_map is not None:return raw_map(z)
        v=np.asarray(subtract_amplitude_coordinates(z,M,False) if changed else z,float).copy();v[free:]*=density_scales
        return v
    def failure(z,mu,state,**diagnostics):
        point=raw(z)
        np.savez_compressed(output/'failed_newton_state.npz',z=z,mu=mu,coefficients=point,initialization_only=True,center_converged=False,M=M,prescription=prescription,numerical_coordinates=numerical_coordinates,**diagnostics,
            **{k:state[k] for k in ('x','y','slack','chi','chiral_slack')})
        write_json(output/'failed_newton_coefficients.json',dict(coordinates=coordinates,prescription=prescription,coefficients=point.tolist(),initialization_only=True,center_converged=False,center_mu=mu))
    def capture(**data):
        state=data.pop('state');g=data['mu_gradient'];data.update(has_newton_gradient=g is not None,mu_gradient=np.zeros_like(data['z']) if g is None else g,newton_step=np.zeros_like(data['z']) if g is None else data['newton_step'],
            **{key:state[key] for key in ('x','y','slack','chi','chiral_slack','density_gradient')},M=M,numerical_coordinates=numerical_coordinates,density_packing='actual_rho',prescription=prescription)
        data['coefficients']=raw(data['z'])
        if data['has_newton_gradient']:np.savez_compressed(output/'center_converged.npz',**data)
        first=output/'center_first.npz'
        if not first.exists():
            np.savez_compressed(first,**data);v=raw(data['z'])
            write_json(output/'center_first_coefficients.json',dict(coordinates=coordinates,prescription=prescription,coefficients=v.tolist(),initialization_only=True,center_mu=data['mu']))
        if data['outer_estimate']<center_best[0]:center_best[0]=data['outer_estimate'];np.savez_compressed(output/'center_best.npz',**data)
        np.savez_compressed(output/'center_last.npz',**data)
        last=raw(data['z'])
        write_json(output/'center_last_coefficients.json',dict(coordinates=coordinates,prescription=prescription,coefficients=last.tolist(),
            initialization_only=True,center_mu=data['mu'],center_converged=data['has_newton_gradient']))
    capture.failure=failure
    return capture
