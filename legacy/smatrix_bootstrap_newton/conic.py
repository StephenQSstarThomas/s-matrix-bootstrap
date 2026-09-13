"""Full-variable cone embedding of the existing finite joint constraints.

This is a feasibility solver, not a hull restriction or representative rule.
Returned points require the original audit, independently of native status.
"""
from pathlib import Path
import time
import numpy as np
from scipy import sparse


def compile_joint_cones(P):
    """Reuse joint maps; lift the actual-density L4 ball with exact SOC identities."""
    import clarabel
    nv=len(P.active);nd=P.p-P.free;total=nv+nd
    fixed=np.asarray(getattr(P,'endpoint_pivots',[]),dtype=int)
    if np.any((fixed>=P.free)&(fixed<P.p)):raise ValueError('Density coordinates must survive the affine slice')
    offset=np.zeros(total);offset[fixed]=np.asarray(P.initial[fixed],float)
    scale=np.ones(total);scale[P.free:P.p]=P.bound/nd**.25
    columns=np.delete(np.arange(total),fixed);D=sparse.diags(scale[columns],format='csc')
    matrices=[];constants=[];cones=[];blocks={};row=0
    def pad(A):
        A=sparse.csc_matrix(A,dtype=float)
        return sparse.hstack((A,sparse.csc_matrix((A.shape[0],total-A.shape[1]))),format='csc')
    def add(name,A,b,types):
        nonlocal row
        A=pad(A);b=np.asarray(b,float).ravel()+A@offset
        matrices.append(-(A[:,columns]@D));constants.append(b);cones.extend(types)
        blocks[name]=[row,row+len(b)];row+=len(b)
    def parabola(name,R,I,t):
        R,I=pad(R),pad(I);n=R.shape[0];root=np.sqrt(2.)
        if t is None:
            A=sparse.vstack((I,root*R,-I),format='csc');width=3;b=[1,0,1]
        else:
            A=sparse.vstack((I,root*R,sparse.diags(np.sqrt(2*np.asarray(t,float)))@I,-I),format='csc');width=4;b=[1,0,0,1]
        order=np.arange(width*n).reshape(width,n).T.ravel()
        add(name,A[order],np.tile(b,n),[clarabel.SecondOrderConeT(width) for _ in range(n)])
    parabola('scattering',P.R,P.I,P.t)
    for j,g in enumerate(P.groups):
        C=np.asarray(P.C[g],float)/P.epsilon
        add(f'chiral_{j}',np.vstack((np.zeros(P.p),C)),np.r_[1.,np.zeros(len(C))],[clarabel.SecondOrderConeT(1+len(C))])
    # (1+w_i,2*r_i,w_i-1) in SOC3 and (1,w) in SOC(nd+1).
    i=np.arange(nd);r=np.r_[3*i,3*i+1,3*i+2];c=np.r_[nv+i,P.free+i,nv+i]
    A=sparse.csc_matrix((np.r_[np.ones(nd),np.full(nd,2/P.bound),np.ones(nd)],(r,c)),shape=(3*nd,total))
    add('density_squares',A,np.tile([1.,0.,-1.],nd),[clarabel.SecondOrderConeT(3) for _ in range(nd)])
    A=sparse.csc_matrix((np.ones(nd),(1+i,nv+i)),shape=(nd+1,total))
    add('density_norm',A,np.r_[1.,np.zeros(nd)],[clarabel.SecondOrderConeT(nd+1)])
    add('gram',P.J.reshape(-1,nv),P.g.ravel(),[clarabel.PSDTriangleConeT(3) for _ in P.keep])
    count=len(P.fc);F=sparse.csc_matrix(P.F.reshape(2*count,nv),dtype=float)
    order=np.arange(3*count).reshape(3,count).T.ravel()
    A=sparse.vstack((sparse.csc_matrix((count,nv)),F[::2],F[1::2]),format='csc')[order]
    add('ff',A,np.column_stack((np.ones(count),P.fc)).ravel(),[clarabel.SecondOrderConeT(3) for _ in range(count)])
    add('moment',np.vstack((-P.W,P.W)),np.r_[1+P.target,1-P.target],[clarabel.NonnegativeConeT(8)])
    if hasattr(P,'asymptotic_R'):parabola('asymptotic',P.asymptotic_R,P.asymptotic_I,None)
    A=sparse.vstack(matrices,format='csc');A.eliminate_zeros();b=np.concatenate(constants)
    if not np.isfinite(A.data).all() or not np.isfinite(b).all():raise ArithmeticError('Nonfinite full conic embedding')
    return dict(A=A,b=b,cones=cones,offset=offset,scale=scale,columns=columns,blocks=blocks,
        joint_variables=nv,auxiliary_density_variables=nd,all_amplitude_directions_retained=True)


def scs_layout(Q):
    """SCS orders linear/SOC/PSD cones and packs PSD lower, not upper, triangles."""
    blocks=Q['blocks'];lo,hi=blocks['moment'];order=list(range(lo,hi));dims=[]
    widths=dict(scattering=4,density_squares=3,ff=3,asymptotic=3)
    for name,(lo,hi) in blocks.items():
        if name in ('moment','gram') or hi==lo:continue
        width=widths.get(name,hi-lo)
        if (hi-lo)%width:raise ValueError('Incomplete SOC blocks')
        dims.extend([width]*((hi-lo)//width));order.extend(range(lo,hi))
    lo,hi=blocks['gram'];count=(hi-lo)//6
    order.extend((lo+6*np.arange(count)[:,None]+[0,1,3,2,4,5]).ravel().tolist())
    if len(order)!=len(Q['b']):raise ValueError('Complete SCS cone ordering required')
    return np.asarray(order,dtype=int),dict(l=8,q=dims,s=[3]*count)


def conic_warm_point(P,Q,point):
    """Map a complete raw candidate to the same affine coordinates; warm data only."""
    from scipy.sparse.linalg import spsolve
    z=np.asarray(point[P.active],float).copy();z[:P.p]/=np.asarray(P.ds,float)
    if P.inverse is not None:z[:P.p]=spsolve(sparse.csc_matrix(P.inverse,dtype=float),z[:P.p])
    if hasattr(P,'endpoint_inverse'):
        z=spsolve(sparse.csc_matrix(P.endpoint_inverse,dtype=float),z)
        z[P.endpoint_pivots]=P.initial[P.endpoint_pivots]
    rho=z[P.free:P.p]/P.bound;w=rho*rho
    w+=max(0.,1-np.linalg.norm(w))/(2*np.sqrt(len(w)))
    return ((np.r_[z,w]-Q['offset'])/Q['scale'])[Q['columns']]


def solve_joint_cones(Q,args,warm,progress=lambda **kw:None):
    """Two standard algorithms for one saved conic matrix; no native certificates."""
    nv=Q['A'].shape[1];started=time.monotonic()
    if args.solver=='scs':
        import scs
        if args.native_backend!='auto':raise ValueError('SCS uses its declared indirect solver; require --native-backend auto')
        order,cones=scs_layout(Q)
        progress(stage='conic_native_setup_start',backend='scs',version=scs.__version__)
        logging={} if not getattr(args,'output',None) else dict(log_csv_filename=str(args.output/'conic_native_trace.csv'))
        solver=scs.SCS(dict(A=Q['A'][order].tocsc(),b=Q['b'][order],c=np.zeros(nv)),cones,
            use_indirect=True,verbose=True,eps_abs=1e-8,eps_rel=1e-8,max_iters=50000,time_limit_secs=args.solver_seconds,**logging)
        progress(stage='conic_native_setup_complete',backend='scs',seconds=time.monotonic()-started)
        warm=np.ascontiguousarray(warm,dtype=float);initial_slack=np.asarray(Q['b']-Q['A']@warm,float)
        sol=solver.solve(warm_start=True,x=warm,y=np.zeros(len(order)),s=initial_slack[order]);info=sol['info'];dual=np.empty(len(order));slack=np.empty(len(order))
        dual[order]=sol['y'];slack[order]=sol['s']
        native=dict(status=info['status'],iterations=info['iter'],seconds=time.monotonic()-started,
            primal_residual=info['res_pri'],dual_residual=info['res_dual'],version=scs.__version__,
            linear_solver='indirect',warm_start_used=True,raw_info=info,primal_infeasibility_certified=False)
        return np.asarray(sol['x']),dual,slack,native
    import clarabel
    settings=clarabel.DefaultSettings();settings.verbose=True;settings.max_threads=8
    settings.direct_solve_method=args.native_backend;settings.time_limit=args.solver_seconds;settings.max_iter=200
    settings.tol_feas=1e-11;settings.tol_gap_abs=settings.tol_gap_rel=1e-10
    progress(stage='conic_native_setup_start',backend='clarabel',version=clarabel.__version__)
    solver=clarabel.DefaultSolver(sparse.csc_matrix((nv,nv)),np.zeros(nv),Q['A'],Q['b'],Q['cones'],settings)
    progress(stage='conic_native_setup_complete',backend='clarabel',seconds=time.monotonic()-started);sol=solver.solve()
    native=dict(status=str(sol.status),iterations=sol.iterations,seconds=sol.solve_time,
        primal_residual=sol.r_prim,dual_residual=sol.r_dual,version=clarabel.__version__,warm_start_used=False,primal_infeasibility_certified=False)
    return np.asarray(sol.x),np.asarray(sol.z),np.asarray(sol.s),native


def initialize_conic(args,H,kap,current,seed,M,L,bound,progress):
    """Solve the original zero-objective joint problem; native status is not proof."""
    from . import read_json,write_json,digest,zero_joint_duals,resolution_report
    from .operators import JointProblem
    from .endpoints import configure_endpoint_coordinates,configure_asymptotic_constraints
    from .certificates import joint_audit
    if args.fixed_x is not None or args.ray:raise ValueError('Full conic feasibility currently requires an unsectioned problem')
    folder=args.coefficients.parent;reference=seed
    inputs=[args.coefficients,folder/'joint.npz'];record=folder/'report.json';cache=folder/'phase_I_state.npz'
    if cache.exists() and record.exists():
        old=read_json(record)
        if old.get('solver',{}).get('phase_I_model')=='current-only-diagonal-v2':
            if any(Path(old['parameters'][k]).resolve()!=getattr(args,k).resolve() for k in ('preparation','current_preparation')):
                raise ValueError('Reference congruence must come from the same joint operators')
            with np.load(cache) as saved:reference=saved['reference_point']
            inputs.extend((cache,record))
    reference=np.asarray(reference,float).copy();p=H.shape[1];bad=np.flatnonzero(reference[p+2*M:]<=0)
    reference[p+2*M+bad]=np.maximum(current['arrays']['k_squared'].ravel()[bad],np.finfo(float).tiny)
    P=JointProblem(H,kap,current,M,L,bound,args.chiral_tolerance,reference,[0,0],absorptive=True,chiral_norm=args.chiral_norm,retain_free_spectra=True)
    if current['metadata'].get('ff_endpoint_order'):configure_endpoint_coordinates(P,current['metadata']['ff_endpoint_order'])
    if current['metadata'].get('asymptotic_unitarity'):configure_asymptotic_constraints(P)
    started=time.monotonic();Q=compile_joint_cones(P);A=Q['A'];nv=A.shape[1]
    sparse.save_npz(args.output/'conic_A.npz',A)
    np.savez_compressed(args.output/'conic_arrays.npz',**{k:Q[k] for k in ('b','offset','scale','columns')})
    meta=dict(shape=A.shape,nonzeros=A.nnz,blocks=Q['blocks'],original_variables=P.p+4*M,
        active_joint_variables=Q['joint_variables'],auxiliary_density_variables=Q['auxiliary_density_variables'],
        fixed_endpoint_coordinates=list(map(int,getattr(P,'endpoint_pivots',[]))),
        density_column_scale=bound/(P.p-P.free)**.25,all_amplitude_directions_retained=True,
        eliminated_high_spectra=False,all_current_Gram_blocks_retained=len(P.keep)==2*M,congruence_reference_spectrum_replacements=bad.tolist(),
        scope='Algebraically equivalent cones; assembled floating coefficients still require original and analytic audits')
    write_json(args.output/'conic_model.json',meta);progress(stage='full_conic_compiled',**meta,seconds=time.monotonic()-started)
    warm=conic_warm_point(P,Q,seed);np.savez_compressed(args.output/'conic_warm.npz',x=warm)
    x,dual,slack,native=solve_joint_cones(Q,args,warm,progress)
    write_json(args.output/'conic_native.json',native)
    np.savez_compressed(args.output/'conic_candidate.npz',x=x,dual=dual,slack=slack)
    if x.shape!=(nv,) or not np.isfinite(x).all():raise ArithmeticError('Conic solver returned no finite full-variable candidate')
    extended=Q['offset'].copy();extended[Q['columns']]+=Q['scale'][Q['columns']]*x
    point=P.raw(np.asarray(extended[:len(P.active)],np.longdouble));duals=zero_joint_duals(M,len(kap),len(P.a['high_energy_indices']))
    point,audit=joint_audit(H,kap,current,point,M,L,bound,args.chiral_tolerance,[0,0],duals,args.bits,chiral_norm=args.chiral_norm)
    method=dict(solver='Full-variable '+args.solver+' conic feasibility',phase_I=True,phase_I_model='full-conic-feasibility-v1',support_search=False,
        all_amplitude_directions_retained=True,original_variables=P.p+4*M,solver_variables=nv,iterations=native['iterations'],
        seconds=time.monotonic()-started,native=native,physical_constraints_relaxed_at_acceptance=False,
        infeasibility_certified=False,representative_center_converged=False,warm_start_used=native['warm_start_used'],
        reference_used_only_for_coordinates_and_congruence=True)
    progress(stage='full_conic_original_audit',native_status=native['status'],original_primal_feasible=audit['primal_feasible'])
    result=resolution_report(args,H,kap,current,point,duals,audit,M,L,bound,method)
    result.update(support_optimality_certified=False,resolution_initializer=True,
        inputs={str(f.resolve()):dict(sha256=digest(f)) for f in inputs},conic_model=meta)
    return result
