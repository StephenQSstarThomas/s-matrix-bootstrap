"""Coordinates and dual supports."""
from flint import arb, ctx
from scipy import sparse
import numpy as np
from . import UVConfig
from .linear import raw_support_residual

def barrier_terms(z, R, I, t, C, epsilon, bound, free, hessian=True, state=None, inverse_kappa=1.,chiral_norm='separate-l2'):
    """SOC/density Schur barrier."""
    from .model import chiral_slices
    x, y = ((R@np.asarray(z, dtype=np.longdouble), I@np.asarray(z, dtype=np.longdouble))
        if state is None else (state['x'],state['y']))
    slack = 2*inverse_kappa*y-x*x-t*y*y if state is None else state['slack']
    rho = np.asarray(z[free:], dtype=np.longdouble)/bound
    a = rho*rho; remaining = 1-np.sum(a*a)
    if not np.isfinite(slack).all() or not np.isfinite(remaining) or min(slack) <= 0 or remaining <= 0: return None
    chi = C@np.asarray(z, dtype=np.longdouble) if epsilon is not None else np.empty(0)
    groups=chiral_slices(chiral_norm) if epsilon is not None else ()
    cs=np.asarray([epsilon**2-np.sum(chi[s]**2) for s in groups],dtype=np.longdouble)
    if state is not None: chi,cs=state['chi'],state['chiral_slack']
    if len(cs)!=len(groups):raise ValueError('Cached chiral grouping mismatch')
    if cs.size and (not np.isfinite(cs).all() or min(cs)<=0):return None
    delta = remaining/(1+len(rho)/2)
    for _ in range(12):
        root = np.sqrt(a*a+2*delta); v = (a+root)/2
        change = (delta+v@v-1)/(1+np.sum(v/root))
        delta -= change
        if abs(change) < 1e-16*delta: break
    v = (a+np.sqrt(a*a+2*delta))/2
    d = delta/(2*v)
    value = -np.log(slack).sum()-np.log(cs).sum()-np.log(d).sum()-np.log(delta)
    result = dict(value=value, x=x, y=y, slack=slack, chi=chi, chiral_slack=cs,
        density_gradient=4*rho*v/delta)
    if not hessian: return result
    S = (2*x/slack)[:,None]*R-(2*(inverse_kappa-t*y)/slack)[:,None]*I
    curved=np.any(R!=0,axis=1); quadratic=t>0
    blocks = [S, np.sqrt(2/slack[curved])[:,None]*R[curved],
              np.sqrt(2*t[quadratic]/slack[quadratic])[:,None]*I[quadratic]]
    rhs = [-np.ones(len(t)), np.zeros(np.sum(curved)+np.sum(quadratic))]
    for g,s in enumerate(groups):
        row=2*chi[s]@C[s]/cs[g]
        blocks.extend([row[None,:],np.sqrt(2/cs[g])*C[s]])
        rhs.extend([-np.ones(1),np.zeros(s.stop-s.start)])
    diagonal = 2/d+8*a/(delta+2*d*d)
    rank = -4*rho*v/(delta+2*d*d)/np.sqrt(1+np.sum(delta/(delta+2*d*d)))
    density_rows = np.zeros((len(rho)+1,len(z)),dtype=np.longdouble)
    density_rows[:-1,free:] = np.diag(np.sqrt(diagonal)/bound)
    density_rows[-1,free:] = rank/bound
    blocks.append(density_rows)
    rhs.append(np.r_[-result['density_gradient']/np.sqrt(diagonal),0])
    result.update(features=np.vstack(blocks), rhs=np.concatenate(rhs))
    return result

def fesr_targets(config=UVConfig(), *, normalized=True):
    from .model import FESR_POWERS
    pi, s0 = np.pi, config.s0
    if config.moment_source == 'printed':
        values = np.array([3.09e-8*(27.38/2+.61), 3.09e-8*27.38/3,
                           4.34e-6*13.26, 4.34e-6*(13.26/2-.41)])
    else:
        g = config.gluon_condensate_gev4/config.pion_mass_gev**4
        j = config.scalar_condensate_gev4/config.pion_mass_gev**4
        scalar = 2*config.mq_squared/(2*pi)**4
        values = np.array([
            scalar*(3*(1+13*config.alpha_s/(3*pi))/(4*pi*(n+2))
                    + (pi*g/4+3*pi*j)/s0**2*(n == 0)) for n in (0, 1)] + [
            (1+config.alpha_s/pi)/(4*pi*(2*pi)**4*(n+2))
            - (pi*g/6+2*pi*j)/(2*(2*pi)**4*s0**2)*(n == 0) for n in (-1, 0)])
    return values if normalized else values*s0**(np.array(FESR_POWERS)+2)

def _fesr_grid(M, config):
    if type(M) is not int or M < 1:
        raise ValueError('A positive integer node count is required')
    step = np.pi/M
    phi = (np.arange(M)+.5)*step
    s = 4/np.cos(phi/2)**2
    jacobian = s*np.tan(phi/2)  # ds/dphi; no Cauchy 1/pi in a FESR.
    if config.cutoff == 'hard-midpoint':
        widths = step*(s <= config.s0)
    else:
        phi0 = 2*np.arccos(2/np.sqrt(config.s0))
        widths = np.clip(phi0-np.arange(M)*step, 0, step)
    return s, widths, jacobian

def joint_hull_candidate(args,H,current,M,L,bound,progress):
    from . import read_cflat, AMPLITUDE_COORDINATES
    from .model import chiral_slices
    import clarabel
    from pathlib import Path
    from . import read_json
    norm=getattr(args,'chiral_norm','separate-l2');groups=chiral_slices(norm);frozen=getattr(args,'fixed_amplitude',False);kind=getattr(args,'resolved_prescription','pv-midpoint')
    if kind not in AMPLITUDE_COORDINATES:raise ValueError('Unknown amplitude prescription')
    expected=dict(M=M,L=L,density_limit=bound,prescription=kind)
    other=getattr(args,'interior_coefficients',None)
    if other is not None:
        fixed=getattr(args,'fixed_x',None);paths=[args.coefficients,other];points=[]
        if fixed is None or any(args.direction):raise ValueError('Joint segment requires fixed x')
        for path in paths:
            r=read_json(path.parent/'report.json')
            if not r.get('joint_feasible') or r.get('status')=='inconclusive':raise ValueError('Require joint-feasible endpoints')
            if r.get('chiral_norm','separate-l2')!=norm or any(r[k]!=v for k,v in (expected|dict(chiral_tolerance=args.chiral_tolerance)).items()) or Path(r['parameters']['current_preparation']).resolve()!=args.current_preparation.resolve():raise ValueError('Same joint model required')
            with np.load(path.parent/'joint.npz') as d:points.append(d['point'].astype(np.longdouble))
        x=np.asarray([H[-2].astype(np.longdouble)@v[:H.shape[1]] for v in points]);weight=(fixed-x[1])/(x[0]-x[1])
        if not 0<weight<1:raise ValueError('Section must lie between endpoints')
        point=np.asarray(weight*points[0]+(1-weight)*points[1],float)
        return point,dict(solver='Convex segment of audited joint points',phase_I=True,support_search=False,
            source_coefficients=[str(p) for p in paths],weights=[float(weight),float(1-weight)],original_variables=len(point))
    data=read_json(args.region_summary);sig=data['model_signature'];a=current['arrays'];m=current['matrices'];p=current['metadata']['amplitude_coefficients']
    if any(sig[k]!=v for k,v in (expected|dict(infinity='free',unitarity_scope='sampled')).items()) or Path(sig['preparation']).resolve()!=args.preparation.resolve():raise ValueError('Scattering model mismatch')
    paths=[str(args.coefficients)] if frozen else list(dict.fromkeys(v['coefficients'] for g in data['regions'] for v in g['points'] if v['feasible']))
    if not paths:raise ValueError('No feasible parent amplitudes')
    ref=getattr(args,'coefficients',None)
    if not frozen and ref is not None and (ref.parent/'report.json').exists():
        source=read_json(ref.parent/'report.json')
        if source.get('joint_feasible') and str(ref.resolve()) not in paths:
            if any(source[k]!=v for k,v in expected.items()):raise ValueError('Joint seed model mismatch')
            paths.append(str(ref.resolve()))
    parents=np.column_stack([read_cflat(path,M,prescription=kind)[0] for path in paths]);nc=len(paths);bounded=np.asarray(m['moment_linear'][:,p+2*M:].getnnz(axis=0))>0;active=np.r_[np.arange(2*M),2*M+np.flatnonzero(bounded)];nv=nc+len(active)
    lam=np.vstack((np.ones(M),np.sqrt((a['nodes']-4)/(np.sqrt(a['nodes'])+2)**2)));scales=np.r_[lam.ravel(),a['k_squared'].ravel()]
    E=sparse.block_diag((sparse.csc_matrix(parents),sparse.diags(scales,format='csc')[:,active]),format='csc');D=np.stack((np.ones((2,M)),1/lam,1/a['k']),axis=-1).reshape(-1,3)[bounded]
    factors=np.array([[d[0]**2,d[0]*d[1],d[1]**2,d[0]*d[2],d[1]*d[2],d[2]**2] for d in D]).ravel();safety=0. if frozen else 1e-6
    rows=(6*np.flatnonzero(bounded)[:,None]+np.arange(6)).ravel();G=sparse.csc_matrix(m['gram_linear'][rows].astype(np.longdouble)@E.astype(np.longdouble),dtype=float);gb=a['gram_constant'][rows]*factors;gb.reshape(-1,6)[:,[0,2,5]]-=safety
    if frozen:
        from flint import arb_mat
        ids=a['native_scattering_indices'].ravel();n=(len(H)-11)//2
        with ctx.workprec(384):
            f=np.array([float(v) for v in (arb_mat(H[np.r_[ids,n+ids]].tolist())*arb_mat(parents.tolist())).entries()]).reshape(2,-1)
        re,im=f*np.tile(np.pi*np.sqrt(1-4/a['nodes']),2);zero=np.zeros(2*M)
        G=G.tolil();G[:,0]=np.stack((-im,np.sqrt(2)*re,im,zero,zero,zero),axis=1).ravel()[rows,None];G=G.tocsc()
    W=m['moment_linear']@E;target=a['moment_targets_raw'];error=a['moment_errors_raw'];ms=np.maximum(abs(target),error)
    high=len(a['high_energy_indices']);fs=np.repeat(1/np.sqrt(a['ff_caps_squared']),3*high);fb=a['ff_cap_constant']*fs;fb[0::3]-=safety
    simplex=sparse.csc_matrix(np.r_[np.ones(nc),np.zeros(len(active))][None,:]);positive=sparse.hstack((-sparse.eye(nc),sparse.csc_matrix((nc,len(active)))))
    A=sparse.vstack((-sparse.diags(factors)@G,sparse.diags(1/ms)@W,-sparse.diags(1/ms)@W,-sparse.diags(fs)@m['ff_cap_linear']@E,simplex,positive),format='csc')
    b=np.r_[gb,(target+error)/ms-safety,(-target+error)/ms-safety,fb,1.,np.zeros(nc)]
    cones=[clarabel.PSDTriangleConeT(3) for _ in range(sum(bounded))]+[clarabel.NonnegativeConeT(8)]+[clarabel.SecondOrderConeT(3) for _ in range(2*high)]+[clarabel.ZeroConeT(1),clarabel.NonnegativeConeT(nc)]
    C=np.asarray(np.asarray(H[-11:-3],np.longdouble)@parents.astype(np.longdouble),float)/args.chiral_tolerance
    chi_rows=[]
    for g in groups:
        size=1+len(C[g]);chi_rows.append(len(b));A=sparse.vstack((A,sparse.hstack((sparse.csc_matrix(np.vstack((np.zeros(nc),-C[g]))),sparse.csc_matrix((size,len(active)))))),format='csc');b=np.r_[b,1-safety,np.zeros(size-1)];cones.append(clarabel.SecondOrderConeT(size))
    cfg=clarabel.DefaultSettings();cfg.verbose=False;cfg.max_threads=8;cfg.direct_solve_method=args.native_backend;cfg.time_limit=args.solver_seconds;cfg.tol_feas=1e-11;cfg.tol_gap_abs=cfg.tol_gap_rel=1e-10
    xy=np.asarray(np.asarray(H[-2:],np.longdouble)@parents.astype(np.longdouble),float);fixed=getattr(args,'fixed_x',None)
    if fixed is not None:A=sparse.vstack((A,sparse.csc_matrix(np.r_[xy[0],np.zeros(len(active))][None,:])),format='csc');b=np.r_[b,fixed];cones.append(clarabel.ZeroConeT(1))
    if getattr(args,'ray',False):
        c0=read_cflat(args.coefficients,M,prescription=kind)[0];target=np.asarray(H[-2:],np.longdouble)@c0.astype(np.longdouble);slope=target[1]/target[0]
        A=sparse.vstack((A,sparse.csc_matrix(np.r_[xy[1]-float(slope)*xy[0],np.zeros(len(active))][None,:])),format='csc');b=np.r_[b,0.];cones.append(clarabel.ZeroConeT(1))
    cost=np.r_[-np.asarray(getattr(args,'direction',[0,0]))@xy,np.zeros(len(active))]
    interior=not frozen and not any(getattr(args,'direction',[0,0]));size=nv+int(interior)
    if interior:
        A=sparse.hstack((A,sparse.csc_matrix((len(b),1))),format='lil');A[chi_rows,nv]=1
        A=sparse.vstack((A,sparse.csc_matrix(([-1.],([0],[nv])),shape=(1,size))),format='csc');b=np.r_[b,0.];cones.append(clarabel.NonnegativeConeT(1));cost=np.r_[cost,-1.]
    if frozen:
        end=6*sum(bounded)+8+6*high;b=b[:end]-A[:end,0].toarray().ravel();A=A[:end,1:];cost=cost[1:];size-=1;cones=cones[:sum(bounded)+1+2*high]
    sol=clarabel.DefaultSolver(sparse.csc_matrix((size,size)),cost,A,b,cones,cfg).solve();x=np.asarray(np.r_[1.,sol.x] if frozen else sol.x,np.longdouble)
    if not np.isfinite(x).all() or not frozen and max(x[:nc])<=0:raise ArithmeticError('No finite joint hull candidate')
    weights=np.ones(1) if frozen else np.maximum(x[:nc],0);weights/=sum(weights);x[:nc]=weights;point=np.asarray(E.astype(np.longdouble)@x[:nv],float)
    if frozen:
        from . import zero_joint_duals
        from .imaginary import symmetric
        d=zero_joint_duals(M,n,high);dual=np.asarray(sol.z);at=6*sum(bounded)
        d['gram'].reshape(2*M,3,3)[bounded]=symmetric(dual[:at].reshape(-1,6));d['gram_congruence'].reshape(2*M,3)[bounded]=D
        d['moment']=(dual[at:at+4]-dual[at+4:at+8])/ms
        d['ff']=-dual[at+8:at+8+6*high].reshape(2,high,3)[:,:,1:]/np.sqrt(a['ff_caps_squared'])[:,None,None]
        np.savez_compressed(args.output/'fiber_duals.npz',**d)
    progress(stage='joint_hull',status=str(sol.status),parents=nc,seconds=sol.solve_time,primal_residual=sol.r_prim)
    return point,dict(solver='Clarabel current SDP',phase_I=True,chiral_interior_fraction=float(x[-1]) if interior else None,native_status=str(sol.status),iterations=sol.iterations,seconds=sol.solve_time,native_primal_residual=sol.r_prim,native_dual_residual=sol.r_dual,numerical_inner_margin=safety,original_variables=p+4*M,solver_variables=size,source_coefficients=paths,weights=weights.astype(float).tolist(),support_search=False,selection='Original constraints audited')

def subtraction_rows(matrix,M):
    from .kernels import density_labels
    source=np.asarray(matrix);free=1+2*M
    if source.ndim!=2 or source.shape[1]!=len(density_labels(M)) or np.iscomplexobj(source):
        raise ValueError('Require full actual-density rows')
    with ctx.workprec(max(128,ctx.prec)):
        angles=[arb(2*j+1)/(4*M) for j in range(M)]
        shifts=[x.sin_pi()/(M*x.cos_pi()) for x in angles]
        b=np.array([np.longdouble(x.mid().str(40,radius=False)) for x in shifts])
    ii,jj=np.triu_indices(M);symmetric_weight=np.where(ii==jj,np.longdouble('.5'),1)
    result=np.empty(source.shape,dtype=np.longdouble if source.dtype==np.longdouble else np.float64)
    for lo in range(0,len(source),128):
        block=np.array(source[lo:lo+128],dtype=np.longdouble,copy=True)
        h0=block[:,0];h1=block[:,1:1+M].copy();h2=block[:,1+M:free].copy()
        block[:,1:1+M]=h1-h0[:,None]*b
        block[:,1+M:free]=h2-2*h0[:,None]*b
        rho1=block[:,free:free+M*M].reshape(-1,M,M)
        rho1+=-2*h1[:,:,None]*b[None,None,:]-b[None,:,None]*h2[:,None,:]+2*h0[:,None,None]*b[:,None]*b[None,:]
        block[:,free+M*M:]+=symmetric_weight*(-h2[:,ii]*b[jj]-h2[:,jj]*b[ii]+2*h0[:,None]*b[ii]*b[jj])
        result[lo:lo+len(block)]=block
    return result

def subtract_amplitude_coordinates(z,M,to_subtracted):
    from .kernels import density_labels
    values=np.asarray(z)
    if values.ndim!=1 or len(values)!=len(density_labels(M)) or np.iscomplexobj(values):
        raise ValueError('Real full actual-density coefficient vector required')
    with ctx.workprec(max(128,ctx.prec)):
        coefficients=[]
        for value in values:
            if isinstance(value,arb):coefficients.append(arb(value))
            else:
                numerator,denominator=np.longdouble(value).as_integer_ratio()
                coefficients.append(arb(numerator)/denominator)
        ii,jj=np.triu_indices(M);off=1+2*M+M*M+np.flatnonzero(ii!=jj)
        for index in off:coefficients[index]*=2
        converted=convert_subtraction(coefficients,M,to_subtracted=to_subtracted)
        for index in off:converted[index]/=2
        return np.array([np.longdouble(value.mid().str(40,radius=False)) for value in converted],dtype=np.longdouble)

def absorptive_change(matrix,point,M,L):
    from .kernels import density_labels
    A=np.asarray(matrix);n=(len(A)-11)//2;p=A.shape[1];free=1+2*M;j=np.arange(M)
    if p!=len(density_labels(M)) or n<3*M*L:raise ValueError('Require full native coordinates')
    anchors=np.r_[n+j*3*L,n+j*3*L+2*L];C=np.asarray(A[anchors],np.longdouble)
    expected=np.zeros((2*M,free));expected[j,1+j]=1.5;expected[j,1+M+j]=1.;expected[M+j,1+M+j]=1.
    if not np.array_equal(C[:,:free],expected):raise ValueError('Native imaginary identities required')
    third=np.longdouble(2)/3;inverse=sparse.eye(p,format='lil',dtype=np.longdouble)
    inverse[1:free,:]=0;inverse[1+j,1+j]=third;inverse[1+j,1+M+j]=-third
    inverse[1+M+j,1+M+j]=1
    inverse[1:1+M,free:]=-third*(C[:M,free:]-C[M:,free:])
    inverse[1+M:free,free:]=-C[M:,free:];inverse=inverse.tocsc()
    out=np.empty_like(A,dtype=np.result_type(A.dtype,np.float64))
    for lo in range(0,len(A),128):out[lo:lo+128]=np.asarray((inverse.T@np.asarray(A[lo:lo+128],np.longdouble).T).T,out.dtype)
    # Exact coordinate identities of the stored operator.
    out[anchors]=0;out[anchors,1+np.arange(2*M)]=1
    z=np.asarray(point,np.longdouble).copy();z[1:free]=C@z
    return out,z,inverse

def convert_subtraction(coefficients, M, *, to_subtracted):
    from .operators import midpoint_grid
    from .kernels import coefficient_blocks
    values,R1,R2=coefficient_blocks(coefficients,M)
    nodes, weights = midpoint_grid(M)
    b = [w/x for x, w in zip(nodes, weights)]
    single1, single2 = values[1:1+M], values[1+M:1+2*M]
    sign = 1 if to_subtracted else -1
    values[0] += sign*sum((b[i]*(single1[i]+2*single2[i]) for i in range(M)), arb(0))
    values[0] += sum(((2*R1[i,j]+R2[i,j])*b[i]*b[j] for i in range(M) for j in range(M)), arb(0))
    values[1:1+M] = [single1[i]+sign*2*sum((R1[i,j]*b[j] for j in range(M)), arb(0)) for i in range(M)]
    values[1+M:1+2*M] = [single2[i]+sign*sum(((R1[j,i]+R2[i,j])*b[j] for j in range(M)), arb(0)) for i in range(M)]
    return values

def normal_dual(H,kappa_rows,M,L,bound,direction,epsilon,coefficients,tail,progress,seconds,*,infinity_zero=False,method="highs-ds",chiral_norm='separate-l2'):
    import time
    from scipy.optimize import linprog
    from .kernels import density_labels
    from .model import chiral_slices
    started=time.monotonic();H=np.asarray(H,float);n=(len(H)-11)//2
    if infinity_zero or tail is not None:raise ValueError('Require sampled/free')
    density=M*M+M*(M+1)//2;full=1+2*M+density;offset=0;free=1+2*M
    kap=np.asarray(kappa_rows,float);kap=kap if len(kap)==n else np.repeat(kap,3*L)
    coeff=np.asarray(coefficients,np.longdouble);direction=np.asarray(direction,float)
    if n<3*M*L or H.shape!=(2*n+11,full) or coeff.shape!=(full,) or kap.shape!=(n,) or direction.shape!=(2,):
        raise ValueError('Complete H/C/kappa/direction required')
    if not all(np.isfinite(v).all() for v in (H,coeff,kap,direction)) or np.any(kap<=0):
        raise ValueError('Require finite inputs, positive kappas')
    if not np.isfinite(bound) or bound<=0 or not np.isfinite(seconds) or seconds<=0:
        raise ValueError('Positive density bound/time required')
    if epsilon is not None and (not np.isfinite(epsilon) or epsilon<=0):raise ValueError('Positive chiral tolerance required')
    density_scales=np.array([2. if f=='rho2' and i!=j else 1. for f,i,j in density_labels(M)[1+2*M:]])
    values=np.empty(len(H),np.longdouble)
    for lo in range(0,len(H),128):values[lo:lo+128]=H[lo:lo+128]@coeff
    nr=np.asarray(kap*values[:n],float);ni=np.asarray(kap*values[n:2*n]-1,float)
    chi=np.asarray(values[2*n:2*n+8],float);groups=chiral_slices(chiral_norm) if epsilon is not None else ();nc=len(groups)
    nt=0;axis_start=n+nc
    axis_rows=np.repeat(np.r_[np.arange(M)*3*L,np.arange(M)*3*L+2*L],2)
    signs=np.tile([1.,-1.],2*M);axis_real=np.zeros(4*M,dtype=bool)
    axis_rows=np.r_[axis_rows,[(M//2)*3*L]*2];signs=np.r_[signs,1.,-1.];axis_real=np.r_[axis_real,True,True]
    count=axis_start+len(axis_rows)
    rows=np.empty((count,full-offset));costs=np.empty(count);scales=np.empty(count);exponents=np.empty(count,dtype=int)
    def store_rows(lo,block,cost):
        block=np.asarray(block,np.longdouble);block[:,free:]*=density_scales
        maximum=np.max(np.abs(block),axis=1)
        exponent=np.where(maximum>0,-np.frexp(maximum)[1],0)
        rows[lo:lo+len(block)]=np.ldexp(block,exponent[:,None])
        costs[lo:lo+len(block)]=np.ldexp(np.asarray(cost,np.longdouble),exponent)
        scales[lo:lo+len(block)]=np.ldexp(np.ones(len(block)),exponent)
        exponents[lo:lo+len(block)]=exponent
    for lo in range(0,n,128):
        hi=min(n,lo+128);a=np.asarray(nr[lo:hi],np.longdouble);b=np.asarray(ni[lo:hi],np.longdouble)
        hyp=np.hypot(a,b);cost=hyp+b
        np.divide(a*a,hyp-b,out=cost,where=b<0)
        block=(a*kap[lo:hi])[:,None]*H[lo:hi,offset:]+(b*kap[lo:hi])[:,None]*H[n+lo:n+hi,offset:]
        store_rows(lo,block,cost)
    for g,s in enumerate(groups):
        block=np.asarray(chi[s],np.longdouble)@H[2*n+s.start:2*n+s.stop,offset:]
        store_rows(n+g,block[None,:],[epsilon*np.linalg.norm(np.asarray(chi[s],np.longdouble))])
    for j,(row,sign,real) in enumerate(zip(axis_rows,signs,axis_real)):
        store_rows(axis_start+j,(np.longdouble(sign)*kap[row]*H[row if real else n+row,offset:])[None,:],
            [1. if real else 1.+sign])
    if not all(np.isfinite(v).all() for v in (rows,costs,scales)):raise ArithmeticError('Nonfinite normal scaling')
    target=(direction@H[-2:])[offset:].copy();target[free:]*=density_scales
    residual=sparse.vstack((sparse.csc_matrix((free,density)),sparse.eye(density,format='csc')),format='csc')
    matrix=sparse.hstack((sparse.csc_matrix(rows.T),residual,-residual),format='csc')
    prices=np.r_[costs,np.full(2*density,bound)]
    progress(stage='normal_dual_lp',physical_normals=n,chiral_normals=nc,tail_normals=nt,axis_normals=len(axis_rows),
        balanced_coefficients=full-offset,retained_density_directions=density,LP_variables=matrix.shape[1])
    lp_started=time.monotonic()
    if method=='highs-primal':
        from . import solve_feasible_basis_lp
        positive=axis_start+2*np.arange(free)
        signed=np.linalg.solve(rows[positive,:free].T,target[:free]);remainder=target[free:]-signed@rows[positive,free:]
        basic=np.r_[positive+(signed<0),count+np.arange(density)+density*(remainder<0)]
        initial=np.zeros(count+2*density);initial[basic]=np.r_[np.abs(signed),np.abs(remainder)]
        result=solve_feasible_basis_lp(prices,matrix,target,basic,initial,seconds)
    else:
        result=linprog(prices,A_eq=matrix,b_eq=target,bounds=(0,None),method=method,options={
            'time_limit':float(seconds),'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    meta=dict(solver=method+' on fixed normals',status=int(result.status),message=result.message,
        LP_seconds=time.monotonic()-lp_started,LP_success=bool(result.success),infinity_zero=bool(infinity_zero),
        physical_normals=n,chiral_normals=nc,tail_normals=nt,axis_normals=len(axis_rows),balanced_coefficients=full-offset,
        retained_density_directions=density,coefficient_directions_truncated=0,chiral_norm=chiral_norm,
        scaling='dyadic normal scaling; actual-rho L1 residual',support_certified=False)
    trial=result.x;fallback=trial is None or not np.isfinite(trial).all()
    if fallback:
        signed=np.linalg.solve(rows[axis_start::2,:free].T,target[:free])
        trial=np.zeros(count+2*density)
        trial[axis_start+2*np.arange(free)+(signed<0)]=np.abs(signed)
        remainder=target[free:]-signed@rows[axis_start::2,free:]
        trial[count:count+density]=np.maximum(remainder,0);trial[count+density:]=np.maximum(-remainder,0)
    gamma=np.maximum(trial[:count],0);weights=gamma*scales
    kr,ki=weights[:n]*nr,weights[:n]*ni;yy=np.zeros(8)
    for g,s in enumerate(groups):yy[s]=weights[n+g]*chi[s]
    ta,tb=np.empty(0),np.empty(0)
    np.add.at(kr,axis_rows[axis_real],weights[axis_start:][axis_real]*signs[axis_real])
    np.add.at(ki,axis_rows[~axis_real],weights[axis_start:][~axis_real]*signs[~axis_real])
    if not all(np.isfinite(v).all() for v in (kr,ki,yy,ta,tb)):raise ArithmeticError('Nonfinite dual weights')
    active=np.flatnonzero(gamma>0)
    meta.update(LP_objective=None if result.fun is None else float(result.fun),candidate_support_cost=float(prices@trial),
        axis_fallback=fallback or getattr(result,'used_initial_fallback',False),active_physical=int(np.sum(active<n)),active_axes=int(np.sum(active>=axis_start)),
        active_chiral=int(np.sum((active>=n)&(active<n+nc))),active_tail=int(np.sum((active>=n+nc)&(active<axis_start))),
        active_density_residual=int(np.sum(trial[count:]>0)),active_normal_rows=active.tolist(),
        active_row_scale_exponents=exponents[active].tolist(),negative_weights_clipped=int(np.sum(trial[:count]<0)),
        LP_equality_residual_max=float(np.max(np.abs(matrix@trial-target))),seconds=time.monotonic()-started)
    meta.update(initial_basis_supplied=method=='highs-primal',simplex_iterations=getattr(result,'nit',None),
        initial_objective=getattr(result,'initial_objective',None),primal_infeasibilities=getattr(result,'primal_infeasibilities',None),solver_version=getattr(result,'solver_version',None))
    progress(stage='normal_dual_candidate',**meta)
    return kr,ki,yy,{'kR':ta.tolist(),'kI':tb.tolist()},meta
