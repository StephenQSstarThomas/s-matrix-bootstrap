"""Convex feasible recovery and current/scattering numerical algebra."""
from flint import arb, arb_mat, ctx
from fractions import Fraction
from pathlib import Path
import numpy as np

def restore_joint_segment(problem,candidate):
    from .operators import JointProblem
    base=JointProblem(problem.H,problem.kap,problem.current,problem.M,problem.L,problem.bound,problem.epsilon,problem.seed,problem.direction,chiral_norm=problem.chiral_norm)
    z0=base.initial;z1=np.asarray(candidate,np.longdouble)[base.active].copy();z1[:base.p]/=base.ds
    v0=base.values(z0);v1=base.values(z1);lo=np.longdouble(0);hi=np.longdouble(1)
    for _ in range(50):
        mid=(lo+hi)/2;v={k:v0[k]+mid*(v1[k]-v0[k]) for k in v0}
        if base.state(z0+mid*(z1-z0),v) is None:hi=mid
        else:lo=mid
    fraction=lo*(1-np.longdouble('1e-7'))
    return np.asarray((1-fraction)*problem.seed.astype(np.longdouble)+fraction*np.asarray(candidate,np.longdouble),float),float(fraction)

def chiral_ratio(energy, channel):
    energy = Fraction(energy)
    if channel == '01':
        return 3*(2*energy-1)/(energy-4)
    if channel == '21':
        return 3*(2-energy)/(energy-4)
    raise ValueError('Chiral channel must be01 or21')

def density_fourth_power(coefficients, M):
    if len(coefficients) != 1+2*M+M*M+M*(M+1)//2:
        raise ValueError('Complete original coefficient vector required')
    at = 1+2*M
    fourth=lambda v:(v*v)*(v*v)
    total = sum(fourth(v) for v in coefficients[at:at+M*M])
    at += M*M
    for i in range(M):
        for j in range(i, M):
            total += fourth(coefficients[at])/(1 if i == j else 16)
            at += 1
    return total

def watsonian_objective(H,kap,current,point,M,L):
    a=current['arrays'];p=H.shape[1];n=len(kap);nodes=a['nodes'];s0=(current['metadata']['uv']['config']['matching_energy_gev']/current['metadata']['uv']['config']['pion_mass_gev'])**2
    low=np.flatnonzero(nodes<=s0);im=np.asarray(point[p:p+2*M]).reshape(2,M)
    F=1+im@a['hilbert_kernel'].T+1j*im
    if np.any(abs(F[:,low])==0):raise ValueError('Watsonian phase undefined at a zero form factor')
    h=(np.asarray(H[:n],np.longdouble)@point[:p]+1j*(np.asarray(H[n:2*n],np.longdouble)@point[:p]))*kap
    weights=np.zeros((n,2),np.longdouble);ceiling=np.longdouble(0)
    for ch in (0,1,2):
        ids=low*3*L+ch*L
        Q=F[ch,low]/F[ch,low].conj() if ch<2 else 1+1j*h[ids]
        scale=(np.sqrt(nodes[low])+2)**2/(nodes[low]-4) if ch==1 else np.ones(len(low))
        weights[ids,0]=scale*Q.imag;weights[ids,1]=-scale*Q.real
        ceiling+=np.sum(scale*(abs(Q)-Q.real))
    count=3*len(low);weights/=count;ceiling/=count
    cost=np.asarray((weights[:,0]*kap)@np.asarray(H[:n],np.longdouble)+(weights[:,1]*kap)@np.asarray(H[n:2*n],np.longdouble),float)
    from . import zero_joint_duals
    duals=zero_joint_duals(M,n,len(a['high_energy_indices']));duals.update(kR=np.asarray(weights[:,0],float),kI=np.asarray(weights[:,1],float))
    return cost,duals,dict(source='2403.10772 Eq.(2.29), GTB_numerics.m 117-120,177-180',channels=['S0','P1','S2'],
        nodes=len(low),normalization=count,weights='Lambda_l(s)^-2; Lambda_0=1, Lambda_1=sqrt(s-4)/(sqrt(s)+2)',
        ideal_disk_ceiling=float(ceiling),initial_value=float(cost@point[:p]),projection_constraints_released=True)

def coherent_boundary(args,progress=lambda **kw:None):
    from .ir import coherent_boundary as run
    return run(args,progress)


def recover_candidate_segment(H,kappa,M,bound,epsilon,candidate,interior,tail,direction,certify,*,infinity_zero=False,bits=384,global_high_spin=False,prescription='pv-midpoint',chiral_norm='combined-l2'):
    from .model import chiral_slices
    from .kernels import density_labels
    if prescription not in ('pv-midpoint','analytic-cardinal') or global_high_spin or infinity_zero or tail is not None:raise ValueError('Active recovery is PV sampled/free')
    groups=chiral_slices(chiral_norm);c,s=[np.asarray(v,float) for v in (candidate,interior)];n=(len(H)-11)//2;free=1+2*M
    original=c.copy();radial=1.;radial_result=None
    if c.shape!=s.shape or c.shape!=(H.shape[1],) or not np.isfinite(c).all() or not np.isfinite(s).all():raise ArithmeticError('Finite complete candidate and interior required')
    def audit(point):
        outer=dict(certify(point),high_spin_necessary=dict(status='not_applicable',violations=[]));passed=outer['primal_feasible']
        if not passed:outer.update(primal_feasible=False,lower=None)
        return passed,outer
    objective=np.asarray(direction,np.longdouble)@np.asarray(H[-2:],np.longdouble)
    def result(point,fraction,outer,coarse):
        return point,outer,dict(candidate_fraction=float(fraction),interior_fraction=float(1-fraction),radial_scale=float(radial),candidate_weight=float(radial*fraction),zero_weight=float((1-radial)*fraction),
            objective_before=float(objective@original),objective_after=float(objective@point),objective_loss=float(objective@(original-point)),repaired=bool(fraction<1 or radial<1),coarse_feasible_fraction=float(coarse),primal_certified=True,global_high_spin_enforced=False,scope='Coefficient mixture on declared constraints; no continuum claim')
    def finish(recovered):
        point,outer,info=recovered;retained=info['candidate_fraction']==0;factor=1-1e-8
        if objective@s>objective@point:
            valid,verified=audit(s)
            if valid:point,outer,info=result(s,0.,verified,0.);retained=True
        contracted=np.asarray(factor*point,float);valid,outer=audit(contracted)
        if not valid:
            valid,verified=audit(s)
            if not valid:raise ArithmeticError('No certified incumbent for final contraction')
            point,_,info=result(s,0.,verified,0.);retained=True;contracted=np.asarray(factor*point,float);valid,outer=audit(contracted)
        if not valid:raise ArithmeticError('Final coefficient contraction failed independent verification')
        info.update(incumbent_retained=retained,incumbent_objective=float(objective@s),final_contraction=factor,candidate_weight=info['candidate_weight']*factor,interior_fraction=info['interior_fraction']*factor,
            objective_after=float(objective@contracted),objective_loss=float(objective@(original-contracted)),repaired=not np.array_equal(original,contracted));info['zero_weight']=1-info['candidate_weight']-info['interior_fraction']
        return contracted,outer,info
    passed,outer=audit(c)
    if passed:return finish(result(c,1.,outer,1.))
    endpoints=np.asarray(H,np.longdouble)@np.asarray([s,c],np.longdouble).T;kap=np.asarray(kappa,np.longdouble)
    if len(kap)!=n:raise ValueError('Segment audit requires one kappa per scattering row')
    density_scales=np.array([2. if f=='rho2' and i!=j else 1. for f,i,j in density_labels(M)[free:]])
    re,im=kap*endpoints[:n,1],kap*endpoints[n:2*n,1];square=re*re+im*im
    if np.all(im>=0):
        limits=[1.,float(np.min(2*im[square>0]/square[square>0])) if np.any(square>0) else 1.];norm=np.linalg.norm(c[free:]/density_scales,ord=4)
        if norm>0:limits.append(bound/norm)
        if epsilon is not None:
            norm=max(np.linalg.norm(endpoints[2*n:2*n+8,1][g]) for g in groups)
            if norm>0:limits.append(epsilon/norm)
        scale=min(limits)*(1-1e-8)
        if np.isfinite(scale) and scale>0:
            radial=float(scale);c=np.asarray(radial*original,float);passed,outer=audit(c)
            if passed:
                radial_result=result(c,1.,outer,1.)
                if radial>=1-1e-5:return finish(radial_result)
                c=original.copy();radial=1.
            endpoints[:,1]=np.asarray(H,np.longdouble)@np.asarray(c,np.longdouble)
    passed,start_outer=audit(s)
    if not passed:
        if radial_result is not None:return finish(radial_result)
        raise ArithmeticError('No certified segment interior is available')
    def better(recovered):return radial_result if radial_result is not None and radial_result[2]['objective_after']>recovered[2]['objective_after'] else recovered
    def point_at(a):return np.asarray((1-np.longdouble(a))*s+np.longdouble(a)*c,float)
    def coarse(a):
        value=(1-a)*endpoints[:,0]+a*endpoints[:,1];re,im=kap*value[:n],kap*value[n:2*n];point=point_at(a)
        return not (np.any(2*im-re*re-im*im<0) or np.linalg.norm(point[free:]/density_scales,ord=4)>bound or epsilon is not None and max(np.linalg.norm(value[2*n:2*n+8][g]) for g in groups)>epsilon)
    lo,hi=0.,1.
    for _ in range(50):
        mid=(lo+hi)/2
        if coarse(mid):lo=mid
        else:hi=mid
    fraction=lo*(1-1e-8)
    for _ in range(12):
        point=point_at(fraction);passed,outer=audit(point)
        if passed:return finish(better(result(point,fraction,outer,lo)))
        fraction*=.5
    return finish(better(result(s,0.,start_outer,lo)))

def symmetric(v):
    out=np.zeros((*v.shape[:-1],3,3),dtype=v.dtype)
    out[...,0,0]=v[...,0];out[...,1,1]=v[...,2];out[...,2,2]=v[...,5]
    for j,(r,s) in zip((1,3,4),((0,1),(0,2),(1,2))):
        out[...,r,s]=out[...,s,r]=v[...,j]/np.sqrt(np.longdouble(2))
    return out

def svec(v):
    root=np.sqrt(np.longdouble(2))
    return np.stack((v[...,0,0],root*v[...,0,1],v[...,1,1],root*v[...,0,2],root*v[...,1,2],v[...,2,2]),axis=-1)

def inverse_cholesky(G):
    if not np.isfinite(G).all() or np.any(G[:,0,0]<=0):return None
    a=np.sqrt(G[:,0,0]);b=G[:,1,0]/a;c=G[:,2,0]/a
    d2=G[:,1,1]-b*b
    if np.any(d2<=0):return None
    d=np.sqrt(d2);e=(G[:,2,1]-b*c)/d;f2=G[:,2,2]-c*c-e*e
    if np.any(f2<=0):return None
    f=np.sqrt(f2);U=np.zeros_like(G)
    U[:,0,0]=1/a;U[:,1,0]=-b/(a*d);U[:,1,1]=1/d
    U[:,2,0]=(b*e-c*d)/(a*d*f);U[:,2,1]=-e/(d*f);U[:,2,2]=1/f
    return U,2*(np.log(a)+np.log(d)+np.log(f)).sum()

def support_outer(H,kappa,kR,kI,y,direction,bound,epsilon=None,bits=256,infinity_zero=False,coefficients=None,tail=None,waves_per_isospin=None,polynomial_dual=None,joint=None,chiral_norm='separate-l2',objective=None):
    from .model import chiral_slices
    from .certificates import verify_scattering_layout
    ctx.prec=bits
    groups=chiral_slices(chiral_norm)
    if len(kappa)==len(kR) and waves_per_isospin is None:
        raise ValueError('Expanded kappa rows require explicit waves_per_isospin')
    n,p=len(kR),H.shape[1]; M=round((np.sqrt(1+24*p)-5)/6); L=waves_per_isospin or n//(3*len(kappa))
    K=[arb(float(v)) for v in (kappa if len(kappa)==n else np.repeat(kappa,3*L))]
    verify_scattering_layout(H,K,M,L)
    kr,ki=[[arb(float(v)) for v in a] for a in (kR,kI)]
    yy=[arb(float(v)) if epsilon is not None else arb(0) for v in y]
    weights=[-a*b for a,b in zip(kr,K)]+[-a*b for a,b in zip(ki,K)]
    goal=None if objective is None else [v if isinstance(v,arb) else arb(float(v)) for v in objective]
    if goal is not None and (len(goal)!=p or any(not v.is_finite() for v in goal)):raise ValueError('Finite complete C_flat objective required')
    weights += [-v for v in yy]+[arb(0)]+([arb(0),arb(0)] if goal is not None else [arb(float(v)) for v in direction])
    residual=arb_mat(1,p) if goal is None else arb_mat(1,p,goal)
    if joint is not None:residual+=joint[0]
    cc=None if coefficients is None else arb_mat(p,1,[arb(float(v)) for v in coefficients]);values=[]
    if tail is not None or polynomial_dual is not None:raise ValueError('Use explicit sampled constraints')
    tail_support,tail_margins=arb(0),[]
    for lo in range(0,len(H),64):
        block=H[lo:lo+64]
        matrix=arb_mat(len(block),p,[arb(float(v)) for v in block.ravel()])
        residual+=arb_mat(1,len(block),weights[lo:lo+len(block)])*matrix
        if cc is not None:values.extend((matrix*cc).entries())
    r=[residual[0,j] for j in range(p)]
    anchor=(M//2)*3*L
    correction=arb(0) if infinity_zero else r[0]/(K[anchor]*arb(float(H[anchor,0])))
    kr[anchor]+=correction
    r=[v-correction*K[anchor]*arb(float(a)) for v,a in zip(r,H[anchor])]
    for j in range(M):
        a,b=j*3*L,j*3*L+2*L
        da=r[1+j]/(arb(3)/2*K[a]); db=r[1+M+j]/K[b]-da
        ki[a]+=da;ki[b]+=db
        for q in range(1+2*M,p):
            r[q]-=da*K[a]*arb(float(H[n+a,q]))+db*K[b]*arb(float(H[n+b,q]))
    disk=arb(0)
    for a,b in zip(kr,ki):
        hyp=(a*a+b*b).sqrt()
        disk+=a*a/(hyp-b) if b<0 else hyp+b
    chi=arb(0) if epsilon is None else arb(float(epsilon))*sum(
        (sum((v*v for v in yy[s]),arb(0)).sqrt() for s in groups),arb(0))
    rho=r[1+2*M:];at=M*M
    for i in range(M):
        for j in range(i,M):
            if i!=j:rho[at]*=2
            at+=1
    norm=sum(((abs(v).upper()**4).root(3) for v in rho),arb(0)).sqrt().sqrt()**3
    density=arb(float(bound))*norm;upper=(disk+chi+density+tail_support+(arb(0) if joint is None else joint[1])).upper()
    result=dict(upper=float(np.nextafter(float(upper),np.inf)),enclosure=upper.str(25),
        disk=float(disk.upper()),chiral=float(chi.upper()),density_residual=float(density.upper()),tail_support=float(tail_support.upper()),
        chiral_norm=chiral_norm,objective_kind='C_flat' if goal is not None else 'target_plane',scope='Exact saved float64 H; angular quadrature error excluded')
    if cc is not None:
        margins=[2*K[i]*values[n+i]-K[i]*K[i]*(values[i]*values[i]+values[n+i]*values[n+i]) for i in range(n)]
        norms=[sum((v*v for v in values[2*n+s.start:2*n+s.stop]),arb(0)).sqrt() for s in groups]
        density_norm=density_fourth_power(cc.entries(),M).root(4)
        passed=all(v>=0 for v in margins+tail_margins) and density_norm<=bound and (not infinity_zero or cc[0,0].is_zero())
        if epsilon is not None:passed=passed and all(v<=arb(float(epsilon)) for v in norms)
        lower=(sum((a*b for a,b in zip(goal,cc.entries())),arb(0)) if goal is not None else sum((arb(float(a))*b for a,b in zip(direction,values[-2:])),arb(0))).lower()
        result.update(primal_feasible=bool(passed),lower=float(np.nextafter(float(lower),-np.inf)) if passed else None,
            minimum_margin_enclosure=min(margins).str(20),density_norm=float(density_norm.upper()),
            targets=[v.str(25) for v in values[-2:]],tail_margin_enclosures=[v.str(20) for v in tail_margins])
        result.update(chiral_norm_enclosures=[] if epsilon is None else [v.str(30) for v in norms],
            chiral_margin_enclosures=[] if epsilon is None else [(arb(float(epsilon))-v).str(30) for v in norms],
            chiral_strict_violation=False if epsilon is None else any(v>arb(float(epsilon)) for v in norms))
    return result
