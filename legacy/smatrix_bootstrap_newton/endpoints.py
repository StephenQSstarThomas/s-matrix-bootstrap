"""Declared affine endpoint completion of the finite analytic family."""
from functools import lru_cache
import numpy as np
from flint import arb,arb_mat,ctx
from scipy import sparse


@lru_cache(maxsize=8)
def endpoint_matrix(M,order,bits):
    from .analytic import cardinal_transform
    if order not in (1,2) or M<order:raise ValueError('Enough modes for the declared endpoint order required')
    with ctx.workprec(bits):
        b=[(arb(2*j+1)/(4*M)).tan_pi()/M for j in range(M)]
        rows=[b]
        if order==2:
            D=cardinal_transform(M)
            rows.append([sum((n*(-1)**(n-1)*D[n-1,j] for n in range(1,M+1)),arb(0)) for j in range(M)])
        E=arb_mat(rows);B=arb_mat([[E[i,j] for j in range(M-order,M)] for i in range(order)])
        if B.det().contains(0):raise ArithmeticError('Endpoint dependent columns not proved independent')
        return E


def complete_imf(values,enabled,order=1):
    """Last order entries are derived algebraic values, never independent floats."""
    v=list(map(arb,values))
    if not enabled:return v
    M=len(v);E=endpoint_matrix(M,order,ctx.prec);free=M-order
    B=arb_mat([[E[i,j] for j in range(free,M)] for i in range(order)])
    rhs=arb_mat(order,1,[arb(int(i==0))-sum((E[i,j]*v[j] for j in range(free)),arb(0)) for i in range(order)])
    return v[:free]+(B.solve(rhs)).entries()


def endpoint_residual(residual,matrix):
    """r.v = constant + reduced.v on E v = (1,0,...)."""
    E=matrix if isinstance(matrix,arb_mat) else arb_mat([matrix]);k,M=E.nrows(),E.ncols();free=M-k
    B=arb_mat([[E[i,j] for j in range(free,M)] for i in range(k)])
    lam=B.transpose().solve(arb_mat(k,1,residual[free:])).entries()
    q=[residual[j]-sum((lam[i]*E[i,j] for i in range(k)),arb(0)) for j in range(free)]+[arb(0)]*k
    return q,lam[0]


def validate_endpoint_model(args,current):
    order=current['metadata'].get('ff_endpoint_order',0);expected=2 if getattr(args,'asymptotic_zeros',False) else 0;enabled=bool(order)
    if order!=expected or (args.infinity=='zero')!=enabled:raise ValueError('Endpoint declaration or order differs from current preparation')
    if enabled and getattr(args,'resolved_prescription',getattr(args,'prescription',None))!='analytic-cardinal':raise ValueError('Endpoint completion requires the analytic family')
    if bool(current['metadata'].get('asymptotic_unitarity'))!=bool(getattr(args,'asymptotic_unitarity',False)):raise ValueError('Asymptotic unitarity declaration differs from current preparation')
    return enabled


def configure_endpoint_coordinates(P,order=2):
    """Invertible coordinate change; five fixed coordinates for double FF zeros."""
    p,M,nv=P.p,P.M,len(P.active);phase=len(P.initial)>nv
    E=endpoint_matrix(M,order,max(256,ctx.prec));F=np.array([[np.longdouble(E[i,j].mid().str(40,radius=False)) for j in range(M)] for i in range(order)])
    constant=np.zeros(nv,np.longdouble)
    constant[:p]=np.eye(1,p,0).ravel() if P.inverse is None else np.asarray(P.inverse.getrow(0).toarray()).ravel()
    eligible=np.arange(P.free)
    if P.inverse is not None:eligible=eligible[(eligible!=0)&(eligible!=P.target_slot)]
    pivot=int(eligible[np.argmax(abs(constant[eligible]))])
    constraints=[constant];targets=[0.];pivots=[pivot]
    for ch in range(2):
        for i in range(order):
            row=np.zeros(nv,np.longdouble);row[p+ch*M:p+(ch+1)*M]=F[i];constraints.append(row);targets.append(float(i==0));pivots.append(p+(ch+1)*M-1-i)
    inverse=sparse.eye(nv,format='csc',dtype=np.longdouble)
    for original,k in zip(constraints,pivots):
        row=np.asarray(original@inverse).ravel();value=row[k]
        if not np.isfinite(value) or value==0:raise ArithmeticError('Endpoint coordinate pivot failed')
        v=-row/value;v[k]=1/value;change=v.copy();change[k]-=1
        P.initial[k]=row@P.initial[:nv]
        for name in ('A','R','I','J','F','W'):
            a=getattr(P,name);width=a.shape[-1]
            if k<width:a+=a[...,k].copy()[...,None]*change[:width]
        P.cost[:nv]+=P.cost[k]*change
        V=sparse.eye(nv,format='lil',dtype=np.longdouble);V[k,:]=v;inverse=inverse@V.tocsc()
    P.C=P.A[2*P.n:2*P.n+8]
    P.endpoint_inverse=inverse;P.endpoint_pivots=np.asarray(pivots);P.endpoint_order=order
    P.initial[P.endpoint_pivots]=targets
    if phase:
        from .imaginary import symmetric
        v=P.values(P.initial);e=np.linalg.eigvalsh(np.asarray(symmetric(v['gram']),float))
        P.initial[-1]=max(0.,float(-e.min()),float(np.max(np.sum(v['ff']**2,axis=1)-1)),float(np.max(abs(v['moment'])-1)))+1


def completed_raw_point(P,z):
    values=P.endpoint_inverse@z if hasattr(P,'endpoint_inverse') else z
    point=np.zeros(P.p+4*P.M,np.longdouble);point[P.active]=values
    point[:P.p]=P.ds*(values[:P.p] if P.inverse is None else P.inverse@values[:P.p])
    result=np.asarray(point,float)
    if hasattr(P,'endpoint_inverse'):
        result[0]=0.
        with ctx.workprec(256):
            for ch in range(2):
                sl=slice(P.p+ch*P.M,P.p+(ch+1)*P.M)
                result[sl]=list(map(float,complete_imf(result[sl],True,P.endpoint_order)))
    return result

def asymptotic_margins(coefficients,M):
    """Necessary leading unitarity coefficients, conditional on exact T0=0."""
    from .kernels import coefficient_blocks
    E=endpoint_matrix(M,2,ctx.prec);a=arb_mat(M,1,[4*E[1,j] for j in range(M)])
    c,R,Q=coefficient_blocks(list(map(arb,coefficients)),M)
    x=(a.transpose()*arb_mat(M,1,c[1:M+1]))[0,0];y=(a.transpose()*arb_mat(M,1,c[M+1:2*M+1]))[0,0]
    r=(a.transpose()*R*a)[0,0];q=(a.transpose()*Q*a)[0,0]
    return [3*x+2*y,y,5*(4*r+q)-arb.pi()*(x+4*y)**2,3*(r-q)-arb.pi()*(x-y)**2,5*(r+q)-arb.pi()*(x+y)**2]

@lru_cache(maxsize=8)
def asymptotic_rows(M,bits):
    """Five exact parabolas 2 I.C - (R.C)^2 >= 0, conditional on T0=0."""
    with ctx.workprec(bits):
        p=1+2*M+M*M+M*(M+1)//2;E=endpoint_matrix(M,2,bits);a=[4*E[1,j] for j in range(M)]
        R=arb_mat(5,p);I=arb_mat(5,p);u=(arb.pi()/5).sqrt();v=(arb.pi()/3).sqrt()
        for j,x in enumerate(a):
            I[0,1+j]=3*x/2;I[0,1+M+j]=x;I[1,1+M+j]=x/2
            R[2,1+j]=u*x;R[2,1+M+j]=4*u*x;R[3,1+j]=v*x;R[3,1+M+j]=-v*x;R[4,1+j]=R[4,1+M+j]=u*x
        at=1+2*M
        for i in range(M):
            for j in range(M):
                q=a[i]*a[j];I[2,at]=2*q;I[3,at]=I[4,at]=q/2;at+=1
        for i in range(M):
            for j in range(i,M):
                q=a[i]*a[j]/2;I[2,at]=I[4,at]=q;I[3,at]=-q;at+=1
        return R,I


def configure_asymptotic_constraints(P):
    if not hasattr(P,'endpoint_inverse'):raise ValueError('Leading high-energy constraints require T0=0 coordinates')
    R,I=asymptotic_rows(P.M,max(256,ctx.prec))
    convert=lambda A:np.array([[np.longdouble(A[i,j].mid().str(40,radius=False)) for j in range(P.p)] for i in range(5)])
    R,I=convert(R)*P.ds,convert(I)*P.ds
    scale=np.maximum(1.,P.bound*np.linalg.norm(I[:,P.free:],ord=4/3,axis=1))
    T=P.endpoint_inverse[:P.p,:P.p]
    if P.inverse is not None:T=P.inverse@T
    P.asymptotic_R=np.asarray(T.T@R.T).T/np.sqrt(scale[:,None]);P.asymptotic_I=np.asarray(T.T@I.T).T/scale[:,None];P.asymptotic_scale=scale


def augment_asymptotic(P,result,z,hessian):
    if not hasattr(P,'asymptotic_R'):return result
    x=P.asymptotic_R@z[:P.p];y=P.asymptotic_I@z[:P.p];slack=2*y-x*x
    if not np.isfinite(slack).all() or np.any(slack<=0):return None
    result['value']-=np.log(slack).sum();result['asymptotic']=dict(x=x,slack=slack)
    if hessian:
        S=(2*x[:,None]*P.asymptotic_R-2*P.asymptotic_I)/slack[:,None]
        F=np.vstack((S,np.sqrt(2/slack[:,None])*P.asymptotic_R))
        result['features']=np.vstack((result['features'],np.pad(F,((0,0),(0,len(z)-P.p)))))
        result['rhs']=np.r_[result['rhs'],-np.ones(5),np.zeros(5)]
    return result


def asymptotic_slope(P,state,dz):
    if 'asymptotic' not in state:return 0.
    t=state['asymptotic'];dx=P.asymptotic_R@dz[:P.p];dy=P.asymptotic_I@dz[:P.p]
    return np.sum((2*t['x']*dx-2*dy)/t['slack'])


def asymptotic_duals(P,state,mu):
    if 'asymptotic' not in state:return {}
    t=state['asymptotic'];w=P.asymptotic_scale
    return dict(asymptotic_kR=np.asarray(2*mu*t['x']/(np.sqrt(w)*t['slack']),float),asymptotic_kI=np.asarray(-2*mu/(w*t['slack']),float))


def asymptotic_audit(coefficients,M,duals):
    R,I=asymptotic_rows(M,ctx.prec);p=len(coefficients);c=arb_mat(p,1,list(map(arb,coefficients)))
    x=(R*c).entries();y=(I*c).entries();u=list(map(arb,duals.get('asymptotic_kR',np.zeros(5))));v=list(map(arb,duals.get('asymptotic_kI',np.zeros(5))))
    if len(u)!=5 or len(v)!=5:raise ValueError('Five complete asymptotic covectors required')
    support=arb(0)
    for a,b in zip(u,v):
        if b<0:support-=a*a/(2*b)
        elif not b.is_zero() or not a.is_zero():raise ValueError('Asymptotic parabola covector has unbounded support')
    residual=-arb_mat(1,5,u)*R-arb_mat(1,5,v)*I
    return residual,support,[2*b-a*a for a,b in zip(x,y)]
