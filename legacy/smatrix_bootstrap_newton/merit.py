"""Self-concordant descent for the joint solver's local cached cone model.

Physical feasibility/support and numerical-center acceptance remain separate.
Only ordinary cones are covered; no Phase I or additional asymptotic barriers.
"""
import numpy as np
from flint import arb, arb_mat, ctx


def self_concordant_decrease(step,descent,curvature):
    """Exact-model lower bound; its floating evaluation is not a certificate."""
    if not np.isfinite([step,descent,curvature]).all() or step<=0 or descent<=0 or curvature<0:return -np.inf
    radius=step*np.sqrt(curvature)
    if radius>=1:return -np.inf
    return step*descent-step*step*curvature/(2*(1-radius))


def _binary(value):
    n,d=np.longdouble(value).as_integer_ratio()
    return arb(n)/arb(d)


def _density_derivatives(rho,velocity):
    if not rho:return arb(0),arb(0)
    a=[v*v for v in rho];remaining=1-sum((v*v for v in a),arb(0));N=len(a)
    if not remaining>0:raise ArithmeticError('Exact local density point is not strict')
    lo=(remaining/(N+1)).lower();hi=(remaining/(1+arb(N)/2)).upper();delta=lo.union(hi)
    for _ in range(16):
        mid=delta.mid();vm=[(v+(v*v+2*mid).sqrt())/2 for v in a]
        f=mid+sum((v*v for v in vm),arb(0))-1
        roots=[(v*v+2*delta).sqrt() for v in a];vals=[(v+r)/2 for v,r in zip(a,roots)]
        derivative=1+sum((v/r for v,r in zip(vals,roots)),arb(0))
        delta=delta.intersection(mid-f/derivative)
        if delta.rad()<abs(delta.mid())*arb(2)**-120:break
    if not delta>0:raise ArithmeticError('Positive density minimizer unresolved')
    vals=[(v+(v*v+2*delta).sqrt())/2 for v in a];d=[delta/(2*v) for v in vals]
    q=[delta+2*v*v for v in d];den=1+sum((delta/v for v in q),arb(0))
    first=sum((4*r*v*u/delta for r,v,u in zip(rho,vals,velocity)),arb(0))
    diagonal=sum(((2/b+8*x/t)*u*u for b,x,t,u in zip(d,a,q,velocity)),arb(0))
    rank=sum((-4*r*v*u/t for r,v,u,t in zip(rho,vals,velocity,q)),arb(0))
    return first,diagonal+rank*rank/den


def certified_local_decrease(P,z,state,dz,direction,mu,step,bits=256):
    """Enclose local directional derivatives without subtracting log values.

    Cached positive disk/chi slacks define their local shifted quadratics.
    Gram entries and long-double inputs are converted as exact binary rationals.
    This certifies the stated local numerical cone model, not an exact physical
    stationary point or the map from accumulated caches to the saved raw H.
    """
    if 'asymptotic' in state or 'moment_slack' in state:
        return dict(certified_descent=False,reason='Unsupported additional or Phase-I barrier')
    try:
        with ctx.workprec(bits):
            A=_binary;v=state['values'];amp=state['amplitude'];prime=arb(0);second=arb(0)
            def quadratic(q,first,curved,weight=arb(1)):
                nonlocal prime,second
                if not q>0 or not curved>=0:raise ArithmeticError('Non-strict local cone')
                prime-=weight*first/q;second+=weight*((first/q)**2+2*curved/q)
            for x,y,dx,dy,t,slack in zip(v['x'],v['y'],direction['x'],direction['y'],P.t,amp['slack']):
                x,y,dx,dy,t,q=map(A,(x,y,dx,dy,t,slack))
                quadratic(q,2*(1-t*y)*dy-2*x*dx,dx*dx+t*dy*dy)
            for group,slack in zip(P.groups,amp['chiral_slack']):
                x=list(map(A,v['chi'][group]));dx=list(map(A,direction['chi'][group]))
                quadratic(A(slack),-2*sum((a*b for a,b in zip(x,dx)),arb(0)),sum((a*a for a in dx),arb(0)),A(P.chiral_weight))
            rho=[A(x)/A(P.bound) for x in z[P.free:P.p]];dr=[A(x)/A(P.bound) for x in dz[P.free:P.p]]
            gd,hd=_density_derivatives(rho,dr);prime+=gd;second+=hd
            root=A(np.sqrt(np.longdouble(2)))
            def symmetric(values):
                x=list(map(A,values))
                return arb_mat([[x[0],x[1]/root,x[3]/root],[x[1]/root,x[2],x[4]/root],[x[3]/root,x[4]/root,x[5]]])
            trace=lambda B:sum((B[j,j] for j in range(3)),arb(0))
            for g,dg in zip(v['gram'],direction['gram']):
                G=symmetric(g)
                if not (G[0,0]>0 and G[0,0]*G[1,1]-G[0,1]**2>0 and G.det()>0):raise ArithmeticError('Local Gram is not strictly PD')
                T=G.inv()*symmetric(dg);prime-=trace(T);second+=trace(T*T)
            for f,df,slack in zip(v['ff'],direction['ff'],state['ffslack']):
                f,df=list(map(A,f)),list(map(A,df))
                quadratic(A(slack),-2*sum((a*b for a,b in zip(f,df)),arb(0)),sum((a*a for a in df),arb(0)))
            for r,dr in zip(v['moment'],direction['moment']):
                r,dr=A(r),A(dr);quadratic(1-r*r,-2*r*dr,dr*dr)
            if not A(mu)>0 or not A(step)>0 or not A(P.chiral_weight)>=1:raise ArithmeticError('Positive step/mu and barrier weights >=1 required')
            descent=sum((A(c)*A(d) for c,d in zip(P.cost,dz)),arb(0))/A(mu)-prime
            if not second>=0:raise ArithmeticError('Nonnegative local curvature unresolved')
            radius=A(step)*second.sqrt()
            lower=(A(step)*descent-A(step)**2*second/(2*(1-radius))).lower() if radius<1 else arb('-inf')
            return dict(certified_descent=bool(descent>0 and radius<1 and lower>0),bits=bits,
                descent=descent.str(30),curvature=second.str(30),radius_upper=radius.upper().str(30),decrease_lower=lower.str(30),
                scope='Exact-binary local cached cone model; independent original-H audit still required')
    except (ArithmeticError,ValueError,ZeroDivisionError,OverflowError) as error:
        return dict(certified_descent=False,reason=str(error),bits=bits)
