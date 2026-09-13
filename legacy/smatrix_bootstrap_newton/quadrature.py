"""Validated angular projections of the one conformal-polynomial family."""
from functools import lru_cache
from flint import arb, acb, arb_mat, ctx
from .kernels import density_row_to_cflat
from .operators import assemble_density_row


@lru_cache(maxsize=4)
def legendre_rule(order,bits):
    return tuple(arb.legendre_p_root(order,j,weight=True) for j in range(order))


def panels(span,gap):
    left=[arb(0)];step=arb(gap)
    while left[-1]+step<span/2:
        left.append(left[-1]+step);step*=2
    left.append(span/2)
    edges=left+[span-x for x in left[-2::-1]]
    return list(zip(edges[:-1],edges[1:]))


def angular_mode_blocks(s,M,spins,order=64):
    """Gauss moments with an explicit Bernstein-ellipse remainder.

    For f analytic on E_rho and |f|<=B, positivity and polynomial exactness
    give |integral f-Q_n f| <= 8 B rho^(-2n)/(1-rho^(-1)).
    Each crossed conformal power has modulus <=1 on the chosen ellipse.
    """
    s=arb(s);spins=sorted(set(spins))
    if type(order) is not int or order<1:raise ValueError('Positive Gauss order required')
    order=max(64,order)
    if not s.is_finite() or not s>0 or not (s<4 or s>4):raise ValueError('Finite positive nonthreshold real s required')
    span=abs(s-4);gap=arb(4) if s>4 else s
    direction=-1 if s>4 else 1
    X=[];Y=[];W={ell:[] for ell in spins};error={ell:arb(0) for ell in spins}
    for lo,hi in panels(span,gap):
        h=(hi-lo)/2;center=(hi+lo)/2
        distance=min((gap+lo).lower(),(gap+span-hi).lower())
        a=(1+distance/(2*h.upper())).lower();rho=(a+(a*a-1).sqrt()).lower()
        R=abs(1-2*center/span)+2*a*h/span
        base=R+(R*R+1).sqrt()
        for ell in spins:
            error[ell]+=4*h/span*base**ell*rho**(-2*order)/(1-1/rho)
        for root,weight in legendre_rule(order,ctx.prec):
            u=center+h*root;mu=1-2*u/span
            t=direction*u;v=direction*(span-u)
            zt=(2-(4-t).sqrt())/(2+(4-t).sqrt())
            zu=(2-(4-v).sqrt())/(2+(4-v).sqrt())
            q=[];r=[];x=y=arb(1)
            for _ in range(M):x*=zt;y*=zu;q.append(x);r.append(y)
            X.append(q);Y.append(r)
            for ell in spins:W[ell].append(h*weight*mu.legendre_p(ell)/(2*span))
    A=arb_mat(X);B=arb_mat(Y);AT=A.transpose();out={}
    for ell in spins:
        J=arb_mat(1,len(X),W[ell])*A
        weighted=arb_mat([[w*y for y in row] for w,row in zip(W[ell],Y)])
        U=AT*weighted;rad=arb(0,error[ell].upper())
        J=[v+rad for v in J.entries()]
        for i in range(M):
            for j in range(i,M):
                v=U[i,j]+rad
                if ell%2 and i==j:
                    if not v.contains(0):raise ArithmeticError('Odd crossed diagonal lost parity')
                    v=arb(0)
                U[i,j]=v;U[j,i]=(-1)**ell*v
        out[ell]=dict(c=arb(1)/2 if ell==0 else arb(0),ell=ell,epsilon=(-1)**ell,
            J=J,U=[[U[i,j] for j in range(M)] for i in range(M)],
            error_bound=error[ell].str(25),order=order,panels=len(panels(span,gap)))
    return out


def polynomial_rows(s,M,waves,node=None,order=64):
    from .basis import mode_row_to_cardinal
    s=arb(s);spins={ell for I,ell in waves}
    if node is not None:
        phi=arb(2*node+1)/(2*M);z=acb(phi.cos_pi(),phi.sin_pi())
    else:
        root=acb(4-s).sqrt()
        if s>4:root=-root
        z=(2-root)/(2+root)
    direct=[z**n for n in range(1,M+1)]
    blocks=angular_mode_blocks(s,M,spins,order);rows={}
    for I,ell in waves:
        raw=density_row_to_cflat(assemble_density_row(direct,blocks[ell],I),M)
        row=mode_row_to_cardinal(raw,M)
        if node is not None:
            for j in range(1+2*M):
                v=arb(0)
                if ell==0 and I in (0,2):
                    if j==1+node and I==0:v=arb(3)/2
                    if j==1+M+node:v=arb(1)
                if not row[j].imag.contains(v):raise ArithmeticError('Native free-column identity failed')
                row[j]=acb(row[j].real,v)
        rows[I,ell]=row
    return rows,{ell:{k:b[k] for k in ('order','panels','error_bound')} for ell,b in blocks.items()}


def reuse_prepared_rows(source,preparation):
    """Same-M row reuse; generate only missing spins in the same analytic family."""
    from pathlib import Path
    import numpy as np
    from . import read_json,digest
    from .kernels import canonical_waves
    path=Path(preparation);meta=read_json(path/'report.json')['amplitude_model']
    M,L=source.M,source.L;oldL=meta['L'];oldn=3*M*oldL;n=source.physical_count;p=source.coefficient_count
    if (meta['prescription'],meta['M'],meta['coordinate_system'])!=('analytic-cardinal',M,source.coordinate_name):
        raise ValueError('Same-M analytic coordinates required for operator reuse')
    if meta.get('source_registry_sha256')!=source.metadata.get('source_registry_sha256') or not meta.get('integral_enclosures'):
        raise ValueError('Matching certified analytic source required')
    with np.load(path/'amplitude.npz') as data:
        H=data['rows'];E=data['radius_upper']
        if not np.array_equal(data['energies'],list(map(float,source.x))) or not np.array_equal(data['waves'],canonical_waves(oldL)) or not np.array_equal(data['coefficient_labels'],np.asarray(source.coefficient_labels,dtype=str)):
            raise ValueError('Native node, wave or coefficient identity mismatch')
    if H.shape!=(2*oldn+11,p) or E.shape!=H.shape or not np.isfinite(H).all() or not np.isfinite(E).all() or np.any(E<0):
        raise ValueError('Complete prior analytic coefficient enclosures required')
    for file in (path/'report.json',path/'amplitude.npz'):source.inputs[str(file.resolve())]=dict(sha256=digest(file))
    source.metadata.update(parent_preparation=str(path.resolve()),parent_arithmetic_bits=meta['bits'],
        reused_rows=6*M*min(L,oldL)+11,generated_rows=6*M*max(0,L-oldL),
        angular_projection='Prior analytic coefficient enclosures; missing spins use validated Gauss projection')
    def oldrow(index):return [arb(float(v),float(e)) for v,e in zip(H[index],E[index])]
    missing=[w for w in source.waves if w[1]//2>=oldL]
    for node in range(M):
        generated={}
        if missing:
            generated,accuracy=polynomial_rows(source.x[node],M,missing,node,source.order)
            source.quadrature_accuracy[str(node)]=accuracy
        for k,(I,ell) in enumerate(source.waves):
            dest=node*3*L+k;entry=dict(node_zero_based=node,isospin=I,ell=ell)
            if (I,ell) in generated:
                row=generated[I,ell]
                yield dest,entry,[v.real for v in row]
                yield n+dest,entry,[v.imag for v in row]
            else:
                src=node*3*oldL+I*oldL+ell//2
                yield dest,entry,oldrow(src)
                yield n+dest,entry,oldrow(oldn+src)
    for j in range(11):yield 2*n+j,dict(kind='retained_analytic_supplement'),oldrow(2*oldn+j)
