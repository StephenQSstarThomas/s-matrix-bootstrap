"""One analytic density interpolation for every Mandelstam amplitude leg.

Sine-cardinal density samples give g_j(z)=sum_n D_nj(z^n-(-1)^n).
This exactly matches K+b+iI at native nodes and has g_j(-1)=0. It is a
declared finite interpolation, not recovery of an unspecified author choice.
"""
from __future__ import annotations

import math
import time

import numpy as np
from flint import arb, acb, arb_mat, acb_mat, ctx

from .precision import PrecisionRows,assemble_rows

_RULES = {}


class SineFamily(PrecisionRows):
    def __init__(self,M,L,dps=40):
        guard=math.ceil(-2*(2*L-1)*math.log10(math.tan(math.pi/(8*M))))
        super().__init__(M,dps+guard+10)
        self.geometry_bits=self.bits
        self.dps,self.L,self.bits=dps,L,math.ceil((dps+15)*math.log2(10))
        self.D=arb_mat([[((arb.pi()*n*(2*j+1)/(2*M)).sin()*(1 if n==M else 2)/M)
                        for j in range(M)] for n in range(1,M+1)])
        self._angular_cache={}

    def cardinal(self,z):
        z=acb(z) if isinstance(z,acb) else arb(z)
        values=[]; power=type(z)(1)
        for n in range(1,self.M+1):
            power*=z;values.append(power-(-1)**n)
        if isinstance(z,acb):
            row=acb_mat([values])*acb_mat(self.D)
        else:
            row=arb_mat([values])*self.D
        return [row[0,j] for j in range(self.M)]

    def cardinal_at(self,s):
        s=arb(s) if isinstance(s,arb) else arb(str(s))
        if s<4:
            r=(4-s).sqrt();return self.cardinal((2-r)/(2+r))
        if s>4:
            r=acb(0,(s-4).sqrt());return self.cardinal((2+r)/(2-r))
        return self.cardinal(arb(1))

    def _values(self,z):
        rows=[]
        for v in z:
            power=arb(1); row=[]
            for n in range(1,self.M+1):
                power*=v;row.append(power-(-1)**n)
            rows.append(row)
        return arb_mat(rows)*self.D

    @staticmethod
    def error_bound(s,ell,n):
        s=arb(s)
        X=abs((s/4).log())/2; C=abs((s+4)/(s-4))
        b=arb.pi()/(4*X);rho=b+(1+b*b).sqrt()
        bound=80*C*X*(C+(C*C+1).sqrt())**ell
        return 8*bound*rho**(-2*n)/(1-1/rho)

    def quadrature(self,s,ell,target,order=None):
        s=arb(s);r=abs((s.sqrt()-2)/(s.sqrt()+2))
        guard=max(0,int(math.ceil(-ell*float(r.log().mid())/math.log(10)))) if ell else 0
        ctx.prec=math.ceil((self.dps+guard+15)*math.log2(10))
        X=abs((s/4).log())/2;C=(s+4)/(s-4)
        if order is None:
            n=16
            while not self.error_bound(s,ell,n)<target:
                n+=16
        else:n=order
        key=n,ctx.prec
        if key not in _RULES:
            _RULES[key]=[arb.legendre_p_root(n,k,weight=True) for k in range(n)]
        roots,weights=zip(*_RULES[key]);xs=[X*v for v in roots]
        mu=np.array([-C*x.tanh() for x in xs],dtype=object)
        jac=np.array([abs(C)*X/(x.cosh()**2)*w for x,w in zip(xs,weights)],dtype=object)
        rt=[((s+4)*x.exp()/(2*x.cosh())).sqrt() for x in xs]
        z=[(2-v)/(2+v) for v in rt]
        gt=self._values(z);gu=arb_mat(gt.tolist()[::-1])
        pol=np.array([v.legendre_p(ell) for v in mu],dtype=object)
        return gt,gu,jac*pol,self.error_bound(s,ell,n),n

    def rows(self,isospin,ell,*,point=None,node=None):
        """Source coefficient rows, with an enclosure for the angular integral error."""
        if (point is None)==(node is None):raise ValueError('Specify a native node or a point')
        s=self.s[node] if node is not None else (arb(point) if isinstance(point,arb) else arb(str(point)))
        key=(node,None if node is not None else s.mid().str(self.dps+10,radius=False),ell)
        if key not in self._angular_cache:
            ctx.prec=self.geometry_bits
            # Account for the centrifugal order before the cone's row scaling.
            centrifugal=abs((s.sqrt()-2)/(s.sqrt()+2))**ell if node is not None else arb(1)
            target=arb(10)**(-(self.dps-2))*centrifugal/(arb.pi()*arb(self.lay.n).sqrt())
            gt,gu,weights,error,n=self.quadrature(s,ell,target)
            weighted=arb_mat([[gu[i,j]*weights[i] for j in range(self.M)] for i in range(n)])
            J=arb_mat([weights.tolist()])*gt
            D=gt.transpose()*weighted
            self._angular_cache[key]=(np.array(J.tolist()[0],dtype=object),
                                      np.array(D.tolist(),dtype=object),error,n,ctx.prec)
        J,D,error,n,bits=self._angular_cache[key];ctx.prec=bits
        if node is not None:
            im=np.zeros(self.M,dtype=object);im[node]=1
            cr=[self.on[node],im]
        else:
            cr=[np.array(self.cardinal_at(s),dtype=object),np.zeros(self.M,dtype=object)]
        rows=assemble_rows(self.lay,cr,J,D,isospin,ell)
        radius=arb(0,error.abs_upper())
        return np.array([[v+radius for v in row] for row in rows],dtype=object)

    def wave(self,c,isospin,ell,*,point=None,node=None,error_target='1e-12',order=None):
        """Direct A(s,t,u) quadrature with its analytic tail and arithmetic enclosures."""
        ctx.prec=self.geometry_bits
        c=[arb(v) if isinstance(v,arb) else arb(float(v)) for v in c]
        if len(c)!=self.lay.n or (point is None)==(node is None):
            raise ValueError('Complete coefficients and exactly one evaluation point required')
        s=self.s[node] if node is not None else arb(str(point))
        norm=sum((abs(v) for v in c),arb(0));target=arb(error_target)/max(norm,arb(1))
        gt,gu,weights,error,n=self.quadrature(s,ell,target,order)
        M=self.M;lay=self.lay
        sig1=arb_mat([[v] for v in c[lay.s1]]);sig2=arb_mat([[v] for v in c[lay.s2]])
        r1=arb_mat(M,M,c[lay.r1]);r2=arb_mat(M,M)
        for (i,j),v in zip(zip(*lay.triu),c[lay.r2]):r2[i,j]=r2[j,i]=v
        gs=acb_mat([self.cardinal_at(s)])
        sigs1=(gs*acb_mat(sig1))[0,0];sigs2=(gs*acb_mat(sig2))[0,0]
        sigt1,sigu1,sigt2,sigu2=gt*sig1,gu*sig1,gt*sig2,gu*sig2
        t1,u1,t2=gt*r1,gu*r1,gt*r2
        s1=gs*acb_mat(r1);s2=gs*acb_mat(r2)
        st1,su1=s1*acb_mat(gt.transpose()),s1*acb_mat(gu.transpose())
        ts1,us1=acb_mat(t1)*gs.transpose(),acb_mat(u1)*gs.transpose()
        su2,st2=s2*acb_mat(gu.transpose()),s2*acb_mat(gt.transpose())
        answer=acb(0)
        for i in range(n):
            tu1=sum((t1[i,j]*gu[i,j] for j in range(M)),arb(0))
            ut1=sum((u1[i,j]*gt[i,j] for j in range(M)),arb(0))
            tu2=sum((t2[i,j]*gu[i,j] for j in range(M)),arb(0))
            A=c[0]+sigs1+sigt2[i,0]+sigu2[i,0]+st1[0,i]+su1[0,i]+tu2
            B=c[0]+sigt1[i,0]+sigs2+sigu2[i,0]+ts1[i,0]+tu1+su2[0,i]
            Cc=c[0]+sigu1[i,0]+sigt2[i,0]+sigs2+ut1+us1[i,0]+st2[0,i]
            T={0:3*A+B+Cc,1:B-Cc,2:B+Cc}[isospin]
            answer+=weights[i]*T/4
        radius=(error*norm).abs_upper()
        answer+=acb(arb(0,radius),arb(0,radius) if s>4 else arb(0))
        return answer,{'quadrature_nodes':n,'working_bits':ctx.prec,
                       'analytic_tail_bound':str(radius.upper().fmpq()),'coefficient_l1':str(norm.upper().fmpq())}


def replay(source_report,outdir,dps=40):
    """Two existing locations in one saved witness, without optimizing new curves."""
    from pathlib import Path
    import json,hashlib
    from .sdpb import write_json
    source=Path(source_report).resolve();old=json.loads(source.read_text());M=old['spec']['M']
    family=SineFamily(M,old['spec']['L'],dps)
    ctx.prec=family.bits;basis=np.load(source.parent/'basis.npy')
    y=[arb(v) for v in (source.parent/'out/y.txt').read_text().splitlines()[1:] if v.strip()]
    av=arb_mat([[v] for v in y[:basis.shape[1]]])
    c=[(arb_mat([[arb(float(x)) for x in row]])*av)[0,0] for row in basis]
    low=max(i for i,s in enumerate(family.s) if s<arb(3600)/49)
    rows=[];start=time.monotonic()
    for I,ell,point,node in ((0,0,'.5',None),(1,1,None,low)):
        f,info=family.wave(c,I,ell,point=point,node=node)
        rows.append({'isospin':I,'ell':ell,'point':point,'node':node,
                     'f_real':f.real.str(40),'f_imag':f.imag.str(40),**info})
    result={'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'prescription':'finite-sine-cardinal candidate','dps':dps,'rows':rows,
            'seconds':time.monotonic()-start,'optimization_performed':False,
            'scope':'same saved coefficients under a different, consistently analytic interpolation; no feasibility transfer'}
    dest=Path(outdir);dest.mkdir(parents=True,exist_ok=True);write_json(dest/'report.json',result)
    return result


def mma_check(outdir):
    """Compare source rows to a separate Mathematica evaluation of the original angular integral."""
    from pathlib import Path
    import subprocess,json,hashlib
    from .sdpb import stop_owned_process,wait_process,write_json
    from .mma import MMA_TIMEOUT_SECONDS
    import os
    root=Path(__file__).resolve().parents[3];dest=Path(outdir).resolve();dest.mkdir(parents=True,exist_ok=True)
    script=root/'scripts/mma/audit_sine_family.wls';cid=dest/'mma.cid'
    args=[str(root/'scripts/mma/wolfram.sh'),str(script),f'rw:{dest}:/audit','--','/audit/mma.json']
    with (dest/'mma.log').open('w') as log:
        proc=subprocess.Popen(args,stdout=log,stderr=subprocess.STDOUT,
                env=dict(os.environ,MMA_CIDFILE=str(cid)),start_new_session=True)
        try:rc=wait_process(proc,MMA_TIMEOUT_SECONDS)
        finally:
            if proc.poll() is None:stop_owned_process(proc,cid,log)
    if rc:
        write_json(dest/'report.json',{'passed':False,'returncode':rc,'stage':'Mathematica launcher',
                    'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest()})
        raise RuntimeError(f'Mathematica audit failed: {rc}')
    ref=json.loads((dest/'mma.json').read_text());p=SineFamily(8,2,dps=40);rows=[]
    for key,I,l,index in [('sigma1_node4_S0',0,0,p.lay.s1.start+4),
                          ('rho1_1_4_P1',1,1,p.lay.r1.start+8+4)]:
        c=[arb(0) for _ in range(p.lay.n)];c[index]=1
        f,info=p.wave(c,I,l,node=4,error_target='1e-30')
        error=max(float(abs(f.real-arb(ref[key][0])).abs_upper()),float(abs(f.imag-arb(ref[key][1])).abs_upper()))
        rows.append({'case':key,'max_absolute_error':error,'passed':error<1e-22,**info})
    result={'passed':all(r['passed'] for r in rows),'rows':rows,'MMA':ref,
            'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),
            'scope':'independent original angular integrals for single and double density controls'}
    write_json(dest/'report.json',result);return result
