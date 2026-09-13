"""Task 5a.10: compare the from-scratch operators against independent targets.

Nothing here is an implementation source.  The repository's historical
``kernels``/``operators``/``model`` modules (Arb arithmetic, a completely
separate code path written for the Newton mainline) are loaded *only* here, so
that the new rows can be diffed against them row by row.

Packing note (the one documented convention difference).  The historical row is
indexed by ``kernels.density_labels`` -- ``T0``, ``sigma1_i``, ``sigma2_i``,
``rho1_{ij}`` (full M x M), ``rho2_{ij}`` for ``i <= j`` -- which is the same
ordering as :class:`smatrix_bootstrap.sdp.projector.Layout`.  It then applies
``kernels.density_row_to_cflat``, halving every off-diagonal ``rho2``
coefficient, because its free variable is "C_flat" = 2 rho2_{ij} for i != j
(see ``kernels.coefficient_blocks``, which divides by two on the way back).
Our packed variable is rho2_{ij} itself, so our coefficient is the sum over the
symmetric pair.  The two describe the same bilinear form; comparison therefore
applies the same halving to the new row.
"""
from __future__ import annotations

import numpy as np


def old_row(M: int, L: int, isospin: int, ell: int, s: float, node: int | None, bits: int = 384):
    """One partial-wave row from the historical Arb implementation."""
    src = _src(M, L, bits)
    if node is not None:
        s = float(src.x[node].str(25, radius=False))   # its own float64 node value
    return np.array(_to_complex(src.row(s, ell, isospin, node=node)))


def to_cflat(row: np.ndarray, lay) -> np.ndarray:
    """Convert a Layout-packed row to the historical C_flat convention."""
    out = row.copy()
    i, j = lay.triu
    blk = out[lay.r2].copy()
    blk[i != j] *= 0.5
    out[lay.r2] = blk
    return out


_SRC: dict = {}


def _src(M, L, bits):
    if (M, L, bits) not in _SRC:
        from ..kernels import PVSourceRows
        _SRC[(M, L, bits)] = PVSourceRows(M=M, L=L, bits=bits, subtracted=False)
    return _SRC[(M, L, bits)]


def _to_complex(row):
    from flint import acb
    out = []
    for v in row:
        a = acb(v)
        out.append(complex(float(a.real.str(20, radius=False)),
                           float(a.imag.str(20, radius=False))))
    return out


def compare_rows(new_op, M: int, L: int, cases) -> list[dict]:
    """Row-by-row diff for a list of ``(isospin, ell, s, node)`` cases."""
    out = []
    for isospin, ell, s, node in cases:
        new = new_op.rows(isospin, ell, s, node)
        new_c = new[0] + 1j * new[1]
        new_c = to_cflat(new_c, new_op.lay)
        old = old_row(M, L, isospin, ell, s, node)
        scale = max(np.abs(old).max(), np.abs(new_c).max(), 1e-300)
        out.append({"isospin": isospin, "ell": ell, "s": float(s), "node": node,
                    "max_abs_diff": float(np.abs(new_c - old).max()),
                    "row_scale": float(scale),
                    "max_rel_diff": float(np.abs(new_c - old).max() / scale)})
    return out


def compare_hilbert_kernel(M: int) -> float:
    """New (3.67) vs the historical ``kernels.pv_matrix`` (Arb)."""
    from ..kernels import pv_matrix
    from .hilbert import hilbert_kernel
    old = pv_matrix(M)
    o = np.array([[float(old[i, j].str(20, radius=False)) for j in range(M)] for i in range(M)])
    return float(np.abs(hilbert_kernel(M) - o).max())


def compare_kinematic_squares(s_values) -> float:
    """New (2.33)^2 vs the historical ``model.current_kinematic_squares``."""
    from ..model import current_kinematic_squares
    from .formfactor import kinematic_factor
    worst = 0.0
    for s in s_values:
        old = current_kinematic_squares(float(s))
        for ell in (0, 1):
            o = float(old[ell].str(25, radius=False))
            n = float(kinematic_factor(ell, np.array([float(s)]))[0] ** 2)
            worst = max(worst, abs(n / o - 1.0))
    return worst


def compare_grid(M: int) -> dict:
    """New (3.60)-(3.61) grid and weights vs ``operators.midpoint_grid`` (Arb)."""
    from ..operators import midpoint_grid
    from .grid import dsdphi, s_nodes
    x, w = midpoint_grid(M)
    xo = np.array([float(v.str(25, radius=False)) for v in x])
    wo = np.array([float(v.str(25, radius=False)) for v in w])
    return {"nodes_max_rel": float(np.abs(s_nodes(M) / xo - 1).max()),
            "weights_max_rel": float(np.abs((dsdphi(M) / M) / wo - 1).max())}


class OriginalAngle:
    """Independent sine polynomials integrated in the original mu variable.

    This reference does not call production cardinal, angular or packing code.
    Gauss order agreement is empirical; no rigorous reference tail is claimed.
    """
    def __init__(self, M, coefficients, dps=140):
        from mpmath import mp
        self.M,self.dps,self.rules=M,dps,{}
        with mp.workdps(dps):
            c=[mp.mpf(v) for v in coefficients]
            if len(c)!=1+2*M+M*M+M*(M+1)//2:raise ValueError('Complete density coefficients required')
            B=mp.matrix(M+1,M)
            for j in range(M):
                phi=mp.pi*(2*j+1)/(2*M)
                for n in range(1,M+1):B[n,j]=mp.mpf(1 if n==M else 2)*mp.sin(n*phi)/M
                B[0,j]=-mp.fsum((-1)**n*B[n,j] for n in range(1,M+1))
            r1=mp.matrix([c[1+2*M+i*M:1+2*M+(i+1)*M] for i in range(M)])
            r2=mp.matrix(M);k=1+2*M+M*M
            for i in range(M):
                for j in range(i,M):r2[i,j]=r2[j,i]=c[k];k+=1
            self.constant=c[0]
            self.single=[B*mp.matrix(c[1:1+M]),B*mp.matrix(c[1+M:1+2*M])]
            self.double=[B*r1*B.T,B*r2*B.T]

    def integrate(self,isospin,ell,s,degree):
        """Order 3*2**(degree-1) Gauss-Legendre on [-1,1], without x/tanh mapping."""
        from mpmath import mp
        with mp.workdps(self.dps):
            s=mp.mpf(s)
            def z(v):
                r=-mp.j*mp.sqrt(v-4) if v>4 else mp.sqrt(4-v)
                return (2-r)/(2+r)
            a,b,c={0:(3,1,1),1:(0,1,-1),2:(0,1,1)}[isospin]
            zs=z(s);powers=mp.matrix([zs**n for n in range(self.M+1)])
            s1,s2=self.single;r1,r2=self.double
            left=(powers.T*r1).T;right=r1*powers;current=r2*powers
            const=(a+b+c)*self.constant+a*(powers.T*s1)[0]+(b+c)*(powers.T*s2)[0]
            pt=list(b*s1+(a+c)*s2+a*left+b*right+c*current)[::-1]
            pu=list(c*s1+(a+b)*s2+a*left+c*right+b*current)[::-1]
            Q=b*r1+c*r1.T+a*r2
            table=[[Q[i,j] for j in range(self.M,-1,-1)] for i in range(self.M,-1,-1)]
            def horner(poly,x):
                value=mp.mpf(0)
                for coefficient in poly:value=value*x+coefficient
                return value
            if degree not in self.rules:
                self.rules[degree]=mp._gauss_legendre.calc_nodes(degree,mp.prec)
            values=[]
            for mu,weight in self.rules[degree]:
                t=-(s-4)*(1-mu)/2;u=4-s-t;zt,zu=z(t),z(u)
                mixed=horner([horner(row,zu) for row in table],zt)
                values.append(weight*mp.legendre(ell,mu)*(const+horner(pt,zt)+horner(pu,zu)+mixed)/4)
            return mp.fsum(values)


def _angle_identity(record,snapshot=None):
    """Authenticate the source family, allowing unrelated sine audit-helper edits."""
    import ast,gzip,hashlib,json
    from pathlib import Path
    old={} if snapshot is None else json.loads(gzip.decompress(Path(snapshot).read_bytes()))
    checked={}
    for name in ('sine.py','precision.py','projector.py','legendreq.py'):
        key='src/smatrix_bootstrap/sdp/'+name;now=Path(__file__).with_name(name).read_text()
        digest=hashlib.sha256(now.encode()).hexdigest();expected=record['source_sha256'].get(key)
        if digest==expected:checked[name]='identical file';continue
        if name!='sine.py' or key not in old or hashlib.sha256(old[key].encode()).hexdigest()!=expected:
            raise ValueError('Source family changed or unauthenticated: '+name)
        def family(text):
            # This documented widening changes no central value or native row.
            text=text.replace('arb(0,radius) if node is not None else arb(0)',
                              'arb(0,radius) if s>4 else arb(0)')
            body=[]
            for node in ast.parse(text).body:
                if isinstance(node,ast.FunctionDef):break
                if not isinstance(node,ast.Expr):body.append(node)
            return ast.dump(ast.Module(body=body,type_ignores=[]),include_attributes=False)
        if family(now)!=family(old[key]):raise ValueError('SineFamily or its numerical imports changed')
        checked[name]='authenticated snapshot; same scientific code, allowing documented physical-point error-radius widening'
    return checked


def audit_angles(source_report,outdir,bits=512,timeout=600,snapshot=None,cases=None,points=None):
    """Three preselected frozen-witness controls, with per-order progress records."""
    import hashlib,importlib.metadata,json,math,platform,signal,time
    from pathlib import Path
    from flint import arb,acb,arb_mat,ctx
    from mpmath import mp
    from .sine import SineFamily
    from .sdpb import write_json
    if not math.isfinite(timeout) or timeout<=0:raise ValueError('Require a positive finite timeout')
    start=time.monotonic();source=Path(source_report).resolve();dest=Path(outdir).resolve()
    if dest.exists() and any(dest.iterdir()):raise FileExistsError('Preserve existing angular audit')
    rec=json.loads(source.read_text())
    if (rec.get('solver')!='SDPB' or rec.get('status') not in ('numerically_accepted','not_accepted')
            or not isinstance(rec.get('verification'),dict) or 'sdpb' not in rec):
        raise ValueError('Select a completed leaf SDPB report with saved verification')
    spec=rec['spec'];M=spec['M'];na=rec['pmp']['n_a']
    if spec.get('scattering_prescription')!='sine-cardinal':raise ValueError('Select a sine-cardinal witness')
    if not spec.get('reduce_basis'):raise ValueError('angle-audit currently requires a reduced-basis witness')
    identity=_angle_identity(rec,snapshot)
    points=points or []
    cases=cases or ([] if points else [(0,0,42),(1,1,33),(1,19,0)])
    if any(I not in (0,1,2) or ell<0 or ell%2!=I%2 or not 0<=k<=M-1 or 2*k!=int(2*k)
           for I,ell,k in cases):
        raise ValueError('Select valid native rows or phi midpoints between adjacent nodes')
    with mp.workdps(140):
        if any(I not in (0,1,2) or ell<0 or ell%2!=I%2 or not mp.isfinite(mp.mpf(s))
               or mp.mpf(s)<=0 or mp.mpf(s)==4 for I,ell,s in points):
            raise ValueError('Explicit points require positive s other than the threshold and valid isospin/parity')
    files=[source,source.parent/'basis.npy',source.parent/'out/y.txt',source.parent/'solution.npz',source.parent/'pmp.json']
    if snapshot:files.append(Path(snapshot).resolve())
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.exists()}
    for name,wanted in [('basis.npy',rec['basis']['sha256']),('solution.npz',rec['solution_sha256']),('pmp.json',rec['pmp']['sha256'])]:
        if hashes.get(str(source.parent/name))!=wanted:raise ValueError('Frozen input hash mismatch: '+name)
    lines=(source.parent/'out/y.txt').read_text().splitlines();texts=[v.strip() for v in lines[1:] if v.strip()]
    try:header=[int(v) for v in lines[0].split()]
    except (ValueError,IndexError) as exc:raise ValueError('Invalid y text header') from exc
    if header!=[rec['pmp']['n_vars']-1,1]:raise ValueError('y text header must be exactly [n_vars-1,1]')
    if len(texts)!=rec['pmp']['n_vars']-1:raise ValueError('Complete y text required')
    if not all(arb(v).is_finite() for v in texts):raise ValueError('Nonfinite saved y value')
    rounded=np.array(texts,dtype=float)
    if not np.all(np.isfinite(rounded)):raise ValueError('Nonfinite binary64 saved y value')
    if rec.get('y_sha256') and rec['y_sha256']!=hashes[str(source.parent/'out/y.txt')]:
        raise ValueError('Recorded full y text hash mismatch')
    with np.load(source.parent/'solution.npz') as sol:
        if not np.array_equal(rounded,sol['y']):raise ValueError('y text differs from authenticated saved y')
    basis=np.load(source.parent/'basis.npy',allow_pickle=False)
    if (basis.dtype!=np.dtype('float64') or basis.shape!=(1+2*M+M*M+M*(M+1)//2,na)
            or not np.all(np.isfinite(basis))):
        raise ValueError('Saved basis must be finite binary64 with the recorded dimensions')
    dps=max(140,math.ceil(bits*math.log10(2)));ctx.prec=math.ceil((dps+25)*math.log2(10))
    report={'source_report':str(source),'input_sha256':hashes,'source_identity':identity,
        'producer_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                           for name in ('crosscheck.py','__main__.py','sine.py','precision.py','projector.py')},
        'source_spec':spec,'reference_dps':dps,'coefficient_bits':ctx.prec,'timeout_seconds':timeout,
        'versions':{'python':platform.python_version(),**{name:importlib.metadata.version(name)
                     for name in ('mpmath','python-flint','numpy')}},
        'selected_rows':[list(v) for v in cases],'results':[],'status':'running','stage':'coefficient reconstruction',
        'selected_points':[list(v) for v in points],
        'coefficient_reconstruction':'exact stored binary64 basis times complete saved y decimal text in Arb',
        'y_integrity':('full text matches the solve-time recorded hash' if rec.get('y_sha256') else
            'fresh full-text hash; binary64 equality to authenticated NPZ y; solve-time full-text hash unavailable'),
        'scope':'Independent original-mu sine-polynomial quadrature; two-order agreement is empirical, not rigorous reference integration, continuum unitarity or optimizer certification.',
        'reference_rule':'mpmath Gauss-Legendre roots/weights; order=3*2**(degree-1); original mu in [-1,1]',
        'passed':False}
    dest.mkdir(parents=True,exist_ok=True)
    def save():
        report['seconds']=time.monotonic()-start;write_json(dest/'progress.tmp',report)
        (dest/'progress.tmp').replace(dest/'report.json')
    def expired(*_):raise TimeoutError('Angular audit runtime budget exhausted')
    save();previous=signal.signal(signal.SIGALRM,expired)
    try:
        signal.setitimer(signal.ITIMER_REAL,max(.000001,timeout-(time.monotonic()-start)))
        y=arb_mat([[arb(v)] for v in texts[:na]])
        coefficients=[(arb_mat([[arb(float(v)) for v in row]])*y)[0,0] for row in basis]
        norm=sum((abs(v) for v in coefficients),arb(0));report['coefficient_l1_ball']=norm.str(40)
        report['stage']='independent polynomial transform';save()
        reference=OriginalAngle(M,[v.mid().str(dps+10,radius=False) for v in coefficients],dps)
        family=SineFamily(M,spec['L'],spec['operator_dps'])
        for I,ell,k,point in [(*v,None) for v in cases]+[(I,ell,None,s) for I,ell,s in points]:
            with mp.workdps(dps):
                s=mp.mpf(point) if point is not None else 8/(1+mp.cos(mp.pi*(2*k+1)/(2*M)))
                r=abs((mp.sqrt(s)-2)/(mp.sqrt(s)+2))
                scale=r**ell if ell>1 else mp.mpf(1);target=scale*mp.mpf('1e-16')
                row={'I':I,'ell':ell,'node_zero_based':k,'s':mp.nstr(s,40),'energy_GeV':mp.nstr(mp.mpf('.14')*mp.sqrt(s),20),
                     'comparison_scale':mp.nstr(scale,40),'absolute_target':mp.nstr(target,40),'reference_orders':[]}
                if k is not None and k!=int(k):
                    row.update(node_zero_based=None,midpoint_after_node=int(k),sampling='phi midpoint; diagnostic only, not an imposed constraint')
                if point is not None:row.update(explicit_point=str(point),sampling='explicit s; diagnostic only')
                report['results'].append(row);values=[]
                for degree in ((5,6) if ell>1 and k==0 else (7,8)):
                    report['stage']=f'original-angle I={I} ell={ell} node={k} order={3*2**(degree-1)}';save()
                    value=reference.integrate(I,ell,s,degree);values.append(value)
                    row['reference_orders'].append({'order':3*2**(degree-1),'real':mp.nstr(value.real,dps),'imag':mp.nstr(value.imag,dps)})
                    save()
                report['stage']=f'production source I={I} ell={ell} node={k}';save()
                position={'node':int(k)} if k is not None and k==int(k) else {'point':mp.nstr(s,dps)}
                f,info=family.wave(coefficients,I,ell,**position,error_target=mp.nstr(target/100,50))
                delta=max(abs(f.real-arb(mp.nstr(values[-1].real,dps))).abs_upper(),abs(f.imag-arb(mp.nstr(values[-1].imag,dps))).abs_upper())
                change=abs(values[1]-values[0])
                row.update(production_real_ball=f.real.str(50),production_imag_ball=f.imag.str(50),production_quadrature=info,
                           two_order_difference=mp.nstr(change,50),source_difference_upper=delta.str(50),
                           source_difference_upper_dyadic=str(delta.upper().fmpq()),
                           passed=bool(change<target and delta<arb(mp.nstr(target,50))))
                if s>4:
                    h=arb.pi()*(1-4/arb(mp.nstr(s,dps))).sqrt()*f
                    S=1+acb(0,1)*h;slack=2*h.imag-(h*h.conjugate()).real
                    from .arbaudit import _enclosure
                    row.update(S_real=_enclosure(S.real),S_imag=_enclosure(S.imag),eta=_enclosure(abs(S)),
                               unitarity_slack=_enclosure(slack),
                               unitarity_verdict='pass' if slack>=0 else 'fail' if slack<0 else 'inconclusive')
                else:row['unitarity_verdict']='not_applicable_below_threshold'
                save()
        report['status']='completed';report['stage']='finished';report['passed']=all(r['passed'] for r in report['results'])
    except (TimeoutError,KeyboardInterrupt) as exc:
        report['status']='partial';report['error']=str(exc) or 'interrupted'
    except Exception as exc:
        report['status']='failed';report['error']=repr(exc);raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,previous);save()
    return report
