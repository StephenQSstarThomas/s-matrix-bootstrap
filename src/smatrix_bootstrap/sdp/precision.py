"""Arbitrary-precision finite 2309 rows and their coordinate projection.

The geometry is rebuilt before projection; no float operator is promoted.
Production Q uses a guarded recurrence, independently audited by the Arb
hypergeometric formula and original angular integrals in the tests.
"""
from __future__ import annotations

import math
import time

import numpy as np
from flint import arb, arb_mat, ctx
from mpmath import mp

from .legendreq import _q_column_mp, _required_dps
from .projector import Layout


class PrecisionRows:
    def __init__(self, M, dps=50):
        self.M, self.dps = M, dps
        self.bits = math.ceil((dps+15)*math.log2(10))
        ctx.prec = self.bits
        self.lay = Layout(M)
        phi = [arb.pi()*(2*j+1)/(2*M) for j in range(M)]
        self.s = np.array([8/(1+p.cos()) for p in phi], dtype=object)
        self.w = np.array([8*p.sin()/(M*(1+p.cos())**2) for p in phi], dtype=object)
        self.kap = np.array([arb.pi()*((x-4)/x).sqrt() for x in self.s], dtype=object)
        def kt(n):
            if n == 0 or n % 2 == 0:
                return arb(0)
            a = arb.pi()*n/(2*M)
            return a.cos()/(M*a.sin())
        self.K = np.array([[kt(i+j-2*M-1)-kt(i-j) for j in range(1,M+1)]
                           for i in range(1,M+1)], dtype=object)
        self.on = self.K + (self.w/self.s)[None, :]
        self._qcache = {}

    def _q(self, z, ell):
        ztext = z.mid().str(self.dps+12, radius=False)
        key = (ztext, ell)
        if key not in self._qcache:
            guard = _required_dps(abs(float(z.mid())), ell)
            with mp.workdps(self.dps + guard + 15):
                zm = mp.mpf(ztext)
                value = _q_column_mp(zm, ell)[ell]
                self._qcache[key] = arb(mp.nstr(value, self.dps+12))
        return self._qcache[key]

    def rows(self, isospin, ell, *, point=None, node=None):
        """Real and imaginary coefficient rows, in the original density packing."""
        ctx.prec = self.bits
        if (point is None) == (node is None):
            raise ValueError("Specify exactly one of subthreshold point or native node")
        s = self.s[node] if node is not None else (arb(point) if isinstance(point,arb) else arb(str(point)))
        if node is None and not (s > 0 and s < 4):
            raise ValueError("Off-node rows are defined only in the subthreshold region")
        M, lay, w = self.M, self.lay, self.w
        d, sign, p0 = s-4, (-1)**ell, 2 if ell == 0 else 0
        q = np.array([4*self._q(1+2*x/d, ell)/d for x in self.s], dtype=object)
        Q = w*q
        D = np.array([[(q[i]+sign*q[j])/(d+self.s[i]+self.s[j]) for j in range(M)]
                      for i in range(M)], dtype=object)*w[:,None]*w[None,:]
        if node is None:
            cr = [w/(self.s-s), np.zeros(M, dtype=object)]
        else:
            im = np.zeros(M, dtype=object); im[node] = 1
            cr = [self.on[node], im]
        return assemble_rows(lay,cr,Q,D,isospin,ell)

    def projected(self, rows, basis):
        ctx.prec = self.bits
        if basis is None:
            return np.asarray(rows, dtype=object)
        matrix = basis if isinstance(basis, arb_mat) else arb_mat(
            [[arb(float(v)) for v in row] for row in basis])
        return np.array((arb_mat(np.asarray(rows, dtype=object).tolist())*matrix).tolist(), dtype=object)

    def uv_data(self, spec):
        """Current coefficients rebuilt from (2.33), (2.56), (3.72–75)."""
        from . import constraints as C
        ctx.prec = self.bits
        pi, s0, M = arb.pi(), arb(3600)/49, self.M
        def kin(e, s):
            return ((6*pi).sqrt()/(16*pi**3)/s.root(4)*((s-4)/4).root(4) if e == 0 else
                    (4*pi/3).sqrt()/(8*pi**3)/s.root(4)*((s-4)/4).root(4)**3)
        g = {e:np.array([kin(e,s) for s in self.s],dtype=object) for e in (0,1)}
        idx = np.array([i for i,s in enumerate(self.s) if s>s0],dtype=int)
        mq = (arb(113)/2800 if spec.m_q == C.M_Q else
              ((arb(4)**2+arb('7.3')**2)/2).sqrt()/140 if spec.m_q == C.M_Q_RMS else arb(str(spec.m_q)))
        eps = arb(str(spec.eps_ff))
        bound = {0:(2*mq*mq*eps).sqrt(),1:(eps/2).sqrt()}
        fk = {e:np.array([kin(e,s0) if spec.ff_frozen_at_s0 else g[e][i] for i in idx],dtype=object)
              for e in (0,1)}
        target, tol, moments = {}, {}, {}
        for e,wave in ((0,'S0'),(1,'P1')):
            for n in C.MOMENTS[e]:
                value = (arb('3.09e-8')*(arb('27.38')/(n+2)+(arb('.61') if n==0 else 0))
                         if e==0 else -arb('4.34e-6')*(-arb('13.26')/(n+2)+(arb('.41') if n==0 else 0)))
                target[wave,n] = value*s0**(n+2)
                tol[wave,n] = (arb('.002') if spec.sr_caliber=='SR-a' else
                                abs(target[wave,n])*arb('.1' if spec.sr_caliber=='SR-b' else '.2'))
                moments[e,n] = np.array([pi*self.w[i]*s**n*g[e][i]**2 if s<=s0 else arb(0)
                                          for i,s in enumerate(self.s)],dtype=object)
        return self.K, idx, bound, fk, target, tol, moments


def benchmark(M, basis_path, outdir, dps=40, digits=30):
    """Small source/coordinate test on the actual M50 basis, without optimization."""
    from pathlib import Path
    import hashlib
    from .arbaudit import ArbAudit
    from .sdpb import write_json
    start = time.monotonic()
    p = PrecisionRows(M, dps)
    basis = np.load(basis_path)
    v = arb_mat([[arb(float(x)) for x in row] for row in basis])
    build_time = time.monotonic()-start
    cases = [(0,0,3,None),(0,0,None,M//2),(1,1,None,M//2),(2,0,None,M//2)]
    projected, timing = [], []
    for I, ell, point, node in cases:
        ts = time.monotonic(); rows = p.rows(I,ell,point=point,node=node)
        source_time = time.monotonic()-ts; ts = time.monotonic()
        result = p.projected(rows,v); projected.append(result)
        timing.append({"I": I, "ell": ell, "point": point, "node": node,
                       "source_seconds": source_time, "projection_seconds": time.monotonic()-ts})
    report = {"M": M, "dps": dps, "working_bits_with_guard": p.bits,
              "basis_shape": list(basis.shape), "basis_build_seconds": build_time,
              "rows": timing, "seconds": time.monotonic()-start,
              "basis_sha256": hashlib.sha256(Path(basis_path).read_bytes()).hexdigest()}
    # Stress exactly the huge coefficients which exposed the old double path.
    yp = Path(basis_path).parent/'out/y.txt'
    if yp.exists():
        ys = [arb(t) for t in yp.read_text().splitlines()[1:] if t.strip()]
        av = arb_mat([[x] for x in ys[:basis.shape[1]]]); cv = v*av
        c = [cv[i,0] for i in range(cv.nrows())]; lay = p.lay
        r2 = [[arb(0) for _ in range(M)] for _ in range(M)]
        for (i,j), val in zip(zip(*lay.triu), c[lay.r2]):
            r2[i][j] = r2[j][i] = val
        rr = c[lay.r1]
        dens = {'T0':c[0], 'sigma1':c[lay.s1], 'sigma2':c[lay.s2],
                'rho1':[rr[i*M:(i+1)*M] for i in range(M)], 'rho2':r2}
        audit = ArbAudit(M,1,384); errors = []
        for case, result in zip(cases,projected):
            I,ell,point,node = case
            ref = audit.partial_wave(dens,I,ell,node or 0,point=point)
            emitted = arb_mat([[arb(x.mid().str(digits,radius=False)) for x in row] for row in result])*av
            errors.append(max(float(abs(emitted[0,0]-ref.real).abs_upper()),
                              float(abs(emitted[1,0]-ref.imag).abs_upper())))
        report['stress_test'] = {'output_digits':digits,'absolute_errors':errors,
            'max_error':max(errors), 'budget':1e-8, 'passed':max(errors)<1e-8,
            'y_sha256':hashlib.sha256(yp.read_bytes()).hexdigest()}
    dest = Path(outdir); dest.mkdir(parents=True, exist_ok=True)
    write_json(dest/'report.json', report)
    write_json(dest/'projected_rows.json', [[[x.mid().str(dps,radius=False) for x in row]
                                            for row in r] for r in projected])
    return report


def assemble_rows(lay, cr, Q, D, isospin, ell):
    """The common isospin/partial-wave algebra, with explicit leg kernels."""
    sign,p0=(-1)**ell,2 if ell==0 else 0
    out = np.zeros((3,2,lay.n), dtype=object)
    for p in range(2):
        c, real = cr[p], p == 0
        out[0,p,0] = p0 if real else 0
        out[0,p,lay.s1] = p0*c
        out[0,p,lay.s2] = (1+sign)*Q if real else 0
        out[0,p,lay.r1] = ((1+sign)*np.outer(c,Q)).ravel()
        out[0,p,lay.r2] = lay.fold_sym(D) if real else 0
        out[1,p,0] = p0 if real else 0
        out[1,p,lay.s1] = Q if real else 0
        out[1,p,lay.s2] = p0*c + (sign*Q if real else 0)
        out[1,p,lay.r1] = (np.outer(Q,c) + (D if real else 0)).ravel()
        out[1,p,lay.r2] = lay.fold_sym(sign*np.outer(c,Q))
        out[2,p,0] = p0 if real else 0
        out[2,p,lay.s1] = sign*Q if real else 0
        out[2,p,lay.s2] = p0*c + (Q if real else 0)
        out[2,p,lay.r1] = (sign*np.outer(Q,c) + (sign*D if real else 0)).ravel()
        out[2,p,lay.r2] = lay.fold_sym(np.outer(c,Q))
    total = {0: 3*out[0]+out[1]+out[2], 1: out[1]-out[2], 2: out[1]+out[2]}[isospin]
    # Avoid Python float division of structural integer entries.
    return np.array([[arb(v)/4 for v in row] for row in total], dtype=object)
