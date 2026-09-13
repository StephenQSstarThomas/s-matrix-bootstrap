"""Assembly and caching of the partial-wave operator rows for one (M, L).

Split out of :mod:`smatrix_bootstrap.sdp.problem` so that operator assembly and
conic-model construction stay separate responsibilities (and each file stays
inside the repository's size convention).
"""
from __future__ import annotations

import numpy as np

from . import constraints as C
from .grid import kappa, lambda_rescale, s_nodes
from .projector import Layout, PartialWaveOperator, ells_for

PROJECTION_S = 3.0            # the plane is (f00(3), f11(3)), section 4.1
BASIS_PACKING = 'T0,sigma1,sigma2,rho1-row-major,rho2-upper-unscaled-v1'


class Operators:
    """Cached partial-wave rows for one (M, L)."""

    _cache: dict = {}

    def __new__(cls, M: int, L: int, prescription='mixed-pv', dps=40):
        key=(M,L,prescription,dps if prescription=='sine-cardinal' else None)
        if key in cls._cache:
            return cls._cache[key]
        self = super().__new__(cls)
        if prescription=='sine-cardinal':self._init_sine(M,L,dps)
        elif prescription=='mixed-pv':self._init(M,L)
        else:raise ValueError(prescription)
        self.prescription=prescription
        cls._cache[key] = self
        return self

    def _init_sine(self,M,L,dps):
        from .sine import SineFamily
        from flint import arb
        import time
        self.M,self.L=M,L
        self.op=self.precise_source=SineFamily(M,L,dps)
        p=self.op;self.lay=p.lay
        self.s=np.array([float(x.mid()) for x in p.s]);self.kap=np.array([float(x.mid()) for x in p.kap])
        self.index=[(I,e) for I in (0,1,2) for e in ells_for(I,L)]
        self.precise_h=np.empty((3*M*L,2,self.lay.n),dtype=object)
        started=time.monotonic()
        for a,(I,ell) in enumerate(self.index):
            for k in range(M):self.precise_h[a*M+k]=p.rows(I,ell,node=k)*p.kap[k]
            if a%3==0 or a+1==len(self.index):
                print(f'analytic source waves {a+1}/{len(self.index)} in {time.monotonic()-started:.1f}s',flush=True)
        self.h_re=np.array([[float(v.mid()) for v in row] for row in self.precise_h[:,0]])
        self.h_im=np.array([[float(v.mid()) for v in row] for row in self.precise_h[:,1]])
        self.nu_measured=np.maximum(abs(self.h_re).max(axis=1),abs(self.h_im).max(axis=1))
        self.row_norm_raw=self.nu_measured.copy()
        self.lambda_analytic=np.concatenate([lambda_rescale(self.s,e)**2 for I,e in self.index])
        self.lambda_agreement=self.nu_measured/self.lambda_analytic
        self.set_cone_scaling('none')
        sub=[p.rows(0,0,point=3)[0],p.rows(1,1,point=3)[0]]
        for sj in (arb(1)/2,arb(1),arb(3)/2,arb(2)):
            f0,f1,f2=(p.rows(I,e,point=sj)[0] for I,e in ((0,0),(1,1),(2,0)))
            sub.extend([f0-3*(2*sj-1)/(sj-4)*f1,f2-3*(2-sj)/(sj-4)*f1])
        self.precise_sub=np.array(sub,dtype=object)
        sub_float=np.array([[float(v.mid()) for v in row] for row in sub])
        self.f_proj={'f00':sub_float[0],'f11':sub_float[1]};self.chi_rows=sub_float[2:]
        self.gram_rows={}
        for e,I in ((0,0),(1,1)):
            start=self.index.index((I,e))*M
            self.gram_rows[e]=self.h_re[start:start+M],self.h_im[start:start+M]
        self.basis=None

    def _init(self, M: int, L: int) -> None:
        self.M, self.L = M, L
        self.op = PartialWaveOperator(M)
        self.lay = Layout(M)
        self.s = s_nodes(M)
        self.kap = kappa(self.s)
        self.index = [(I, e) for I in (0, 1, 2) for e in ells_for(I, L)]
        n = self.lay.n
        nw = len(self.index)
        re = np.empty((nw, M, n))
        im = np.empty((nw, M, n))
        lam = np.empty((nw, M))
        for a, (I, e) in enumerate(self.index):
            lam[a] = lambda_rescale(self.s, e)
            for k in range(M):
                r = self.op.rows(I, e, self.s[k], k)
                re[a, k], im[a, k] = r[0], r[1]
        # Fold kappa in.  The cone may then be rescaled by any positive
        # per-row Lambda -- |h|^2 <= 2 Im h  <=>  |h/L|^2 <= 2 (Im h)/L^2, just
        # multiply through by L^2 > 0 -- which is exact for every choice.
        #
        # Store the original rows independently of solver conditioning. The
        # SDPB ModelSpec selects rownorm by default; any positive scale is a
        # congruence and changes no exact inequality. Tiny raw rows cannot be
        # assumed inactive merely from their coefficient magnitude.
        hre = re.reshape(-1, n) * np.repeat(self.kap[None, :], nw, 0).reshape(-1, 1)
        him = im.reshape(-1, n) * np.repeat(self.kap[None, :], nw, 0).reshape(-1, 1)
        nu = np.maximum(np.abs(hre).max(axis=1), np.abs(him).max(axis=1))
        self.lambda_analytic = lam.reshape(-1) ** 2
        self.nu_measured = nu
        ok = self.lambda_analytic > 1e-7
        self.lambda_agreement = nu[ok] / self.lambda_analytic[ok]
        self.h_re, self.h_im = hre, him
        self.row_norm_raw = np.maximum(np.abs(hre), np.abs(him)).max(axis=1)
        self.set_cone_scaling("none", 0.0)
        # unscaled S0 / P1 rows, needed by the Gram blocks
        self.gram_rows = {}
        for ell, I in ((0, 0), (1, 1)):
            a = self.index.index((I, ell))
            self.gram_rows[ell] = (re[a] * self.kap[:, None], im[a] * self.kap[:, None])
        # subthreshold rows
        self.f_proj = {"f00": self.op.rows(0, 0, PROJECTION_S)[0],
                       "f11": self.op.rows(1, 1, PROJECTION_S)[0]}
        self.chi_rows = []
        for sj in C.CHIRAL_POINTS:
            r01, r21 = C.chiral_ratios(sj)
            f00 = self.op.rows(0, 0, sj)[0]
            f11 = self.op.rows(1, 1, sj)[0]
            f20 = self.op.rows(2, 0, sj)[0]
            self.chi_rows.append(f00 - r01 * f11)
            self.chi_rows.append(f20 - r21 * f11)
        self.chi_rows = np.array(self.chi_rows)
        self.basis = None

    def build_basis(self, tol: float = 1e-12, mask=None, include_gram=False) -> np.ndarray:
        """SVD of canonical unscaled effective rows; a numerical truncation.

        Equal effective masks give identical SVD inputs across cone scaling and
        replacement of primary disks by Grams. A small projection residual is
        not an exact rank certificate or a bound on arbitrary omitted amplitudes.
        """
        R = self.basis_rows(mask,include_gram)
        _, sv, Vt = np.linalg.svd(R, full_matrices=False)
        k = int((sv > tol * sv[0]).sum())
        self.basis_singular_values = sv
        self.basis_tol = tol
        self.basis = Vt[:k].T
        self.basis_projection_error = float(np.max(np.abs(
            R - (R @ self.basis) @ self.basis.T)))
        return self.basis

    def basis_rows(self,mask=None,include_gram=False):
        """Canonical native order, retaining precisely disk-mask union Gram rows."""
        keep=np.ones(len(self.index)*self.M,dtype=bool) if mask is None else np.asarray(mask,dtype=bool).copy()
        if keep.shape!=(len(self.index)*self.M,):raise ValueError('Basis disk mask dimensions mismatch')
        if include_gram:
            for wave in ((0,0),(1,1)):
                start=self.index.index(wave)*self.M;keep[start:start+self.M]=True
        R=np.vstack([self.h_re[keep],self.h_im[keep],self.f_proj['f00'][None,:],
                     self.f_proj['f11'][None,:],self.chi_rows])
        return R/np.maximum(np.abs(R).max(axis=1,keepdims=True),1e-300)

    def adopt_basis(self,basis,mask=None,include_gram=False):
        """Use the exact saved coordinates; measure their action on fresh rows."""
        R=self.basis_rows(mask,include_gram);self.basis=basis
        self.basis_projection_error=float(np.max(np.abs(R-(R@basis)@basis.T)))
        return basis

    def set_cone_scaling(self, mode: str, drop_tiny: float = 0.0) -> None:
        """Choose the (exact) per-row rescaling of the unitarity cone."""
        if mode == "none":
            lam2 = np.ones_like(self.nu_measured)
        elif mode == "centrifugal":                  # Lambda_ell(s)^2 of (task 5)
            lam2 = self.lambda_analytic.copy()
            if np.any(lam2 <= 0):
                raise ValueError("analytic centrifugal scale underflow; increase assembly precision")
        elif mode == "rownorm":                      # measured row magnitude
            lam2 = np.where(self.nu_measured > 0.0, self.nu_measured, 1.0)
        else:
            raise ValueError(mode)
        lam1 = np.sqrt(lam2)
        self.cone_scaling = mode
        self.P_re = self.h_re / lam1[:, None]
        self.P_im = self.h_im / lam1[:, None]
        self.R_im = self.h_im / lam2[:, None]
        if drop_tiny > 0.0:
            # Entries this far below their row's scale cannot move a constraint
            # by more than ``drop_tiny`` relative, i.e. far below tol_feas; they
            # only destroy sparsity and conditioning of the KKT system.  The
            # exact operator is kept in ``h_re``/``h_im`` so every solution is
            # re-verified against the unsparsified constraints.
            for Mx in (self.P_re, self.P_im, self.R_im):
                sc = np.abs(Mx).max(axis=1, keepdims=True)
                Mx[np.abs(Mx) < drop_tiny * np.maximum(sc, 1e-300)] = 0.0
        self.sparsify = drop_tiny
        self.row_norm_scaled = np.maximum(np.abs(self.P_re), np.abs(self.R_im)).max(axis=1)
        self.row_scale = lam2


def load_saved_basis(spec,source_report):
    """Authenticate coordinates only; no constraints, solution or feasibility transfer."""
    import hashlib,io,json
    from pathlib import Path
    if not spec.reduce_basis:raise ValueError('--basis-source-report requires --reduce-basis')
    path=Path(source_report).resolve();record=json.loads(path.read_text())
    if (record.get('solver')!='SDPB' or record.get('status')!='numerically_accepted'
            or record.get('accepted') is not True
            or record.get('verification',{}).get('primal_feasible') is not True
            or record.get('convergence',{}).get('solver_optimal') is not True):
        raise ValueError('Basis source must be a completed accepted SDPB leaf')
    old=record.get('spec',{});pmp=record.get('pmp',{});meta=record.get('basis',{})
    for key in ('M','scattering_prescription'):
        if old.get(key)!=getattr(spec,key):raise ValueError('Basis source differs in '+key)
    if record['verification'].get('source_audit',{}).get('scattering_prescription')!=spec.scattering_prescription:
        raise ValueError('Basis source lacks matching source acceptance')
    projector='src/smatrix_bootstrap/sdp/projector.py'
    packing_sha=hashlib.sha256(Path(__file__).with_name('projector.py').read_bytes()).hexdigest()
    if (record.get('source_sha256',{}).get(projector)!=packing_sha
            or meta.get('packing',BASIS_PACKING)!=BASIS_PACKING):
        raise ValueError('Basis source amplitude packing is not authenticated')
    if record.get('source_sha256',{}).get('src/smatrix_bootstrap/sdp/grid.py')!=hashlib.sha256(Path(__file__).with_name('grid.py').read_bytes()).hexdigest():
        raise ValueError('Basis source grid and packing helpers are not authenticated')
    n=1+2*spec.M+spec.M**2+spec.M*(spec.M+1)//2;na=pmp.get('n_a')
    if (not old.get('reduce_basis') or not pmp.get('reduce_basis') or type(na) is not int
            or not 0<na<=n or pmp.get('n_vars')!=1+na+(4*spec.M if old.get('uv') else 0)):
        raise ValueError('Basis source variable layout mismatch')
    basis_path=path.parent/'basis.npy';raw=basis_path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
    if digest!=meta.get('sha256'):raise ValueError('Basis source hash mismatch')
    basis=np.load(io.BytesIO(raw),allow_pickle=False)
    if basis.dtype!=np.dtype('float64') or basis.shape!=(n,na) or not np.all(np.isfinite(basis)):
        raise ValueError('Basis source must be finite binary64 with the recorded layout')
    info={'source_report':str(path),'source_report_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
          'source_basis_path':str(basis_path),'source_basis_sha256':digest,'packing':BASIS_PACKING,
          'source_L':old['L'],'target_L':spec.L,'cross_L':old['L']!=spec.L,
          'source_tolerance':meta.get('tolerance'),'target_svd_tolerance_applied':False,
          'feasibility_transferred':False,'constraints_reassembled':True,
          'scope':'exact saved coordinate matrix; fresh source rows, PMP and verification required'}
    return basis,info


def recorded_basis_key(record):
    """An absent reduced-basis identity never identifies a shared finite problem."""
    if not record.get('spec',{}).get('reduce_basis'):return 'unreduced'
    meta=record.get('basis') or {};identity=record.get('basis_identity')
    if identity is not None and not identity.get('verified'):return None
    key=meta.get('sha256')
    if (not isinstance(key,str) or len(key)!=64 or any(v not in '0123456789abcdef' for v in key)
            or meta.get('packing',BASIS_PACKING)!=BASIS_PACKING):return None
    return key


def basis_identity(record,report_path,point=None):
    """Inspect the actual saved basis, including leaf paths in aggregate reports."""
    import hashlib
    from pathlib import Path
    spec=record.get('spec',{});key=recorded_basis_key(record)
    if not spec.get('reduce_basis'):
        return {'verified':True,'key':'unreduced','scope':'full original density coordinates'}
    info={'verified':False,'recorded_sha256':key,'key':None}
    if key is None:return dict(info,reason='reduced basis identity missing or invalid')
    root=Path(report_path).resolve().parent;candidates=[]
    if point:candidates.append(root/point/'basis.npy')
    candidates.append(root/'basis.npy')
    pmp=record.get('pmp',{})
    if pmp.get('path'):candidates.append(Path(pmp['path']).parent/'basis.npy')
    path=next((p for p in candidates if p.is_file()),None)
    if path is None:return dict(info,reason='saved basis file unavailable')
    actual=hashlib.sha256(path.read_bytes()).hexdigest();info.update(path=str(path.resolve()),actual_sha256=actual)
    if actual!=key:return dict(info,reason='saved basis hash differs from accepted report')
    basis=np.load(path,mmap_mode='r',allow_pickle=False);M=spec.get('M',0)
    expected=(1+2*M+M*M+M*(M+1)//2,pmp.get('n_a'))
    if basis.dtype!=np.dtype('float64') or basis.shape!=expected or not np.all(np.isfinite(basis)):
        return dict(info,reason='saved basis layout or values are invalid')
    packing=record.get('source_sha256',{}).get('src/smatrix_bootstrap/sdp/projector.py')
    if packing!=hashlib.sha256(Path(__file__).with_name('projector.py').read_bytes()).hexdigest():
        return dict(info,reason='amplitude packing is not authenticated')
    if record.get('source_sha256',{}).get('src/smatrix_bootstrap/sdp/grid.py')!=hashlib.sha256(Path(__file__).with_name('grid.py').read_bytes()).hexdigest():
        return dict(info,reason='grid and packing helpers are not authenticated')
    return dict(info,verified=True,key=key,shape=list(basis.shape),packing=BASIS_PACKING,
                scope='exact saved matrix identity; no omitted-direction or rank certificate')
