"""One conformal-polynomial amplitude and authenticated analytic row reuse.

The stored source basis consists of z^n, z_s^n(z_t^m+z_u^m), and
(z_t^n z_u^m+z_t^m z_u^n)/2. No source producer is loaded or executed.
"""
from pathlib import Path
from fractions import Fraction
from flint import arb, acb, arb_mat, acb_mat, ctx
import numpy as np
from functools import lru_cache
from . import read_json, write_json, workspace, digest, decode_real_ball, encode_real_ball
from .kernels import PVSourceRows, density_row_to_cflat, canonical_waves
from .analytic import cardinal_transform


def stored_path(value):
    path=Path(value)
    if str(path).startswith('/tmp/'):
        path=workspace()/'results/analytic_source_storage'/path.relative_to('/tmp')
    elif not path.is_absolute():
        path=workspace()/path
    return path


@lru_cache(maxsize=4)
def transform_data(M,bits):
    from .operators import midpoint_grid
    D=acb_mat(cardinal_transform(M));x,w=midpoint_grid(M)
    b=acb_mat(M,1,[wi/xi for xi,wi in zip(x,w)])
    return D,D.transpose(),b,b*b.transpose()


def mode_row_to_cardinal(row,M):
    """Dual of the full DST and exact subtraction map, retaining every column."""
    D,DT,b,bb=transform_data(M,ctx.prec)
    g=list(map(acb,row));free=1+2*M
    if len(g)!=free+M*M+M*(M+1)//2:raise ValueError('Complete conformal row required')
    h1=acb_mat(1,M,g[1:M+1])*D;h2=acb_mat(1,M,g[M+1:free])*D
    R=DT*acb_mat(M,M,g[free:free+M*M])*D
    Q=acb_mat(M,M);at=free+M*M
    for i in range(M):
        for j in range(i,M):
            Q[i,j]=Q[j,i]=g[at];at+=1
    Q=DT*Q*D
    R+=2*h1.transpose()*b.transpose()+b*h2+2*g[0]*bb
    Q+=(h2.transpose()*b.transpose()+b*h2)/2+g[0]*bb
    return [g[0]]+(h1+g[0]*b.transpose()).entries()+(
        h2+2*g[0]*b.transpose()).entries()+R.entries()+[
            Q[i,j] for i in range(M) for j in range(i,M)]


class CardinalSourceRows(PVSourceRows):
    """Exact analytic rows, in unsubtracted nodal-density C_flat coordinates."""
    def __init__(self,M=50,L=10,bits=384,order=24,subtracted=False,registry=None):
        super().__init__(M,L,bits,order,subtracted)
        if subtracted:raise ValueError('Use unsubtracted analytic-cardinal coordinates')
        self.registry_path=Path(registry) if registry else (workspace()/'results/analytic_source_storage/analytic_source_registry_stage/registry.json' if M==50 else None)
        self.registry=read_json(self.registry_path) if self.registry_path else None
        if self.registry is not None:
            r=self.registry
            if r.get('schema')!='smatrix-complete-analytic-source-registry-v1' or r.get('scattering_rank_normalization')!='f partial waves; G=pi*beta*f is an exact positive row scaling at each physical node':
                raise ValueError('Recognized analytic registry and f normalization required')
            if r['M']!=M or r['N']!=M:raise ValueError('Registry does not cover this M')
            if r['canonical_wave_order']!=list(map(list,canonical_waves(len(r['canonical_wave_order'])//3))) or len(r['physical_nodes'])!=M or r['shape']!=[2*M*len(r['canonical_wave_order'])+11,self.coefficient_count]:
                raise ValueError('Complete canonical analytic registry required')
        self.quadrature_accuracy={};self.continuous_cache={}
        self.coordinate_name='unsubtracted analytic-cardinal C_flat';self.cache={};self.supplement=None
        self.inputs={str(self.registry_path.resolve()):dict(sha256=digest(self.registry_path))} if self.registry_path else {}
        self.metadata.update(prescription='analytic-cardinal',coordinate_system=self.coordinate_name,
            analytic_function_family_defined=True,nyquist_retained=True,integral_enclosures=True,
            analyticity_status='One conformal polynomial from all M sine-cardinal density modes',
            angular_projection='Authenticated rows where available; otherwise Arb Gauss-Legendre with ellipse remainder',angular_order=max(64,order),
            spectral_quadrature='Exact Cauchy transform of finite sine density',
            physical_evaluation='One analytic family at positive real energies',
            source_registry=str(self.registry_path.resolve()) if self.registry_path else None,source_registry_sha256=digest(self.registry_path) if self.registry_path else None,
            infinity_row_meaning='Unsubtracted constant; recorded but no endpoint equality imposed')
        self.metadata.pop('residue_identity',None)


    def attach_preparation(self,preparation):
        path=Path(preparation);meta=read_json(path/'report.json')['amplitude_model']
        if (meta['prescription'],meta['M'],meta['L'],meta['coordinate_system'])!=('analytic-cardinal',self.M,self.L,self.coordinate_name):
            raise ValueError('Prepared analytic rows do not match this provider')
        if meta.get('source_registry_sha256')!=self.metadata.get('source_registry_sha256'):
            raise ValueError('Prepared analytic registry differs')
        with np.load(path/'amplitude.npz') as data:
            H=data['rows'];E=data['radius_upper'];n=(len(H)-11)//2
            if n!=self.physical_count:
                from .sampling import prepared_sampling
                prepared_sampling(data,self.M,self.L);ids=np.r_[np.arange(self.physical_count),n+np.arange(self.physical_count),np.arange(2*n,2*n+11)];H,E=H[ids],E[ids]
        if H.shape!=self.shape or E.shape!=self.shape or not np.isfinite(H).all() or not np.isfinite(E).all() or np.any(E<0):
            raise ValueError('Complete finite analytic coefficient boxes required')
        self.prepared=(H,E)
        for file in (path/'report.json',path/'amplitude.npz'):
            self.inputs[str(file.resolve())]=dict(sha256=digest(file))
        self.metadata['evaluation_source']='Prepared analytic coefficient enclosures'

    def exact_energy(self,s):
        """Resolve only CLI float aliases; exact/interval inputs keep their meaning."""
        if type(s) is float:return next((x for x in self.x if float(x)==s),arb(s))
        return arb(s.numerator)/s.denominator if isinstance(s,Fraction) else s

    def _read(self,item,prefix='rows'):
        path=stored_path(item[prefix+'_path']);sha=item[prefix+'_sha256']
        data=read_json(path,sha);self.inputs[str(path)]=dict(sha256=sha)
        labels=([['constant',0,0]]+[[name,n,0] for name in ('single_s','single_tu') for n in range(1,self.M+1)]
            +[['double_s_tu',n,m] for n in range(1,self.M+1) for m in range(1,self.M+1)]
            +[['double_tu',n,m] for n in range(1,self.M+1) for m in range(n,self.M+1)])
        if data.get('N')!=self.M or data.get('source_basis_labels')!=labels:
            raise ValueError('Conformal basis labels/order mismatch')
        return data

    def _node(self,node):
        if node in self.cache:return self.cache[node]
        out={}
        if self.registry is not None:
            data=self._read(self.registry['physical_nodes'][node])
            if data['node_zero_based']!=node or len(data['source_basis_labels'])!=self.coefficient_count:
                raise ValueError('Analytic source identity mismatch')
            for record in data['rows']:
                I,ell=record['isospin'],record['ell']
                if (I,ell) not in self.waves:continue
                if (I,ell) in out or self.registry['canonical_wave_order'][record['canonical_wave_index']]!=[I,ell]:
                    raise ValueError('Duplicate or mislabeled analytic wave')
                raw=[acb(decode_real_ball(v[0]),decode_real_ball(v[1])) for v in record['coefficients']]
                row=mode_row_to_cardinal(raw,self.M)
                # These columns follow exactly from orthogonality and the nodal jump.
                for j in range(1+2*self.M):
                    value=arb(0)
                    if ell==0 and I in (0,2):
                        if j==1+node and I==0:value=arb(3)/2
                        if j==1+self.M+node:value=arb(1)
                    if not row[j].imag.contains(value):raise ArithmeticError('Analytic free-column identity failed')
                    row[j]=acb(row[j].real,value)
                out[I,ell]=row
        missing=[w for w in self.waves if w not in out]
        if missing:
            from .quadrature import polynomial_rows
            generated,accuracy=polynomial_rows(self.x[node],self.M,missing,node,self.order)
            out.update(generated);self.quadrature_accuracy[str(node)]=accuracy
        if set(out)!=set(self.waves):raise ValueError('Incomplete retained analytic wave coverage')
        self.cache={node:out}
        return out

    def row(self,s,ell,isospin,node=None,spins=None):
        ctx.prec=self.working_bits
        native=next((j for j,x in enumerate(self.x) if s is x),None)
        explicit_alias=type(s) is float and node is not None and 0<=node<self.M and s==float(self.x[node])
        s=arb(s.numerator)/s.denominator if isinstance(s,Fraction) else arb(s)
        if not s.is_finite() or not s>0 or not (s<4 or s==4 or s>4):raise ValueError('Finite positive s not crossing threshold required')
        if type(ell) is not int or ell<0 or isospin not in (0,1,2) or ell%2!=int(isospin==1):raise ValueError('Allowed nonnegative partial wave required')
        if node is None:node=native
        if node is not None:
            if type(node) is not int or not 0<=node<self.M or not (native==node or explicit_alias):raise ValueError('Explicit native energy object or float alias required')
            if getattr(self,'prepared',None) is not None:
                if not 0<=isospin<3 or ell%2!=int(isospin==1) or not 0<=ell//2<self.L:
                    raise ValueError('Retained canonical wave required')
                H,E=self.prepared;k=node*3*self.L+isospin*self.L+ell//2;n=self.physical_count
                return [acb(arb(float(a),float(e)),arb(float(b),float(f)))
                    for a,b,e,f in zip(H[k],H[n+k],E[k],E[n+k])]
            return self._node(node)[isospin,ell]
        if s==4:
            raw=[arb(0)]*self.coefficient_count
            if ell==0 and isospin in (0,2):
                raw[0]=arb(5)/2 if isospin==0 else arb(1)
                raw[1:self.M+1]=[arb(3)/2 if isospin==0 else arb(0)]*self.M
                raw[self.M+1:1+2*self.M]=[arb(1)]*self.M
            return mode_row_to_cardinal(raw,self.M)
        item=next((v for v in (self.registry['subthreshold_inputs'] if self.registry else []) if s==arb(Fraction(v['exact_energy']).numerator)/Fraction(v['exact_energy']).denominator),None)
        if item is not None:
            for record in self._read(item)['rows']:
                if (record['isospin'],record['ell'])==(isospin,ell):
                    return mode_row_to_cardinal([acb(decode_real_ball(v[0]),decode_real_ball(v[1])) for v in record['coefficients']],self.M)
        key=(repr(encode_real_ball(s)),isospin,ell)
        if key not in self.continuous_cache:
            from .quadrature import polynomial_rows
            rows,_=polynomial_rows(s,self.M,[(isospin,ell)],order=self.order)
            self.continuous_cache[key]=rows[isospin,ell]
        return self.continuous_cache[key]

    def subthreshold_row(self,s,ell,isospin):
        row=self.row(s,ell,isospin)
        if any(not v.imag.contains(0) for v in row):raise ArithmeticError('Subthreshold analyticity/reality failed')
        return [v.real for v in row]

    def iter_rows(self,canonical_ids=None):
        ids=set(range(self.shape[0])) if canonical_ids is None else set(canonical_ids);n=self.physical_count
        if any(type(i) is not int or not 0<=i<self.shape[0] for i in ids):
            raise ValueError('Canonical analytic row outside retained range')
        if self.registry is None:
            yield from super().iter_rows(ids)
            return
        for node in range(self.M):
            if not any(i in ids for i in list(range(node*3*self.L,(node+1)*3*self.L))+list(range(n+node*3*self.L,n+(node+1)*3*self.L))):continue
            rows=self._node(node)
            for k,(I,ell) in enumerate(self.waves):
                for offset,part in ((0,'real'),(n,'imag')):
                    index=offset+node*3*self.L+k
                    if index in ids:yield index,dict(kind='scattering_'+part,node_zero_based=node,isospin=I,ell=ell),[getattr(v,part) for v in rows[I,ell]]
        item=dict(rows_path=self.registry['supplement_path'],rows_sha256=self.registry['supplement_sha256'])
        if any(i>=2*n for i in ids):
            data=self._read(item);original_n=self.registry['shape'][0]-11
            for record in data['rows']:
                index=2*n+record['canonical_row']-original_n
                if index not in ids:continue
                if index==2*n+8:row=[arb(1)]+[arb(0)]*(self.coefficient_count-1)
                else:row=[v.real for v in mode_row_to_cardinal(list(map(decode_real_ball,record['coefficients'])),self.M)]
                yield index,dict(kind=record['kind']),row

def prepare_cardinal(args,progress):
    if getattr(args,'unitarity_energies',[]):
        from .sampling import prepare_additional
        return prepare_additional(args,progress)
    from .operators import midpoint_grid
    from .certificates import verify_scattering_layout
    source=CardinalSourceRows(args.nodes,args.waves,args.bits,args.angular_order,args.subtracted,args.source_registry)
    H=np.full(source.shape,np.nan);errors=np.full(source.shape,np.nan);n=source.physical_count;seen=set()
    from .quadrature import reuse_prepared_rows
    rows=source.iter_rows() if getattr(args,'preparation',None) is None else reuse_prepared_rows(source,args.preparation)
    for index,meta,row in rows:
        if index in seen or not 0<=index<source.shape[0] or len(row)!=source.shape[1]:
            raise ValueError('Duplicate, invalid or incomplete analytic row')
        seen.add(index)
        mid=np.asarray([float(v.mid()) for v in row]);H[index]=mid
        errors[index]=[float(np.nextafter(float(abs(v-arb(float(x))).upper()),np.inf)) if not v.is_zero() else 0. for v,x in zip(row,mid)]
        if index<n and index%(3*source.L)==0:progress(stage='analytic_source_reuse',node=meta['node_zero_based'],M=source.M)
    if seen!=set(range(source.shape[0])):raise ValueError('Incomplete analytic row coverage')
    k=np.repeat([arb.pi()*(1-4/x).sqrt() for x in source.x],3*source.L)
    verify_scattering_layout(H,k,source.M,source.L)
    path=args.output/'amplitude.npz'
    np.savez_compressed(path,rows=H,radius_upper=errors,energies=list(map(float,source.x)),
        waves=source.waves,coefficient_labels=np.asarray(source.coefficient_labels,dtype=str))
    bound=args.density_limit or 100*(source.coefficient_count-1-2*source.M)
    return dict(status='prepared',amplitude_model=source.metadata,inputs=source.inputs,
        restrictions=dict(infinity='free',density_limit=bound,density_rule='declared 100*double_density_count',density_norm='actual-rho L4'),
        operator=dict(path=str(path),shape=list(H.shape),storage='Float64 midpoints plus outward bounds around exact analytic coefficient enclosures'),
        optimization_performed=False,old_feasibility_transferred=False,
        source_registry=str(source.registry_path) if source.registry_path else None,quadrature_accuracy=source.quadrature_accuracy,maximum_operator_radius_upper=float(errors.max()))


def audit_cardinal(args,progress):
    if args.profile_kind=='analytic-window':return audit_direct_samples(args,progress)
    if args.current_preparation is not None:
        from .certificates import analytic_joint_audit,analytic_fiber_audit
        fixed=read_json(args.coefficients.parent/'report.json').get('fixed_amplitude')
        return (analytic_fiber_audit if fixed else analytic_joint_audit)(args,progress)
    if args.preparation is not None:
        from .certificates import analytic_scattering_audit
        return analytic_scattering_audit(args,progress)
    from .analytic import CardinalAmplitude,_complex_record
    if args.coefficients is None:raise ValueError('Complete frozen coefficients required')
    from . import read_cflat
    M=args.nodes;values,input_meta=read_cflat(args.coefficients,M)
    source=CardinalSourceRows(M,args.waves,args.bits,registry=args.source_registry)
    amp=CardinalAmplitude(values,M);records=[]
    for node in sorted(set((0,M//2,M-1))):
        for I,ell in ((0,0),(2,0),(1,1)):
            row=source.row(source.x[node],ell,I,node=node)
            value=sum((v*c for v,c in zip(row,amp.raw)),acb(0))
            independent=amp.partial_wave(source.x[node],node,I,ell)
            match=(value-independent).contains(0)
            records.append(dict(node=node,isospin=I,ell=ell,source=_complex_record(value),
                independent=_complex_record(independent),match=match))
            progress(stage='source_angular_check',node=node,isospin=I,match=match)
    inputs=dict(source.inputs)
    if source.registry_path:inputs[str(source.registry_path)]=dict(sha256=digest(source.registry_path))
    inputs[str(args.coefficients)]=dict(sha256=digest(args.coefficients))
    result=dict(status='analytic_source_verified' if all(v['match'] for v in records) else 'analytic_source_mismatch',
        amplitude_model=source.metadata,inputs=inputs,checks=records,
        source_basis='z^n; symmetrized tu products divided by two',input_prescription=input_meta['prescription'],
        old_feasibility_transferred=False,optimization_performed=False)
    write_json(args.output/'source_audit.json',result)
    return result


def read_cflat(path,M,key='coefficients',*,energies=None,waves=None,prescription=None):
    from . import AMPLITUDE_COORDINATES
    labels=dict(AMPLITUDE_COORDINATES,**{'finite-sine-cardinal':'unsubtracted sine-cardinal C_flat'})
    path=Path(path);metadata=dict(path=str(path),key=key)
    if path.suffix=='.npz' or key!='coefficients':raise ValueError('CG retired; use C_flat JSON')
    record=read_json(path);vector=np.asarray(record['coefficients'],float);metadata['format']='declared_JSON_C_flat'
    kind=next((k for k,v in labels.items() if v==record.get('coordinates')),None)
    if kind is None or record.get('prescription',kind)!=kind or prescription is not None and kind!=prescription:raise ValueError('Invalid declared coordinates/prescription')
    metadata.update(coordinates=labels[kind],prescription=kind,sha256=digest(path),legacy_data_only=kind=='finite-sine-cardinal')
    if vector.shape!=(1+2*M+M*M+M*(M+1)//2,) or not np.isfinite(vector).all():
        raise ValueError('Finite complete C_flat required')
    return vector,metadata


def scattering_input(args,M,target):
    """Known nodal data may seed a newly audited analytic problem; no old provider."""
    from . import AMPLITUDE_COORDINATES
    if target not in AMPLITUDE_COORDINATES:raise ValueError('Unknown active target prescription')
    c,record=read_cflat(args.coefficients,M,getattr(args,'coefficient_key','coefficients'))
    kind=record['prescription'];legacy=kind=='finite-sine-cardinal'
    if kind!=target and not (legacy and target=='analytic-cardinal' and args.prescription=='analytic-cardinal'):
        raise ValueError('Cross-model input requires explicit analytic-cardinal for legacy sine data; PV transfer forbidden')
    record.update(target_prescription=target,coordinate_transform='identity',explicit_legacy_replay=legacy,old_feasibility_or_bound_transferred=False,
        purpose='Numerical C only; recompute feasibility and support on the prepared target H')
    return c,record


def candidate_weights(path,n,record,coefficients):
    """Read arbitrary finite covectors; old values, margins and bounds are ignored."""
    path=Path(path);kind=record['prescription'];legacy=kind=='finite-sine-cardinal'
    if kind!=record['target_prescription'] and not (legacy and record['target_prescription']=='analytic-cardinal' and record.get('explicit_legacy_replay')):
        raise ValueError('Candidate target/source prescription mismatch')
    with np.load(path,allow_pickle=False) as old:
        tagged='prescription' in old;declared=str(old['prescription']) if tagged else kind if legacy else None
        if declared!=kind or 'coordinates' in old and str(old['coordinates'])!=record['coordinates']:
            raise ValueError('Candidate covectors use another input prescription')
        if any(bool(old.get(k,False)) for k in ('polynomial_dual_present','global_high_spin_enforced','tail_conditions_applied')) or any(np.any(old.get(k,[])) for k in ('tail_kR','tail_kI','polynomial_dual')):
            raise ValueError('Tail/polynomial duals cannot enter a sampled recheck')
        raw=[old[k] for k in ('kR','kI','y')]
        if any(np.iscomplexobj(v) for v in raw):raise ValueError('Real candidate covectors required')
        weights=[np.asarray(v,float).copy() for v in raw]
        if any(v.shape!=(size,) or not np.isfinite(v).all() for v,size in zip(weights,(n,n,8))):
            raise ValueError('Complete finite kR/kI/y on every target H row required')
        if 'coefficients' in old and not np.array_equal(old['coefficients'],coefficients):raise ValueError('Candidate/JSON C mismatch')
        if legacy and 'coefficients' not in old:raise ValueError('Legacy candidate must retain its complete C')
    meta=dict(path=str(path),sha256=digest(path),input_prescription=kind,target_prescription=record['target_prescription'],
        scattering_rows=n,identity_source='checkpoint tag' if tagged else 'legacy C_flat JSON',old_feasibility_or_bound_transferred=False,old_values_and_margins_used=False,
        purpose='Arbitrary candidate covectors; fresh target-H residual and support calculation')
    return (*weights,meta)


def audit_direct_samples(args,progress):
    """Independent rational-circle integrals for requested direct evaluations."""
    from .analytic import CardinalAmplitude,_complex_record
    from .model import unitarity_margin
    ctx.prec=args.bits;records=[];inputs={}
    for folder in args.profile_runs:
        file=Path(folder)/'evaluation.json';d=read_json(file);meta=d['amplitude_model'];path=Path(d['coefficients'])
        c,_=read_cflat(path,meta['M'],prescription='analytic-cardinal');amp=CardinalAmplitude(c,meta['M'])
        for value in args.energies:
            j=d['energies'].index(value);s=arb(value)
            for k,(I,ell) in enumerate(d['waves']):
                old=d['f_enclosures'][j][k];saved=acb(decode_real_ball(old[0]),decode_real_ball(old[1]))
                direct=amp.partial_wave(s,None,I,ell);margin=unitarity_margin(direct,s)
                match=(direct-saved).contains(0)
                records.append(dict(role=d['selection'].get('role'),energy=float(value),isospin=I,ell=ell,match=match,
                    direct=_complex_record(direct),stored=_complex_record(saved),unitarity_margin=margin.str(30),strict_violation=bool(margin<0)))
                progress(stage='independent_off_node_integral',energy=float(value),isospin=I,match=match)
        inputs.update({str(p.resolve()):dict(sha256=digest(p)) for p in (file,path)})
    if not records:raise ValueError('At least one independent sample required')
    result=dict(status='independent_analytic_samples_verified' if all(r['match'] for r in records) else 'independent_analytic_mismatch',
        checks=records,inputs=inputs,optimization_performed=False,coefficients_modified=False,
        scope='Independent full-amplitude integral versus Gauss source rows; only requested samples')
    write_json(args.output/'source_audit.json',result)
    return result
