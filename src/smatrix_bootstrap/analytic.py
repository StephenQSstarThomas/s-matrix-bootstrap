"""Analytic cardinal reconstruction and audits of the native PV approximation.

This is an independently derived comparison, not a continuation of PV values
or a transfer of historical sine feasibility/rank results. Angular integrals
use a rational circle parameter so Arb can certify holomorphy directly.
"""
from pathlib import Path
from flint import arb, arb_mat, acb, acb_mat, ctx
import numpy as np
from . import read_json, write_json, digest

def upper_float(x):
    return float(np.nextafter(float(x.upper()),np.inf))


def cardinal_transform(M):
    if type(M) is not int or M < 1:
        raise ValueError('Positive integer M required')
    return arb_mat([[arb(1 if n == M else 2)/M *
        (arb(n*(2*j+1))/(2*M)).sin_pi() for j in range(M)]
        for n in range(1,M+1)])


def modes(z,M,folded=False):
    z=acb(z);powers=[];v=acb(1)
    for _ in range(2*M if folded else M):
        v*=z;powers.append(v)
    if not folded:
        return powers
    denominator=1+powers[-1]
    return [(powers[n-1]+powers[2*M-n-1])/denominator
        for n in range(1,M)]+[2*powers[M-1]/denominator]


def reconstruction_identity(M,s):
    """Exact folded DST representation of the midpoint Cauchy quadrature."""
    from .operators import midpoint_grid
    x,w=midpoint_grid(M);s=arb(s)
    if not s < 4:
        raise ValueError('This identity compares real off-cut values s<4')
    z=(2-(4-s).sqrt())/(2+(4-s).sqrt());D=cardinal_transform(M)
    folded=acb_mat(1,M,modes(z,M,True))*acb_mat(D)
    polynomial=acb_mat(1,M,modes(z,M))*acb_mat(D)
    rational=[wi*s/(xi*(xi-s)) for xi,wi in zip(x,w)]
    errors=[folded[0,j]-rational[j] for j in range(M)]
    return dict(M=M,s=float(s),folded_identity_enclosed=all(v.contains(0) for v in errors),
        maximum_identity_error_upper=max(upper_float(abs(v)) for v in errors),
        polynomial_difference_upper=max(upper_float(abs(polynomial[0,j]-rational[j])) for j in range(M)))


class CardinalAmplitude:
    """One entire finite polynomial in each conformal variable, with full C."""
    def __init__(self,coefficients,M):
        from .quotient import convert_subtraction
        from .kernels import coefficient_blocks
        self.M=M;self.raw=list(map(arb,coefficients));D=cardinal_transform(M)
        converted=convert_subtraction(self.raw,M,to_subtracted=True)
        c,R,Q=coefficient_blocks(converted,M)
        self.T=acb(c[0]);self.T_unsub=self.raw[0]
        self.a=acb_mat((D*arb_mat(M,1,c[1:M+1])).transpose())
        self.b=acb_mat((D*arb_mat(M,1,c[M+1:2*M+1])).transpose())
        self.R=acb_mat(D*R*D.transpose());self.Q=acb_mat(D*Q*D.transpose())

    def amplitude(self,zs,zt,zu):
        p,q,r=[acb_mat(self.M,1,modes(z,self.M)) for z in (zs,zt,zu)]
        return (self.T+(self.a*p)[0,0]+(self.b*(q+r))[0,0]+
            (p.transpose()*self.R*(q+r))[0,0]+(q.transpose()*self.Q*r)[0,0])

    def uniform_error_bound(self,s0,quadrature_nodes=None):
        """All physical s<=s0: polynomial versus folded crossed kernels.

        Both compared boundary prescriptions use the same direct polynomial.
        This is not a definition or feasibility certificate for off-node PV.
        """
        r=((arb(s0).sqrt()-2)/(arb(s0).sqrt()+2)).upper();M=self.M
        N=M if quadrature_nodes is None else quadrature_nodes
        if type(N) is not int or N<M:
            raise ValueError('Quadrature nodes must be an integer at least the mode count')
        if not 0<r<1:
            raise ValueError('A finite matching scale above threshold is required')
        powers=[r**n for n in range(1,M+1)]
        error=[(r**(2*N-n)+r**(2*N+n))/(1-r**(2*N)) for n in range(1,M+1)]
        R,RT,Q=self.R,self.R.transpose(),self.Q;result={}
        for I,B,J,K in ((0,self.a+4*self.b,3*R+RT+Q,R+RT+3*Q),
                        (2,self.a+self.b,RT+Q,R+RT),
                        (1,self.a-self.b,RT-Q,R-RT)):
            linear=[abs(B[0,j]).upper()+sum(abs(J[i,j]).upper() for i in range(M)) for j in range(M)]
            bound=sum(linear[j]*error[j] for j in range(M))
            bound+=sum(abs(K[i,j]).upper()*(powers[i]*error[j]+error[i]*powers[j]+error[i]*error[j])/2
                       for i in range(M) for j in range(M))
            result[I]=bound.upper()
        return result

    def partial_wave(self,s,node,isospin,ell,folded=False):
        """Validated full-amplitude integral; node=None uses the supplied energy."""
        if isospin not in (0,1,2) or type(ell) is not int or ell<0 or ell%2!=int(isospin==1):
            raise ValueError('Allowed isospin/parity and nonnegative spin required')
        s=arb(s)
        if not s.is_finite() or not s>4 or node is not None and (type(node) is not int or not 0<=node<self.M):raise ValueError('Physical energy or native node required')
        if node is None:
            root=-acb(4-s).sqrt();zs=(2-root)/(2+root)
        else:
            angle=arb(2*node+1)/(2*self.M);native=4/(angle/2).cos_pi()**2
            if not s.overlaps(native):raise ValueError('Energy enclosure must cover the indexed native point')
            s=native;zs=acb(angle.cos_pi(),angle.sin_pi())
        p=acb_mat(1,self.M,modes(zs,self.M))
        R,RT,Q=self.R,self.R.transpose(),self.Q
        if isospin==0:
            constant=5*self.T+((3*self.a+2*self.b)*p.transpose())[0,0]
            linear=self.a+4*self.b+3*p*R+p*RT+p*Q;quadratic=R+RT+3*Q
        elif isospin==2:
            constant=2*self.T+(2*self.b*p.transpose())[0,0]
            linear=self.a+self.b+p*RT+p*Q;quadratic=R+RT
        else:
            constant=acb(0);linear=self.a-self.b+p*RT-p*Q;quadratic=R-RT
        radius=(s+4).sqrt();lo=2/(radius+s.sqrt());hi=s.sqrt()/(radius+2)
        def integrand(q,analytic):
            # Entire rational expression in q; poles give nonfinite balls.
            x=radius*(1-q*q)/(1+q*q);y=2*radius*q/(1+q*q)
            t=4-x*x;mu=1+2*t/(s-4)
            a=acb_mat(self.M,1,modes((2-x)/(2+x),self.M,folded))
            b=acb_mat(self.M,1,modes((2-y)/(2+y),self.M,folded))
            value=constant+(linear*(a-b if isospin==1 else a+b))[0,0]+(a.transpose()*quadratic*b)[0,0]
            jacobian=4*(s+4)*q*(1-q*q)/((s-4)*(1+q*q)**3)
            return value*jacobian*mu.legendre_p(ell)
        value=acb.integral(integrand,lo,hi,rel_tol=arb(2)**-65,
            abs_tol=arb(2)**(-85 if ell<=1 else -min(300,ctx.prec-32)),eval_limit=20000,depth_limit=40)
        if not value.is_finite():
            raise ArithmeticError('Validated angular integration did not resolve')
        return value


def _complex_record(z):
    z=acb(z)
    midpoint=[float(z.real.mid()),float(z.imag.mid())]
    return dict(real=z.real.str(30),imaginary=z.imag.str(30),
        midpoint=midpoint,
        radius_upper=max(upper_float(abs(z.real-midpoint[0])),upper_float(abs(z.imag-midpoint[1]))),
        arithmetic_radius_upper=max(upper_float(z.real.rad()),upper_float(z.imag.rad())))


def audit_reconstruction(args,progress):
    from .run import support_data
    from .kernels import PVSourceRows
    from .model import scattering_matrix,unitarity_margin
    H,kap,ss,ww,M,L,_,_=support_data(args,progress)
    if args.comparison_manifest is not None:
        data=read_json(args.comparison_manifest)
        if data.get('status')!='analytic_reconstruction_audited' or (data['M'],data['L'])!=(M,L):
            raise ValueError('Matching completed analytic audit required')
        for file,record in data['inputs'].items():
            if digest(file)!=record['sha256']:raise ValueError('Frozen analytic-audit input changed')
        with ctx.workprec(args.bits):
            from .operators import midpoint_grid
            native,weights=midpoint_grid(M)
            for record in data['profiles']:
                amplitude=CardinalAmplitude(read_json(record['coefficients'])['coefficients'],M)
                bounds=amplitude.uniform_error_bound((arb('1.2')/arb('.14'))**2)
                record['uniform_below_matching']={str(I):dict(f_difference_bound=v.str(30),S_difference_bound=(arb.pi()*v).str(30)) for I,v in bounds.items()}
                record['quadrature_repair_bound']=quadrature_repair_bound(amplitude,native[-1])
                for sample in record['samples']:
                    residue=-weights[sample['node']]*arb(sample['PV']['imaginary'])
                    sample.update(PV_offcut_continuation_pole_residue=residue.str(30),nonzero_pole_certified=bool(not residue.contains(0)))
                    diff=sample['difference'];delta=acb(arb(diff['real']),arb(diff['imaginary']))
                    sample['delta_S_absolute_upper']=upper_float(arb.pi()*(1-4/native[sample['node']]).sqrt()*abs(delta))
                    for name in ('PV','PV_integral','analytic_cardinal','difference'):
                        value=sample[name];value['radius_upper']=max(upper_float(abs(arb(value[part])-value['midpoint'][j])) for j,part in enumerate(('real','imaginary')))
        data.update(physical_input_audit=input_consistency(M),angular_integrals_reused_from=str(args.comparison_manifest),uniform_bound_scope='Same direct boundary modes; crossed-kernel error for 4<s<=s0, all ell; no unitarity transfer')
        data['local_logic_counterexamples']=local_logic_counterexamples()
        plot_audit(args.output,data)
        write_json(args.output/'analytic_audit.json',data)
        return data
    if not args.profile_runs:
        raise ValueError('Explicit unchanged coefficient directories required')
    source=PVSourceRows(M,L,args.bits);nodes=source.x
    low=[j for j,s in enumerate(nodes) if s<=(arb('1.2')/arb('.14'))**2]
    selected=sorted(set([0,M//2,low[-1],M-1]));records=[];inputs={}
    with ctx.workprec(args.bits):
        identities=[reconstruction_identity(M,s) for s in (3,-4,-100)]
        for path in args.profile_runs:
            path=Path(path);file=path/'coefficients.json';data=read_json(file)
            if data['coordinates']!='unsubtracted PV-midpoint C_flat':
                raise ValueError('Complete unsubtracted PV coordinates required')
            inputs[str(file.resolve())]=dict(sha256=digest(file))
            amplitude=CardinalAmplitude(data['coefficients'],M)
            row=dict(coefficients=str(file.resolve()),T0=amplitude.T_unsub.str(30),samples=[])
            row['analytic_limit_unitarity_defect']={
                'S0':(-arb.pi()**2*(arb(5)/2*amplitude.T_unsub)**2).str(30),
                'S2':(-arb.pi()**2*amplitude.T_unsub**2).str(30)}
            for node in selected:
                s=nodes[node]
                for I,ell in ((0,0),(2,0),(1,1)):
                    direct=sum((a*b for a,b in zip(source.row(s,ell,I,node=node),amplitude.raw)),acb(0))
                    integral=amplitude.partial_wave(s,node,I,ell,True)
                    consistent=amplitude.partial_wave(s,node,I,ell,False)
                    S=scattering_matrix(direct,s);Snew=scattering_matrix(consistent,s)
                    margin=unitarity_margin(consistent,s)
                    item=dict(node=node,energy_gev=float(arb('.14')*s.sqrt()),isospin=I,ell=ell,
                        PV=_complex_record(direct),PV_integral=_complex_record(integral),
                        PV_integral_matches_closed_projection=(direct-integral).contains(0),
                        analytic_cardinal=_complex_record(consistent),difference=_complex_record(consistent-direct),
                        delta_S_absolute_upper=upper_float(abs(Snew-S)),analytic_eta=str(abs(Snew)),
                        analytic_unitarity_defect=margin.str(30),analytic_disk='violated' if margin<0 else 'passed' if margin>=0 else 'unresolved')
                    row['samples'].append(item)
                    progress(stage='analytic_comparison',coefficients=str(file),node=node,isospin=I,ell=ell,
                        PV_integral_match=item['PV_integral_matches_closed_projection'],delta_S=item['delta_S_absolute_upper'])
                    write_json(args.output/'partial.json',dict(completed=records+[row]))
            current=path/'current.json'
            if current.is_file():
                inputs[str(current.resolve())]=dict(sha256=digest(current));im=read_json(current)['ImF']
                D=cardinal_transform(M);signs=arb_mat(1,M,[arb((-1)**n) for n in range(1,M+1)])
                row['form_factor_limits']=[(1+(signs*D*arb_mat(M,1,list(map(arb,v))))[0,0]).str(30) for v in im]
            records.append(row)
        result=dict(status='analytic_reconstruction_audited',M=M,L=L,inputs=inputs,identities=identities,
            native_node_rule='first, phi midpoint, last below s0, last',native_nodes=selected,profiles=records,
            angular_integrals='Validated Arb integrals of rational circle parametrization',
            original_model_changed=False,old_feasibility_transferred=False,coefficients_refitted=False,
            paper_curve_inputs=False,continuum_unitarity_certified=False,
            scope='Independent analytic reconstruction compared with the native PV approximation; no new optimization')
        write_json(args.output/'analytic_audit.json',result)
        return result


def plot_audit(out,data):
    import matplotlib.pyplot as plt
    from . import save_figure
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    for p,style in zip(data['profiles'],('o-','s--')):
        role=Path(p['coefficients']).parent.name
        for I,color,name in ((0,'tab:red','S0'),(2,'tab:green','S2'),(1,'tab:blue','P1')):
            values=[v for v in p['samples'] if v['isospin']==I]
            axes[0].loglog([v['energy_gev'] for v in values],[v['delta_S_absolute_upper'] for v in values],style,color=color,label=role+' '+name,ms=4)
            high=values[-1];eta=arb(high['analytic_eta'])
            axes[1].scatter(name+' / '+role,float(eta.mid()),color=color)
    axes[0].axvline(1.2,color='.5',ls=':');axes[0].set(xlabel='E [GeV]',ylabel='Upper bound on |S_cardinal - S_PV|',title='Frozen coefficients, validated integrals');axes[0].legend(fontsize=7)
    axes[1].axhline(1,color='.5',ls=':');axes[1].set(yscale='log',ylabel='|S_cardinal| at last native node',title='Analytic lift fails the original disk test');axes[1].tick_params(axis='x',rotation=60)
    save_figure(fig,Path(out)/'analytic_comparison')


def input_consistency(M):
    """Audit stated inputs, not phase curves; a moment-only counterexample."""
    from . import UVConfig
    from .model import fesr_data
    base=fesr_data(M,UVConfig());other={}
    for rule in ('arithmetic-mean','rms'):
        data=fesr_data(M,UVConfig(moment_source='eq250',mq_rule=rule))
        other[rule]=dict(targets=data['targets_raw'],relative_to_printed=[a/b-1 for a,b in zip(data['targets_raw'],base['targets_raw'])],
            scope='Same Eq2.50 formula with a declared mass rule; RMS is not proven to describe the same scalar current')
    W=[list(map(arb,row)) for row in base['weights_raw'][2:]];s=list(map(arb,base['nodes']));target=list(map(arb,base['targets_raw'][2:]))
    matrix=arb_mat([[sum(w),sum(x*y for x,y in zip(w,s))] for w in W]);ab=matrix.solve(arb_mat(2,1,target))
    spectrum=[ab[0,0]+ab[1,0]*x for x in s];residual=arb_mat(W)*arb_mat(M,1,spectrum)-arb_mat(2,1,target)
    ratio=arb('.0001')/(arb('.14')**2*arb('.092')**2)
    return dict(printed=base,eq250=other,GMOR=dict(ratio=ratio.str(30),fpi_from_stated_condensate_gev=(arb('.0001')/arb('.14')**2).sqrt().str(30),
        scope='Stated central values do not satisfy exact leading-order GMOR; no change to inputs'),
        monotone_P1_moment_counterexample=dict(a=ab[0,0].str(30),b=ab[1,0].str(30),positive=all(v>0 for v in spectrum),
            increasing=bool(ab[1,0]>0),moment_residuals=[v.str(30) for v in residual.entries()],
            determinant=matrix.det().str(30),
            exact_center_match=all(v.contains(0) for v in residual.entries()),
            scope='R_i=a+b*s_i matches both finite P1 moments; not a counterexample to all joint constraints'))


def quadrature_repair_bound(amplitude,smax,tolerance='1e-6'):
    """Certify a quadrature budget for this fixed analytic polynomial.

    N changes quadrature work, not modes, variables, constraints or coefficients.
    This bounds polynomial versus N-point folded crossed modes at every physical
    s<=smax and every ell. It does not make the old PV optimum analytically feasible.
    """
    M=amplitude.M;tol=arb(tolerance);cache={}
    if not tol>0 or not tol.is_finite():
        raise ValueError('A finite positive absolute S error required')
    def bound(N):
        if N not in cache:
            cache[N]={I:arb.pi()*v for I,v in amplitude.uniform_error_bound(smax,N).items()}
        return max(v.upper() for v in cache[N].values())
    lo=M-1;hi=M
    while not bound(hi)<tol:
        lo=hi;hi*=2
        if hi>64*M:
            raise ArithmeticError('Bounded quadrature budget exhausted')
    while hi-lo>1:
        mid=(hi+lo)//2
        if bound(mid)<tol:hi=mid
        else:lo=mid
    bound(hi)
    return dict(mode_count=M,maximum_s=arb(smax).str(30),
        maximum_energy_gev=(arb('.14')*arb(smax).sqrt()).str(30),
        absolute_S_tolerance=tol.str(30),quadrature_nodes_sufficient=hi,
        S_bounds={str(I):v.str(30) for I,v in cache[hi].items()},
        previous_bound=None if lo<M else bound(lo).str(30),
        native_quadrature_bound=bound(M).str(30),evaluations=len(cache),
        minimality='Smallest N satisfying this sufficient majorant, not a necessary quadrature order',
        scope='Fixed coefficients, all 4<s<=smax and all spins; same direct polynomial',
        model_changed=False,old_feasibility_transferred=False)


def local_logic_counterexamples():
    """Exact small examples refute implications, not the full QCD bootstrap."""
    from flint import fmpq_mat
    B=fmpq_mat([[1,0,1],[0,1,1]]);G=B.transpose()*B
    initial=fmpq_mat([['1','1/2','1'],['1/2','1','1'],['1','1','4/3']])
    return dict(
        PSD_boundary_not_rank_one=dict(gram=[[str(G[i,j]) for j in range(3)] for i in range(3)],
            factor_B=[[str(B[i,j]) for j in range(3)] for i in range(2)],
            determinant=str(G.det()),rank=G.rank(),eigenvalues=[0,1,3],
            spectral_excess=1,inelasticity_defect=1,
            proof='G=B^T B, rank(B)=2. det(G)=0 while R-|F|^2=1 and |S|=0.',
            scope='Local Gram cone; not a crossing/FESR counterexample'),
        Watson_functionals_not_identical=dict(
            old_gram=[[str(initial[i,j]) for j in range(3)] for i in range(3)],
            old_determinant=str(initial.det()),old_rank=initial.rank(),
            old_S='1/2',old_F='1',old_h='i/2',
            coefficients_on_Re_h_Im_h=dict(eq_28=['1/2','-1'],eq_215_v2=['0','-1']),
            proof='v2=(-Q,1,0), Q=F/F*. v2^dag G v2=2-2 Re(Q* S).',
            scope='2505v3 Eq2.8 and Eq2.15 differ for a valid old local Gram; no full model changed'),
        rank_one_fixed_points=dict(parameter='0<alpha<pi',
            S='exp(2 i alpha)',F='exp(i alpha)',R=1,
            h='2 sin(alpha) exp(i alpha)',
            proof='Both proposed directions fix every member on the local disk boundary.',
            scope='No generic uniqueness theorem from this update alone; not multiple full GTB solutions'))
