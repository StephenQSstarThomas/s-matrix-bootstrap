from fractions import Fraction
import pytest
flint=pytest.importorskip('flint')
import mpmath as mp
from flint import arb,acb,ctx
from smatrix_bootstrap_newton.kernels import density_labels, pv_matrix, physical_node_row, offcut_row
from smatrix_bootstrap_newton.operators import midpoint_grid, angular_kernels
def as_arb(value):
    value=Fraction(value);return arb(value.numerator)/value.denominator
def inside(answer,reference):
    if isinstance(answer,acb):
        return answer.real.contains(arb(mp.nstr(mp.re(reference),110)).mid()) and answer.imag.contains(arb(mp.nstr(mp.im(reference),110)).mid())
    return answer.contains(arb(mp.nstr(reference,110)).mid())
@pytest.mark.parametrize('s,ell',[(Fraction(1,2),0),(1,1),(Fraction(3,2),2),(2,3),(3,4),(9,0)])
def test_closed_subtracted_angular_kernels_against_independent_direct_integrals(s,ell):
    with ctx.workprec(256),mp.workdps(120):
        ss=mp.mpf(s.numerator)/s.denominator if isinstance(s,Fraction) else mp.mpf(s)
        nodes=[arb(5),arb(13)];weights=[arb(2),arb(3)]
        kernels=angular_kernels(as_arb(s),nodes,weights,ell)
        def h(i,v):return [2,3][i]*(1/([mp.mpf(5),mp.mpf(13)][i]-v)-1/[mp.mpf(5),mp.mpf(13)][i])
        t=lambda mu:-(ss-4)*(1-mu)/2
        u=lambda mu:4-ss-t(mu)
        for i in range(2):
            truth=mp.quad(lambda mu:mp.legendre(ell,mu)*h(i,t(mu))/4,[-1,0,1])
            assert inside(kernels['J'][i],truth)
        for i,j in ((0,0),(0,1),(1,1)):
            truth=mp.quad(lambda mu:mp.legendre(ell,mu)*h(i,t(mu))*h(j,u(mu))/4,[-1,0,1])
            if ell%2 and i==j:
                assert kernels['U'][i][j].is_zero() and abs(truth)<mp.mpf('1e-110')
            else:assert inside(kernels['U'][i][j],truth)
def test_sine_hilbert_map_and_unit_imaginary_jump_retain_Nyquist():
    with ctx.workprec(256):
        M=4;nodes,weights=midpoint_grid(M);K=pv_matrix(M)
        for n in range(1,M+1):
            for k in range(M):
                value=sum((K[k,i]*(arb((2*i+1)*n)/(2*M)).sin_pi() for i in range(M)),arb(0))
                exact=(arb((2*k+1)*n)/(2*M)).cos_pi()
                assert (value-exact).contains(0)
        row=physical_node_row(2,nodes,weights,K,0,0)
        assert all(isinstance(value,acb) for value in row)
        for i in range(M):
            assert row[1+i].imag== (arb(3)/2 if i==2 else 0)
            assert row[1+M+i].imag== (arb(1) if i==2 else 0)
        assert row[0].real==arb(5)/2 and row[0].imag.is_zero()
        assert len(density_labels(50))==3876
def direct_amplitude(hs,ht,hu,coefficients,M):
    at=0;value=coefficients[at];at+=1
    for i in range(M):value+=coefficients[at]*hs[i];at+=1
    for i in range(M):value+=coefficients[at]*(ht[i]+hu[i]);at+=1
    for i in range(M):
        for j in range(M):value+=coefficients[at]*hs[i]*(ht[j]+hu[j]);at+=1
    for i in range(M):
        for j in range(i,M):
            value+=coefficients[at]*(ht[i]*hu[i] if i==j else ht[i]*hu[j]+ht[j]*hu[i]);at+=1
    return value
@pytest.mark.parametrize('isospin,ell',[(0,0),(0,2),(1,1),(1,3),(2,0),(2,2)])
@pytest.mark.parametrize('physical',[False,True])
def test_rational_calibration_row_against_raw_isospin_angular_amplitude(isospin,ell,physical):
    with ctx.workprec(256):
        M=2;nodes,weights=midpoint_grid(M);K=pv_matrix(M)
        with mp.workdps(120):
            phis=[mp.pi*(j+mp.mpf('.5'))/M for j in range(M)]
            xx=[4/mp.cos(phi/2)**2 for phi in phis];ww=[x*mp.tan(phi/2)/M for x,phi in zip(xx,phis)]
            index=1;ss=xx[index] if physical else mp.mpf(3)
            coeff=[mp.mpf((-1)**k*(k+1))/(k+2) for k in range(len(density_labels(M)))]
            cb=[arb((-1)**k*(k+1))/(k+2) for k in range(len(coeff))]
            def h(v):return [wi*(1/(xi-v)-1/xi) for xi,wi in zip(xx,ww)]
            hs=([sum(2*mp.cos(n*phis[index])*mp.sin(n*phis[i])/M for n in range(1,M))+mp.j*int(index==i) for i in range(M)] if physical else h(ss))
            def integrand(mu):
                tt=-(ss-4)*(1-mu)/2;uu=4-ss-tt;ht=h(tt);hu=h(uu)
                aa=direct_amplitude(hs,ht,hu,coeff,M);bb=direct_amplitude(ht,hs,hu,coeff,M)
                cc=direct_amplitude(hu,ht,hs,coeff,M)
                amplitude=(3*aa+bb+cc if isospin==0 else bb-cc if isospin==1 else bb+cc)
                return mp.legendre(ell,mu)*amplitude/4
            truth=mp.quad(integrand,[-1,0,1])
            row=(physical_node_row(index,nodes,weights,K,ell,isospin) if physical else offcut_row(arb(3),nodes,weights,ell,isospin))
            answer=sum((a*b for a,b in zip(row,cb)),acb(0) if physical else arb(0))
            assert inside(answer,truth)
def test_exact_threshold_limits_and_no_off_node_physical_pole_evaluation():
    with ctx.workprec(128):
        nodes=[arb(5),arb(9)];weights=[arb(1),arb(2)]
        kernels=angular_kernels(arb(4),nodes,weights,0)
        assert all(value.is_zero() for value in kernels['J'])
        assert all(value.is_zero() for row in kernels['U'] for value in row)
        with pytest.raises(ValueError,match='sub-cut|off-cut|0<s<4'):
            offcut_row(arb(5),nodes,weights,0,0)
def test_rational_subtraction_matches_full_density_coordinate_change():
    with ctx.workprec(256):
        M=2;nodes,weights=midpoint_grid(M);K=pv_matrix(M)
        a=[arb(2)/3,-arb(3)/5];b=[arb(1)/7,arb(4)/9];T=arb(5)/11
        R1=[[arb(1),arb(2)/3],[-arb(3)/4,arb(5)/6]]
        R2=[[arb(2)/5,-arb(1)/3],[-arb(1)/3,arb(7)/8]]
        shift=[w/x for x,w in zip(nodes,weights)]
        anew=[a[i]+2*sum(R1[i][j]*shift[j] for j in range(M)) for i in range(M)]
        bnew=[b[i]+sum((R1[j][i]+R2[i][j])*shift[j] for j in range(M)) for i in range(M)]
        Tnew=T+sum((a[i]+2*b[i])*shift[i] for i in range(M))+sum((2*R1[i][j]+R2[i][j])*shift[i]*shift[j] for i in range(M) for j in range(M))
        tail=[R1[i][j] for i in range(M) for j in range(M)]+[R2[i][j] for i in range(M) for j in range(i,M)]
        original=[T]+a+b+tail;transformed=[Tnew]+anew+bnew+tail
        for I,ell in [(0,0),(1,1),(2,0)]:
            for physical in [False,True]:
                rows=[(physical_node_row(0,nodes,weights,K,ell,I,subtracted=sub) if physical
                    else offcut_row(arb(3),nodes,weights,ell,I,subtracted=sub)) for sub in (False,True)]
                first=sum(x*y for x,y in zip(rows[0],original));second=sum(x*y for x,y in zip(rows[1],transformed))
                difference=first-second
                if physical:assert difference.real.contains(0) and difference.imag.contains(0)
                else:assert difference.contains(0)
def test_pv_source_free_constant_threshold_and_original_coordinate_change():
    from smatrix_bootstrap_newton.kernels import PVSourceRows
    with ctx.workprec(ctx.prec):
        original = PVSourceRows(2, 2, bits=256, order=24)
        subtracted = PVSourceRows(2, 2, bits=256, order=24, subtracted=True)
        # A constant amplitude remains a genuine variable; its high-energy value is recorded.
        coeff = [arb(1)]+[arb(0)]*(original.coefficient_count-1)
        values = original.evaluate(coeff, [3, 4, float(original.x[0])], [(0, 0), (1, 1), (2, 0)])
        for row in values:
            assert row == [acb(arb(5)/2), acb(0), acb(1)]
        assert original.infinity_row() == coeff
        assert not original.metadata['infinity_equality_imposed']
        # H=q+b; even a pure unsubtracted single density shifts T0 in q coordinates.
        coeff[1], coeff[0] = arb(1), arb(0)
        transformed = list(coeff); transformed[0] = original.b[0]
        a = original.evaluate(coeff, [3, 4, float(original.x[0])], [(0, 0), (1, 1), (2, 0)])
        b = subtracted.evaluate(transformed, [3, 4, float(original.x[0])], [(0, 0), (1, 1), (2, 0)])
        assert max(float(abs(x-y).upper()) for row, other in zip(a, b) for x, y in zip(row, other)) < 1e-65
        assert max(float(abs(v).upper()) for v in original.row(4, 3, 1)) == 0
def test_full_pv_density_coordinate_change_preserves_amplitude_and_density_norm():
    import numpy as np
    from smatrix_bootstrap_newton.quotient import convert_subtraction
    from smatrix_bootstrap_newton.kernels import PVSourceRows; from smatrix_bootstrap_newton.quotient import subtraction_rows, subtract_amplitude_coordinates
    from smatrix_bootstrap_newton.imaginary import density_fourth_power
    with ctx.workprec(ctx.prec):
        original,other=[PVSourceRows(3,2,bits=256,order=24,subtracted=sub) for sub in (False,True)]
        c=[arb((-1)**i)/2**(i+1) for i in range(original.coefficient_count)]
        d=convert_subtraction(c,3,to_subtracted=True);recovered=convert_subtraction(d,3,to_subtracted=False)
        assert all((a-b).contains(0) for a,b in zip(c,recovered)) and d[7:]==c[7:]
        assert density_fourth_power(c,3)==density_fourth_power(d,3)
        for s in [3,4,float(original.x[0])]:
            a,b=[source.evaluate(v,[s])[0] for source,v in ((original,c),(other,d))]
            assert max(float(abs(x-y).upper()) for x,y in zip(a,b))<1e-65
        asymptotic=sum((a*b for a,b in zip(other.infinity_row(),d)),arb(0));assert (asymptotic-c[0]).contains(0)
        scale=np.array([2 if family=='C_flat' and i!=j else 1 for family,i,j in original.coefficient_labels])
        rows=[np.tile(np.array([[float(v.mid()) for v in row] for _,_,row in source.iter_rows()]),(3,1))*scale for source in (original,other)]
        np.testing.assert_allclose(subtraction_rows(rows[0],3),rows[1],rtol=1e-12,atol=1e-14)
        mid=lambda values:np.array([np.longdouble(v.mid().str(40,radius=False)) for v in values])
        z=mid(c)/scale;zs=subtract_amplitude_coordinates(z,3,True);back=subtract_amplitude_coordinates(zs,3,False)
        assert np.max(abs(zs-mid(d)/scale))<np.longdouble('5e-19') and np.max(abs(back-z))<np.longdouble('5e-19')
        assert np.array_equal(zs[7:],z[7:])
        zero=z.copy();zero[0]=0;shifted=subtract_amplitude_coordinates(zero,3,True)
        assert shifted[0]!=0 and abs(subtract_amplitude_coordinates(shifted,3,False)[0])<np.longdouble('5e-19')
def test_eliminated_density_barrier_gradient_and_hessian_by_variation():
    import numpy as np
    from smatrix_bootstrap_newton.quotient import barrier_terms
    z = np.array([.2,.4,.1,-.2])
    R=np.array([[1.,0,0,0]]); I=np.array([[0.,1,0,0]])
    t=np.ones(1); C=np.zeros((8,4))
    args=(R,I,t,C,None,1.,2)
    result=barrier_terms(z,*args)
    gradient=-result['features'].T@result['rhs']
    hessian=result['features'].T@result['features']
    numerical=np.zeros(4); curvature=np.zeros((4,4)); step=1e-5
    for j in range(4):
        delta=np.eye(4)[j]*step
        plus,minus=barrier_terms(z+delta,*args),barrier_terms(z-delta,*args)
        numerical[j]=float((plus['value']-minus['value'])/(2*step))
        curvature[:,j]=(-plus['features'].T@plus['rhs']+minus['features'].T@minus['rhs'])/(2*step)
    np.testing.assert_allclose(gradient,numerical,rtol=1e-7,atol=1e-8)
    np.testing.assert_allclose(hessian,curvature,rtol=1e-7,atol=1e-8)
@pytest.mark.parametrize('eta', ['1', '.6', '1.1'])
def test_unitarity_from_independent_elastic_inelastic_and_violating_phase(eta):
    from smatrix_bootstrap_newton.model import scattering_matrix, unitarity_margin, scattering_gram
    with ctx.workprec(256):
        eta = arb(eta); phase = arb(2)/3; s = arb(9)
        expected_S = eta*acb(phase.cos(), phase.sin())
        f = (expected_S-1)/acb(0, arb.pi()*((s-4)/s).sqrt())
        margin, gram = unitarity_margin(f, s), scattering_gram(f, s)
        assert (scattering_matrix(f, s)-expected_S).contains(0)
        assert (margin-(1-eta**2)).contains(0)
        assert gram[0, 0] == gram[1, 1] == 1
        assert (gram[1, 0]-gram[0, 1].conjugate()).contains(0)
        assert (gram.det()-margin).contains(0)
        assert f.imag > 0  # Also true for eta>1: positivity alone is insufficient.
        if eta < 1: assert margin > 0 and gram.det().real > 0
        elif eta > 1: assert margin < 0 and gram.det().real < 0
        else: assert margin.contains(0) and gram.det().contains(0)
def test_unitarity_finite_partial_wave_threshold_limit_and_domain():
    from smatrix_bootstrap_newton.model import phase_space, scattering_matrix, unitarity_margin, scattering_gram
    f = acb(2, 3)
    assert phase_space(4).is_zero() and scattering_matrix(f, 4) == 1
    assert unitarity_margin(f, 4).is_zero() and scattering_gram(f, 4).det().is_zero()
    assert unitarity_margin(f, arb(4)+arb(2)**-30) > 0
    for s in [3, arb('4 +/- .1'), arb('inf')]:
        with pytest.raises(ValueError, match='Physical'): unitarity_margin(f, s)
    for check in [scattering_matrix, unitarity_margin, scattering_gram]:
        with pytest.raises(ValueError, match='finite partial wave'): check(acb('inf'), 4)
def test_pv_provider_is_independent_and_preparation_keeps_the_native_rows(tmp_path):
    import numpy as np
    from types import SimpleNamespace
    from smatrix_bootstrap_newton.kernels import PVSourceRows
    from smatrix_bootstrap_newton.operators import prepare_amplitude
    assert PVSourceRows.__bases__==(object,)
    args=SimpleNamespace(nodes=2,waves=1,bits=256,angular_order=24,subtracted=False,prescription='pv-midpoint',sampling_factor=1,unitarity_energies=[],density_limit=100.,
        moment_source='printed',sr_error='raw-absolute',fesr_cutoff='hard-midpoint',mq_rule='arithmetic-mean',processes=1,output=tmp_path,infinity='free',chiral_tolerance=.002)
    report=prepare_amplitude(args);source=PVSourceRows(2,1,256);H=np.load(tmp_path/'amplitude.npz')['rows']
    expected=np.empty(source.shape)
    for index,_,row in source.iter_rows():expected[index]=[float(v.mid()) for v in row]
    np.testing.assert_array_equal(H,expected)
    assert H.shape==(23,12) and report['prescription']=='pv-midpoint' and report['amplitude_model']['physical_evaluation'].startswith('native nodes')
    with pytest.raises(ValueError,match='native nodes'):source.row(9,0,0)
def test_cardinal_reconstruction_from_continuous_dispersion_and_midpoint_aliases():
    from smatrix_bootstrap_newton.analytic import cardinal_transform,modes,reconstruction_identity,_complex_record
    from flint import arb_mat,acb_mat
    with ctx.workprec(256),mp.workdps(90):
        M=3;z=arb(1)/3;D=cardinal_transform(M);x,w=midpoint_grid(M)
        for n in range(1,M+1):
            density=arb_mat(M,1,[(arb(n*(2*j+1))/(2*M)).sin_pi() for j in range(M)])
            poly=(acb_mat(1,M,modes(z,M))*acb_mat(D*density))[0,0]
            exact=mp.quad(lambda phi:2*mp.mpf(1)/3*mp.sin(phi)*mp.sin(n*phi)/(1-2*mp.cos(phi)/3+mp.mpf(1)/9)/mp.pi,[0,mp.pi])
            assert inside(poly,exact) and (poly-z**n).contains(0)
            record=_complex_record(poly)
            assert arb(record['midpoint'][0],record['radius_upper']).contains(poly.real)
            alias=sum(wi*3/(xi*(xi-3))*density[j,0] for j,(xi,wi) in enumerate(zip(x,w)))
            expected=(2*z**M if n==M else z**n+z**(2*M-n))/(1+z**(2*M))
            assert (alias-expected).contains(0) and alias>poly.real
        assert reconstruction_identity(M,3)['folded_identity_enclosed']
        limit=(arb_mat(1,M,[arb((-1)**n) for n in range(1,M+1)])*D).entries()
        assert all((a+wi/xi).contains(0) for a,xi,wi in zip(limit,x,w))
        # Independently sum the Cauchy quadrature with N distinct from degree M.
        for N in (M+1,2*M):
            xx,ww=midpoint_grid(N);zz=arb(-7)/8;ss=16*zz/(1+zz)**2
            for n in range(1,M+1):
                sampled=[(arb(n*(2*j+1))/(2*N)).sin_pi() for j in range(N)]
                value=sum(wj*ss/(xj*(xj-ss))*v for xj,wj,v in zip(xx,ww,sampled))
                expected=(zz**n+zz**(2*N-n))/(1+zz**(2*N))
                error=(abs(zz)**(2*N-n)+abs(zz)**(2*N+n))/(1-abs(zz)**(2*N))
                assert (value-expected).contains(0) and abs(value-zz**n)<error
def test_rational_circle_projection_against_direct_angle_and_native_PV():
    from smatrix_bootstrap_newton.analytic import CardinalAmplitude
    from smatrix_bootstrap_newton.kernels import PVSourceRows
    with ctx.workprec(192):
        M=2;source=PVSourceRows(M,1,192);coeff=[arb((-1)**j)/(j+2) for j in range(source.coefficient_count)]
        A=CardinalAmplitude(coeff,M);s=source.x[0];angle=arb(1)/(2*M);zs=acb(angle.cos_pi(),angle.sin_pi())
        for I,ell in ((0,0),(2,0),(1,1),(0,2),(1,3)):
            def integrand(mu,analytic):
                t=-(s-4)*(1-mu)/2;u=4-s-t
                rt=(4-t).sqrt(analytic=analytic);ru=(4-u).sqrt(analytic=analytic)
                zt,zu=(2-rt)/(2+rt),(2-ru)/(2+ru)
                a,b,c=A.amplitude(zs,zt,zu),A.amplitude(zt,zs,zu),A.amplitude(zu,zt,zs)
                T=3*a+b+c if I==0 else b+c if I==2 else b-c
                return T*mu.legendre_p(ell)/4
            direct=acb.integral(integrand,-1,1,rel_tol=arb(2)**-65,abs_tol=arb(2)**-85)
            continuous=A.partial_wave(s,0,I,ell)
            assert (continuous-direct).contains(0)
            pv=sum((x*y for x,y in zip(source.row(s,ell,I,node=0),coeff)),acb(0))
            assert (A.partial_wave(s,0,I,ell,True)-pv).contains(0)
            assert abs(continuous-pv)<A.uniform_error_bound(9)[I]
            from smatrix_bootstrap_newton.quadrature import polynomial_rows
            rows,_=polynomial_rows(s,M,[(I,ell)],node=0)
            quadrature=sum((x*y for x,y in zip(rows[I,ell],coeff)),acb(0))
            assert (quadrature-continuous).contains(0)
        from smatrix_bootstrap_newton.analytic import quadrature_repair_bound
        budget=quadrature_repair_bound(A,9)
        assert all(arb(v)<arb('1e-6') for v in budget['S_bounds'].values())
        with pytest.raises(ValueError,match='Quadrature'):
            A.uniform_error_bound(9,M-1)
def test_local_gram_boundary_and_Watson_non_equivalence():
    from smatrix_bootstrap_newton.analytic import local_logic_counterexamples
    from flint import fmpq_mat
    data=local_logic_counterexamples();item=data['PSD_boundary_not_rank_one']
    G=fmpq_mat(item['gram']);B=fmpq_mat(item['factor_B'])
    assert G==B.transpose()*B and G.rank()==2 and G.det()==0
    old=fmpq_mat(data['Watson_functionals_not_identical']['old_gram'])
    assert old.det()==0
    # All principal minors, not just the determinant, certify this old Gram.
    for i,j in ((0,1),(0,2),(1,2)):
        assert old[i,i]*old[j,j]-old[i,j]**2>0
    assert all(old[i,i]>0 for i in range(3))
    h=acb(0,arb(1)/2);F=acb(1);Q=F/F.conjugate()
    replaced=abs(h)*F/abs(F)
    assert replaced.real==arb(1)/2 and Q.imag==0
def test_complete_cardinal_dual_transform_on_independent_monomials(tmp_path):
    from smatrix_bootstrap_newton.basis import mode_row_to_cardinal
    from smatrix_bootstrap_newton.analytic import CardinalAmplitude
    with ctx.workprec(256):
        M=3;zs,zt,zu=map(acb,('0.2','-0.1','0.3'))
        p,q,r=([z**n for n in range(1,M+1)] for z in (zs,zt,zu))
        raw=[acb(1)]+p+[a+b for a,b in zip(q,r)]
        raw += [a*(b+c) for a in p for b,c in zip(q,r)]
        raw += [(q[i]*r[j]+q[j]*r[i])/2 for i in range(M) for j in range(i,M)]
        row=mode_row_to_cardinal(raw,M)
        c=[arb((-1)**j)/(j+2) for j in range(len(row))]
        exact=CardinalAmplitude(c,M).amplitude(zs,zt,zu)
        assert (sum((a*b for a,b in zip(row,c)),acb(0))-exact).contains(0)
        from types import SimpleNamespace
        from smatrix_bootstrap_newton.basis import prepare_cardinal,CardinalSourceRows;from smatrix_bootstrap_newton import write_json
        args=SimpleNamespace(nodes=M,waves=1,bits=256,angular_order=64,subtracted=False,source_registry=None,output=tmp_path,density_limit=1.)
        report=prepare_cardinal(args,lambda **kw:None);write_json(tmp_path/'report.json',report)
        args.preparation=tmp_path;args.output=tmp_path/'extended';args.output.mkdir();args.waves=2
        report=prepare_cardinal(args,lambda **kw:None);write_json(args.output/'report.json',report)
        assert report['amplitude_model']['reused_rows']==29 and report['amplitude_model']['generated_rows']==18
        src=CardinalSourceRows(M,2,256);src.attach_preparation(args.output);src.waves=[(0,0),(2,0),(1,1)]
        value=sum((a*b for a,b in zip(src.row(src.x[0],0,2,node=0),c)),acb(0))
        assert (value-CardinalAmplitude(c,M).partial_wave(src.x[0],0,2,0)).contains(0)
        value=sum((a*b for a,b in zip(src.row(src.x[0],3,1,node=0),c)),acb(0))
        assert (value-CardinalAmplitude(c,M).partial_wave(src.x[0],0,1,3)).contains(0)
        wide=arb(src.x[0],1e-10);assert src.exact_energy(wide) is wide
        with pytest.raises(ValueError,match='native'):src.row(wide,0,2,node=0)
        assert max(v.imag.rad() for v in src.row(wide,0,2))>arb('1e-12')
        from smatrix_bootstrap_newton.quadrature import angular_mode_blocks
        for bad in (arb('inf'),arb(4,.01)):
            with pytest.raises(ValueError,match='Finite'):angular_mode_blocks(bad,M,[0])
@pytest.mark.parametrize('indices',[[0],[0,0]])
def test_analytic_prepare_rejects_missing_or_duplicate_rows(tmp_path,monkeypatch,indices):
    from types import SimpleNamespace
    import smatrix_bootstrap_newton.basis as basis
    fake=SimpleNamespace(M=1,L=1,physical_count=3,shape=(17,5),
        iter_rows=lambda:iter((j,dict(node_zero_based=0),[arb(0)]*5) for j in indices))
    monkeypatch.setattr(basis,'CardinalSourceRows',lambda *a,**kw:fake)
    args=SimpleNamespace(nodes=1,waves=1,bits=128,angular_order=24,
        subtracted=False,source_registry=None,output=tmp_path,density_limit=100)
    with pytest.raises(ValueError,match='[Dd]uplicate|Incomplete'):
        basis.prepare_cardinal(args,lambda **kw:None)
def test_watson_replay_keeps_fixed_goal_current_center_and_direct_inputs(tmp_path):
    import json,numpy as np;from types import SimpleNamespace;from smatrix_bootstrap_newton import write_json,read_json,digest,zero_joint_duals,watson_lineage
    from smatrix_bootstrap_newton.sampling import replay_support_data;from smatrix_bootstrap_newton.run import calculation_inputs
    seed,src,out=[tmp_path/name for name in ('seed','src','out')]
    for path in (seed,src,out):path.mkdir()
    args=SimpleNamespace(preparation=tmp_path/'H',current_preparation=tmp_path/'J',coefficients=src/'center_last_coefficients.json',interior_coefficients=None,output=out,density_limit=1.,chiral_tolerance=.002,chiral_norm='separate-l2')
    sig=dict(M=1,preparation=str(args.preparation),current_preparation=str(args.current_preparation),density_limit=1.,chiral_tolerance=.002,chiral_norm='separate-l2')
    for path in (seed,src,out):
        write_json(path/'coefficients.json',dict(coefficients=[1.,2.]));write_json(path/'center_last_coefficients.json',dict(coefficients=[3.,4.]));write_json(path/'current.json',{});write_json(path/'selection.json',{});write_json(path/'report.json',{});np.savez(path/'joint.npz',point=np.zeros(9))
    write_json(seed/'selection.json',dict(role='tip',model_signature=sig,target=[0.,0.],representative_center_converged=True,coefficient_sha256='old',watson_iteration=4))
    fixed=dict(source_coefficients=str(seed/'center_last_coefficients.json'),identity='old_Q');goal=np.array([2.,3.]);duals=zero_joint_duals(1,3,0);duals['kR'][:]=7
    np.savez(src/'objective.npz',coefficients=goal,metadata=json.dumps(fixed),kR=np.full(3,-99.));np.savez(src/'joint.npz',point=np.zeros(9),**duals);np.savez(src/'barrier_state.npz',mu=1e-5)
    prior=dict(objective_kind='watson',functional=fixed,support_optimality_certified=True,parameters=sig,support_direction=None);retained,normal,latest,cost,meta=replay_support_data(args,prior,1,3,0);assert retained and normal==[0,0] and np.array_equal(cost,goal) and meta==fixed and np.array_equal(latest['kR'],duals['kR'])
    part=replay_support_data(args,dict(prior,support_optimality_certified=False),1,3,0);assert not part[0] and np.array_equal(part[2]['kR'],duals['kR']) and np.array_equal(part[3],goal) and part[4]==fixed
    result=dict(solver=dict(representative_center_converged=False),joint_feasible=True,support_optimality_certified=False);watson_lineage(args,np.array([.1,.2]),fixed,write_json,read_json,result);record=read_json(out/'selection.json');assert record['watson_iteration']==5 and not record['representative_center_converged'] and record['coefficient_sha256']==digest(out/'coefficients.json')
    result['solver']['representative_center_converged']=True;result['support_optimality_certified']=True;watson_lineage(args,np.array([.1,.2]),fixed,write_json,read_json,result);assert read_json(out/'selection.json')['watson_iteration']==5
    snapshot=calculation_inputs(args)
    for path in [args.coefficients,seed/'center_last_coefficients.json',src/'barrier_state.npz',src/'objective.npz',src/'selection.json',src/'current.json',*(seed/name for name in ('coefficients.json','joint.npz','current.json','selection.json','report.json'))]:assert snapshot[str(path.resolve())]['sha256']==digest(path)
    write_json(tmp_path/'evaluation.json',dict(coefficients=str(args.coefficients)));profile_args=SimpleNamespace(**(vars(args)|dict(coefficients=None,current_preparation=None,profile_runs=[tmp_path])))
    assert set(snapshot)<=set(calculation_inputs(profile_args))
    from smatrix_bootstrap_newton import _phase_data;from smatrix_bootstrap_newton.spectra import direct_profiles,plot_profiles;import matplotlib.pyplot as plt
    ev=dict(coefficients=str(out/'coefficients.json'),coefficient_sha256=digest(out/'coefficients.json'),selection=read_json(out/'selection.json'),energies=[4.,5.,6.],waves=[[0,0],[2,0],[1,1]],unitarity=dict(status='sampled_passed'),f=np.zeros((3,3,2)));write_json(tmp_path/'evaluation.json',ev)
    consumers=[_phase_data,lambda p:direct_profiles(SimpleNamespace(profile_runs=[p],output=p)),lambda p:plot_profiles(SimpleNamespace(profile_runs=[p],profile_kind='phase'))]
    for name in ('coefficients.json','current.json','joint.npz'):
        file=out/name;original=file.read_bytes();file.write_bytes(original+b' ')
        for consumer in consumers:
            with pytest.raises(ValueError,match='identity|joint point'):consumer(tmp_path)
        file.write_bytes(original)
    plt.close('all')
