import json;from pathlib import Path
import pytest;flint=pytest.importorskip('flint')
from flint import arb,ctx,fmpq;import numpy as np
from types import SimpleNamespace
from smatrix_bootstrap_newton.imaginary import support_outer, density_fourth_power; from smatrix_bootstrap_newton.operators import midpoint_grid
from smatrix_bootstrap_newton.kernels import density_labels, density_row_to_cflat
from smatrix_bootstrap_newton.model import unitarity_margin, weinberg_waves, current_operators, current_gram; from smatrix_bootstrap_newton.analysis import phase_shifts; from smatrix_bootstrap_newton import UVConfig
from smatrix_bootstrap_newton.linear import barrier_support
def test_exact_half_conversion_affects_only_rho2_offdiagonals():
    row=[fmpq(i+1) for i in range(22)];converted=density_row_to_cflat(row,3); assert all(converted[i]==(row[i]/2 if i in (17,18,20) else row[i]) for i in range(22))
def test_density_norm_counts_actual_upper_density_entries_once():
    coefficients = [fmpq(0)]*5 + list(map(fmpq, [1, 2, 3, 4, 5, 12, 7])); assert density_fourth_power(coefficients, 2) == sum(i**4 for i in range(1, 8))
def _disk_support_control():  # Independent native disks: max Re f00=1 for kappa=1.
    H=np.zeros((17,5)); H[0,0]=2.5;H[2,0]=1.; H[3,1:3]=[1.5,1.];H[4,4]=1.;H[5,2]=1.
    H[14,0]=1.;H[15]=H[0]
    return H
@pytest.mark.parametrize('packed',[False,True])
def test_parallel_extended_products_preserve_cancellation_and_adjoint(packed):
    from concurrent.futures import ThreadPoolExecutor
    from smatrix_bootstrap_newton import extended_product, extended_sparse_blocks
    A=np.array([[2**54,1,-2**54],[2**53,-1,-2**53]],dtype=np.longdouble);A=extended_sparse_blocks(A) if packed else A
    with ThreadPoolExecutor(2) as pool:
        assert np.array_equal(extended_product(A,np.ones(3,np.longdouble),pool=pool),[1,-1])
        assert np.array_equal(extended_product(A,np.array([1,-2],np.longdouble),True,pool),[0,3,0])
def test_zero_native_candidate_keeps_the_verified_nontrivial_incumbent():
    from smatrix_bootstrap_newton.imaginary import recover_candidate_segment
    H=_disk_support_control();initial=np.array([.1/2.5,1/3,.5,.1,.01])
    def certify(c):
        return support_outer(H,[1.],np.zeros(3),np.zeros(3),np.zeros(8),(1.,0.),.25, coefficients=c)
    c,outer,recovery=recover_candidate_segment(H,np.ones(3),1,.25,None,np.zeros(5), initial,None,(1.,0.),certify,infinity_zero=False)
    assert recovery['incumbent_retained'] and recovery['candidate_weight']==0
    assert recovery['interior_fraction']==1-1e-8 and recovery['final_contraction']==1-1e-8
    assert np.array_equal(c,(1-1e-8)*initial) and outer['primal_feasible']
    assert recovery['objective_before']==0 and recovery['objective_after']>.09999999
    assert recovery['objective_loss']<0  # An inferior native return did not erase the incumbent.
def test_full_rows_enforce_the_tighter_disk_with_original_dual():
    H=_disk_support_control();full=np.vstack((H[:3],2*H[0:1],H[3:6],H[5:6],H[6:]))
    seed=np.array([.05,1/3,.5,.1,.01]);solve=barrier_support
    c,kr,ki,y,report=solve(full,np.ones(4),1,1,.25,(1.,0.),None,seed,lambda **kw:None,seconds=10.,raw=True,start_mu=1e-5,gap=1e-5)
    w=report['tail_weights'];outer=support_outer(full,np.ones(4),kr,ki,y,(1.,0.),.25,coefficients=c,waves_per_isospin=1,polynomial_dual=report.get('polynomial_dual'))
    assert outer['primal_feasible'] and .49999<outer['lower']<=.5<=outer['upper']<.50001
def test_weinberg_projection_and_threshold_anchored_phases():
    z,w=np.polynomial.legendre.leggauss(4);F=92/140
    for s in [.5,1.,3.,4.,9.]:
        t=-(s-4)*(1-z)/2;u=4-s-t;A=lambda v:(v-1)/(8*np.pi**2*F**2)
        direct=np.array([w@(3*A(s)+A(t)+A(u)),w@(A(t)+A(u)),(w*z)@(A(t)-A(u))])/4
        assert np.allclose(weinberg_waves(s),direct,rtol=1e-13,atol=1e-15)
    s=np.array([4.,5.,9.,16.,25.]);degrees=np.array([0.,40.,80.,100.,130.]);eta=np.array([1.,.9,.8,.9,.7])
    S=eta*np.exp(2j*np.deg2rad(degrees));f=np.zeros((5,1),complex);f[1:,0]=(S[1:]-1)/(1j*np.pi*np.sqrt(1-4/s[1:]))
    delta,inelas=phase_shifts(f,s)
    assert np.allclose(np.rad2deg(delta[:,0]),degrees) and np.allclose(inelas[:,0],eta)
    undefined=f.copy();undefined[2,0]=complex('nan')
    with pytest.raises(ValueError,match='nonzero S'):phase_shifts(undefined,s)
@pytest.mark.parametrize('zero_imaginary',[False,True])
def test_current_operator_connects_complex_gram_moments_and_high_FF(zero_imaginary):
    M=4;phi=(np.arange(M)+.5)*np.pi/M;nodes=4/np.cos(phi/2)**2
    ss=np.repeat(nodes,3)[::-1];ww=np.tile([[0,0],[1,1],[2,0]],(M,1))[::-1]
    N=len(ss);p=len(density_labels(M));rng=np.random.default_rng(51);kap=np.pi*np.sqrt(1-4/ss)
    H=rng.normal(size=(2*N+11,p))*.01;w=rng.normal(size=p+4*M);cfg=UVConfig()
    if zero_imaginary:w[p:p+2*M]=0
    op=current_operators(H,kap,ss,ww,M,cfg,256);a=op['arrays'];mat=op['matrices']
    tilde=lambda m:0. if m%2==0 else 1/(M*np.tan(m*np.pi/(2*M)))
    K=np.array([[tilde(i+j-2*M-1)-tilde(i-j) for j in range(1,M+1)] for i in range(1,M+1)])
    np.testing.assert_allclose(a['hilbert_kernel'],K,rtol=1e-14,atol=1e-15)
    v=w[p:p+2*M].reshape(2,M);F=1+v@K.T+1j*v;rho=w[p+2*M:].reshape(2,M)
    if zero_imaginary:np.testing.assert_array_equal(F,np.ones((2,M)))
    k=np.array([np.sqrt(6*np.pi)/(16*np.pi**3)*nodes**(-.25)*((nodes-4)/4)**.25, np.sqrt(4*np.pi/3)/(8*np.pi**3)*nodes**(-.25)*((nodes-4)/4)**.75])
    np.testing.assert_allclose(a['k_squared'],k*k,rtol=1e-14,atol=0)
    U=np.array([[1,1j,0],[1,-1j,0],[0,0,np.sqrt(2)]])/np.sqrt(2)
    packed=(a['gram_constant']+mat['gram_linear']@w).reshape(2,M,6)
    for ch in range(2):
        for j,s in enumerate(nodes):
            i=np.flatnonzero((ss==s)&(ww[:,0]==ch)&(ww[:,1]==ch))[0]
            S=1+1j*kap[i]*((H[i]+1j*H[N+i])@w[:p]);f=k[ch,j]*F[ch,j];r=rho[ch,j]
            B=np.array([[1,S,f],[S.conjugate(),1,f.conjugate()],[f.conjugate(),f,r]])
            z=packed[ch,j];G=np.array([[z[0],z[1]/np.sqrt(2),z[3]/np.sqrt(2)], [z[1]/np.sqrt(2),z[2],z[4]/np.sqrt(2)],[z[3]/np.sqrt(2),z[4]/np.sqrt(2),z[5]]])
            np.testing.assert_allclose(G,U.conj().T@B@U,rtol=1e-13,atol=2e-15)
            C=current_gram(S,F[ch,j],r,k[ch,j]**2);C=np.array([[complex(C[a,b]) for b in range(3)] for a in range(3)])
            D=np.diag([1,1,k[ch,j]]);np.testing.assert_allclose(D@C@D,B,rtol=1e-13,atol=2e-15)
    weights=np.pi/M*nodes*np.tan(phi/2)*(nodes<=cfg.s0);powers=np.array([0,1,-1,0])
    expected=[(weights*nodes**n)@rho[ch] for ch,n in [(0,0),(0,1),(1,-1),(1,0)]]
    np.testing.assert_allclose(mat['moment_linear']@w,expected,rtol=1e-13,atol=1e-13)
    targets=np.array([3.09e-8*(27.38/2+.61),3.09e-8*27.38/3,4.34e-6*13.26,4.34e-6*(13.26/2-.41)])
    np.testing.assert_allclose(a['moment_targets_raw'],targets*cfg.s0**(powers+2));np.testing.assert_array_equal(a['moment_errors_raw'],[.002]*4)
    high=np.flatnonzero(nodes>cfg.s0);assert high.tolist()==[3]
    caps=np.array([2*((.004+.0073)/(.140*2))**2*6e-5,3e-5]);slack=(a['ff_cap_constant']+mat['ff_cap_linear']@w).reshape(2,len(high),3)
    expected=np.stack((np.broadcast_to(np.sqrt(caps)[:,None],(2,len(high))),k[:,high]*F[:,high].real,k[:,high]*F[:,high].imag),axis=-1)
    np.testing.assert_allclose(slack,expected,rtol=1e-13,atol=2e-15)
    assert op['metadata']['total_variables']==p+4*M and op['metadata']['gram_block_sizes']==[3]*(2*M)
@pytest.mark.parametrize("extra",[False,True])
def test_joint_hull_imposes_chiral_balls_on_a_common_full_coefficient_mix(tmp_path,extra):
    """Non-QCD toy: ordinary-feasible parents individually violate both chi balls."""
    from smatrix_bootstrap_newton.quotient import joint_hull_candidate
    from smatrix_bootstrap_newton.kernels import joint_audit;from smatrix_bootstrap_newton import write_json
    H=_disk_support_control();H[6,0]=1;H[10,0]=2
    if extra:H=np.vstack((H[:3],H[:3]/2,H[3:6],H[3:6]/2,H[6:]))
    ss=np.repeat([8.,16.] if extra else [8.],3);n=len(ss);kap=np.pi*np.sqrt(1-4/ss)
    parents=np.array([[.02,.1,.2,.1,.25],[-.02,.1,.2,.1,.25]]).T;eps=.001
    current=current_operators(H,kap,ss,np.tile([[0,0],[1,1],[2,0]],(n//3,1)),1)
    a=current['arrays'];seed=np.r_[parents.mean(axis=1),[0,0],3*a['k_squared'].ravel()]
    a['moment_targets_raw']=current['matrices']['moment_linear']@seed;a['moment_errors_raw']=np.full(4,.001)
    points=[]
    for j,c in enumerate(parents.T):
        path=tmp_path/f'parent{j}.json';write_json(path,dict(coefficients=c.tolist(),prescription='pv-midpoint',coordinates='unsubtracted PV-midpoint C_flat'));points.append(dict(coefficients=str(path),feasible=True))
        plain=support_outer(H,kap,np.zeros(n),np.zeros(n),np.zeros(8),[0,0],1.,coefficients=c,waves_per_isospin=1)
        assert plain['primal_feasible'] and np.all(np.linalg.norm((H[2*n:2*n+8]@c).reshape(2,4),axis=1)>eps)
    args=SimpleNamespace(region_summary=tmp_path/'regions.json',preparation=tmp_path/'prepare',chiral_tolerance=eps,native_backend='faer',solver_seconds=5.)
    sig=dict(M=1,L=1,density_limit=1.,prescription='pv-midpoint',infinity='free',unitarity_scope='sampled',preparation=str(args.preparation))
    write_json(args.region_summary,dict(model_signature=sig,regions=[dict(points=points)]))
    candidate,meta=joint_hull_candidate(args,H,current,1,1,1.,lambda **kw:None)
    complete,audit=joint_audit(H,kap,current,candidate,1,1,1.,eps,[0,0]);weights=np.asarray(meta['weights'])
    assert len(complete)==9 and meta['original_variables']==9 and meta['support_search'] is False
    assert np.all(weights>=0) and sum(weights)==pytest.approx(1.) and len(weights)==2
    np.testing.assert_allclose(candidate[:5],parents@weights,rtol=1e-14,atol=1e-16)
    np.testing.assert_array_equal(complete[:5],candidate[:5])
    assert audit['joint_primal_feasible'] and audit['amplitude_primal_feasible']
    assert np.all(np.linalg.norm((H[2*n:2*n+8]@complete[:5]).reshape(2,4),axis=1)<=eps)
@pytest.mark.parametrize('chiral_norm',['separate-l2','combined-l2'])
@pytest.mark.parametrize('chiral_weight',[1.,7.,750.])
@pytest.mark.parametrize('phase_one',[False,True])
def test_joint_barrier_gradient_and_curvature_by_independent_variation(tmp_path,chiral_weight,chiral_norm,phase_one):
    from smatrix_bootstrap_newton.operators import JointProblem;M=4;p=len(density_labels(M));n=3*M;free=1+2*M
    if phase_one:from smatrix_bootstrap_newton.operators import PhaseOneProblem as JointProblem
    with ctx.workprec(256):ss=np.repeat([float(s) for s in midpoint_grid(M)[0]],3)
    kap=np.pi*np.sqrt(1-4/ss);waves=np.tile([[0,0],[1,1],[2,0]],(M,1));H=np.zeros((2*n+11,p));c=np.zeros(p)
    v=np.sin((np.arange(M)+.5)*np.pi/M);c[1:free]=np.tile(v/8,2);c[free:free+M*M]=np.outer(v,v).ravel()
    c[free+M*M:]=[v[i]*v[j]/(2 if i==j else 1) for i in range(M) for j in range(i,M)]
    for j in range(n):H[n+j,free+j]=.7/kap[j]
    rng=np.random.default_rng(51);H[2*n:2*n+8]=rng.normal(size=(8,p))*.0001;H[-2,free]=1
    current=current_operators(H,kap,ss,waves,M);a=current['arrays'];point=np.r_[c,np.zeros(2*M),3*a['k_squared'].ravel()]
    a['moment_targets_raw']=current['matrices']['moment_linear']@point;a['moment_errors_raw']=np.full(4,.01);a['ff_caps_squared']=4*a['k_squared'][:,-1]
    if phase_one and chiral_weight==7:point[-2*M]=-a['k_squared'][0,0];point[-M]=0
    from smatrix_bootstrap_newton.operators import positive_spectral_reference;prior=point.copy();point,replaced=positive_spectral_reference(point,a['k_squared']);assert np.array_equal(point[:p+2*M],prior[:p+2*M]) and np.all(point[-2*M:]>0) and np.array_equal(point[-2*M:][prior[-2*M:]>0],prior[-2*M:][prior[-2*M:]>0])
    problem=JointProblem(H,kap,current,M,1,10.,.002,point,[1,0],chiral_norm=chiral_norm,**({'retain_free_spectra':True} if chiral_weight==750 else {}));problem.chiral_weight=chiral_weight;z=problem.initial;s=problem.state(z,hessian=True)
    if phase_one and chiral_weight==750:
        from smatrix_bootstrap_newton.endpoints import configure_endpoint_coordinates,configure_asymptotic_constraints
        configure_endpoint_coordinates(problem);configure_asymptotic_constraints(problem);z=problem.initial;s=problem.state(z,hessian=True)
    if chiral_weight==750 and chiral_norm=='separate-l2':
        from smatrix_bootstrap_newton.conic import compile_joint_cones,scs_layout,solve_joint_cones
        Q=compile_joint_cones(problem);assert sorted(scs_layout(Q)[0])==list(range(len(Q['b'])));nv=len(problem.active);assert nv==p+4*M;rho=z[problem.free:problem.p]/problem.bound;ext=np.r_[z[:nv],rho*rho];u=((ext-Q['offset'])/Q['scale'])[Q['columns']]
        slack=Q['b']-Q['A']@u;block={k:slack[slice(*v)] for k,v in Q['blocks'].items()};disk=block['scattering'].reshape(-1,4);values=problem.values(z)
        np.testing.assert_allclose(disk[:,0]**2-np.sum(disk[:,1:]**2,axis=1),2*s['amplitude']['slack'],rtol=1e-10,atol=1e-12)
        np.testing.assert_allclose(slack[scs_layout(Q)[0]][-6*len(problem.keep):].reshape(-1,6)[:,[0,1,3,2,4,5]],values['gram'],rtol=1e-11,atol=1e-11);np.testing.assert_allclose(block['moment'],np.r_[1-values['moment'],1+values['moment']],atol=1e-11)
        np.testing.assert_allclose(block['ff'].reshape(-1,3)[:,1:],values['ff'],atol=1e-11)
        for j,g in enumerate(problem.groups):np.testing.assert_allclose(block[f'chiral_{j}'][1:],values['chi'][g]/problem.epsilon,atol=1e-11)
        sq=block['density_squares'].reshape(-1,3);np.testing.assert_allclose(sq[:,0]**2-np.sum(sq[:,1:]**2,axis=1),0,atol=1e-14);assert block['density_norm'][0]**2-np.sum(block['density_norm'][1:]**2)==pytest.approx(float(1-np.sum(rho**4)),abs=1e-14)
        if not phase_one:assert solve_joint_cones(Q,SimpleNamespace(solver='scs',native_backend='auto',solver_seconds=2.,output=tmp_path),u)[-1]['primal_residual']<1e-7
        if 'asymptotic' in block:
            tail=block['asymptotic'].reshape(-1,3);np.testing.assert_allclose(tail[:,0]**2-np.sum(tail[:,1:]**2,axis=1),2*s['asymptotic']['slack'],atol=1e-11)
    d=rng.normal(size=len(z))*np.maximum(abs(z),1e-5);h=np.longdouble('0.00001');d[getattr(problem,'endpoint_pivots',[])]=0
    plus=problem.state(z+h*d);minus=problem.state(z-h*d);g=-s['features'].T@s['rhs']
    shifted=problem.advance(z,s,d,problem.delta(d),h)
    for key in ('slack','chiral_slack'):np.testing.assert_allclose(shifted['amplitude'][key],plus['amplitude'][key],rtol=1e-12,atol=1e-18)
    slope=problem.slope(s,problem.delta(d),d);quadratic=np.sum((s['features']@d)**2)
    assert float(slope)==pytest.approx(float(g@d),rel=1e-10,abs=1e-10)
    assert float((plus['value']-minus['value'])/(2*h))==pytest.approx(float(g@d),rel=2e-6,abs=1e-7)
    assert float((plus['value']+minus['value']-2*s['value'])/h**2)==pytest.approx(float(quadratic),rel=2e-5,abs=1e-5)
    if not phase_one and 'asymptotic' not in s:
        from smatrix_bootstrap_newton.merit import certified_local_decrease
        cost=problem.cost.copy();problem.cost[:]=0;sd=-d if plus['value']>minus['value'] else d;proof=certified_local_decrease(problem,z,s,sd,problem.delta(sd),1.,h);problem.cost[:]=cost
        assert proof['certified_descent'] and 0<float(arb(proof['decrease_lower']))<=float(s['value']-min(plus['value'],minus['value']))+1e-8
@pytest.mark.parametrize("section",[False,True])
def test_stiff_coordinates_and_inexact_descent_do_not_weaken_center_acceptance(section):
    from smatrix_bootstrap_newton import radial_coordinates;from smatrix_bootstrap_newton.linear import joint_direction_status; rng=np.random.default_rng(72);F=rng.normal(size=(25,12)).astype(np.longdouble);F[[3,8]]*=np.array([1e12,1e8])[:,None];c=rng.normal(size=12).astype(np.longdouble)
    full,fullc=F.copy(),c.copy();v=rng.normal(size=11).astype(np.longdouble);F,c=(F[:,1:]+F[:,0,None]*v,c[1:]+c[0]*v) if section else (F,c);free=5-int(section)
    A,q,restore=radial_coordinates(F,c,[3,8],free);w=rng.normal(size=len(c)).astype(np.longdouble);z=restore(w);dz=np.r_[v@z,z] if section else z;assert np.linalg.norm(full@dz-A@w)<1e-16*np.linalg.norm(full@dz) and abs(fullc@dz-q@w)<1e-15 and np.array_equal(z[free:],w[free:])
    H=np.diag([1.,1e-8]);g=np.array([1.,1e-4]);d=np.array([.99999,1e4]);gd=g@d;curvature=d@H@d;step=gd/curvature;status=joint_direction_status(1e-5,gd,gd,curvature);assert status['usable'] and not status['centered'] and .5*step**2*curvature-step*gd<0
    assert not joint_direction_status(1e-5,1e-12,1e-12,1e-12)['centered'] and not joint_direction_status(1e-5,1.,-1.,1.)['usable'] and joint_direction_status(0.,0.,0.,0.)['centered'] and not joint_direction_status(1e-7,3.34e-11,3.34e-11,3.34e-11)['centered'] and joint_direction_status(1e-7,1e-17,1e-17,1e-17)['centered']
    from smatrix_bootstrap_newton.linear import self_concordant_decrease; e=np.longdouble('1e-7'); lower=self_concordant_decrease(1,e*e,e*e)
    with ctx.workprec(192): truth=(1+arb('1e-7')).log()-(1-arb('1e-7'))*arb('1e-7'); assert truth>float(lower)>0 and np.longdouble(1e6)-np.longdouble(np.longdouble(1e6)-np.longdouble(float(truth)))==0
    from smatrix_bootstrap_newton.merit import certified_local_decrease
    xx=np.longdouble('.25');dd=np.longdouble('1e-9');cc=np.longdouble(8)/15+np.longdouble('1e-8');P=SimpleNamespace(free=1,p=1,bound=1.,groups=(),chiral_weight=1.,t=np.ones(1),cost=np.array([cc]))
    vals=dict(x=np.array([xx]),y=np.ones(1),chi=np.empty(0),gram=np.empty((0,6)),ff=np.empty((0,2)),moment=np.empty(0))
    state=dict(values=vals,amplitude=dict(slack=np.array([1-xx*xx]),chiral_slack=np.empty(0)),ffslack=np.empty(0))
    dv=dict(x=np.array([dd]),y=np.zeros(1),chi=np.empty(0),gram=np.empty((0,6)),ff=np.empty((0,2)),moment=np.empty(0));proof=certified_local_decrease(P,np.array([xx]),state,np.array([dd]),dv,1.,1.)
    with ctx.workprec(256):
        A=lambda v:arb(v.as_integer_ratio()[0])/arb(v.as_integer_ratio()[1]);X,D,C=map(A,(xx,dd,cc));exact=C*D+((-2*X*D-D*D)/(1-X*X)).log1p();assert proof['certified_descent'] and 0<arb(proof['decrease_lower'])<exact
    dv['x']*=-1;assert not certified_local_decrease(P,np.array([xx]),state,np.array([-dd]),dv,1.,1.)['certified_descent']
def test_watsonian_functional_matches_physical_complex_S_identity():
    from smatrix_bootstrap_newton.imaginary import watsonian_objective; rng=np.random.default_rng(54);M=4;p=len(density_labels(M));n=3*M;H=rng.normal(size=(2*n+11,p))*.01
    nodes=np.array([float(x) for x in midpoint_grid(M)[0]]);ss=np.repeat(nodes,3);kap=np.pi*np.sqrt(1-4/ss);ww=np.tile([[0,0],[1,1],[2,0]],(M,1))
    current=current_operators(H,kap,ss,ww,M);old=rng.normal(size=p+4*M)*.1;new=rng.normal(size=p+4*M)*.1
    c,duals,meta=watsonian_objective(H,kap,current,old,M,1);im=old[p:p+2*M].reshape(2,M);F=1+im@current['arrays']['hilbert_kernel'].T+1j*im
    Sold=1+1j*kap*((H[:n]+1j*H[n:2*n])@old[:p]);Snew=1+1j*kap*((H[:n]+1j*H[n:2*n])@new[:p]);truth=0.;upper=0.;count=0
    for i in np.flatnonzero(nodes<=(1.2/.14)**2):
        for ch in range(3):
            Q=F[ch,i]/F[ch,i].conjugate() if ch<2 else Sold[3*i+ch];w=(np.sqrt(nodes[i])+2)/(np.sqrt(nodes[i])-2) if ch==1 else 1
            truth+=w*np.real(Q.conjugate()*(Snew[3*i+ch]-1));upper+=w*(abs(Q)-Q.real);count+=1
    assert c@new[:p]==pytest.approx(truth/count,abs=1e-12) and meta['ideal_disk_ceiling']==pytest.approx(upper/count,abs=1e-12)
@pytest.mark.parametrize('lo,hi,verdict',[(-.875,.5,True),(-.5,.875,False),(-.75,.75,None)])
def test_certified_reference_section_uses_feasible_segments_and_true_supports(lo,hi,verdict):
    from smatrix_bootstrap_newton.kernels import certify_reference_section;from smatrix_bootstrap_newton.run import encode_real_ball
    def box(low,high,label):
        vertices=[(-1,low),(1,low),(1,high),(-1,high)];directions=[(-1,0),(0,-1),(1,0),(0,1)]
        return [dict(target_enclosures=[encode_real_ball(arb(v)) for v in p],direction=d,upper=float(max(np.dot(d,q) for q in vertices)),coefficients=label+str(i),report=label+str(i)) for i,(p,d) in enumerate(zip(vertices,directions))]
    result=certify_reference_section(box(-1,1,'IR'),box(lo,hi,'UV'),0.)
    for name,value in [('upper_shrink',1-hi),('lower_rise',lo+1),('upper_minus_lower',-hi-lo),('upper_to_lower_ratio',(1-hi)/(lo+1))]:
        a,b=result[name];assert a<=value<=b and b-a<1e-14
    assert result['upper_change_larger'] is verdict
def test_ff_guided_lift_preserves_S_and_records_native_aliasing():
    from smatrix_bootstrap_newton.analysis import ff_guided_phase
    d=np.array([.1,1.9,2.8]);S=np.array([.9,.8,.7])*np.exp(2j*d);F=np.exp(1j*d)
    for sign in (1.,-1.):
        lift=ff_guided_phase(S,sign*F,sign);np.testing.assert_allclose(lift,d,atol=1e-14);np.testing.assert_allclose(abs(S)*np.exp(2j*lift),S,atol=1e-14)
def test_original_zero_objective_farkas_bound_from_nonnegative_current_spectrum():
    """PSD rho0>=0 contradicts a negative upper bound on its positive weighted moment."""
    from smatrix_bootstrap_newton.kernels import joint_audit;H=_disk_support_control();kap=np.full(3,np.pi/np.sqrt(2))
    current=current_operators(H,kap,np.full(3,8.),np.array([[0,0],[1,1],[2,0]]),1);a=current['arrays']
    a['moment_targets_raw'][0]=-1.;a['moment_errors_raw'][0]=.25
    weight=float(current['matrices']['moment_linear'][0,7]);assert weight>0
    Z=np.zeros((2,1,3,3));Z[0,0,2,2]=weight
    point=np.r_[[0.,.1,.2,.1,.25],np.zeros(4)];duals=dict(moment=[1,0,0,0],gram=Z,gram_congruence=np.ones((2,1,3)))
    _,audit=joint_audit(H,kap,current,point,1,1,1.,None,[0,0],duals,objective=None)
    assert not audit['primal_feasible'] and audit['lower'] is None
    assert -.75<=audit['upper']<0 and arb(audit['enclosure'])<0
def test_pv_candidate_recovery_keeps_one_exact_full_coefficient_segment():
    from smatrix_bootstrap_newton.imaginary import recover_candidate_segment
    H=_disk_support_control();initial=np.array([.04,1/3,.5,.1,.01]);candidate=np.array([.8,1/3,.5,.1,.01])
    def certify(c):return support_outer(H,[1.],np.zeros(3),np.zeros(3),np.zeros(8),(1.,0.),1.,coefficients=c)
    c,outer,meta=recover_candidate_segment(H,np.ones(3),1,1.,None,candidate,initial,None,(1.,0.),certify)
    assert outer['primal_feasible'] and abs(meta['candidate_fraction']-9/19)<2e-8
    np.testing.assert_allclose(c,meta['candidate_weight']*candidate+meta['interior_fraction']*initial,atol=2e-16)
    assert .9999999<outer['lower']<=1<=outer['upper']<1.00000001
def test_full_C_json_keeps_data_identity_and_rejects_old_CG(tmp_path):
    from smatrix_bootstrap_newton import read_cflat, write_json
    c=np.arange(22.);path=tmp_path/'coefficients.json'
    write_json(path,dict(coordinates='unsubtracted PV-midpoint C_flat',prescription='pv-midpoint',coefficients=c))
    result,record=read_cflat(path,3);np.testing.assert_array_equal(result,c);assert record['prescription']=='pv-midpoint'
    with pytest.raises(ValueError,match='complete C_flat'):read_cflat(path,4)
    write_json(path,dict(coordinates='unsubtracted sine-cardinal C_flat',coefficients=c))
    result,record=read_cflat(path,3);assert record['prescription']=='finite-sine-cardinal' and np.array_equal(result,c)
    with pytest.raises(ValueError,match='CG'):read_cflat(tmp_path/'old.npz',3,'raw_candidate')
def test_joint_report_helpers_keep_complete_dimensions_and_phase_one_scope():
    from smatrix_bootstrap_newton import center_report,phase_one_report
    P=SimpleNamespace(p=3876,M=50,chiral_weight=750)
    center=center_report(P,np.array([2.,500.]),.07,np.longdouble(0),False,[],3,2,7,1.5)
    phase=phase_one_report(P,np.r_[np.zeros(4062),.1],4,2.,.001,False)
    assert center['original_variables']==phase['original_variables']==4076
    assert center['support_direction']==[2.,500.] and center['fixed_x']==.07 and center['chiral_barrier_weight']==750
    assert phase['solver_variables']==4063 and phase['phase_I'] and phase['tau']==.1 and phase['chiral_barrier_weight']==750
    assert phase['scattering_chiral_density_always_hard'] and not phase['physical_constraints_relaxed_at_acceptance']
    from smatrix_bootstrap_newton.certificates import model_signature
    r=dict(M=1,L=1,density_limit=1.,infinity='zero',unitarity_scope='sampled',prescription='analytic-cardinal',chiral_tolerance=.002,scattering_samples=3,parameters=dict(preparation='H',current_preparation='K'),current_model=dict(uv=dict(config={}),ff_endpoint_order=2,asymptotic_unitarity=True));sig=model_signature(r)
    assert sig['tail_conditions_applied'] and sig['ff_endpoint_order']==2
def test_spectral_floor_from_independent_complex_schur_complement():
    from smatrix_bootstrap_newton.model import spectral_budget
    S=np.array([.3+.4j,-.8+.1j]);F=np.array([1.+2j,-.4+.7j]);k2=np.array([.01,.03]);floor=[]
    for s,f,k in zip(S,F,k2):
        A=np.array([[1,s],[s.conjugate(),1]]);v=np.array([f,f.conjugate()])
        floor.append(k*np.vdot(v,np.linalg.solve(A,v)).real)
    rho=np.array(floor)+.05;b=spectral_budget(S,F,rho,k2)
    np.testing.assert_allclose(b['rho_minimum'],floor,rtol=1e-14)
    np.testing.assert_allclose(b['spectral_excess'],.05,rtol=1e-13)
    for j,(s,f,r,k) in enumerate(zip(S,F,rho,k2)):
        G=np.array([[1,s,np.sqrt(k)*f],[s.conjugate(),1,np.sqrt(k)*f.conjugate()], [np.sqrt(k)*f.conjugate(),np.sqrt(k)*f,r]])
        assert np.linalg.eigvalsh(G).min()>0
        assert b['determinant_slack'][j]==pytest.approx(np.linalg.det(G).real,rel=1e-13)
        G[2,2]=floor[j]-.01;assert np.linalg.eigvalsh(G).min()<0
def test_spectral_diagnostic_does_not_divide_through_elastic_zero_or_violation():
    from smatrix_bootstrap_newton.model import spectral_budget
    b=spectral_budget(np.array([1.,1.+1e-13,0.]),np.array([1.,1.,0.]),np.ones(3),np.ones(3))
    assert b['phase_penalty']==[None,None,0.] and b['rho_minimum']==[None,None,0.]
    assert b['determinant_slack'][1]<0 and b['spectral_excess'][2]==1.
def test_native_peak_reporting_distinguishes_interior_peak_and_endpoint_maximum():
    from smatrix_bootstrap_newton.model import native_peaks;from smatrix_bootstrap_newton.certificates import native_peak_bounds,native_rho_certificate;x=np.array([.3,.5,.7,.9,1.1]);r=native_peaks(x,[0.,.8,.2,.9,.3])
    assert [v['index'] for v in r['local']]==[1,3] and not r['maximum']['endpoint'];assert len(native_peak_bounds([1,2,3],[arb(0),arb(1,.1),arb(0)])['strict_peaks'])==1 and native_peak_bounds([1,2,3],[arb(0),arb(0,1),arb(0)])['strict_peaks']==[]
    r=native_peaks(x,[0.,.1,.2,.3,.4]);assert r['local']==[] and r['maximum']['endpoint'] and not r['pole_determined'];nodes=midpoint_grid(2)[0];d=dict(amplitude_model=dict(prescription='pv-midpoint'),selection=dict(model_signature=dict(M=2)),energies=[4,float(nodes[0])*(1+1e-14),float(nodes[1])],waves=[[0,0],[2,0],[1,1]],f_enclosures=[[[None,None]]*3]*3)
    with pytest.raises(ValueError,match='native physical nodes'):native_rho_certificate(d)
def test_current_residual_bound_dominates_independent_exact_ellipsoid_support():
    from smatrix_bootstrap_newton.model import affine_norm_support
    from flint import arb_mat
    K=np.array([[.3,-.2],[.5,.1]]);d=np.array([.2,.7]);r=np.array([.4,-.8]);D=np.diag(d)
    C=K.T@D@K+D;g=K.T@d;center=-np.linalg.solve(C,g)
    radius2=2-d.sum()+g@np.linalg.solve(C,g)
    exact=r@center+np.sqrt(radius2*(r@np.linalg.solve(C,r)))
    with ctx.workprec(256):
        bound=affine_norm_support(arb_mat(K.tolist()),list(map(arb,d)),list(map(arb,r)))
        assert bound>arb(float(exact)) and float(bound.upper())<3.
@pytest.mark.parametrize("extra",[False,True])
@pytest.mark.parametrize('excluded',[False,True])
def test_fixed_current_fiber_certificate_excludes_only_the_frozen_amplitude(tmp_path,excluded,extra):
    """Independent toy: rho=3k² permits S0>0, but cannot accommodate S0=-.9."""
    from smatrix_bootstrap_newton.quotient import joint_hull_candidate;from smatrix_bootstrap_newton.kernels import joint_audit;from smatrix_bootstrap_newton import write_json
    H=_disk_support_control();c=np.array([0.,.1,.2,.1,.25])
    if extra:H=np.vstack((H[:3],H[:3]/2,H[3:6],H[3:6]/2,H[6:]))
    ss=np.repeat([8.,16.] if extra else [8.],3);n=len(ss);kap=np.pi*np.sqrt(1-4/ss)
    current=current_operators(H,kap,ss,np.tile([[0,0],[1,1],[2,0]],(n//3,1)),1);a=current['arrays']
    witness=np.r_[c,[0.,0.],3*a['k_squared'].ravel()]
    a['moment_targets_raw']=current['matrices']['moment_linear']@witness;a['moment_errors_raw']=np.full(4,1e-9)
    _,full=joint_audit(H,kap,current,witness,1,1,1.,.002,[0,0]);assert full['primal_feasible']
    if excluded:c[1]=(1.9/kap[0]-c[2])/1.5
    path=tmp_path/'coefficients.json';write_json(path,dict(coefficients=c,prescription='pv-midpoint',coordinates='unsubtracted PV-midpoint C_flat'))
    args=SimpleNamespace(region_summary=tmp_path/'regions.json',preparation=tmp_path/'H',coefficients=path,output=tmp_path, chiral_tolerance=.002,chiral_norm='combined-l2',fixed_amplitude=True,direction=[0,0],native_backend='faer',solver_seconds=5.)
    sig=dict(M=1,L=1,density_limit=1.,prescription='pv-midpoint',infinity='free',unitarity_scope='sampled',preparation=str(args.preparation))
    write_json(args.region_summary,dict(model_signature=sig,regions=[]))
    point,_=joint_hull_candidate(args,H,current,1,1,1.,lambda **kw:None)
    np.testing.assert_array_equal(point[:5],c)
    with np.load(tmp_path/'fiber_duals.npz') as saved:duals={k:saved[k] for k in saved.files}
    assert duals['kR'].shape==duals['kI'].shape==(n,)
    _,audit=joint_audit(H,kap,current,point,1,1,1.,.002,[0,0],duals,fixed_amplitude=True)
    assert audit['amplitude_primal_feasible'] and audit['fixed_amplitude'] and (audit['upper']<0)==excluded
    if not excluded:assert audit['primal_feasible']
def test_gram_dual_repair_respects_unequal_diagonal_scales():
    from smatrix_bootstrap_newton.kernels import joint_audit;H=_disk_support_control();kap=np.full(3,np.pi/np.sqrt(2))
    current=current_operators(H,kap,np.full(3,8.),np.array([[0,0],[1,1],[2,0]]),1)
    point=np.r_[[0.,.1,.2,.1,.25],np.zeros(2),3*current['arrays']['k_squared'].ravel()]
    v=np.array([1e8,1.,1e-4]);B=np.outer(v,v)-1e-8*np.diag(v*v);Z=np.array([[B],[B]])
    _,audit=joint_audit(H,kap,current,point,1,1,1.,.002,[0,0],dict(gram=Z))
    repaired=np.array(audit['repaired_normalized_gram_dual'])[0];change=repaired-B
    assert change[0,0]>0 and 0<change[2,2]<change[0,0]*1e-20
    scaled=repaired/v[:,None]/v
    assert np.linalg.eigvalsh(scaled).min()>-1e-14
def test_dual_rejects_a_layout_with_a_truly_unbounded_free_objective():
    H=_disk_support_control();H[3,1]=0;H[-2]=0;H[-2,1]=1
    c=np.array([0.,1e30,0.,0.,0.])
    np.testing.assert_array_equal(H[:6]@c,np.zeros(6))
    with pytest.raises(ValueError,match='native imaginary'):
        support_outer(H,[1.],np.zeros(3),np.zeros(3),np.zeros(8),[1,0],1.,coefficients=c)
def test_paper_marker_outputs_cannot_drive_new_scientific_selection():
    from smatrix_bootstrap_newton.analysis import select_gauge
    with pytest.raises(ValueError,match='comparison outputs'):
        select_gauge(SimpleNamespace(selection_reference=Path('unread_paper_markers.csv')))
@pytest.mark.parametrize('epsilon,expected',[(None,(1+np.sqrt(.75))/1.5),(1.,1.)])
def test_fixed_section_extrema(epsilon,expected):
    from smatrix_bootstrap_newton.scattering import barrier_support
    H=np.zeros((17,5));H[0,0]=H[1,4]=H[2,4]=1;H[3,1:3]=[1.5,1];H[4,3]=H[5,2]=1
    H[-2,0]=H[-1,1]=1;H[6,1]=1
    z,kr,ki,y,r=barrier_support(H,np.ones(3),1,1,1.,[0,1],epsilon,np.array([.5,.1,.5,.1,.1]),lambda **kw:None,raw=True,fixed_x=.5,require_center=True,chiral_weight=750,seconds=3,gap=1e-6)
    assert r['representative_center_converged'] and abs(z[0]-.5)<1e-12
    assert z[1]==pytest.approx(expected,abs=2e-6)
    a=support_outer(H,np.ones(3),kr,ki,y,r['support_direction'],1.,epsilon,coefficients=z,waves_per_isospin=1)
    assert a['primal_feasible'] and a['upper']-a['lower']<2e-6
def test_joint_resume_uses_checkpoint_mu_only_for_same_problem(tmp_path):
    from smatrix_bootstrap_newton import resolve_start_mu,write_json
    base=dict(command='boundary',mode='gauge',objective='projection',start_mu=None, coefficients=tmp_path/'coefficients.json',interior_coefficients=None,joint_feasibility=False,resume_objective=False,
        chiral_norm='separate-l2',chiral_tolerance=.002,infinity='free',unitarity_scope='sampled',fixed_x=None,ray=False, preparation=tmp_path/'H',current_preparation=tmp_path/'J',direction=[1.,0.],density_limit=None,chiral_barrier_weight=None)
    old={k:str(v) if isinstance(v,Path) else v for k,v in base.items()}
    old['start_mu']=.01
    write_json(tmp_path/'report.json',dict(M=3,density_limit=1500,objective_kind='projection',parameters=old,solver={}))
    np.savez(tmp_path/'barrier_state.npz',mu=2.5e-7)
    args=SimpleNamespace(**base);resolve_start_mu(args)
    assert args.start_mu==2.5e-7 and args.start_mu_provenance['kind']=='matching_joint_checkpoint'
    for change in (dict(direction=[0.,1.]),dict(preparation=tmp_path/'other_H'),dict(chiral_norm='combined-l2')):
        args=SimpleNamespace(**(base|change));resolve_start_mu(args);assert args.start_mu==.01
    args=SimpleNamespace(**(base|dict(start_mu=.001)));resolve_start_mu(args);assert args.start_mu==.001
def test_joint_audit_rejects_incomplete_energy_metadata_before_using_current_rows():
    from smatrix_bootstrap_newton.certificates import joint_audit
    current={'arrays':{},'metadata':{}}
    with pytest.raises(ValueError,match='scattering metadata'):
        joint_audit(np.zeros((19,5)),np.ones(3),current,np.zeros(9),1,1,1.,.002,[0,0])
