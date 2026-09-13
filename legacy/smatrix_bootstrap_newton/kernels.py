"""PV amplitude kernels and original-operator audits."""
from flint import acb, arb, arb_mat, ctx
from fractions import Fraction
from math import factorial
from pathlib import Path
from scipy.spatial import ConvexHull, QhullError
import numpy as np
from . import Model, decode_real_ball, fraction_ball as _fraction

class PVSourceRows:
    """Native PV collocation and rational off-cut rows; no exact finite analytic family."""
    def __init__(self,M=50,L=10,bits=384,order=24,subtracted=False):
        from .operators import midpoint_grid
        Model(M,L,bits);ctx.prec=bits
        self.M,self.L,self.working_bits,self.order,self.subtracted=M,L,bits,order,subtracted
        self.x,self.w=midpoint_grid(M);self.b=[w/x for x,w in zip(self.x,self.w)];self.K=pv_matrix(M)
        self.coefficient_labels,self.raw_density_labels=cflat_labels(M),density_labels(M)
        self.coefficient_count=len(self.coefficient_labels);self.waves=canonical_waves(L)
        self.physical_count=3*M*L;self.shape=(2*self.physical_count+11,self.coefficient_count)
        self.coordinate_name=('subtracted' if subtracted else 'unsubtracted')+' PV-midpoint C_flat'
        self.metadata=dict(prescription='pv-midpoint',M=M,L=L,bits=bits,coordinate_system=self.coordinate_name,
            coefficient_count=self.coefficient_count,physical_samples=self.physical_count,nu0=0,nyquist_retained=None,
            infinity_equality_imposed=False,density_limit_imposed=False,quadrature_order_per_panel=None,
            arithmetic='Arb',integral_enclosures=False,authors_implementation_recovered=False,analytic_function_family_defined=False,
            angular_projection='analytic exterior Legendre Q',spectral_quadrature='midpoint in phi; native PV cot map and nodal jump',
            physical_evaluation='native nodes only; no off-grid continuation specified',
            analyticity_status='Mixed collocation approximation; native jumps and rational off-cut rows are not one exact analytic family',
            residue_identity='Res[f_offcut,s_i] = -w_i Im[f_PV(s_i)]',
            infinity_row_meaning='unsubtracted constant coordinate; not a certified high-energy limit')

    def row(self,s,ell,isospin,node=None,spins=None):
        from .operators import angular_kernels, assemble_density_row
        ctx.prec=self.working_bits;s=_fraction(s) if isinstance(s,Fraction) else arb(s)
        if not s.is_finite() or not s>0:raise ValueError('Finite positive s required')
        if node is None and s>4:
            node=next((j for j,x in enumerate(self.x) if s==x or s.is_exact() and float(s)==float(x)),None)
        if node is not None:
            if type(node) is not int or not 0<=node<self.M or float(s)!=float(self.x[node]):raise ValueError('Matching native PV node required')
            row=physical_node_row(node,self.x,self.w,self.K,ell,isospin,subtracted=self.subtracted)
        elif 0<s<4:row=offcut_row(s,self.x,self.w,ell,isospin,subtracted=self.subtracted)
        elif s==4:
            blocks=angular_kernels(s,self.x,self.w,ell,subtracted=self.subtracted)
            direct=[w/(x-s)-(b if self.subtracted else 0) for x,w,b in zip(self.x,self.w,self.b)]
            row=assemble_density_row(direct,blocks,isospin)
        else:raise ValueError('PV physical evaluation requires native nodes; off-grid is unspecified')
        return density_row_to_cflat(row,self.M)

    def evaluate(self,coefficients,energies,waves=None):
        ctx.prec=self.working_bits;values=list(map(arb,coefficients));waves=self.waves if waves is None else list(waves)
        if len(values)!=self.coefficient_count or any(not v.is_finite() for v in values):raise ValueError('Finite complete PV C_flat coefficients required')
        return [[sum((a*b for a,b in zip(self.row(s,ell,I),values)),acb(0)) for I,ell in waves] for s in energies]

    def iter_rows(self,canonical_ids=None):
        requested=set(range(self.shape[0])) if canonical_ids is None else set(canonical_ids)
        if any(type(i) is not int or not 0<=i<self.shape[0] for i in requested):raise ValueError('Canonical source row outside retained range')
        yielded=set();offset=self.physical_count
        coordinates=self.coordinate_name
        for node in range(self.M):
            for wave_index,(I,ell) in enumerate(self.waves):
                real_index=node*len(self.waves)+wave_index;imag_index=offset+real_index
                if real_index not in requested and imag_index not in requested:continue
                row=self.physical_row(node,ell,I)
                for index,part in ((real_index,'real'),(imag_index,'imaginary')):
                    if index not in requested:continue
                    values=[value.real if part=='real' else value.imag for value in row]
                    yielded.add(index)
                    yield index,dict(kind='scattering_'+part,node_zero_based=node,isospin=I,ell=ell,
                        normalization='f',amplitude_coordinates=coordinates),values
        base=2*offset;energies=[Fraction(1,2),Fraction(1),Fraction(3,2),Fraction(2)]
        for channel in ('01','21'):
            for position,s in enumerate(energies):
                index=base+position+(0 if channel=='01' else 4)
                if index not in requested:continue
                ratio=3*(2*s-1)/(s-4) if channel=='01' else 3*(2-s)/(s-4)
                first=self.subthreshold_row(s,0,0 if channel=='01' else 2);pwave=self.subthreshold_row(s,1,1)
                row=[a-_fraction(ratio)*b for a,b in zip(first,pwave)];yielded.add(index)
                yield index,dict(kind='chiral_'+channel,exact_energy=str(s),P1_ratio=str(ratio),amplitude_coordinates=coordinates),row
        if base+8 in requested:
            yielded.add(base+8);yield base+8,dict(kind='unsubtracted_constant_coordinate',amplitude_coordinates=coordinates),self.infinity_row()
        for extra,I,ell,label in ((9,0,0,'target_f00_at_3'),(10,1,1,'target_f11_at_3')):
            if base+extra in requested:
                yielded.add(base+extra);yield base+extra,dict(kind=label,exact_energy='3',amplitude_coordinates=coordinates),self.subthreshold_row(Fraction(3),ell,I)
        if yielded!=requested:raise ValueError('Canonical PV row coverage incomplete')

    def physical_row(self, node, ell, isospin):
        if type(node) is not int or not 0 <= node < self.M:
            raise ValueError('Node outside the physical grid')
        return self.row(self.x[node], ell, isospin, node=node)

    def subthreshold_row(self, s, ell, isospin):
        return self.row(s, ell, isospin)

    def infinity_row(self):
        if not self.subtracted:return [arb(1)]+[arb(0)]*(self.coefficient_count-1)
        return ([arb(1)]+[-v for v in self.b]+[-2*v for v in self.b]
            +[2*a*b for a in self.b for b in self.b]
            +[self.b[i]*self.b[j] for i in range(self.M) for j in range(i,self.M)])

def density_labels(M):
    if type(M) is not int or M<1:raise ValueError('A positive integer node count is required')
    return ([('T0',0,0)]+[('sigma1',i,0) for i in range(M)]+[('sigma2',i,0) for i in range(M)]
        +[('rho1',i,j) for i in range(M) for j in range(M)]
        +[('rho2',i,j) for i in range(M) for j in range(i,M)])

def pv_matrix(M):
    density_labels(M);values=[]
    for m in range(2*M):
        if m%2==0:values.append(arb(0))
        else:
            angle=arb(m)/(2*M)
            values.append(angle.cos_pi()/(M*angle.sin_pi()))
    return arb_mat(M,M,[-values[(k-i)%(2*M)]+values[(k+i+1)%(2*M)]
        for k in range(M) for i in range(M)])

def exterior_legendre_q(ell,z):
    if type(ell) is not int or ell<0:raise ValueError('Nonnegative integer spin required')
    z=arb(z)
    if z>1:sign=1
    elif z<-1:sign=(-1)**(ell+1);z=-z
    else:raise ValueError('A separated real exterior Legendre-Q argument is required')
    coefficient=arb(2**ell*factorial(ell)**2)/factorial(2*ell+1)
    argument=1/(z*z)
    hyper=argument.hypgeom_2f1(arb(ell+1)/2,arb(ell+2)/2,arb(2*ell+3)/2)
    answer=sign*coefficient*z**(-ell-1)*hyper
    if not answer.is_finite():raise ArithmeticError('Exterior Legendre-Q enclosure is unresolved')
    return answer

def offcut_row(s,nodes,weights,ell,isospin,*,subtracted=True):
    from .operators import angular_kernels, assemble_density_row
    s=arb(s)
    if not 0<s<4:raise ValueError('This off-cut API is restricted to0<s<4')
    blocks=angular_kernels(s,nodes,weights,ell,subtracted=subtracted)
    direct=[w*(1/(x-s)-(1/x if subtracted else 0)) for x,w in zip(blocks['nodes'],blocks['weights'])]
    return assemble_density_row(direct,blocks,isospin)

def physical_node_row(node,nodes,weights,pv,ell,isospin,*,subtracted=True):
    from .operators import angular_kernels, assemble_density_row
    M=len(nodes)
    if type(node) is not int or not 0<=node<M or pv.nrows()!=M or pv.ncols()!=M:
        raise ValueError('Matching node and full M-square PV matrix required')
    blocks=angular_kernels(nodes[node],nodes,weights,ell,subtracted=subtracted)
    direct=[acb(pv[node,i]+(0 if subtracted else blocks['b'][i]),int(node==i)) for i in range(M)]
    return assemble_density_row(direct,blocks,isospin)

def cflat_labels(M):
    return [('C_flat' if family=='rho2' else family,i,j) for family,i,j in density_labels(M)]

def density_row_to_cflat(row,M):
    labels=density_labels(M)
    if len(row)!=len(labels):raise ValueError('Full density row required')
    return [value/2 if family=='rho2' and i!=j else value for value,(family,i,j) in zip(row,labels)]

def coefficient_blocks(coefficients,M):
    if len(coefficients)!=len(density_labels(M)):raise ValueError('Complete C_flat coefficients required')
    values=list(map(arb,coefficients)); free=1+2*M
    R1=arb_mat(M,M,values[free:free+M*M]); R2=arb_mat(M,M); at=free+M*M
    for i in range(M):
        for j in range(i,M):
            R2[i,j]=R2[j,i]=values[at]/(1 if i==j else 2); at+=1
    return values,R1,R2

def canonical_waves(L):
    return [(I,2*k+int(I==1)) for I in range(3) for k in range(L)]

def source_rows(M=50,L=10,bits=384,order=24,subtracted=False,*,prescription='pv-midpoint',registry=None):
    if prescription=='analytic-cardinal':
        from .basis import CardinalSourceRows
        return CardinalSourceRows(M,L,bits,order,subtracted,registry)
    if prescription!='pv-midpoint' or registry is not None:raise ValueError('Unknown provider or registry incompatible with PV')
    return PVSourceRows(M,L,bits,order,subtracted)

def region_geometry(region):
    from .analysis import published_inner_vertices
    from functools import cmp_to_key
    from .run import encode_real_ball
    with ctx.workprec(384):
        up=lambda x:float(np.nextafter(float(x.upper()),np.inf))
        dot=lambda a,b:sum((x*y for x,y in zip(a,b)),arb(0))
        norm=lambda a:dot(a,a).sqrt()
        cross=lambda a,b:a[0]*b[1]-a[1]*b[0]
        epsilon=region.get('epsilon');a,b,c=(arb(1),arb(0),arb(1)) if epsilon is None else (16*arb.pi()**2*(arb(92)/140)**2/5,1/(15*arb(epsilon)),1/arb(epsilon))
        metric=lambda p:[a*p[0],b*p[0]+c*p[1]]
        dual=lambda d:[(d[0]-b*d[1]/c)/a,d[1]/c]
        points=region['points'];D=[list(map(arb,p['direction'])) for p in points];U=[arb(p['upper']) for p in points]
        P=[[decode_real_ball(v) for v in p['target_enclosures']] for p in points];mixtures=None
        if region.get('published_positive_x',False):
            P,mixtures=published_inner_vertices(P);D.append([arb(-1),arb(0)]);U.append(arb(0))
        Q=list(map(metric,P))
        mids=np.array([[float(v.mid()) for v in q] for q in Q]);directions=np.array([[float(v) for v in d] for d in D])
        angles=np.mod(np.arctan2(directions[:,1],directions[:,0]),2*np.pi)
        half=lambda d:int(d[1]<0 or d[1].is_zero() and d[0]<0)
        def compare(i,j):
            h=half(D[i])-half(D[j]);z=cross(D[i],D[j])
            return h if h else -1 if z>0 else 1 if z<0 else 0
        order=sorted(range(len(D)),key=cmp_to_key(compare))
        turns=[(cross(D[i],D[j]),dot(D[i],D[j])) for i,j in zip(order,np.roll(order,-1))]
        span=bool(any(x>0 for x,y in turns) and all(x>0 or x.is_zero() and y>0 for x,y in turns))
        unique=np.unique(angles);gaps=np.diff(np.r_[unique,unique[0]+2*np.pi]);angle=unique[np.argmax(gaps)]+max(gaps)/2
        raw=[(arb(p['upper'])-arb(p['lower']))/norm(d) for p,d in zip(points,D)]
        metric_gaps=[(arb(p['upper'])-arb(p['lower']))/norm(dual(d)) for p,d in zip(points,D)]
        vertices=[];uncertain=False
        if span:
            for i in range(len(D)):
                for j in range(i):
                    det=cross(D[i],D[j])
                    if det.is_zero():continue
                    if det.contains(0):uncertain=True;continue
                    v=[(U[i]*D[j][1]-D[i][1]*U[j])/det,(D[i][0]*U[j]-U[i]*D[j][0])/det]
                    if not any(dot(d,v)>u for d,u in zip(D,U)):vertices.append(v)
        indices=np.arange(len(P))
        try:indices=ConvexHull(mids).vertices if len(P)>=3 else indices
        except QhullError:pass
        edges=list(zip(indices,np.roll(indices,-1)));distances=[];normals=[];witnesses=[]
        for v in vertices:
            q=metric(v);vm=np.array([float(x.mid()) for x in q]);choices=[]
            for i,j in edges:
                w=mids[j]-mids[i];alpha=float(np.clip((vm-mids[i])@w/(w@w),0,1)) if w@w else 0.
                choices.append((np.linalg.norm(vm-mids[i]-alpha*w),i,j,alpha))
            _,i,j,alpha=min(choices);t=arb(alpha);z=[(1-t)*x+t*y for x,y in zip(Q[i],Q[j])]
            delta=[x-y for x,y in zip(q,z)];distances.append(norm([abs(x).upper() for x in delta]))
            witnesses.append(dict(endpoints=[int(i),int(j)],alpha=alpha,distance_upper=up(distances[-1])))
            normals.append(vm-((1-alpha)*mids[i]+alpha*mids[j]))
        closed=bool(span and vertices and not uncertain);distance=max(map(up,distances),default=None) if closed else None
        ready=bool(closed and all(x<=arb(1)/100 for x in distances))
        d=np.array([np.cos(angle),np.sin(angle)]);action='add_direction'
        if distances and np.linalg.norm(normals[np.argmax(list(map(up,distances)))]):
            n=normals[np.argmax(list(map(up,distances)))];n=n/np.linalg.norm(n);d=np.array([float((a*n[0]+b*n[1]).mid()),float((c*n[1]).mid())])
        d/=float(norm(dual(list(map(arb,d)))).mid());budget=((arb(1)/400)*norm(dual(list(map(arb,d))))).lower()
        return dict(geometry_ready=ready,outer_closed=closed,max_distance_upper=distance,
            published_positive_x=mixtures is not None,inner_vertices=[[float(v.mid()) for v in p] for p in P],
            inner_vertex_enclosures=[[encode_real_ball(v) for v in p] for p in P],inner_mixtures=mixtures,
            metric='q=(x,y)' if epsilon is None else 'q=(x/(5/(16*pi^2*(92/140)^2)),(y+x/15)/epsilon)',epsilon=epsilon,
            budgets=dict(distance=.01,metric_support=.0025),distance_witnesses=witnesses,
            next_support=None if ready else dict(action=action,direction=d.tolist(),gap=float(np.nextafter(float(budget.lower()),-np.inf)),normal_units='unit normal in declared region metric'),
            metric_support_gap_upper=max(map(up,metric_gaps)),raw_support_gap_upper=max(map(up,raw)),
            outer_vertices=[[float(x.mid()) for x in v] for v in vertices],outer_vertex_enclosures=[[encode_real_ball(x) for x in v] for v in vertices],
            direction_coverage=dict(unique_normals=len(unique),angles_degrees=np.degrees(unique).tolist(),largest_gap_degrees=float(np.degrees(max(gaps))),positive_span=span,domain_half_planes=int(mixtures is not None)),
            geometry_scope='Arb bounds; float hull only selects segments. Finite H, excluding quadrature/continuum errors.')

def certify_reference_section(ir,uv,x):
    from . import read_json
    with ctx.workprec(384):
        xx=arb(float(x));sources={}
        def section(points,planes,name):
            P=[list(map(decode_real_ball,p['target_enclosures'])) if 'target_enclosures' in p else
               list(map(arb,read_json(Path(p['report']))['outer']['targets'])) for p in points]
            feasible=[]
            for i,a in enumerate(P):
                if (a[0]-xx).is_zero():feasible.append((a[1],[i,i],arb(0)))
                for j,b in enumerate(P[:i]):
                    den=b[0]-a[0]
                    if den.contains(0):continue
                    t=(xx-a[0])/den
                    if t>=0 and t<=1:feasible.append(((1-t)*a[1]+t*b[1],[i,j],t))
            lower=min(feasible,key=lambda v:float(v[0].upper()));upper=max(feasible,key=lambda v:float(v[0].lower()))
            support=[((arb(p['upper'])-arb(p['direction'][0])*xx)/arb(p['direction'][1]),p) for p in planes if p['direction'][1]!=0]
            lo=max((v for v in support if v[1]['direction'][1]<0),key=lambda v:float(v[0].lower()))
            hi=min((v for v in support if v[1]['direction'][1]>0),key=lambda v:float(v[0].upper()))
            sources[name]=dict(inner=[dict(coefficients=[points[i]['coefficients'] for i in v[1]],second_weight=v[2].str(25)) for v in (lower,upper)],
                outer=[v[1]['report'] for v in (lo,hi)])
            return [lo[0].lower(),lower[0].upper()],[upper[0].lower(),hi[0].upper()]
        il,iu=section(ir,ir,'IR');ul,uu=section(uv,ir+uv,'UV')
        top=[iu[0]-uu[1],iu[1]-uu[0]];bottom=[ul[0]-il[1],ul[1]-il[0]];difference=[top[0]-bottom[1],top[1]-bottom[0]]
        interval=lambda v:[float(np.nextafter(float(v[0].lower()),-np.inf)),float(np.nextafter(float(v[1].upper()),np.inf))]
        ratio=[top[0]/bottom[1],top[1]/bottom[0]] if bottom[0]>0 else None
        return dict(x=x,IR_lower=interval(il),IR_upper=interval(iu),UV_lower=interval(ul),UV_upper=interval(uu),
            upper_shrink=interval(top),lower_rise=interval(bottom),upper_minus_lower=interval(difference),upper_to_lower_ratio=None if ratio is None else interval(ratio),
            upper_change_larger=True if difference[0]>0 else False if difference[1]<0 else None,sources=sources,
            scope='Arb at saved x: feasible segments and global supports; no dense-region or continuum claim')

def cardinal_projection(old,new):
    if old==new:return np.eye(old)
    n=np.arange(1,min(old,new)+1);phi=(np.arange(old)+.5)*np.pi/old
    W=2/old*np.sin(n[:,None]*phi);W[n==old]/=2
    return np.sin((np.arange(new)[:,None]+.5)*np.pi/new*n)@W


from .certificates import joint_audit
