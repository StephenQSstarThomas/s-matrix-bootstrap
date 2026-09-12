"""Single finite-problem calculation entry and worker supervision."""
import sys

# Route SDP before importing the historical optimization implementation.
if __name__ == '__main__' and sys.argv[1:2] == ['sdp']:
    from .sdp.__main__ import main as sdp_main
    raise SystemExit(sdp_main(sys.argv[2:]))

from math import isfinite
from pathlib import Path
import argparse
import json
import numpy as np
import os
import signal
import subprocess
import sys
import time
from . import encode_real_ball, read_json, workspace, write_json

def resolve_prescription(args,record):
    coords=record.get('coordinates',record.get('coordinate_system',''));name=record.get('prescription','pv-midpoint' if 'PV-midpoint' in coords else None)
    from . import AMPLITUDE_COORDINATES
    if name not in AMPLITUDE_COORDINATES or getattr(args,'prescription','auto') not in ('auto',name):raise ValueError('Amplitude prescription mismatch')
    args.resolved_prescription=name
    return name

def source_contract(args):
    from . import AMPLITUDE_COORDINATES
    name=getattr(args,'resolved_prescription','pv-midpoint')
    if name=='pv-midpoint' and args.command in ('boundary','dual') and (
        getattr(args,'unitarity_scope','sampled')!='sampled' or args.infinity!='free'):
        raise ValueError('PV requires sampled unitarity and free infinity')
    return name,AMPLITUDE_COORDINATES[name]

def support_data(args,progress):
    parent=read_json(args.preparation/'report.json');meta=parent['amplitude_model']
    scheme=resolve_prescription(args,meta);_,coords=source_contract(args)
    if meta['coordinate_system']!=coords:raise ValueError('Require unsubtracted source coordinates')
    if getattr(args,'additional_constraints',None) is not None:raise ValueError('Use prepare --unitarity-energies for authenticated additional analytic rows')
    data=np.load(args.preparation/'amplitude.npz');H=data['rows'];M=meta['M'];L=meta['L']
    from .sampling import prepared_sampling
    energies,waves=prepared_sampling(data,M,L);base_count=3*M*L
    high_spin=np.empty((0,2))
    from .certificates import verify_scattering_layout
    from flint import arb
    verify_scattering_layout(H,[arb(float(v)) for v in np.pi*np.sqrt(1-4/energies)],M,L)
    return H,np.pi*np.sqrt(1-4/energies),energies,waves,M,L,high_spin,base_count

def calculation_inputs(args):
    from .io import calculation_inputs as collect
    return collect(args)


def worker(args):
    from . import json_scalar, save_evaluation; from .analysis import gauge_phases, gauge_regions, select_gauge; from .imaginary import coherent_boundary; from .operators import prepare_amplitude
    import flint, numpy, scipy
    from flint import ctx
    ctx.prec = 2048
    start = time.monotonic()
    folder = args.output;input_snapshot=calculation_inputs(args)
    report = dict(command=args.command, status='running', artifacts={}, inputs={},
        parameters={key: str(value) if isinstance(value, Path) else value
            for key, value in vars(args).items() if key != 'worker'},
        environment=dict(python=sys.version, flint=flint.__version__, numpy=numpy.__version__, scipy=scipy.__version__,blas_threads=os.getenv('OPENBLAS_NUM_THREADS')),
        unitarity_scope=getattr(args,'unitarity_scope','sampled'),physical_unitarity_accepted=False,
        scope='Finite model; regions/phases require verification.')
    write_json(folder/'worker.json', report)

    def progress(**data):
        item = dict(elapsed_seconds=time.monotonic()-start, **data)
        with (folder/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps(item,default=json_scalar)+'\n')
        if data.get('stage') in ('conic_optimization','native_conic'):
            report['solver_attempt'] = item
            write_json(folder/'worker.json', report)

    try:
        if args.command in ('boundary','dual'):
            report.update(coherent_boundary(args, progress))
        elif args.command=='support-probe':
            from .probes import support_probe
            report.update(support_probe(args,progress))
        elif args.command=='regions':
            from .io import summarize_regions
            report.update(summarize_regions(folder,args.support_runs,args.reference_directory,args.nodes,args.chiral_norm))
        elif args.command=='resolution':
            from .operators import resolution_run as run
            report.update(run(args,progress))
        elif args.command=='compare':
            from .model import compare_selected_profiles
            report.update(compare_selected_profiles(args))
        elif args.command=='analytic-audit':
            from .analytic import audit_reconstruction
            report.update(audit_reconstruction(args,progress))
        elif args.command=='source-audit':
            from .basis import audit_cardinal
            report.update(audit_cardinal(args,progress))
        elif args.command in ('gauge-regions','select-gauge','gauge-phases'):
            report.update({'gauge-regions':gauge_regions,'select-gauge':select_gauge,'gauge-phases':gauge_phases}[args.command](args))
        elif args.command=='select':
            from .io import select_chiral_amplitude
            report.update(select_chiral_amplitude(args))
        elif args.command=='profiles':
            from .io import plot_profiles
            from .spectra import direct_profiles
            report.update(direct_profiles(args) if args.profile_kind=='analytic-window' else plot_profiles(args))
            signature=report['profiles'][0]['selection']['model_signature']
            if any(p['selection']['model_signature']!=signature for p in report['profiles']):raise ValueError('Profiles require one declared model')
            report.update(model_signature=signature,scope='Saved profiles; no continuum certificate.')
            if args.profile_kind=='selection':
                if len(report['profiles'])!=1:raise ValueError('One representative required')
                write_json(folder/'coefficients.json',read_json(report['profiles'][0]['coefficients']))
                write_json(folder/'selection.json',report['profiles'][0]['selection'])
        elif args.command == 'prepare':
            report.update(prepare_amplitude(args, progress))
        elif args.command == 'prepare-current':
            from scipy.sparse import save_npz
            from .model import current_operators
            from . import UVConfig
            H,kap,ss,ww,M,L,_,_=support_data(args,progress)
            uv=UVConfig(moment_source=args.moment_source,sr_error=args.sr_error,cutoff=args.fesr_cutoff,mq_rule=args.mq_rule)
            current=current_operators(H,kap,ss,ww,M,uv,args.bits)
            current['metadata']['ff_endpoint_order']=2 if args.asymptotic_zeros else 0
            current['metadata']['asymptotic_unitarity']=args.asymptotic_unitarity
            for name,matrix in current['matrices'].items():save_npz(folder/(name+'.npz'),matrix)
            numpy.savez_compressed(folder/'current_data.npz',**current['arrays'])
            report.update(status='prepared_current',M=M,L=L,current_model=current['metadata'],optimization_performed=False,
                matrix_shapes={k:list(v.shape) for k,v in current['matrices'].items()},scattering_samples=len(kap))
        elif args.command == 'evaluate':
            from .kernels import source_rows
            if args.coefficients is None:
                raise ValueError('Require complete --coefficients JSON')
            from . import digest
            coefficient_hash=digest(args.coefficients);report['inputs'][str(args.coefficients.resolve())]=dict(sha256=coefficient_hash)
            record=read_json(args.coefficients);scheme=resolve_prescription(args,record);bits=args.bits
            if args.scan_nodes:raise ValueError('The declared provider requires native physical nodes')
            registry=getattr(args,'source_registry',None)
            if scheme=='analytic-cardinal':
                from . import digest
                if args.preparation is None:raise ValueError('Analytic evaluation requires its preparation')
                preparation=_path(args.preparation)/'report.json';meta=read_json(preparation)['amplitude_model']
                if (meta['prescription'],meta['M'],meta['L'],meta['coordinate_system'])!=(scheme,args.nodes,args.waves,record.get('coordinates')):
                    raise ValueError('Analytic evaluation/preparation identity mismatch')
                saved_registry=_path(meta['source_registry']) if meta.get('source_registry') else None;expected=meta.get('source_registry_sha256')
                if registry is not None and _path(registry)!=saved_registry:raise ValueError('Registry override differs from preparation')
                if saved_registry and digest(saved_registry)!=expected:raise ValueError('Prepared registry hash changed')
                registry=saved_registry;report['inputs'][str(preparation)]=dict(sha256=digest(preparation))
            source=source_rows(args.nodes,args.waves,bits,args.angular_order,args.subtracted,prescription=scheme,registry=registry)
            if scheme=='analytic-cardinal' and source.metadata['source_registry_sha256']!=expected:raise ValueError('Registry changed during source construction')
            if scheme=='analytic-cardinal':source.attach_preparation(args.preparation)
            if record.get('coordinates') != source.coordinate_name:
                raise ValueError('Coefficient coordinates mismatch')
            if args.primary_waves:source.waves=[(0,0),(2,0),(1,1)]
            kinematic_energies=[source.exact_energy(s) for s in args.energies] if scheme=='analytic-cardinal' else args.energies
            values = source.evaluate(record['coefficients'], kinematic_energies)
            report['inputs'].update(getattr(source,'inputs',{}))
            evaluation = dict(energies=args.energies, waves=source.waves,
                f=[[[float(v.real.mid()), float(v.imag.mid())] for v in row] for row in values],
                f_enclosures=[[[encode_real_ball(v.real),encode_real_ball(v.imag)] for v in row] for row in values],
                amplitude_model=source.metadata, constraints_checked=['unitarity at supplied physical energy/wave pairs'])
            evaluation['coefficients']=str(args.coefficients.resolve());evaluation['coefficient_sha256']=coefficient_hash
            selected=args.coefficients.parent/'selection.json'
            if selected.exists():
                evaluation['selection']=read_json(selected);report['inputs'][str(selected.resolve())]=dict(sha256=digest(selected))
                if evaluation['selection'].get('coefficient_sha256',coefficient_hash)!=coefficient_hash:raise ValueError('Selected amplitude changed before evaluation')
                for name,suffix in (('current','json'),('joint','npz')):
                    expected=evaluation['selection'].get(name+'_sha256')
                    if expected and digest(args.coefficients.parent/(name+'.'+suffix))!=expected:raise ValueError('Selected joint point changed before evaluation')
            from .model import scattering_matrix, unitarity_margin
            evaluation['eta'] = [[float(abs(scattering_matrix(v,s)).upper()) for v in row]
                if s>=4 else None for s,row in zip(kinematic_energies,values)]
            evaluation['maximum_eta'] = max((max(row) for row in evaluation['eta'] if row is not None),default=None)
            margins = [[unitarity_margin(v,s) for v in row] if s>=4 else None
                for s,row in zip(kinematic_energies,values)]
            physical = [v for row in margins if row is not None for v in row]
            evaluation['unitarity'] = dict(
                criterion='Eq.(2.12): 1-|S|^2 = kappa*(2 Im f-kappa*|f|^2) >= 0',
                samples_checked=len(physical),
                status=('violated' if any(v<0 for v in physical) else
                    'sampled_passed' if physical and all(v>=0 for v in physical) else 'unresolved'),
                margin_intervals=[[[float(np.nextafter(float(v.lower()),-np.inf)),float(np.nextafter(float(v.upper()),np.inf))] for v in row] if row is not None else None for row in margins],
                minimum_margin=float(min(v.lower() for v in physical)) if physical else None,
                scope='Declared sampled pairs; analytic source enclosures included' if scheme=='analytic-cardinal' else 'Declared sampled pairs; discretization errors excluded',
                continuum_certified=False, omitted_spins_checked=False, high_energy_tail_certified=False)
            evaluation['unitarity'].update(analytic_tail=dict(status='not_applicable'),high_spin_necessary=dict(status='not_applicable',violations=[]))
            violations=[dict(energy=s,isospin=I,ell=ell,margin=v.str(25)) for s,row in zip(args.energies,margins) if row is not None
                for (I,ell),v in zip(source.waves,row) if v<0]
            write_json(folder/'violations.json',dict(violations=violations,high_spin=[],scope='Supplied physical energy/wave pairs'))
            write_json(folder/'constraints.json',dict(violations=violations,high_spin=[]))
            save_evaluation(folder,evaluation)
            report.update(status='evaluated', amplitude_model=source.metadata, unitarity=evaluation['unitarity'])
        else:
            raise ValueError('Unsupported calculation command: '+args.command)
    except Exception as error:
        report.update(status='inconclusive', error=str(error), exception_type=type(error).__name__)
    report.setdefault('inputs',{}).update(input_snapshot)
    origin=args.start_mu_provenance
    if origin.get('checkpoint'):
        report.setdefault('inputs',{}).update({origin['checkpoint']:dict(sha256=origin['sha256']),origin['report']:dict(sha256=origin['report_sha256'])})
    model=report.get('amplitude_model') or report.get('model_signature') or {}
    report.update(elapsed_seconds=time.monotonic()-start,prescription=model.get('prescription',report.get('prescription',getattr(args,'resolved_prescription','pv-midpoint'))),unitarity_scope=model.get('unitarity_scope',report['unitarity_scope']))
    write_json(folder/'worker.json', report)

def main():
    from . import check_core_unchanged, digest, resolve_start_mu, snapshot_core
    parser = argparse.ArgumentParser(description=__doc__)
    add=parser.add_argument
    add('command', choices=['status','prepare','prepare-current','evaluate','boundary','dual','regions','select','profiles','gauge-regions','select-gauge','gauge-phases','resolution','compare','analytic-audit','source-audit','support-probe'])
    add('--comparison-manifest',type=Path)
    add('--region-summary',type=Path,default=Path('results/runs/stage_B_delivery_20260907/regions/regions.json'))
    add('--primary-waves',action='store_true',help='S0/S2/P1 with full coefficients')
    add('--profile-runs',nargs='+',default=[])
    add('--profile-kind',choices=['subthreshold','selection','phase','analytic-window'],default='subthreshold')
    add('--phase-reference',type=Path)
    add('--baseline-profile',type=Path,help='IR evaluation for UV comparison')
    add('--selection-reference',type=Path)
    add('--support-runs',nargs='+',default=[])
    add('--reference-directory',type=Path,default=Path('results/runs/stage_B_reference_20260906'))
    add('--output',type=Path)
    add('--density-limit', type=float)
    add('--preparation', type=Path)
    add('--prescription',choices=['auto','pv-midpoint','analytic-cardinal'],default='auto')
    add('--source-registry',type=Path)
    add('--solver',choices=['barrier','clarabel','centered','scs'],default='barrier')
    add('--unitarity-scope',choices=['sampled'],default='sampled')
    add('--native-backend',choices=['auto','qdldl','faer'],default='auto')
    add('--native-coordinates',choices=['unsubtracted','subtracted','absorptive'],default='unsubtracted')
    add('--current-preparation',type=Path)
    add('--joint-feasibility',action='store_true')
    add('--fixed-amplitude',action='store_true',help='Test the current fiber of the unchanged C_flat')
    add('--solver-seconds', type=float, default=300.)
    add('--dual-seconds',type=float,default=120.)
    add('--dual-method',choices=['highs-ds','highs-ipm','highs-primal'],default='highs-ds')
    add('--additional-constraints',type=Path,help='Physical rows JSON/NPZ')
    add('--start-mu', type=float)
    add('--unit-newton-step',action='store_true',help='Joint Newton: minimize on the feasible segment with step <=1; numerical path only')
    add('--gap', type=float, default=1e-4)
    add('--nodes', type=int, default=50)
    add('--processes',type=int,default=1)
    add('--waves', type=int, default=10, help='Partial waves per isospin')
    add('--bits', type=int, default=384)
    add('--angular-order', type=int, default=24, help='Gauss order per panel')
    add('--sampling-factor',type=int,default=1,help='Odd refinement')
    add('--unitarity-energies',type=float,nargs='*',default=[],help='Extra energies')
    add('--subtracted', action='store_true', help='Use q=H-H(0) coordinates')
    add('--infinity', choices=['free'], default='free', help='Asymptotic constant')
    add('--asymptotic-zeros',action='store_true',help='Analytic slice: T0=0; joint currents additionally have F0/F1 zeros of order at least two at z=-1')
    add('--asymptotic-unitarity',action='store_true',help='Also enforce five derived leading high-energy unitarity inequalities')
    add('--moment-source', choices=['printed', 'eq250'], default='printed')
    add('--sr-error', choices=['raw-absolute', 'normalized-absolute'], default='raw-absolute')
    add('--fesr-cutoff', choices=['hard-midpoint', 'clipped-phi'], default='hard-midpoint')
    add('--mq-rule', choices=['arithmetic-mean', 'rms'], default='arithmetic-mean')
    add('--coefficients', type=Path, help='C_flat JSON')
    add('--coefficient-key',default='coefficients',choices=['coefficients'],help='Saved vector key')
    add('--interior-coefficients',type=Path)
    add('--energies', type=float, nargs='+', default=[.5, 1., 1.5, 2., 3., 4., 9., 25.])
    add('--scan-nodes',type=int,default=0)
    add('--mode', choices=['pure', 'chiral','gauge'], default='chiral')
    add('--chiral-tolerance', type=float, default=.002)
    add('--chiral-norm',choices=['separate-l2','combined-l2'])
    add('--chiral-barrier-weight',type=float)
    add('--objective',choices=['projection','watson'],default='projection')
    add('--resume-objective',action='store_true')
    add('--probe-anchor',type=Path)
    add('--probe-node',type=int)
    add('--probe-wave',choices=['S0','S2','P1'])
    add('--probe-component',choices=['real','imag'])
    add('--probe-sign',type=int,choices=[-1,1],default=1)
    add('--probe-budget',type=float,default=.02)
    add('--probe-audit-only',action='store_true')
    add('--probe-xref',type=float,help='Audit a section inner bound from two --support-runs')
    add('--require-center',action='store_true',help='Select an audited converged center for phase comparisons')
    add('--direction', type=float, nargs=2, default=[1., 0.])
    add('--fixed-x', type=float)
    add('--ray',action='store_true',help='Saved joint projection ray')
    add('--seconds', type=float)
    add('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.unit_newton_step and not (args.solver in ('barrier','centered') and (args.command=='resolution' and args.joint_feasibility and not args.profile_runs or args.command=='boundary' and args.mode=='gauge' and not args.joint_feasibility)):
        parser.error('--unit-newton-step applies to full joint Newton calculations')
    if args.asymptotic_unitarity and not args.asymptotic_zeros:parser.error('Leading asymptotic unitarity requires --asymptotic-zeros')
    if args.asymptotic_zeros:
        if args.prescription!='analytic-cardinal' or not (args.command in ('prepare-current','resolution','source-audit','evaluate') or args.command in ('boundary','dual')):parser.error('Endpoint slice requires an explicit analytic calculation')
        args.infinity='zero'
    from .model import declared_chiral_norm
    try:args.chiral_norm=declared_chiral_norm(args)
    except ValueError as error:parser.error(str(error))
    needs_preparation=args.command in ('boundary','dual','prepare-current','support-probe') or args.command=='resolution' and not args.profile_runs or args.current_preparation is not None
    if needs_preparation and args.preparation is None:parser.error('This calculation requires an explicit --preparation')
    if args.fixed_amplitude and not (args.command in ('boundary','dual') and args.mode=='gauge' and args.joint_feasibility and args.coefficients is not None and args.fixed_x is None and not args.ray and not any(args.direction)):parser.error('Use fixed amplitude with gauge joint-feasibility, coefficients, zero direction')
    if args.processes<1:parser.error('--processes must be positive')
    os.chdir(workspace())
    resolve_start_mu(args)
    if args.command == 'status':
        print((workspace()/'STATUS.md').read_text())
        return
    if args.output is None:
        parser.error('--output is required for a calculation')
    args.output = args.output.resolve()
    run_root = (workspace()/'results/runs').resolve()
    if args.output == run_root or not args.output.is_relative_to(run_root):
        parser.error('--output must be a new subdirectory of results/runs')
    if args.command in ('regions','gauge-regions') and not args.support_runs:parser.error(args.command+' needs --support-runs')
    cap=args.seconds if args.seconds is not None else {'boundary':900,'dual':900,'prepare':900,'prepare-current':300,'evaluate':300,'regions':60,'select':120,'profiles':60,'gauge-regions':60,'select-gauge':60,'gauge-phases':60,'resolution':300,'compare':60,'analytic-audit':600,'source-audit':300,'support-probe':900}[args.command]
    if not isfinite(cap) or cap<=0:parser.error('--seconds must be finite and positive')
    if args.worker:
        worker(args)
        return
    sources={} if args.command=='regions' else snapshot_core(args.output)
    if args.command=='regions':args.output.mkdir(parents=True,exist_ok=False)
    command = [sys.executable, '-m', 'smatrix_bootstrap.run', *sys.argv[1:], '--worker']
    start = time.monotonic()
    timed_out = False
    threads=os.getenv('SMATRIX_BLAS_THREADS','8')
    with (args.output/'stdout.log').open('w') as stdout, (args.output/'stderr.log').open('w') as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr, env=dict(os.environ,OPENBLAS_NUM_THREADS=threads,OMP_NUM_THREADS=threads,MKL_NUM_THREADS=threads),start_new_session=True)
        print(json.dumps(dict(pid=process.pid, output=str(args.output), cap_seconds=cap)), flush=True)
        try:
            code = process.wait(timeout=cap)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid,signal.SIGKILL)
            code = process.wait()
    path = args.output/'worker.json'
    report = read_json(path) if path.exists() else dict(status='inconclusive')
    changes=[] if args.command=='regions' else check_core_unchanged(sources)
    inputs_changed = [name for name, item in ({} if args.command=='regions' else report.get('inputs', {})).items()
        if not Path(name).is_file() or digest(name) != item['sha256']]
    report.update(core_sources=sources, source_changes=changes, input_changes=inputs_changed,
        exit_code=code, timeout=timed_out, cap_seconds=cap, supervisor_seconds=time.monotonic()-start)
    if code != 0 or timed_out or changes or inputs_changed:
        report.update(status='inconclusive',support_optimality_certified=False,physical_unitarity_accepted=False)
    write_json(args.output/'report.json', report)
    print(json.dumps({k: report.get(k) for k in ('status','objective','sampled_feasible','joint_feasible','support_optimality_certified','supervisor_seconds')}))

from .gauge import gauge_support

def _path(value):
    return (workspace()/Path(value)).resolve()

from .certificates import model_signature as _signature

if __name__=='__main__':main()
