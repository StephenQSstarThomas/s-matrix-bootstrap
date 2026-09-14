"""Public SDPB-only calculation entry, checks and figure delivery."""
from __future__ import annotations

import argparse
import json
import os
import sys

from . import constraints as C


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] in (["solve"], ["sdpb"]):
        from .cli import main as sdpb_main
        return sdpb_main(argv[1:])
    p = argparse.ArgumentParser("smatrix_bootstrap.sdp")
    p.add_argument("command", choices=["selfcheck", "prereg", "figures", "ir-figures", "mma-audit", "audit-run", "angle-audit", "currents", "subthreshold", "amplitude-export", "precision-probe", "contrast", "ir-select", "uv-plan", "uv-select", "support", "refine", "sine-replay", "sine-check"])
    p.add_argument("--out")
    p.add_argument("--M", type=int, default=50)
    p.add_argument("--source-report")
    p.add_argument("--compare-report")
    p.add_argument('--ref-report',help='uv-select: completed upper support at physical x_ref')
    p.add_argument('--near-report',help='uv-select: completed upper support at the frozen nearby x')
    p.add_argument("--source-snapshot", help="authenticated non-executable source JSON.gz for historical family checks")
    p.add_argument("--angle-row", nargs=3, type=int, action="append", metavar=("I","ELL","NODE"))
    p.add_argument("--angle-midpoint", nargs=3, type=int, action="append", metavar=("I","ELL","LEFT_NODE"))
    p.add_argument("--angle-point", nargs=3, action="append", metavar=("I","ELL","S"))
    p.add_argument("--bits", type=int, default=384)
    p.add_argument('--precision',type=int,help='refine: SDPB arithmetic bits; input coefficients are unchanged')
    p.add_argument("--error-target", default="1e-12", help="currents: absolute angular tail target")
    p.add_argument("--timeout", type=float, default=600, help="diagnostic runtime budget; explicit support values override its saved solver budget")
    p.add_argument('--resume',action='store_true',help='refine: resume a terminal checkpoint with exactly the same PMP, precision and local rank layout')
    p.add_argument('--warm-start',action='store_true',help='support: initialize from the accepted source checkpoint; same block structure, precision and local rank layout required')
    p.add_argument("--basis")
    p.add_argument("--operator-dps", type=int, default=40)
    p.add_argument("--digits", type=int, default=30)
    p.add_argument("--ir-report")
    p.add_argument("--uv-report")
    p.add_argument('--ir-selection-report',help='contrast: revalidate a geometry-only IR endpoint selection receipt')
    p.add_argument('--subthreshold-replay',action='append',help='figures: verified post-solve evaluation overlay; repeat for multiple source amplitudes')
    p.add_argument("--reference-only",action="store_true",help="contrast: compare accepted ref leaves; registered three-point chain remains incomplete")
    p.add_argument("--point",choices=['tip','ref','mid'],help='support: registered point; ref by default')
    p.add_argument('--direction',type=float,nargs=2,help='support: arbitrary nonzero finite linear direction')
    p.add_argument('--fix-f00',type=float,help='support: optional fixed f00(3), requires --direction')
    p.add_argument('--face-margin',type=float,help='support: face diagnostic; slab d.(f00,f11) >= source value - margin, keeping the source direction and section')
    p.add_argument('--functional',nargs=4,metavar=('KIND','WAVE','NODE','SENSE'),
                   help='support: secondary objective on the face, e.g. ImKH P1 38 max (kinds ImKH ImS ImF rho; waves S0 P1)')
    p.add_argument("--duality-gap",type=float)
    a = p.parse_args(argv)
    if a.resume and a.command!='refine':p.error('--resume is supported only by refine')
    if a.warm_start and a.command!='support':p.error('--warm-start is supported only by support')
    if a.command!='uv-select' and (a.ref_report or a.near_report):
        p.error('--ref-report and --near-report are supported only by uv-select')
    if a.command!='support' and (a.direction is not None or a.fix_f00 is not None):
        p.error('--direction and --fix-f00 are supported only by support')
    if a.command!='support' and (a.face_margin is not None or a.functional is not None):
        p.error('--face-margin and --functional are supported only by support')
    if a.command not in ('contrast','ir-figures') and a.ir_selection_report is not None:
        p.error('--ir-selection-report is supported only by contrast or ir-figures')
    if a.command!='figures' and a.subthreshold_replay:
        p.error('--subthreshold-replay is supported only by figures')
    if a.command == "selfcheck":
        import subprocess
        here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        try:threads=int(os.environ.get('SDP_TEST_THREADS','1'))
        except ValueError:p.error('SDP_TEST_THREADS must be a positive integer')
        if threads<1:p.error('SDP_TEST_THREADS must be a positive integer')
        env=dict(os.environ,**{name:str(threads) for name in
                 ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')})
        return subprocess.call([sys.executable, "-m", "pytest", "-q",
                                os.path.join(here, "tests", "sdp")],env=env)
    if not a.out:
        p.error("--out is required")
    os.environ.setdefault("SDP_CACHE", os.path.join(a.out, "operator-cache"))
    if a.command=='ir-figures':
        if not a.source_report or not a.ir_selection_report:p.error('ir-figures requires --source-report and --ir-selection-report')
        if a.source_snapshot or any(t.split('=')[0] in ('--bits','--precision','--digits','--operator-dps') for t in argv):
            p.error('ir-figures uses the authenticated saved amplitude; precision overrides do not apply')
        from .delivery import deliver_ir
        r=deliver_ir(a.source_report,a.ir_selection_report,a.out)
        print(json.dumps({'report':os.path.join(a.out,'report.json'),'C4_source':r['C4_source']['verdict'],
                          'selected':r['selected'],'complete_C3':False}))
        return 0
    if a.command in ('uv-plan','uv-select'):
        if not a.source_report:p.error('--source-report must select the completed UV tip or frozen UV plan')
        if a.source_snapshot or any(t.split('=')[0] in ('--bits','--precision','--digits','--operator-dps') for t in argv):
            p.error('UV selection inherits authenticated source data; precision overrides do not apply')
        from .uv_selection import plan_uv,select_uv
        if a.command=='uv-plan':
            r=plan_uv(a.source_report,a.out)
            summary={k:r['plan'][k] for k in ('x_ref','x_near','x_tip','rule')}
        else:
            if not a.ref_report or not a.near_report:p.error('uv-select requires --ref-report and --near-report')
            r=select_uv(a.source_report,a.ref_report,a.near_report,a.out)
            summary={'role_mapping':r['selection']['role_mapping'],'legacy_registered_chain_complete':False}
        print(json.dumps({'report':os.path.join(a.out,'report.json'),'schema':r['schema'],**summary}))
        return 0
    if a.command == "mma-audit":
        from .mma import audit
        return 0 if audit(a.M, a.out)["passed"] else 1
    if a.command=='refine':
        if not a.source_report or a.precision is None or a.precision<64:
            p.error('refine requires --source-report and --precision >=64')
        from pathlib import Path
        from .sdpb import run_once,Settings
        from .spec import ModelSpec
        source=json.loads(Path(a.source_report).read_text())
        if source.get('status') not in ('numerically_accepted','not_accepted','solver_failed','readback_failed','preprocess_failed') or 'pmp' not in source:
            p.error('Select a terminal leaf with a complete saved PMP; do not duplicate a live run')
        if a.resume and a.precision!=source['settings']['precision']:p.error('resume cannot change checkpoint precision')
        cfg=Settings(**source['settings']);cfg.precision=a.precision;cfg.checkpoint_interval=300
        if any(t.split('=')[0]=='--timeout' for t in argv):
            if a.timeout<=0 or not a.timeout.is_integer():p.error('refine --timeout must be a positive whole number of seconds')
            cfg.timeout=int(a.timeout)
        if a.duality_gap is not None:cfg.duality_gap=a.duality_gap
        r,_,_=run_once(ModelSpec(**source['spec']),a.out,source['direction'],source.get('fix_f00'),cfg,
                       prepared=a.source_report,exact_input=True,**({'resume':True} if a.resume else {}))
        return 0 if r['accepted'] else 1
    if a.command == "support":
        if not a.source_report:p.error('--source-report is required')
        if a.point is not None and (a.direction is not None or a.fix_f00 is not None):
            p.error('--point conflicts with --direction/--fix-f00')
        if a.direction is None and a.fix_f00 is not None:p.error('--fix-f00 requires --direction')
        if a.face_margin is not None and a.functional is None:p.error('--face-margin needs --functional')
        if a.functional is not None and (a.point is not None or a.direction is not None or a.fix_f00 is not None):
            p.error('a functional run keeps the source direction and section; drop --point/--direction/--fix-f00')
        from .sdpb import support_from_saved
        extra={'warm_start':True} if a.warm_start else {}
        if a.functional is not None:
            kind,wave,node,sense=a.functional
            if not node.lstrip('-').isdigit():p.error('--functional NODE must be an integer node index')
            extra.update(functional={'kind':kind,'wave':wave,'node':int(node),'sense':sense})
            if a.face_margin is not None:extra['face_margin']=a.face_margin
        if any(t.split('=')[0]=='--timeout' for t in argv):
            if a.timeout<=0 or not a.timeout.is_integer():p.error('support --timeout must be a positive whole number of seconds')
            extra['timeout']=int(a.timeout)
        r=support_from_saved(a.source_report,a.out,a.point,a.duality_gap,direction=a.direction,fix_f00=a.fix_f00,**extra)
        return 0 if r['accepted'] else 1
    if a.command == 'sine-replay':
        if not a.source_report:p.error('--source-report is required')
        from .sine import replay
        print(json.dumps(replay(a.source_report,a.out,a.operator_dps),indent=2))
        return 0
    if a.command == 'sine-check':
        from .sine import mma_check
        r=mma_check(a.out);print(json.dumps(r,indent=2));return 0 if r['passed'] else 1
    if a.command == "contrast":
        if not a.ir_report or not a.uv_report:
            p.error("--ir-report and --uv-report must contain accepted chains or --reference-only ref leaves")
        from .delivery import deliver
        print(json.dumps(deliver(a.ir_report,a.uv_report,a.out,reference_only=a.reference_only,ir_selection_report=a.ir_selection_report),indent=2))
        return 0
    if a.command=='ir-select':
        if not a.source_report or not a.compare_report:p.error('ir-select requires upper --source-report and lower --compare-report')
        from .ir_selection import select_ir
        receipt=select_ir(a.source_report,a.compare_report,a.out);selection=receipt['selection']
        print(json.dumps({'receipt':os.path.join(a.out,'report.json'),'status':selection['status'],
                          'selected':selection['selected'],'scope':selection['rule']}))
        return 0 if selection['status']=='selected' else 1
    if a.command == "precision-probe":
        if not a.basis:
            p.error("--basis must select a saved coordinate basis")
        from .precision import benchmark
        print(json.dumps(benchmark(a.M, a.basis, a.out, a.operator_dps, a.digits), indent=2))
        return 0
    if a.command == "audit-run":
        if not a.source_report:
            p.error("--source-report must select a completed leaf SDPB report")
        from .accuracy import audit_saved_run
        audit_saved_run(a.source_report, a.out, a.bits, a.compare_report)
        return 0
    if a.command == 'amplitude-export':
        if not a.source_report:p.error('--source-report is required')
        if a.source_snapshot or any(t.split('=')[0] in ('--bits','--precision','--digits','--operator-dps') for t in argv):
            p.error('amplitude-export converts saved coefficients exactly; precision overrides do not apply')
        from .amplitude_export import export_amplitude
        r=export_amplitude(a.source_report,a.out,a.timeout)
        print(json.dumps({'report':os.path.join(a.out,'report.json'),'status':r['status'],
                          'n_c':r['n_c'],'seconds':r['seconds'],'representation_error':r['representation_error']}))
        return 0
    if a.command == 'subthreshold':
        if not a.source_report:p.error('--source-report is required')
        if a.source_snapshot or any(t.split('=')[0] in ('--bits','--precision','--digits','--operator-dps') for t in argv):
            p.error('subthreshold requires matching source kernels and inherits the recorded source precision')
        from .evaluation import evaluate_subthreshold
        r=evaluate_subthreshold(a.source_report,a.out,a.timeout)
        print(json.dumps({'report':os.path.join(a.out,'report.json'),'status':r['status'],
                          'points_per_wave':len(r['subthreshold']['s']),'seconds':r['seconds'],
                          'optimization_performed':False}))
        return 0
    if a.command == "angle-audit":
        if not a.source_report:p.error('--source-report is required')
        from .crosscheck import audit_angles
        rows=(a.angle_row or [])+[(I,ell,k+.5) for I,ell,k in (a.angle_midpoint or [])]
        points=[(int(I),int(ell),s) for I,ell,s in (a.angle_point or [])]
        r=audit_angles(a.source_report,a.out,a.bits,a.timeout,a.source_snapshot,rows or None,points)
        print(json.dumps(_angle_summary(r,a.out)))
        return 0 if r['passed'] else 1
    if a.command == "currents":
        if not a.source_report:
            p.error("--source-report must select a completed precise sine-cardinal UV leaf")
        from .accuracy import audit_saved_currents
        r = audit_saved_currents(a.source_report, a.out, a.bits, a.error_target, a.timeout,snapshot=a.source_snapshot)
        print(json.dumps({"report":os.path.join(a.out, "report.json"), "nodes":len(r["nodes"]),
                          "uv_ok":r["uv_checks"]["uv_ok"], "scope":r["scope"]}, indent=2))
        return 0
    if a.command == "prereg":
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, "preregistration.json"), "x") as fh:
            json.dump(preregistration(), fh, indent=1)
        return 0
    from .figures import build_all
    build_all(a.out,subthreshold_reports=a.subthreshold_replay or ())
    return 0


def _angle_summary(report,outdir):
    counts={name:0 for name in ('pass','fail','inconclusive','not_applicable')}
    for row in report.get('results',[]):
        verdict=row.get('unitarity_verdict')
        if isinstance(verdict,str) and verdict.startswith('not_applicable'):verdict='not_applicable'
        if verdict in counts:counts[verdict]+=1
    selected=len(report.get('selected_rows',[]))+len(report.get('selected_points',[]))
    return {'report':os.path.join(outdir,'report.json'),'status':report['status'],
            'integrators_agree':report['passed'],'selected_unitarity_counts':counts,
            'selected_unitarity_not_evaluated':max(0,selected-sum(counts.values())),
            'scope':'Exit status reports numerical integration agreement; selected-row unitarity is separate and does not establish full physical acceptance.'}


def preregistration() -> dict:
    """The three under-determined calibers and the acceptance thresholds.

    Frozen before P3; nothing here may be added to, removed or re-tuned
    afterwards (task section 4 and section 7).
    """
    return {
        "written_before": "P3",
        "paper": "arXiv:2309.12402v3",
        "calibers": {
            "chiral_norm": {
                "source": "(3.64) says 'with some norm'",
                "chi-a": "literal L-infinity box: 8 separate |r| <= eps_chi",
                "chi-b": "single 8-dimensional L2 ball (the authors' 2403 code caliber)",
                "chi-c": "two separate 4-dimensional L2 balls",
            },
            "fesr_tolerance": {
                "source": "(3.74) and (2.56) do not close dimensionally",
                "SR-a": "raw absolute: |M_n - T_n| <= 2e-3 per moment",
                "SR-b": "relative 10%: |M_n - T_n| <= 0.10 |T_n|",
                "SR-c": "relative 20%",
                "SR-d": "raw per-wave L2: ||(M_n - T_n)_n||_2 <= 2e-3 per wave",
            },
            "density_regularisation_B": {
                "source": "absent from the paper; present in the authors' 2403 code as "
                          "norm(rho/Mrho,4) <= 1e2, i.e. ||rho||_4 <= 377500",
                "default": "absent",
                "grid_if_needed": [1e4, 1e5, 377500.0, 1e6],
                "must_report": "whether B is active (||rho||_4 >= 0.9 B)",
            },
        },
        "primary_caliber": {"chiral_norm": "chi-b", "fesr_tolerance": "SR-b"},
        "grid": "3 x 3 calibers on the x_ref section and the three representative points",
        "fixed_inputs": {"M": 50, "L": 10, "s0": C.S0, "alpha_s": 0.4,
                         "m_q_mean": C.M_Q, "m_q_rms": C.M_Q_RMS,
                         "eps_chi_main": C.EPS_CHI_MAIN, "eps_chi_grid": list(C.EPS_CHI_GRID),
                         "eps_sr": C.EPS_SR, "eps_ff": C.EPS_FF,
                         "projection_plane": "(f00(3), f11(3))",
                         "chiral_points": list(C.CHIRAL_POINTS)},
        "representative_points": {
            "tip": "argmax f00(3)",
            "ref": "max f11(3) on the section f00(3) = x_ref = 0.0733214",
            "mid": "max f11(3) on the section f00(3) = (x_tip + x_ref)/2",
        },
        "acceptance": {
            "C1": "f00 extrema in [-2.902, 2.233] +-2%; f11 in [-0.734, 0.079] +-2%",
            "C2": "eps=2e-3 +x end 0.0826 +-5%; x_ref section width 0.00076 +-20%; "
                  "monotone in eps",
            "C3": "subthreshold rms <= 8% of f00(3); S0 zero near 0.425 / 0.305 / none",
            "C4": "delta00(0.9 GeV) in [85,110] deg; delta11(1.2 GeV) <= 25 deg",
            "C5": "x_ref section: upper shrinks 2.3e-4 +-25%, lower rises <= 1.0e-4, "
                  "ratio >= 4; UV +x end 0.0811 +-5%",
            "C6": "90 deg crossing in [795, 845] MeV for all three points, spread "
                  "<= 20 MeV, min eta(P1) >= 0.9",
            "C7": "delta00(1.196 GeV) in [85,110] deg; delta20 in [-40,-15] deg; "
                  "node rms <= 10 deg",
            "C8": "L in {8,10,12}: rho spread <= 20 MeV; M in {45,50,60}: spread "
                  "40-70 MeV with M60 < M45 < M50; five delta00(1.0 GeV) in [85,105]",
        },
        "prediction": "under SR-b/c the ref/mid rho lands at 822 +- 40 MeV, "
                      "min eta(P1) >= 0.9, and the S0 spectrum returns to ~950 MeV",
    }


if __name__ == "__main__":
    raise SystemExit(main())
