"""Retired pre-SDPB report API; original producer bytes remain authenticated.

The old renderer hard-coded obsolete solver, B-regulator and rank conclusions.
It is deliberately unavailable for new scientific reports. Historical sources:
results/runs/sdpb_mainline_20260912/legacy_report_retirement_20260913/.
"""


def legacy_report_unavailable(*args, **kwargs):
    raise RuntimeError(
        "The pre-SDPB report generator is retired because it hard-coded obsolete "
        "physical and solver conclusions. Use python -m smatrix_bootstrap.run sdp "
        "figures or contrast with explicit completed scientific inputs; see "
        "HANDOFF_PHYSICS_AUDIT_ZH.md. No report has been written."
    )


collect = render = build = legacy_report_unavailable
