#!/usr/bin/env python
"""Retired pre-SDPB report command; use the single scientific CLI.

Authenticated original producer bytes and the reason for retirement are kept in
results/runs/sdpb_mainline_20260912/legacy_report_retirement_20260913/.
"""
from smatrix_bootstrap.sdp.report import legacy_report_unavailable


if __name__ == "__main__":
    legacy_report_unavailable()
