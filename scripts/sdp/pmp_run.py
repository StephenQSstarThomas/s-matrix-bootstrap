"""Compatibility wrapper for the self-contained SDPB CLI. No other solver is loaded."""
from smatrix_bootstrap.sdp.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
