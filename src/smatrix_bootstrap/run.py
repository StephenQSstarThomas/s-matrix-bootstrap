"""Calculation entry: ``python -m smatrix_bootstrap.run sdp <subcommand> ...``.

Only the SDPB route is dispatched here.  The retired Newton mainline (``prepare``, ``boundary``, ``select``,
``evaluate``, ...) now lives in ``legacy/smatrix_bootstrap_newton`` and is run as
``PYTHONPATH=legacy python -m smatrix_bootstrap_newton.run ...``.
"""
import sys


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["sdp"]:
        from .sdp.__main__ import main as sdp_main
        return sdp_main(argv[1:])
    sys.stderr.write("usage: python -m smatrix_bootstrap.run sdp <subcommand> ...\n"
                     "The Newton mainline was retired to legacy/smatrix_bootstrap_newton (see legacy/README.md).\n")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
