"""Make the retired Newton package importable as ``smatrix_bootstrap_newton`` for its own tests."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
