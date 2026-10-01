# Provenance of the authors' code (third_party/gauge-theory-bootstrap)

- Repository: https://github.com/hyfysics/gauge-theory-bootstrap (Yifei He, Martin Kruczenski)
- Pinned commit: `801684d9ece3de2918a20145178f569459081098` (committed 2025-10-08T11:01:27+02:00, branch main)
- Contents used: the arXiv:2403.10772 snapshot (`theories/qcd/papers/arxiv-2403.10772/`, identical in substance to the
  v3/v4 arXiv ancillary files and to GitHub commit c840933) and the arXiv:2505.19332 snapshot.
- Licence: the repository has no licence file (checked 2026-10-01 through the GitHub API: `license: null`). It is
  therefore included here **by reference** as a git submodule, not copied or redistributed. Nothing under this
  directory is executed by our pipeline except `mma/gtb2403/authors_reference.wls`, which evaluates the 2403 notebook
  unchanged (only its `SetDirectory` and `Export` cells are skipped) to produce a reference export for comparison.
- Role: reference only. Our kernels are derived and implemented independently in `mma/gtb2403/GTBKernels.wl`
  (derivation in `docs/gtb2403/DERIVATION_ZH.md`); the solver is SDPB, not the MATLAB/CVX/MOSEK programs below.

| file (relative to the submodule) | sha256 |
|---|---|
| `theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb` | 9a3a81d7e99ef35355aa5a2bbf10ce9a21e6af929e0cb8371c4afaff133922bd |
| `theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m` | b12bf411407fc4cb8ae56aea41c88aa336e0f73aa7ae44ebef2b1bca3fc08ad2 |
| `theories/qcd/papers/arxiv-2505.19332/src/mathematica/gtb_qcd_01_generate_matrices.nb` | 1a3a018683797cd85ee67e4fe8e37d7b44d89e5776cd5a95c01ce63ccc8f97e8 |
| `theories/qcd/papers/arxiv-2505.19332/src/matlab/gtb_qcd_02_optimize.m` | dcb0ee7c2c555393e8e4bcd424fe0fb7aac82dfca3432810b890ae4a7c4b588c |
| `theories/qcd/papers/arxiv-2505.19332/src/matlab/optimize_core.m` | 0e7c883cf7d79c8fd118d9875640414cac4fd110a47ee5da06267348a19a78fb |
| `theories/qcd/papers/arxiv-2505.19332/src/mathematica/gtb_qcd_03_analyze_plot.nb` | f4fe8ba1010c2115e0213bd206a7156676b26138bc8cfbf0b1c7ac2667555ec5 |
