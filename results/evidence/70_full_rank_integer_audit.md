# Independent integer audit of the full analytic rank certificate

An exact integer calculation independently confirms the decisive 1510 by
1510 residual bound in [chapter 69](69_full_analytic_rank_and_finite_target_plane.md).
It uses the **same original interval minor and the same dyadic inverse
candidate**, with no new inverse, selection, physical masters or Arb matrix
products. The [audit index](../remaining_rank_exact_integer_audit.json)
retains every input hash and links the complete integer/rational certificate.

The resulting exact rational bound is approximately

\[
\boxed{\delta_{\mathbb Z}
\simeq2.2069515031620150242231683812890449\times10^{-68}<1.}
\]

Every one of the 1510 exact integer row bounds is below its previously
saved Arb bound. Together with the established operator interpretation,
the same rank chain follows: selected minor rank1510, complete quotient
rank1511, and **analytic augmented rank3011**.

## Exact integer computation

The authenticated midpoint/radius encodings give the exact lifts

\[
R_{Q,0}=2^{-1019}R,\qquad
\operatorname{mid}B_Q=2^{-2072}U,\qquad
\operatorname{rad}B_Q=2^{-861}V,
\]

where R and U are signed integer matrices and V is nonnegative. The
candidate radii are all exactly zero. Maximum integer widths are 1220,
2072 and 454 bits respectively.

One native integer product computes RU. For each output row i define

\[
b_i=\sum_{j=1}^{1510}
\left|2^{3091}\delta_{ij}-(RU)_{ij}\right|
+2^{1211}\sum_{k=1}^{1510}|R_{ik}|
\sum_{j=1}^{1510}V_{kj}.
\]

Then, uniformly for every point matrix in the original B enclosure,

\[
\|I-R_{Q,0}B_Q\|_\infty
\le2^{-3091}\max_i b_i=\delta_{\mathbb Z}.
\]

The radius contribution uses one nonnegative integer matrix-vector
product after summing the radius rows. All additions, shifts, products,
absolute values and comparisons are exact. No precision parameter enters
this residual calculation.

The certificate saves all b_i in hexadecimal and rational form, together
with their separate midpoint-residual and uncertainty contributions. Root
independently rechecked all 1510 decompositions, integer-to-rational
normalizations, the strict maximum below one, and comparisons with the
earlier Arb bounds. This final check does not claim a third large product.

## Results and resource record

| Quantity | Result |
|---|---:|
| Entries in each original matrix | 2,280,100 |
| Integer square products / radius matrix-vector products | 1 / 1 |
| Maximum midpoint-only residual row bound | about 1.2387885717e-237 |
| Maximum uncertainty contribution | about 2.2069515032e-68 |
| Exact total residual upper bound | about 2.2069515032e-68 |
| Integer square-product time | 134.65 seconds |
| Total successful supervisor time | 182.69 seconds |
| Predeclared whole-worker cap | 1800 seconds |
| Peak resident memory | 8.70 GiB |
| Timeouts / retries / input or source changes | none |

The uncertainty term dominates; the inverse candidate's midpoint residual
is already much smaller. Nevertheless the total bound is far below one,
so no first-node refinement, Neumann correction or new candidate is needed
for the full-rank conclusion. The independently materialized
[first-node interval intersection](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/67_first_source_interval_intersection.md`)
was not used in either successful full-rank certificate.

The worst integer residual row is zero-based index1185 and corresponds to
selected quotient coordinate1642. These are **selected-column domain**
indices for the left residual \(I-R_{Q,0}B_Q\), not labels of physical
source partial waves.

The prior [preflight](../remaining_rank_integer_preflight.json)
authenticated both matrices and all row/column/scaling bindings, and
measured integer widths and memory before any large product. It ran no
matrix multiplication. The later audit retained its frozen preflight,
used one declared attempt, and terminated successfully.

## Reproduction and scope

The full integer row certificate is stored at
`results/analytic_source_storage/source_remaining_integer_audit_attempt01/integer_row_certificate.json`.
Its 5,741,573 bytes include the exact row arrays, dyadic plan and authenticated
input ledger. The run embeds its driver and helper sources.

The reproducible driver is
[audit_remaining_rank_integer.py](../runs/repo_reorganization_20260906/recovery.json) (historical member: `scripts/audit_remaining_rank_integer.py`),
using the already independently tested
[dyadic_inverse_residual.py](../runs/repo_reorganization_20260906/recovery.json) (historical member: `src/smatrix_bootstrap/dyadic_inverse_residual.py`).
An existing output directory is rejected:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
PYTHONPATH=/tmp/collocation_arb:src python scripts/audit_remaining_rank_integer.py \
  --directory NEW_IMMUTABLE_AUDIT_DIRECTORY
```

This strengthens the arithmetic audit of the existing analytic rank
theorem. It adds no model restriction or bootstrap optimization. The
finite target-plane consequence and the remaining published-region,
phase and continuum questions retain exactly the scope of chapter69.
