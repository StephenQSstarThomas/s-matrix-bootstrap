# Superseded numerical diagnostics

These files retain the initial investigation, including failures. They are not
accepted reference results. The sequence was ordinary angular quadrature and
unitarity-disk cones; logarithmic angular coordinates and invertible SVD
preconditioning; scaled rotated cones and the constant endpoint condition;
then analytic threshold projection. Earlier runs can have spurious tiny
high-spin contributions from angular cancellation, inaccurate optimizer
termination, and substantial independent unitarity violations.

`pure_threshold_N4.json` and `pure_threshold_N6.json` here precede the addition
of solver-gap and asymptotic diagnostics; their current counterparts are at
the parent level. Original output hashes document the source state then used,
but these archived development stages are not all exposed as selectable modes
of the current solver. The maintained, rerunnable calculations and their exact
commands are specified in `../../docs/04_numerical_investigation.md`.
