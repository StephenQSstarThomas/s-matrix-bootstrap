# Original v3 settings needed to resolve the remaining comparison

This is a local draft; no message has been sent.

We are reproducing Fig.3–4 of arXiv:2309.12402v3. We use the stated M=50 grid and ten allowed partial waves per isospin, retain all 3876 scattering coefficients, and use f=(1/4) integral P_l T and S=1+i pi sqrt(1-4/s) f. The current explicit comparison uses the later authors' Legendre-Q/PV-midpoint kernels, a free constant, and only the 1500 native scattering disks. It does not add our separate finite-sine tail or high-spin conditions.

The four original numerical details that would resolve the ambiguity are:

1. The amplitude-density regularizer used for Fig.3 and Fig.4: norm, variable packing, any normalization by M or the number of density coefficients, bound, and whether the two figures used the same value.
2. The precise norm and grouping in Eq.(3.64), including any normalization of the residual vectors. We currently implement two separate four-component Euclidean balls; later public code uses a combined vector in a different setup.
3. The original scattering matrices or interpolation/subtraction prescription used in v3, including treatment of the free constant. The later public implementation is useful evidence but does not establish identity of the original figure settings.
4. Fig.4's boundary-generation inputs or numerical point/coefficient files, including its scan grid and solution tolerances. The displayed green branches are nearly uniform in x; the general radial example alone does not specify the actual scan.

For context, with the later-code interpretation B=377500 and our stated chiral norm, independently bounded +x supports are approximately [.10588510,.10594018] at epsilon=.002 and [.05446749,.05452032] at epsilon=.0002. The displayed v3 rightmost samples are approximately .082573 and .022363. We do not treat these plotted sample maxima as rigorous upper bounds on the authors' regions, or infer an error in the paper from a comparison with incompletely recovered settings. Original code, matrices, or the relevant parameter lines would be sufficient to settle the differences.
