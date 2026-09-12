"""SDP-route reproduction of He-Kruczenski arXiv:2309.12402v3.

Self-contained re-derivation of every operator from the paper's equations
(TeX source ``references/2309.12402v3-source/prd_submission_2.tex``).  Nothing
in this subpackage imports the repository's historical Newton-mainline modules
(``kernels``, ``operators``, ``model``, ``linear``, ``scattering``, ...); those
are only ever loaded by :mod:`smatrix_bootstrap.sdp.crosscheck` as an
*independent* comparison target.

Equation numbers quoted in docstrings refer to the published paper
(arXiv:2309.12402v3), e.g. (2.7)=Adef, (2.9)=fdef, (3.58)=zmap, (3.67)=h37.
"""

__all__ = ["grid", "legendreq", "hilbert", "projector", "formfactor",
           "constraints", "problem", "observables"]
