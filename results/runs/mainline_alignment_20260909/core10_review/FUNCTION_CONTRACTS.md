# Helper contracts

Short helper documentation moved from code during consolidation; all executable expressions remain in the ten modules.

- `figures._clip`: Floating display polygon; authoritative evidence remains its support lines.
- `figures._interval_verdict`: Strict interval sign; an unresolved comparison is not a negative result.
- `kernels.midpoint_grid`: Exact-source Arb nodes and w=Delta_phi*x'/pi=x'/M, center nu0=0.
- `kernels.pv_matrix`: Integrated subtracted PV map: sin(nphi)->cos(nphi), including Nyquist.
- `kernels.exterior_legendre_q`: Real integral-defined Q outside[-1,1], explicitly continued to z<-1.
- `kernels.assemble_density_row`: rho2 upper entries are actual symmetric density values, not half-basis coefficients.
- `kernels.offcut_row`: PV row below threshold, separate from the sine amplitude.
- `kernels.physical_node_row`: PV native nodal jump and rational crossed projection.
- `kernels.density_row_to_cflat`: Exact dual-coordinate conversion; no coefficient or row is discarded.
- `linear.symmetric`: Inverse of the orthonormal upper-triangle svec convention.
- `linear.inverse_cholesky`: Batched 3x3 Cholesky in longdouble, avoiding a float64 eigensolver.
- `linear.restore_joint_segment`: Locate a complete joint feasible segment before one independent rounded audit.
- `linear.radial_center`: Minimize the same barrier on the feasible positive ray through z.
- `model.chiral_slices`: Explicit groups of the same eight physical chiral residuals.
- `model.weinberg_waves`: Eq.(2.18), columns S0, S2, P1, in pion-mass units.
- `model.phase_shifts`: Threshold-anchored delta=unwrap(arg S)/2; no elastic saturation imposed.
- `model.ff_guided_phase`: A conditional mod-pi lift of the same S; native samples do not prove continuity.
- `model.phase_space`: Physical beta, including its exact threshold value; not a sub-cut test.
- `model.scattering_matrix`: v3 Eq.(2.11); at s=4 the value S=1 assumes a finite partial wave f.
- `model.unitarity_margin`: Stable 1-|S|^2: nonnegative iff the elastic-channel S is contractive.
- `model.scattering_gram`: v3 Eq.(3.71): eigenvalues 1 +/- |S|; PSD iff |S| <= 1.
- `model.current_gram`: The exact positive-k congruence of the source's current Gram matrix.
- `model.fesr_targets`: Four J_n=s0^(-n-2) integral rho(x)x^n dx, or the raw integrals.
- `model.fesr_weights`: Four nodal FESR rows; clipped-phi retains partial cell overlap, not exact integrals.
- `model.fesr_data`: JSON-ready A2 contract; D2 can impose four independent scalar errors.
- `model.density_fourth_power`: Ordinary norm of rho1 plus the upper triangle of actual rho2.
- `model.chiral_bound`: One mode/tolerance contract shared by construction and independent repair.
- `model.convert_subtraction`: Exact H=q+b coordinate change in C_flat; the two densities are unchanged.
- `gauge.joint_observable_change`: Native S/F changes between complete joint points, independent of projection gap.
- `gauge.watsonian_objective`: 2403 Eq.(2.29): frozen form-factor phases, author Lambda^-2 weights.
- `gauge.advance`: Exact quadratic slack updates, retaining relative precision near active cones.
- `selection.watson_lineage`: Track the fixed initial role and the author-defined change of amplitude.
- `selection.plot_continuation`: Show coordinate motion and physical saturation without moving phase branches.
- `operators.raw_support_residual`: Remove unbounded raw sigma/T0 residuals using the original nodal rows.
- `operators.energy_rows`: Independent energy worker; native nodes keep the exact nodal jump identity.
- `operators.coefficient_blocks`: Unpack full C_flat coefficients into the two actual density matrices.
- `resolution.cardinal_projection`: Sine coefficients truncated/padded, not a PV amplitude identity.
- `resolution.project_coefficients`: Galerkin densities in subtracted coordinates; only a numerical initial guess.
- `resolution.resolution_report`: Selection binds the actual new model; old feasibility/role labels never transfer.
- `resolution.transfer`: Rebuild each finite problem and re-audit a numerical transfer, including dual bounds.
- `resolution.original_duals`: Covectors for the original zero-objective problem, not the relaxed one.
- `__init__.solve_feasible_basis_lp`: Public HiGHS adapter for min cost*x, matrix*x=target, x>=0 with a feasible basis.
- `__init__.extended_triangular`: Back substitution retains longdouble; scipy's triangular solver casts it.
- `__init__.extended_sparse_blocks`: Cache exact nonzeros and transpose views in eight blocks balanced by work.
- `__init__.extended_product`: Same extended-precision products, partitioned by independent row blocks.
- `__init__.support_gradient_width`: Bound the gradient-error functional's width using disks and the actual L4 ball.
- `__init__.mq_squared`: Effective m_q^2 in pion units, also used by the scalar FF cap.
- `quotient.subtraction_rows`: Rewrite full actual-density rows in subtracted coordinates; never fix T0s.
- `quotient.subtract_amplitude_coordinates`: Invertible full actual-density coordinate change, returned as longdouble midpoints.
- `quotient.absorptive_change`: Full triangular coordinates (T0_sub, native Im f00, Im f20, actual rho).
- `run.support_data`: Load scattering rows.
- `io.decode_real_ball`: Decode authoritative midpoint-radius enclosures.
- `io.save_evaluation`: Attach current observables to a fresh, selected primary-wave evaluation.
- `io.resolve_start_mu`: Continue a saved fixed objective at its checkpoint mu; explicit values take precedence.
- `io.published_inner_vertices`: Feasible mixtures in x>=0; zero is the known zero scattering amplitude.
- `io.joint_result`: One original-variable serialization and acceptance rule shared by D, E and F.
- `io.save_figure`: Identical PDF/PNG products for every figure command.

Additional contracts duplicated by saved report fields:

- `certify_reference_section`: Arb section bounds from certified feasible segments and global support planes.
- `barrier_support`: Feasible inexact Newton path; float QR preconditions longdouble CG.
- `center_joint_support`: Feasible Newton line minimization with fully centered half-mu continuation.
- `current_operators`: D1 affine current blocks in original C_flat + ImF0,ImF1,rho0,rho1.
- `watson_diagnostics`: Physical S/FF consistency and two-pion spectral fraction on evaluated nodes.
- `prepare_amplitude`: A1–A3 amplitude and UV operators.
- `initialize`: Full current Phase I with every original scattering direction and constraint.
- `recover_candidate_segment`: One complete coefficient segment, followed by the independent original-row audit.
- `joint_hull_candidate`: Construct a joint witness in saved feasible amplitudes' convex hull, not a support.
- `UVConfig`: v3 inputs in GeV, with explicit conditional choices for missing numerics.
- `normal_dual`: Select fixed normals and nodal disk axes; independent outer checks remain required.
- `support_outer`: Original-row support and primal certificate for the complete PV finite problem.
- `joint_audit`: Independent finite H/K/k current audit; only unconstrained high spectra are restored.
- `gauge_support`: Joint PV scattering/current SDP, with original-variable primal and dual audit.
- `gauge_regions`: Compare UV supports with the same-model chiral region; no hull-only completion.
- `gauge_phases`: Use the already frozen three amplitudes for Fig.9 P1 and Fig.10 S0/S2.
- `select_gauge`: Freeze three audited joint boundary approximations using only geometry.
- `compare_phase`: Descriptive native-sample errors against an extracted PDF curve; no fit.
- `select_chiral_amplitude`: A feasible upper-boundary interpolation at a fixed IR coupling, before phases.
- `compare_selected_profiles`: Compare saved curves without modifying amplitudes, parameters or phase branches.
- `resolution_compare`: Fig.11 native curves; interpolation is only a display/comparison diagnostic.

Additional contracts duplicated by saved report fields:



Additional contracts duplicated by saved report fields:



Additional contracts duplicated by saved report fields:



Additional contracts duplicated by saved report fields:



Additional contracts duplicated by saved report fields:



Additional contracts duplicated by saved report fields:



Additional contracts duplicated by saved report fields:


