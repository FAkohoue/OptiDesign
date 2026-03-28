# Package-level setup helpers
utils::globalVariables("mod")

#' OptiDesign: Optimized Experimental Field Design for Plant Breeding
#'
#' @description
#' `OptiDesign` provides tools for constructing, evaluating, and optimizing
#' experimental field designs for plant breeding and related agricultural
#' applications. The package integrates field structure, treatment grouping,
#' genetic relationship information, spatial modelling, and mixed-model
#' efficiency evaluation directly into the design stage.
#'
#' @details
#' ## Design constructors
#'
#' The package provides two main design construction functions:
#'
#' **`prep_famoptg()`** constructs repeated-check block designs with flexible
#' replication, supporting three design classes: augmented repeated-check
#' designs (checks repeated, all entries unreplicated), partially replicated
#' (p-rep) repeated-check designs (mixture of replicated and unreplicated
#' entries), and RCBD-type repeated-check designs (all entries equally
#' replicated). The p-rep constraint - that no replicated entry appears twice
#' in the same block - is enforced by construction. Supports optional grouping
#' from family labels or relationship matrices and optional genetic dispersion
#' optimisation. Construction only; evaluation is performed separately by
#' `evaluate_famoptg_efficiency()`.
#'
#' **`alpha_rc_stream()`** constructs fixed-grid alpha row-column stream
#' designs where the field is treated as a global traversal stream. Replicates
#' are contiguous stream segments, incomplete blocks may be unequal in size,
#' checks are included in every block, and unused cells appear only at the end
#' of the stream. Construction only; evaluation is performed separately by
#' `evaluate_alpha_efficiency()`.
#'
#' ## Design evaluation
#'
#' **`evaluate_famoptg_efficiency()`** evaluates the statistical efficiency of
#' a design produced by `prep_famoptg()`. The mixed model contains Block + Row
#' + Column random effects (no replicate or incomplete-block nesting). Variance
#' component `sigma_b2` controls block variance. Supported criteria:
#'
#' - **A-criterion**: mean pairwise contrast variance (fixed) or mean PEV
#'   (random). Lower is better.
#' - **D-criterion**: geometric mean of contrast covariance eigenvalues
#'   (fixed only). Lower is better.
#' - **CDmean**: mean coefficient of determination for GEBV prediction
#'   (Rincent et al. 2012). Higher is better. Requires random treatment
#'   effects and a genomic prediction model.
#'
#' **`evaluate_alpha_efficiency()`** evaluates the statistical efficiency of
#' a design produced by `alpha_rc_stream()`. The mixed model contains
#' Rep + IBlock(Rep) + Row + Column random effects. Variance components
#' `sigma_rep2` and `sigma_ib2` control replicate and incomplete-block
#' variance respectively. Supports the same A, D, and CDmean criteria as
#' `evaluate_famoptg_efficiency()`.
#'
#' Both evaluation functions are fully decoupled from construction so the same
#' field book can be evaluated multiple times under different model assumptions
#' without rebuilding the layout.
#'
#' ## Design optimisation
#'
#' **`optimize_famoptg()`** wraps `prep_famoptg()` and
#' `evaluate_famoptg_efficiency()` in a **Random Restart (RS)** optimisation
#' loop. RS is used exclusively because the p-rep constraint (no replicated
#' treatment in the same block twice) is enforced by construction at every
#' `prep_famoptg()` call - permutation-based methods such as SA and GA would
#' require block-aware swap logic to preserve this constraint. Every candidate
#' is a valid design with all constraints satisfied. Supports A, D, both, and
#' CDmean criteria.
#'
#' **`optimize_alpha_rc()`** wraps `alpha_rc_stream()` and
#' `evaluate_alpha_efficiency()` in an optimisation loop with three search
#' strategies:
#'
#' - **RS** (Random Restart): generate many independent designs, return the
#'   best by criterion. Simple and fully safe by construction.
#' - **SA** (Simulated Annealing): iteratively propose entry permutation swaps,
#'   accepting improvements always and degradations with a temperature-governed
#'   probability. Repeated across multiple restarts.
#' - **GA** (Genetic Algorithm): maintain a population of entry permutations
#'   and evolve it using Order Crossover (OX1), swap mutation, and elitism.
#'
#' All methods in both optimisers preserve structural constraints by
#' construction. Both implement a four-point integrity checking strategy:
#' post-construction validation, re-check before storing as best, final check
#' before return, and emergency fallback if no valid design is found.
#'
#' ## Internal helpers
#'
#' **`alpha_rc_helpers.R`** provides internal functions shared across all six
#' exported functions: sparse incidence matrix construction, AR1 precision
#' matrix generation, the mixed model solver, the Hutchinson stochastic trace
#' estimator, Chebyshev neighbourhood enumeration, and the dispersion scoring
#' function. These functions are prefixed with `.` and are not exported.
#'
#' @section Function overview:
#'
#' | Function | Role | Design family |
#' |---|---|---|
#' | `prep_famoptg()` | Construction | Repeated-check block |
#' | `evaluate_famoptg_efficiency()` | Evaluation | Repeated-check block |
#' | `optimize_famoptg()` | Optimisation (RS) | Repeated-check block |
#' | `alpha_rc_stream()` | Construction | Alpha row-column stream |
#' | `evaluate_alpha_efficiency()` | Evaluation | Alpha row-column stream |
#' | `optimize_alpha_rc()` | Optimisation (RS/SA/GA) | Alpha row-column stream |
#'
#' @section Typical workflow - repeated-check block design:
#'
#' ```r
#' # 1. Construct
#' design <- prep_famoptg(
#'   check_treatments        = checks,
#'   check_families          = check_fam,
#'   p_rep_treatments        = prep_trts,
#'   p_rep_reps              = rep(2L, length(prep_trts)),
#'   p_rep_families          = prep_fam,
#'   unreplicated_treatments = unrep_trts,
#'   unreplicated_families   = unrep_fam,
#'   n_blocks = 5, n_rows = 15, n_cols = 20
#' )
#'
#' # 2. Evaluate
#' eff <- evaluate_famoptg_efficiency(
#'   field_book         = design$field_book,
#'   n_rows             = 15,
#'   n_cols             = 20,
#'   check_treatments   = checks,
#'   treatment_effect   = "fixed",
#'   residual_structure = "AR1xAR1",
#'   rho_row            = 0.10,
#'   rho_col            = 0.10
#' )
#'
#' # 3. Or optimise directly
#' opt <- optimize_famoptg(
#'   check_treatments        = checks,
#'   check_families          = check_fam,
#'   p_rep_treatments        = prep_trts,
#'   p_rep_reps              = rep(2L, length(prep_trts)),
#'   p_rep_families          = prep_fam,
#'   unreplicated_treatments = unrep_trts,
#'   unreplicated_families   = unrep_fam,
#'   n_blocks           = 5,
#'   n_rows             = 15,
#'   n_cols             = 20,
#'   treatment_effect   = "fixed",
#'   residual_structure = "AR1xAR1",
#'   rho_row            = 0.10,
#'   rho_col            = 0.10,
#'   criterion          = "A",
#'   n_restarts         = 50
#' )
#' ```
#'
#' @section Typical workflow - alpha row-column stream design:
#'
#' ```r
#' # 1. Construct
#' design <- alpha_rc_stream(
#'   check_treatments = checks,
#'   check_families   = check_fam,
#'   entry_treatments = entries,
#'   entry_families   = entry_fam,
#'   n_reps           = 3,
#'   n_rows           = 30,
#'   n_cols           = 20,
#'   min_block_size   = 19,
#'   max_block_size   = 20
#' )
#'
#' # 2. Evaluate
#' eff <- evaluate_alpha_efficiency(
#'   field_book         = design$field_book,
#'   n_rows             = 30,
#'   n_cols             = 20,
#'   check_treatments   = checks,
#'   treatment_effect   = "fixed",
#'   residual_structure = "AR1xAR1",
#'   rho_row            = 0.10,
#'   rho_col            = 0.10
#' )
#'
#' # 3. Or optimise directly
#' opt <- optimize_alpha_rc(
#'   check_treatments   = checks,
#'   check_families     = check_fam,
#'   entry_treatments   = entries,
#'   entry_families     = entry_fam,
#'   n_reps             = 3,
#'   n_rows             = 30,
#'   n_cols             = 20,
#'   min_block_size     = 19,
#'   max_block_size     = 20,
#'   treatment_effect   = "fixed",
#'   residual_structure = "AR1xAR1",
#'   rho_row            = 0.10,
#'   rho_col            = 0.10,
#'   method             = "RS",
#'   criterion          = "A",
#'   n_restarts         = 50
#' )
#' ```
#'
#' @section Key differences between the two design families:
#'
#' | Feature | `prep_famoptg` family | `alpha_rc_stream` family |
#' |---|---|---|
#' | Blocking structure | Flat blocks | Replicates -> incomplete blocks |
#' | Replication | Flexible per-entry | Uniform across entries |
#' | Design types | Augmented, p-rep, RCBD-type | Alpha-lattice |
#' | Block variance component | `sigma_b2` | `sigma_rep2` + `sigma_ib2` |
#' | Optimisation methods | RS only | RS, SA, GA |
#' | P-rep constraint enforced | Yes | Not applicable |
#'
#' @section Relationship matrices and grouping:
#'
#' Both design families accept one of three grouping sources for adjacency
#' scoring and clustering within blocks:
#'
#' - Direct family labels (`cluster_source = "Family"`)
#' - A genomic relationship matrix (`cluster_source = "GRM"`, argument `GRM`)
#' - A pedigree relationship matrix (`cluster_source = "A"`, argument `A`)
#'
#' A separate relationship matrix `K` can be provided independently for:
#'
#' - GBLUP/PBLUP prediction efficiency (`prediction_type = "GBLUP"`)
#' - CDmean computation for genomic selection training population optimisation
#' - Spatial dispersion optimisation (`dispersion_source = "K"`)
#'
#' When treatment names do not match matrix rownames, mapping tables `id_map`
#' and `line_id_map` reconcile identifiers.
#'
#' @section Included example dataset:
#'
#' The package ships with `OptiDesign_example_data`, containing synthetic
#' objects for treatment and family lists, example relationship matrices, and
#' example argument lists for both design families. Load with:
#'
#' ```r
#' data("OptiDesign_example_data", package = "OptiDesign")
#' ```
#'
#' @section Dependencies:
#'
#' - **Matrix**: sparse matrix operations used throughout efficiency evaluation,
#'   AR1 precision matrix construction, and the mixed model coefficient matrix.
#'   Required by all six exported functions.
#' - **pracma**: `mod()` is imported for serpentine parity handling in
#'   `prep_famoptg()` (alternating row/column traversal direction).
#'
#' @section References:
#'
#' Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla, P.,
#' ... & Moreau, L. (2012). Maximizing the reliability of genomic selection by
#' optimizing the calibration set of reference individuals: comparison of
#' methods in two diverse groups of maize inbreds. *Genetics*, 192(2), 715-728.
#'
#' Jones, B., Allen-Moyer, K., & Goos, P. (2021). A-optimal versus D-optimal
#' design of screening experiments. *Journal of Quality Technology*, 53(4),
#' 369-382.
#'
#' Hutchinson, M.F. (1990). A stochastic estimator of the trace of the
#' influence matrix for Laplacian smoothing splines. *Communications in
#' Statistics - Simulation and Computation*, 19(2), 433-450.
#'
#' Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983). Optimization by
#' simulated annealing. *Science*, 220(4598), 671-680.
#'
#' Holland, J. H. (1992). *Adaptation in Natural and Artificial Systems*.
#' MIT Press.
#'
#' @keywords internal
"_PACKAGE"
