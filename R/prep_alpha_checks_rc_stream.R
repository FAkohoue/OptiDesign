#' Create a fixed-grid alpha row-column design with unequal incomplete blocks and checks in every block
#'
#' Construct a fixed-field row-column design in which the full `n_rows × n_cols`
#' grid is first converted into a single ordered stream of positions, then split
#' into contiguous replicate segments, then subdivided into incomplete blocks of
#' possibly unequal size. Within each replicate, every entry appears exactly once,
#' every incomplete block contains all checks, and any unused field cells are left
#' as trailing `NA` positions at the end of the global stream.
#'
#' @description
#' `prep_alpha_checks_rc_stream()` is intended for practical field situations where:
#'
#' - the overall field size is fixed in advance,
#' - replicate boundaries are determined by field-book order rather than by rigid
#'   rectangular subfields,
#' - replicate 2 begins exactly where replicate 1 ends in the traversal stream,
#' - incomplete blocks may differ in size,
#' - checks must be repeated in every incomplete block,
#' - all entries must appear exactly once in each replicate,
#' - any leftover field cells should remain unused and appear at the end of the stream.
#'
#' The function:
#'
#' 1. Builds the full field stream using `order` and `serpentine`.
#' 2. Determines how many incomplete blocks can be supported in each replicate.
#' 3. Divides each replicate into incomplete blocks.
#' 4. Allocates checks and entries subject to the block plan.
#' 5. Arranges entries heuristically to reduce clustering of similar materials.
#' 6. Leaves any surplus cells as trailing `NA`.
#' 7. Optionally performs dispersion optimization.
#' 8. Optionally computes mixed-model efficiency diagnostics.
#'
#' @section Conceptual design:
#' The key feature of this function is that it is **stream-based** rather than
#' **rectangle-based**.
#'
#' Instead of cutting the field into rectangular replicate areas, the function:
#'
#' - lists all grid cells in a single global order,
#' - cuts that order into `n_reps` contiguous replicate segments,
#' - then cuts each replicate segment into incomplete blocks.
#'
#' This is useful when:
#'
#' - planting follows a continuous field-book order,
#' - operational field movement is row-wise or column-wise,
#' - replicate boundaries are administrative or logistical rather than geometric,
#' - the user wants all unused cells to appear only at the tail of the field stream.
#'
#' For each replicate:
#'
#' - each entry is used exactly once,
#' - each incomplete block receives every check,
#' - the remaining block capacity is filled with entries,
#' - allocation is then improved heuristically to reduce local grouping conflicts.
#'
#' @section What this function is most useful for:
#' This function is particularly useful when:
#'
#' - the field geometry is fixed and should never be auto-resized,
#' - the design requires repeated checks in every incomplete block,
#' - the user wants a row-column design without forcing rectangular replicate areas,
#' - practical planting or harvesting follows a stream order,
#' - family, pedigree, or genomic structure should influence entry placement,
#' - the user wants to evaluate the resulting design using a mixed-model framework.
#'
#' @section Dependency guide:
#' Many arguments are active only under particular modes.
#'
#' **Grouping source**
#'
#' - If `cluster_source = "Family"`:
#'   - grouping comes directly from `check_families` and `entry_families`,
#'   - `GRM`, `A`, `id_map`, `cluster_method`, `cluster_seed`,
#'     `cluster_attempts`, and `n_pcs_use` are ignored.
#'
#' - If `cluster_source = "GRM"`:
#'   - `GRM` is required,
#'   - `id_map` is only needed if entry names differ from `rownames(GRM)`,
#'   - grouping labels are derived by PCA followed by clustering.
#'
#' - If `cluster_source = "A"`:
#'   - `A` is required,
#'   - `id_map` is only needed if entry names differ from `rownames(A)`,
#'   - grouping labels are derived by PCA followed by clustering.
#'
#' **Derived block structure**
#'
#' - `min_entry_slots_per_block` controls the minimum number of non-check entry
#'   slots permitted in an incomplete block.
#' - `max_blocks_per_rep` optionally limits the number of blocks per replicate.
#' - If `max_blocks_per_rep = NULL`, block count is derived from capacity and
#'   `min_entry_slots_per_block`.
#'
#' **Efficiency evaluation**
#'
#' - If `eval_efficiency = FALSE`, all efficiency arguments are ignored.
#'
#' - If `eval_efficiency = TRUE` and `treatment_effect = "fixed"`:
#'   - fixed-treatment contrast precision is computed,
#'   - `prediction_type`, `K`, and `line_id_map` are ignored.
#'
#' - If `eval_efficiency = TRUE` and `treatment_effect = "random"`:
#'   - `prediction_type` becomes active.
#'
#' - If `prediction_type = "none"`:
#'   - random-effect efficiency is skipped.
#'
#' - If `prediction_type = "IID"`:
#'   - entry random effects are treated as independent,
#'   - `K` and `line_id_map` are ignored.
#'
#' - If `prediction_type %in% c("GBLUP", "PBLUP")`:
#'   - `K` is required,
#'   - `line_id_map` may be required if treatment IDs differ from `rownames(K)`.
#'
#' **Residual structure**
#'
#' - If `residual_structure = "IID"`:
#'   - `rho_row` and `rho_col` are ignored.
#'
#' - If `residual_structure = "AR1"`:
#'   - only `rho_row` is used.
#'
#' - If `residual_structure = "AR1xAR1"`:
#'   - both `rho_row` and `rho_col` are used.
#'
#' **Dispersion optimization**
#'
#' - If `use_dispersion = FALSE`, all dispersion arguments are ignored.
#'
#' - If `use_dispersion = TRUE`, `dispersion_source` selects which matrix is used.
#'
#' - If `dispersion_source = "K"`:
#'   - `K` is required,
#'   - `line_id_map` may be needed if treatment IDs differ from matrix row names.
#'
#' - If `dispersion_source = "A"`:
#'   - `A` is required.
#'
#' - If `dispersion_source = "GRM"`:
#'   - `GRM` is required.
#'
#' **Check placement**
#'
#' - `check_placement = "systematic"` spreads checks approximately evenly within
#'   each block stream.
#' - `check_placement = "random"` randomizes check positions.
#' - `check_position_pattern` is retained for interface continuity but is not
#'   used by the current stream-based placement logic.
#'
#' @details
#' Let:
#'
#' - `E` = number of entries,
#' - `C` = number of checks,
#' - `R` = number of replicates,
#' - `b` = number of incomplete blocks per replicate.
#'
#' Then the number of **used plots per replicate** is:
#'
#' `E + b * C`
#'
#' and the **total number of used plots** is:
#'
#' `R * (E + b * C)`
#'
#' This total must fit within the fixed field capacity:
#'
#' `n_rows * n_cols`
#'
#' The function chooses a block plan that respects:
#'
#' - the number of entries,
#' - the number of checks,
#' - the minimum entry capacity per block,
#' - any user-supplied cap on blocks per replicate,
#' - the total fixed field capacity.
#'
#' In efficiency evaluation, the function uses only the observed non-`NA` plots.
#' Depending on `treatment_effect`, the returned metrics represent either:
#'
#' - precision of fixed treatment contrasts, or
#' - average prediction error variance of random treatment effects.
#'
#' The arguments `warn_and_correct`, `fix_rows`, `spatial_engine`,
#' `dense_max_n`, and `check_position_pattern` are retained for interface
#' continuity, even where they do not materially alter the stream-based design logic.
#'
#' @param check_treatments Character vector of check IDs.
#'
#' Every check is placed exactly once in every incomplete block.
#'
#' This argument is always active and is central to the design because repeated
#' checks define the benchmark structure across blocks and replicates.
#'
#' Use this when the trial requires standard controls or benchmark entries in
#' every incomplete block.
#'
#' Example use case:
#' a trial with 3 commercial checks that must appear in all incomplete blocks.
#'
#' @param check_families Character vector of the same length as `check_treatments`.
#'
#' Family labels for checks.
#'
#' These are used directly when `cluster_source = "Family"` and are also stored in
#' the output field book under all settings.
#'
#' This argument must align exactly with `check_treatments`.
#'
#' Use it when checks belong to known families or control groups, or when the user
#' wants family labels preserved in the output.
#'
#' @param entry_treatments Character vector of distinct entries.
#'
#' Each entry appears exactly once **per replicate**.
#'
#' This is the main set of test entries being evaluated.
#'
#' Use this for breeding lines, hybrids, progeny, or other candidate materials
#' that should be represented once in each replicate.
#'
#' Entries must not overlap with `check_treatments`.
#'
#' @param entry_families Character vector of the same length as `entry_treatments`.
#'
#' Family labels for entries.
#'
#' These are used directly when `cluster_source = "Family"`. When
#' `cluster_source %in% c("GRM", "A")`, they still help determine the target number
#' of clusters among entries.
#'
#' Use this when family structure matters for spacing, grouping, or interpretation.
#'
#' @param n_reps Integer giving the number of replicates.
#'
#' Each replicate is defined as a contiguous segment of the global stream.
#'
#' Use more replicates when stronger replication is needed for entry precision,
#' provided the field can still support repeated checks in every block.
#'
#' This argument affects:
#' - total number of used plots,
#' - the length of each replicate segment,
#' - the total number of times entries appear.
#'
#' @param n_rows Integer giving the number of field rows.
#'
#' The field size is fixed in this function and is never altered.
#'
#' Use this to reflect the true row dimension of the physical field.
#'
#' It directly affects:
#' - total field capacity,
#' - the global traversal stream,
#' - spatial coordinates in the returned `layout_matrix` and `field_book`.
#'
#' @param n_cols Integer giving the number of field columns.
#'
#' The field size is fixed in this function and is never altered.
#'
#' Use this to reflect the true column dimension of the physical field.
#'
#' Together with `n_rows`, it determines field capacity and stream structure.
#'
#' @param order Character specifying the global traversal order.
#'
#' Allowed values:
#' - `"column"`: column-major traversal,
#' - `"row"`: row-major traversal.
#'
#' This determines:
#' - the order of the global field stream,
#' - where replicate boundaries occur,
#' - where incomplete block boundaries occur,
#' - where trailing `NA` cells are placed.
#'
#' Use `"row"` when operational movement follows rows.
#' Use `"column"` when operational movement follows columns.
#'
#' This argument interacts directly with `serpentine`.
#'
#' @param serpentine Logical indicating whether alternate rows or columns reverse
#' direction during stream generation.
#'
#' If `TRUE`, traversal alternates direction:
#' - by row when `order = "row"`,
#' - by column when `order = "column"`.
#'
#' Use this when field-book order should mimic serpentine movement in the field.
#'
#' This changes the stream order, which in turn changes replicate and incomplete
#' block boundaries.
#'
#' Depends on:
#' - `order`
#'
#' @param seed Optional integer seed for reproducibility.
#'
#' If `NULL`, a seed is generated internally and returned as `seed_used`.
#'
#' Use this when the layout should be reproducible.
#'
#' It affects:
#' - entry allocation,
#' - within-block arrangement,
#' - random check placement,
#' - optional dispersion optimization (unless a separate dispersion seed is supplied).
#'
#' @param attempts Integer giving the number of swap proposals used when improving
#' entry allocation across blocks within each replicate.
#'
#' Larger values increase search effort but do not change the design algorithm itself.
#'
#' Use larger values when:
#' - family structure is strongly imbalanced,
#' - the user wants better separation of similar materials,
#' - the number of entries per replicate is large.
#'
#' @param warn_and_correct Logical retained for interface continuity.
#'
#' In this function, the field dimensions are fixed and are **not** altered.
#'
#' This argument is included so the interface stays similar to related design
#' functions, but it does not drive field resizing here.
#'
#' @param fix_rows Logical retained for interface continuity.
#'
#' In this function, the field dimensions are fixed and are **not** altered.
#'
#' Included mainly for consistency with related interfaces.
#'
#' @param cluster_source Character specifying the grouping source used during
#' allocation and arrangement.
#'
#' Allowed values:
#' - `"Family"`
#' - `"GRM"`
#' - `"A"`
#'
#' Use `"Family"` when user-supplied family labels are adequate and interpretable.
#'
#' Use `"GRM"` when genomic structure should define similarity more precisely.
#'
#' Use `"A"` when pedigree relationships are available and should define grouping.
#'
#' This argument determines whether `GRM`, `A`, `id_map`, `cluster_method`,
#' `cluster_seed`, `cluster_attempts`, and `n_pcs_use` become active.
#'
#' @param GRM Optional genomic relationship matrix.
#'
#' Required when:
#' - `cluster_source = "GRM"`, or
#' - `use_dispersion = TRUE` and `dispersion_source = "GRM"`.
#'
#' Ignored otherwise.
#'
#' Use this when genomic relatedness should influence grouping or dispersion.
#'
#' Matrix row and column names must match entry IDs, or be reachable through `id_map`.
#'
#' @param A Optional pedigree relationship matrix.
#'
#' Required when:
#' - `cluster_source = "A"`, or
#' - `use_dispersion = TRUE` and `dispersion_source = "A"`.
#'
#' Ignored otherwise.
#'
#' Use this when pedigree relatedness should influence grouping or dispersion.
#'
#' @param id_map Optional `data.frame` with columns `Treatment` and `LineID`.
#'
#' Used only when `cluster_source %in% c("GRM", "A")` and treatment IDs do not
#' already match matrix row/column names.
#'
#' Use this when field-book treatment labels are different from matrix IDs.
#'
#' Example:
#' treatment labels may be breeder-friendly names while `GRM` or `A` uses internal IDs.
#'
#' @param cluster_method Character specifying the clustering method applied after PCA
#' in matrix-based grouping.
#'
#' Allowed values:
#' - `"kmeans"`
#' - `"hclust"`
#'
#' Use `"kmeans"` when reproducible compact clustering is preferred.
#'
#' Use `"hclust"` when a hierarchical approach is preferred.
#'
#' Active only when `cluster_source %in% c("GRM", "A")`.
#'
#' @param cluster_seed Integer seed for k-means initialization.
#'
#' Active only when:
#' - `cluster_source %in% c("GRM", "A")`, and
#' - `cluster_method = "kmeans"`.
#'
#' Use this to make matrix-based grouping reproducible.
#'
#' @param cluster_attempts Integer number of random starts for k-means clustering.
#'
#' Active only when:
#' - `cluster_source %in% c("GRM", "A")`, and
#' - `cluster_method = "kmeans"`.
#'
#' Use larger values when clustering stability is important.
#'
#' @param n_pcs_use Integer or `Inf` giving the number of principal components used
#' for matrix-based clustering.
#'
#' Ignored when `cluster_source = "Family"`.
#'
#' Use smaller values when only broad structure should drive grouping.
#' Use `Inf` when the function should use as many informative components as possible.
#'
#' @param min_entry_slots_per_block Integer giving the minimum number of entry slots
#' allowed in any incomplete block.
#'
#' This argument constrains the automatically chosen number of blocks per replicate.
#'
#' Use smaller values when more blocks per replicate are acceptable.
#' Use larger values when blocks should not become too small once checks are inserted.
#'
#' This is particularly important when checks are numerous, because repeated checks
#' consume capacity in every incomplete block.
#'
#' @param max_blocks_per_rep Optional integer giving an upper bound on the number of
#' incomplete blocks per replicate.
#'
#' If `NULL`, the number of blocks per replicate is derived from capacity constraints.
#'
#' Use this when the user wants to prevent the design from creating too many
#' incomplete blocks, even if more could fit mathematically.
#'
#' Useful when operational simplicity or analysis preference favors fewer blocks.
#'
#' @param eval_efficiency Logical indicating whether efficiency diagnostics should be computed.
#'
#' Use `TRUE` when the user wants to compare design quality under a mixed model.
#'
#' Use `FALSE` when only the layout itself is needed.
#'
#' If `FALSE`, all efficiency-related arguments are ignored.
#'
#' @param treatment_effect Character indicating how entries are treated in efficiency evaluation.
#'
#' Allowed values:
#' - `"random"`
#' - `"fixed"`
#'
#' Active only when `eval_efficiency = TRUE`.
#'
#' Use `"fixed"` when interest lies in contrast precision among entry means.
#'
#' Use `"random"` when interest lies in BLUP-style prediction quality.
#'
#' @param prediction_type Character controlling the random-effect efficiency model.
#'
#' Allowed values:
#' - `"none"`
#' - `"IID"`
#' - `"GBLUP"`
#' - `"PBLUP"`
#'
#' Active only when:
#' - `eval_efficiency = TRUE`, and
#' - `treatment_effect = "random"`.
#'
#' Use `"none"` to skip random-effect efficiency.
#'
#' Use `"IID"` when entry random effects are assumed independent.
#'
#' Use `"GBLUP"` or `"PBLUP"` when prediction should use a relationship matrix `K`.
#'
#' @param K Optional relationship matrix used in two contexts:
#'
#' 1. random-effect efficiency when:
#'    - `eval_efficiency = TRUE`,
#'    - `treatment_effect = "random"`,
#'    - `prediction_type %in% c("GBLUP", "PBLUP")`;
#'
#' 2. dispersion optimization when:
#'    - `use_dispersion = TRUE`,
#'    - `dispersion_source = "K"`.
#'
#' Ignored otherwise.
#'
#' Use this when prediction or dispersion should be based on an explicit
#' relationship matrix.
#'
#' @param line_id_map Optional `data.frame` with columns `Treatment` and `LineID`.
#'
#' Used only when `K` is active and treatment labels differ from `rownames(K)`.
#'
#' This serves the same role for `K` that `id_map` serves for `GRM` or `A`.
#'
#' @param varcomp Named list of variance components used only in efficiency evaluation.
#'
#' Must contain:
#' - `sigma_e2`
#' - `sigma_g2`
#' - `sigma_rep2`
#' - `sigma_ib2`
#' - `sigma_r2`
#' - `sigma_c2`
#'
#' Use realistic values when efficiency should reflect a plausible trial model.
#'
#' Example:
#' use larger `sigma_ib2` when incomplete blocks are expected to differ strongly.
#'
#' @param check_as_fixed Logical indicating whether checks are included as fixed
#' indicators during efficiency evaluation.
#'
#' Active only when `eval_efficiency = TRUE`.
#'
#' Use `TRUE` when checks should be explicitly modeled as benchmark fixed effects.
#'
#' @param residual_structure Character specifying the residual correlation model.
#'
#' Allowed values:
#' - `"IID"`
#' - `"AR1"`
#' - `"AR1xAR1"`
#'
#' Active only when `eval_efficiency = TRUE`.
#'
#' Use `"IID"` when no spatial residual correlation is assumed.
#'
#' Use `"AR1"` when row-wise spatial correlation is expected.
#'
#' Use `"AR1xAR1"` when both row and column spatial correlation are expected.
#'
#' @param rho_row Numeric AR1 row parameter.
#'
#' Active only when:
#' - `eval_efficiency = TRUE`, and
#' - `residual_structure %in% c("AR1", "AR1xAR1")`.
#'
#' Use values near 0 for weak row correlation and larger absolute values for stronger correlation.
#'
#' @param rho_col Numeric AR1 column parameter.
#'
#' Active only when:
#' - `eval_efficiency = TRUE`, and
#' - `residual_structure = "AR1xAR1"`.
#'
#' Use this when column-wise correlation is expected in addition to row-wise correlation.
#'
#' @param spatial_engine Character retained for interface compatibility.
#'
#' Allowed values:
#' - `"auto"`
#' - `"sparse"`
#' - `"dense"`
#'
#' In this function the efficiency code primarily uses sparse Matrix operations,
#' but the argument is retained so the interface remains aligned with related functions.
#'
#' @param dense_max_n Integer retained for interface compatibility.
#'
#' This argument is included for consistency with related functions.
#'
#' @param eff_trace_samples Integer number of Hutchinson trace samples used when
#' approximate efficiency is needed because the target treatment dimension exceeds
#' `eff_full_max`.
#'
#' Larger values improve approximation stability but increase runtime.
#'
#' @param eff_full_max Integer maximum target dimension for exact efficiency extraction.
#'
#' Above this threshold, the function switches to an approximate trace-based method.
#'
#' Use larger values for more exact computation when memory and runtime allow.
#'
#' @param check_placement Character specifying how check positions are chosen within blocks.
#'
#' Allowed values:
#' - `"systematic"`
#' - `"random"`
#'
#' Use `"systematic"` when checks should be approximately evenly spread in each block.
#'
#' Use `"random"` when check positions may vary freely.
#'
#' @param check_position_pattern Argument retained for interface compatibility.
#'
#' It is not used by the current stream-based placement logic.
#'
#' @param use_dispersion Logical; if `TRUE`, apply post-hoc dispersion optimization
#' among non-check treatments.
#'
#' Use this when genetically or pedigree-similar entries should be less clustered
#' within replicates.
#'
#' If `FALSE`, all dispersion-related arguments are ignored.
#'
#' @param dispersion_source Character selecting which matrix is used for dispersion scoring.
#'
#' Allowed values:
#' - `"K"`
#' - `"A"`
#' - `"GRM"`
#'
#' Active only when `use_dispersion = TRUE`.
#'
#' Use `"K"` when a dedicated relationship matrix is preferred.
#' Use `"A"` for pedigree-based dispersion.
#' Use `"GRM"` for genomic-based dispersion.
#'
#' @param dispersion_radius Integer neighborhood radius for dispersion scoring.
#'
#' Two plots are considered neighbors when:
#'
#' `max(|dr|, |dc|) <= dispersion_radius`
#'
#' Use `1` for immediate neighbors only.
#' Use larger values when broader local dispersion matters.
#'
#' @param dispersion_iters Integer number of swap proposals used in the dispersion search.
#'
#' Larger values increase optimization effort and runtime.
#'
#' @param dispersion_seed Optional integer seed used in dispersion search.
#'
#' Use this when dispersion optimization should be reproducible independently of
#' the initial design seed.
#'
#' If left at its default, the function uses the supplied value in the interface.
#'
#' @param verbose Logical; if `TRUE`, print the derived replicate and block structure.
#'
#' Use this for diagnostics, debugging, or understanding how the function chose
#' replicate sizes and incomplete block sizes.
#'
#' @return A list with:
#' \describe{
#'   \item{layout_matrix}{A character matrix of size `n_rows × n_cols`. Unused cells are `NA`.}
#'   \item{field_book}{A plot-level data frame with one row per field cell, including unused cells.}
#'   \item{efficiency}{`NULL` or a list of efficiency summaries computed from non-`NA` plots.}
#'   \item{seed_used}{The actual random seed used internally.}
#'   \item{design_info}{A list summarizing replicate sizes, incomplete block sizes, used and unused cells, and key settings.}
#' }
#'
#' @examples
#' data("OptiDesign_example_data", package = "OptiDesign")
#' x <- OptiDesign_example_data
#'
#' ## Family-based row-column stream design
#' out_alpha_family <- do.call(
#'   prep_alpha_checks_rc_stream,
#'   c(
#'     x$OptiDesign_alpha_example,
#'     x$OptiDesign_alpha_args_family
#'   )
#' )
#'
#' dim(out_alpha_family$layout_matrix)
#' head(out_alpha_family$field_book)
#' out_alpha_family$design_info$n_blocks_per_rep
#'
#' \dontrun{
#' ## GRM-based row-column stream design with dispersion
#' out_alpha_grm <- do.call(
#'   prep_alpha_checks_rc_stream,
#'   c(
#'     x$OptiDesign_alpha_example,
#'     x$OptiDesign_alpha_args_grm
#'   )
#' )
#'
#' dim(out_alpha_grm$layout_matrix)
#' head(out_alpha_grm$field_book)
#' out_alpha_grm$efficiency
#' }
#'
#' @importFrom stats runif setNames
#' @importFrom utils head
#' @export
#' 
prep_alpha_checks_rc_stream <- function(
    check_treatments,
    check_families,
    entry_treatments,
    entry_families,
    n_reps,
    n_rows,
    n_cols,
    order = "column",
    serpentine = FALSE,
    seed = NULL,
    attempts = 5000,
    warn_and_correct = TRUE,
    fix_rows = TRUE,
    cluster_source = c("Family", "GRM", "A"),
    GRM = NULL,
    A = NULL,
    id_map = NULL,
    cluster_method = c("kmeans", "hclust"),
    cluster_seed = 1,
    cluster_attempts = 25,
    n_pcs_use = Inf,
    min_entry_slots_per_block = 8,
    max_blocks_per_rep = NULL,
    eval_efficiency = FALSE,
    treatment_effect = c("random", "fixed"),
    prediction_type = c("none", "IID", "GBLUP", "PBLUP"),
    K = NULL,
    line_id_map = NULL,
    varcomp = list(
      sigma_e2 = 1,
      sigma_g2 = 1,
      sigma_rep2 = 1,
      sigma_ib2 = 1,
      sigma_r2 = 1,
      sigma_c2 = 1
    ),
    check_as_fixed = TRUE,
    residual_structure = c("IID", "AR1", "AR1xAR1"),
    rho_row = 0,
    rho_col = 0,
    spatial_engine = c("auto", "sparse", "dense"),
    dense_max_n = 5000,
    eff_trace_samples = 80,
    eff_full_max = 400,
    check_placement = c("systematic", "random"),
    check_position_pattern = c("spread", "corners_first"),
    use_dispersion = FALSE,
    dispersion_source = c("K", "A", "GRM"),
    dispersion_radius = 1,
    dispersion_iters = 2000,
    dispersion_seed = 1,
    verbose = TRUE
) {
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  
  # ============================================================
  # 0. RNG
  # ============================================================
  seed_used <- seed
  if (is.null(seed_used)) seed_used <- sample.int(.Machine$integer.max, 1)
  set.seed(seed_used)
  
  .with_local_seed <- function(seed_local, expr) {
    old_seed <- NULL
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    set.seed(seed_local)
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    force(expr)
  }
  
  # ============================================================
  # 1. VALIDATION
  # ============================================================
  cluster_source <- match.arg(cluster_source)
  cluster_method <- match.arg(cluster_method)
  treatment_effect <- match.arg(treatment_effect)
  prediction_type <- match.arg(prediction_type)
  residual_structure <- match.arg(residual_structure)
  spatial_engine <- match.arg(spatial_engine)
  check_placement <- match.arg(check_placement)
  check_position_pattern <- match.arg(check_position_pattern)
  dispersion_source <- match.arg(dispersion_source)
  
  if (!order %in% c("row", "column")) {
    stop("Invalid 'order'. Use 'row' or 'column'.")
  }
  if (!is.logical(serpentine) || length(serpentine) != 1) {
    stop("serpentine must be TRUE or FALSE.")
  }
  if (n_reps < 1) stop("n_reps must be >= 1.")
  if (n_rows < 1 || n_cols < 1) stop("n_rows and n_cols must be >= 1.")
  if (length(check_treatments) != length(check_families)) {
    stop("Length of check_families must match length of check_treatments.")
  }
  if (length(entry_treatments) != length(entry_families)) {
    stop("Length of entry_families must match length of entry_treatments.")
  }
  if (anyDuplicated(check_treatments)) stop("Duplicate check_treatments found.")
  if (anyDuplicated(entry_treatments)) stop("Duplicate entry_treatments found.")
  if (length(intersect(check_treatments, entry_treatments)) > 0) {
    stop("A treatment cannot be both a check and an entry.")
  }
  if (!is.list(varcomp) ||
      !all(c("sigma_e2", "sigma_g2", "sigma_rep2", "sigma_ib2", "sigma_r2", "sigma_c2") %in% names(varcomp))) {
    stop("varcomp must contain sigma_e2, sigma_g2, sigma_rep2, sigma_ib2, sigma_r2, sigma_c2.")
  }
  if (abs(rho_row) >= 1) stop("rho_row must satisfy |rho_row| < 1.")
  if (abs(rho_col) >= 1) stop("rho_col must satisfy |rho_col| < 1.")
  if (!is.null(max_blocks_per_rep) && max_blocks_per_rep < 1) stop("max_blocks_per_rep must be >= 1.")
  if (min_entry_slots_per_block < 1) stop("min_entry_slots_per_block must be >= 1.")
  if (dispersion_radius < 1) stop("dispersion_radius must be >= 1.")
  if (dispersion_iters < 0) stop("dispersion_iters must be >= 0.")
  
  total_plots <- n_rows * n_cols
  n_checks <- length(check_treatments)
  n_entries <- length(entry_treatments)
  
  if (n_checks == 0) {
    warning("No checks were provided. Blocks will contain entries only.")
  }
  
  # ============================================================
  # 2. GLOBAL LOOKUPS
  # ============================================================
  family_lookup <- setNames(
    c(check_families, entry_families),
    c(check_treatments, entry_treatments)
  )
  
  gcluster_lookup <- setNames(
    rep(NA_character_, length(c(check_treatments, entry_treatments))),
    c(check_treatments, entry_treatments)
  )
  
  # ============================================================
  # 3. CLUSTERING FROM GRM / A (OPTIONAL)
  # ============================================================
  if (cluster_source %in% c("GRM", "A")) {
    Kc <- if (cluster_source == "GRM") GRM else A
    if (is.null(Kc)) stop(paste0("cluster_source='", cluster_source, "' selected but matrix is NULL."))
    if (is.null(rownames(Kc)) || is.null(colnames(Kc))) stop("GRM/A must have rownames and colnames.")
    
    k_clusters <- length(unique(entry_families))
    if (k_clusters < 2) stop("Need at least 2 unique families for matrix-derived clustering.")
    
    if (is.null(id_map)) {
      line_ids <- setNames(entry_treatments, entry_treatments)
    } else {
      if (!is.data.frame(id_map) || !all(c("Treatment", "LineID") %in% names(id_map))) {
        stop("id_map must contain columns Treatment and LineID.")
      }
      line_ids <- setNames(id_map$LineID, id_map$Treatment)
    }
    
    ids <- unname(line_ids[entry_treatments])
    missing_ids <- setdiff(ids, rownames(Kc))
    if (length(missing_ids) > 0) stop("Some entry LineIDs are missing in GRM/A.")
    
    Ksub <- Kc[ids, ids, drop = FALSE]
    eg <- eigen(Ksub, symmetric = TRUE)
    pos <- which(eg$values > 1e-10)
    if (length(pos) < 2) stop("GRM/A has too few positive eigenvalues for clustering.")
    
    max_possible_pcs <- min(length(pos), nrow(Ksub) - 1)
    n_pcs <- if (is.infinite(n_pcs_use)) max_possible_pcs else min(as.integer(n_pcs_use), max_possible_pcs)
    if (n_pcs < 2) stop("Fewer than 2 PCs available for clustering.")
    
    pcs <- eg$vectors[, pos[seq_len(n_pcs)], drop = FALSE]
    pcs <- sweep(pcs, 2, sqrt(eg$values[pos[seq_len(n_pcs)]]), `*`)
    
    if (cluster_method == "kmeans") {
      clust <- .with_local_seed(cluster_seed, {
        stats::kmeans(pcs, centers = k_clusters, nstart = cluster_attempts)$cluster
      })
    } else {
      hc <- stats::hclust(stats::dist(pcs), method = "ward.D2")
      clust <- stats::cutree(hc, k = k_clusters)
    }
    
    prefix <- if (cluster_source == "GRM") "G" else "A"
    gcluster_lookup[entry_treatments] <- paste0(prefix, clust)
  }
  
  get_adj_group <- function(trt) {
    if (cluster_source == "Family") {
      family_lookup[trt]
    } else if (trt %in% check_treatments) {
      family_lookup[trt]
    } else {
      gcluster_lookup[trt]
    }
  }
  
  # ============================================================
  # 4. HELPERS
  # ============================================================
  build_positions <- function(n_rows, n_cols, order = "column", serpentine = FALSE) {
    positions <- vector("list", n_rows * n_cols)
    k <- 1
    
    if (order == "row") {
      for (r in seq_len(n_rows)) {
        cols <- seq_len(n_cols)
        if (serpentine && r %% 2 == 0) cols <- rev(cols)
        for (c in cols) {
          positions[[k]] <- c(Row = r, Column = c)
          k <- k + 1
        }
      }
    } else {
      for (c in seq_len(n_cols)) {
        rows <- seq_len(n_rows)
        if (serpentine && c %% 2 == 0) rows <- rev(rows)
        for (r in rows) {
          positions[[k]] <- c(Row = r, Column = c)
          k <- k + 1
        }
      }
    }
    
    pos <- as.data.frame(do.call(rbind, positions))
    pos$Row <- as.integer(pos$Row)
    pos$Column <- as.integer(pos$Column)
    pos$PlotStream <- seq_len(nrow(pos))
    pos
  }
  
  split_integer <- function(total, parts) {
    base <- total %/% parts
    rem <- total %% parts
    c(rep(base + 1, rem), rep(base, parts - rem))
  }
  
  choose_n_blocks_per_rep <- function(total_plots,
                                      n_reps,
                                      n_entries,
                                      n_checks,
                                      min_entry_slots_per_block = 8,
                                      max_blocks_per_rep = NULL) {
    max_b_by_entries <- floor(n_entries / min_entry_slots_per_block)
    if (max_b_by_entries < 1) max_b_by_entries <- 1L
    
    if (n_checks > 0) {
      max_b_by_field <- floor((total_plots / n_reps - n_entries) / n_checks)
    } else {
      max_b_by_field <- max_b_by_entries
    }
    
    max_b <- min(max_b_by_entries, max_b_by_field)
    
    if (!is.null(max_blocks_per_rep)) {
      max_b <- min(max_b, max_blocks_per_rep)
    }
    
    if (is.na(max_b) || max_b < 1) {
      stop(
        paste0(
          "Field is too small to allocate ", n_entries,
          " entries across ", n_reps,
          " replicates with repeated checks in every block."
        )
      )
    }
    
    as.integer(if (max_b >= 2) max_b else 1L)
  }
  
  make_block_plan <- function(n_reps, n_entries, n_checks, n_blocks_per_rep) {
    entry_counts <- split_integer(n_entries, n_blocks_per_rep)
    block_sizes <- entry_counts + n_checks
    
    plans <- vector("list", n_reps)
    for (rr in seq_len(n_reps)) {
      plans[[rr]] <- data.frame(
        Rep = rr,
        BlockInRep = seq_len(n_blocks_per_rep),
        EntryCount = entry_counts,
        BlockSize = block_sizes,
        EntryCapacity = entry_counts,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, plans)
  }
  
  allocate_entries_to_blocks <- function(entries, block_entry_caps, family_lookup, get_adj_group, attempts = 1000) {
    b <- length(block_entry_caps)
    alloc <- vector("list", b)
    for (i in seq_len(b)) alloc[[i]] <- character(0)
    
    counts <- integer(b)
    names(counts) <- seq_len(b)
    
    entries_shuffled <- sample(entries)
    
    for (trt in entries_shuffled) {
      eligible <- which(counts < block_entry_caps)
      if (length(eligible) == 0) {
        stop("Internal error: insufficient block entry capacity.")
      }
      
      scores <- sapply(eligible, function(bb) {
        cur <- alloc[[bb]]
        fam_pen <- if (length(cur) == 0) 0 else sum(family_lookup[cur] == family_lookup[trt])
        grp_pen <- if (length(cur) == 0) 0 else sum(vapply(cur, function(x) get_adj_group(x) == get_adj_group(trt), logical(1)))
        load_pen <- counts[bb] / max(block_entry_caps[bb], 1)
        5 * grp_pen + 2 * fam_pen + load_pen + runif(1, 0, 1e-6)
      })
      
      bb_pick <- eligible[which.min(scores)]
      alloc[[bb_pick]] <- c(alloc[[bb_pick]], trt)
      counts[bb_pick] <- counts[bb_pick] + 1
    }
    
    alloc
  }
  
  improve_block_allocation <- function(alloc, block_entry_caps, max_iter = attempts) {
    score_alloc <- function(lst) {
      sc <- 0
      for (bb in seq_along(lst)) {
        x <- lst[[bb]]
        if (length(x) > 1) {
          fams <- family_lookup[x]
          sc <- sc + sum(outer(fams, fams, "==")) - length(x)
          grps <- vapply(x, get_adj_group, character(1))
          sc <- sc + 2 * (sum(outer(grps, grps, "==")) - length(x))
        }
      }
      sc
    }
    
    current <- alloc
    best <- current
    best_score <- score_alloc(best)
    
    if (length(current) < 2) return(best)
    
    for (it in seq_len(max_iter)) {
      b1 <- sample(seq_along(current), 1)
      b2 <- sample(setdiff(seq_along(current), b1), 1)
      if (length(current[[b1]]) == 0 || length(current[[b2]]) == 0) next
      
      i1 <- sample(seq_along(current[[b1]]), 1)
      i2 <- sample(seq_along(current[[b2]]), 1)
      
      cand <- current
      tmp <- cand[[b1]][i1]
      cand[[b1]][i1] <- cand[[b2]][i2]
      cand[[b2]][i2] <- tmp
      
      if (length(cand[[b1]]) > block_entry_caps[b1] || length(cand[[b2]]) > block_entry_caps[b2]) next
      
      sc <- score_alloc(cand)
      if (sc < best_score) {
        best <- cand
        best_score <- sc
        current <- cand
      }
    }
    
    best
  }
  
  get_check_positions_stream <- function(n_cells, n_checks, mode = c("systematic", "random")) {
    mode <- match.arg(mode)
    
    if (n_checks <= 0) return(integer(0))
    if (n_checks > n_cells) stop("n_checks cannot exceed number of cells in the block.")
    
    if (mode == "random") {
      return(sort(sample.int(n_cells, n_checks, replace = FALSE)))
    }
    
    pos <- round(seq(1, n_cells, length.out = n_checks + 2))[2:(n_checks + 1)]
    pos <- pmin(pmax(pos, 1), n_cells)
    
    pos <- unique(pos)
    if (length(pos) < n_checks) {
      candidates <- setdiff(seq_len(n_cells), pos)
      pos <- sort(c(pos, head(candidates, n_checks - length(pos))))
    }
    pos
  }
  
  arrange_block <- function(block_coords,
                            entries_block,
                            checks,
                            check_placement = c("systematic", "random"),
                            check_position_pattern = c("spread", "corners_first"),
                            Klocal = NULL,
                            family_lookup,
                            get_adj_group) {
    check_placement <- match.arg(check_placement)
    check_position_pattern <- match.arg(check_position_pattern)
    
    n <- nrow(block_coords)
    trt_vec <- rep(NA_character_, n)
    
    if (length(checks) > 0) {
      if (length(checks) > n) stop("Block has fewer cells than number of checks.")
      
      check_pos <- get_check_positions_stream(
        n_cells = n,
        n_checks = length(checks),
        mode = check_placement
      )
      
      chk_use <- sample(checks, length(checks), replace = FALSE)
      trt_vec[check_pos] <- chk_use
    }
    
    remaining_pos <- which(is.na(trt_vec))
    remaining_entries <- entries_block
    
    if (length(remaining_entries) > 0) {
      for (pp in remaining_pos) {
        if (length(remaining_entries) == 0) break
        
        neigh_idx <- which(
          pmax(abs(block_coords$Row - block_coords$Row[pp]),
               abs(block_coords$Column - block_coords$Column[pp])) <= 1 &
            seq_len(n) != pp
        )
        neigh_filled <- neigh_idx[!is.na(trt_vec[neigh_idx])]
        neigh_trt <- trt_vec[neigh_filled]
        
        scores <- sapply(seq_along(remaining_entries), function(j) {
          trt <- remaining_entries[j]
          sc <- 0
          if (length(neigh_trt) > 0) {
            sc <- sc + 10 * sum(vapply(neigh_trt, function(x) get_adj_group(x) == get_adj_group(trt), logical(1)))
            sc <- sc + 2 * sum(vapply(neigh_trt, function(x) family_lookup[x] == family_lookup[trt], logical(1)))
          }
          if (!is.null(Klocal) && length(neigh_trt) > 0) {
            neigh_noncheck <- neigh_trt[neigh_trt %in% colnames(Klocal)]
            if (length(neigh_noncheck) > 0 && trt %in% colnames(Klocal)) {
              sc <- sc + sum(Klocal[trt, neigh_noncheck, drop = TRUE])
            }
          }
          sc + runif(1, 0, 1e-6)
        })
        
        pick <- which.min(scores)
        trt_vec[pp] <- remaining_entries[pick]
        remaining_entries <- remaining_entries[-pick]
      }
    }
    
    trt_vec
  }
  
  make_sparse_incidence <- function(levels_vec) {
    lv <- unique(levels_vec[!is.na(levels_vec)])
    nn <- length(levels_vec)
    if (length(lv) == 0) {
      return(list(M = Matrix::Matrix(0, nrow = nn, ncol = 0, sparse = TRUE), levels = character(0)))
    }
    j <- match(levels_vec, lv)
    ok <- !is.na(j)
    M <- Matrix::sparseMatrix(i = which(ok), j = j[ok], x = 1, dims = c(nn, length(lv)))
    colnames(M) <- lv
    list(M = M, levels = lv)
  }
  
  pinv_sym_dense <- function(A, tol = 1e-10) {
    eg <- eigen(A, symmetric = TRUE)
    keep <- eg$values > tol
    if (!any(keep)) return(matrix(0, nrow(A), ncol(A)))
    eg$vectors[, keep, drop = FALSE] %*%
      diag(1 / eg$values[keep]) %*%
      t(eg$vectors[, keep, drop = FALSE])
  }
  
  safe_logdet_psd_dense <- function(A, tol = 1e-10) {
    eg <- eigen(A, symmetric = TRUE)
    lam <- eg$values[eg$values > tol]
    if (length(lam) == 0) return(-Inf)
    sum(log(lam))
  }
  
  pairwise_diff_mean_var <- function(V) {
    p <- nrow(V)
    if (p < 2) return(NA_real_)
    one <- rep(1, p)
    num <- p * sum(diag(V)) - as.numeric(t(one) %*% V %*% one)
    (2 * num) / (p * (p - 1))
  }
  
  ar1_precision_sparse <- function(nn, rho) {
    if (nn <= 0) stop("nn must be >= 1.")
    if (nn == 1) return(Matrix::Diagonal(1, 1))
    a <- 1 / (1 - rho^2)
    d <- rep((1 + rho^2) * a, nn)
    d[1] <- a
    d[nn] <- a
    o <- rep(-rho * a, nn - 1)
    Matrix::sparseMatrix(
      i = c(seq_len(nn), seq_len(nn - 1), 2:nn),
      j = c(seq_len(nn), 2:nn, seq_len(nn - 1)),
      x = c(d, o, o),
      dims = c(nn, nn)
    )
  }
  
  solve_C <- function(Cmat, B) {
    out <- try({
      fac <- Matrix::Cholesky(Cmat, LDL = TRUE, Imult = 0)
      Matrix::solve(fac, B)
    }, silent = TRUE)
    if (inherits(out, "try-error")) out <- Matrix::solve(Cmat, B)
    out
  }
  
  trace_subinv_est <- function(Cmat, idx, m = 80, seed_local = 1) {
    .with_local_seed(seed_local, {
      p <- length(idx)
      if (p == 0) return(NA_real_)
      nn <- nrow(Cmat)
      acc <- 0
      for (k in seq_len(m)) {
        z <- sample(c(-1, 1), p, replace = TRUE)
        u <- Matrix::sparseVector(i = idx, x = z, length = nn)
        x <- solve_C(Cmat, u)
        acc <- acc + as.numeric(Matrix::crossprod(u, x))
      }
      acc / m
    })
  }
  
  build_neighbor_pairs <- function(row, col, radius = 1) {
    n <- length(row)
    if (n < 2) return(matrix(integer(0), ncol = 2, dimnames = list(NULL, c("i", "j"))))
    out <- vector("list", 0)
    kk <- 1
    for (i in seq_len(n - 1)) {
      dr <- abs(row[i] - row[(i + 1):n])
      dc <- abs(col[i] - col[(i + 1):n])
      ok <- pmax(dr, dc) <= radius
      if (any(ok)) {
        jj <- (i + 1):n
        jj <- jj[ok]
        out[[kk]] <- cbind(i = rep(i, length(jj)), j = jj)
        kk <- kk + 1
      }
    }
    if (length(out) == 0) return(matrix(integer(0), ncol = 2, dimnames = list(NULL, c("i", "j"))))
    ans <- do.call(rbind, out)
    colnames(ans) <- c("i", "j")
    ans
  }
  
  score_dispersion <- function(trt_vec, is_check, Ksub, pairs) {
    if (nrow(pairs) == 0) return(0)
    line_levels <- colnames(Ksub)
    line_idx <- rep(NA_integer_, length(trt_vec))
    line_idx[!is_check] <- match(trt_vec[!is_check], line_levels)
    ii <- pairs[, "i"]
    jj <- pairs[, "j"]
    li <- line_idx[ii]
    lj <- line_idx[jj]
    ok <- !is.na(li) & !is.na(lj)
    if (!any(ok)) return(0)
    sum(Ksub[cbind(li[ok], lj[ok])])
  }
  
  # ============================================================
  # 5. GLOBAL STREAM -> REPLICATES -> BLOCKS
  # ============================================================
  pos_all <- build_positions(n_rows, n_cols, order = order, serpentine = serpentine)
  
  n_blocks_per_rep <- choose_n_blocks_per_rep(
    total_plots = total_plots,
    n_reps = n_reps,
    n_entries = n_entries,
    n_checks = n_checks,
    min_entry_slots_per_block = min_entry_slots_per_block,
    max_blocks_per_rep = max_blocks_per_rep
  )
  
  block_plan <- make_block_plan(
    n_reps = n_reps,
    n_entries = n_entries,
    n_checks = n_checks,
    n_blocks_per_rep = n_blocks_per_rep
  )
  
  rep_used_size <- sum(block_plan$BlockSize[block_plan$Rep == 1])
  rep_sizes <- rep(rep_used_size, n_reps)
  total_used_plots <- sum(rep_sizes)
  
  if (total_used_plots > total_plots) {
    stop(
      paste0(
        "The fixed field (", n_rows, " x ", n_cols, " = ", total_plots,
        " plots) is too small. Required used plots = ", total_used_plots, "."
      )
    )
  }
  
  rep_starts <- cumsum(c(1, head(rep_sizes, -1)))
  rep_ends <- cumsum(rep_sizes)
  
  if (verbose) {
    message(
      paste0(
        "Fixed field = ", n_rows, " x ", n_cols,
        "; used plots = ", total_used_plots,
        "; trailing NA plots = ", total_plots - total_used_plots,
        "; replicate used sizes = {", paste(rep_sizes, collapse = ", "), "}",
        "; blocks/rep = ", n_blocks_per_rep,
        "; block sizes in rep 1 = {",
        paste(block_plan$BlockSize[block_plan$Rep == 1], collapse = ", "),
        "}"
      )
    )
  }
  
  # ============================================================
  # 6. BUILD FIELD BOOK TEMPLATE (INCLUDING TRAILING NA CELLS)
  # ============================================================
  field_book <- pos_all
  field_book$Rep <- NA_integer_
  field_book$IBlock <- NA_integer_
  field_book$BlockInRep <- NA_integer_
  field_book$Treatment <- NA_character_
  field_book$Family <- NA_character_
  field_book$Gcluster <- NA_character_
  field_book$Check <- FALSE
  
  for (rr in seq_len(n_reps)) {
    idx <- rep_starts[rr]:rep_ends[rr]
    field_book$Rep[idx] <- rr
  }
  
  global_block_counter <- 1L
  block_meta_list <- list()
  
  for (rr in seq_len(n_reps)) {
    rep_idx <- rep_starts[rr]:rep_ends[rr]
    rep_block_sizes <- block_plan$BlockSize[block_plan$Rep == rr]
    block_starts <- cumsum(c(1, head(rep_block_sizes, -1)))
    block_ends <- cumsum(rep_block_sizes)
    
    for (bb in seq_along(rep_block_sizes)) {
      local_idx <- block_starts[bb]:block_ends[bb]
      global_idx <- rep_idx[local_idx]
      field_book$IBlock[global_idx] <- global_block_counter
      field_book$BlockInRep[global_idx] <- bb
      
      block_meta_list[[length(block_meta_list) + 1]] <- data.frame(
        Rep = rr,
        IBlock = global_block_counter,
        BlockInRep = bb,
        BlockSize = rep_block_sizes[bb],
        EntryCapacity = block_plan$EntryCapacity[block_plan$Rep == rr][bb],
        StartStream = min(global_idx),
        EndStream = max(global_idx),
        stringsAsFactors = FALSE
      )
      
      global_block_counter <- global_block_counter + 1L
    }
  }
  
  block_meta <- do.call(rbind, block_meta_list)
  
  # ============================================================
  # 7. ALLOCATE ENTRIES TO BLOCKS PER REPLICATE
  # ============================================================
  rep_allocations <- vector("list", n_reps)
  names(rep_allocations) <- paste0("Rep", seq_len(n_reps))
  
  for (rr in seq_len(n_reps)) {
    caps <- block_meta$EntryCapacity[block_meta$Rep == rr]
    alloc_rr <- allocate_entries_to_blocks(
      entries = entry_treatments,
      block_entry_caps = caps,
      family_lookup = family_lookup,
      get_adj_group = get_adj_group,
      attempts = attempts
    )
    alloc_rr <- improve_block_allocation(
      alloc = alloc_rr,
      block_entry_caps = caps,
      max_iter = attempts
    )
    rep_allocations[[rr]] <- alloc_rr
  }
  
  # ============================================================
  # 8. PLACE CHECKS AND ENTRIES BLOCK BY BLOCK
  # ============================================================
  K_arrange <- NULL
  if (cluster_source %in% c("GRM", "A")) {
    K_arrange <- if (cluster_source == "GRM") GRM else A
  }
  
  for (rr in seq_len(n_reps)) {
    rep_blocks <- block_meta[block_meta$Rep == rr, , drop = FALSE]
    
    for (bb in seq_len(nrow(rep_blocks))) {
      block_id <- rep_blocks$IBlock[bb]
      block_idx <- which(field_book$IBlock == block_id)
      block_coords <- field_book[block_idx, c("Row", "Column", "PlotStream"), drop = FALSE]
      entries_block <- rep_allocations[[rr]][[bb]]
      
      Klocal <- NULL
      if (!is.null(K_arrange) && length(entries_block) > 1) {
        ids_here <- entries_block
        if (!is.null(id_map)) {
          mapv <- setNames(id_map$LineID, id_map$Treatment)
          ids_here <- unname(mapv[entries_block])
        }
        ok <- !is.na(ids_here) & ids_here %in% rownames(K_arrange)
        if (sum(ok) > 1) {
          Ktmp <- K_arrange[ids_here[ok], ids_here[ok], drop = FALSE]
          colnames(Ktmp) <- entries_block[ok]
          rownames(Ktmp) <- entries_block[ok]
          Klocal <- Ktmp
        }
      }
      
      trt_vec <- arrange_block(
        block_coords = block_coords,
        entries_block = entries_block,
        checks = check_treatments,
        check_placement = check_placement,
        check_position_pattern = check_position_pattern,
        Klocal = Klocal,
        family_lookup = family_lookup,
        get_adj_group = get_adj_group
      )
      
      field_book$Treatment[block_idx] <- trt_vec
      not_na <- !is.na(trt_vec)
      field_book$Family[block_idx[not_na]] <- family_lookup[trt_vec[not_na]]
      field_book$Gcluster[block_idx[not_na]] <- gcluster_lookup[trt_vec[not_na]]
      field_book$Check[block_idx] <- !is.na(trt_vec) & trt_vec %in% check_treatments
    }
  }
  
  # ============================================================
  # 9. OPTIONAL DISPERSION (WITHIN REPLICATE ONLY)
  # ============================================================
  if (isTRUE(use_dispersion)) {
    Kdisp <- switch(dispersion_source, "K" = K, "A" = A, "GRM" = GRM)
    if (is.null(Kdisp)) stop("use_dispersion=TRUE but selected dispersion matrix is NULL.")
    if (is.null(rownames(Kdisp)) || is.null(colnames(Kdisp))) stop("Dispersion matrix must have rownames and colnames.")
    
    .with_local_seed(if (is.null(dispersion_seed)) seed_used else dispersion_seed, {
      trt <- as.character(field_book$Treatment)
      is_check <- !is.na(trt) & trt %in% check_treatments
      movable <- which(!is.na(trt) & !is_check)
      
      if (length(movable) >= 2 && dispersion_iters > 0) {
        non_trt <- unique(trt[!is.na(trt) & !is_check])
        
        if (is.null(line_id_map)) {
          line_ids <- setNames(non_trt, non_trt)
        } else {
          if (!is.data.frame(line_id_map) || !all(c("Treatment", "LineID") %in% names(line_id_map))) {
            stop("line_id_map must contain Treatment and LineID.")
          }
          line_ids <- setNames(line_id_map$LineID, line_id_map$Treatment)
        }
        
        ids <- unname(line_ids[non_trt])
        miss <- setdiff(ids, rownames(Kdisp))
        if (length(miss) > 0) stop("Some LineIDs are missing in selected dispersion matrix.")
        
        Ksub <- Kdisp[ids, ids, drop = FALSE]
        rownames(Ksub) <- colnames(Ksub) <- non_trt
        
        pairs <- build_neighbor_pairs(field_book$Row[!is.na(trt)], field_book$Column[!is.na(trt)], radius = dispersion_radius)
        trt_obs <- trt[!is.na(trt)]
        is_check_obs <- is_check[!is.na(trt)]
        obs_index <- which(!is.na(trt))
        
        best_trt <- trt_obs
        best_score <- score_dispersion(best_trt, is_check_obs, Ksub, pairs)
        
        movable_obs <- which(!is_check_obs)
        
        for (it in seq_len(dispersion_iters)) {
          ij <- sample(movable_obs, 2, replace = FALSE)
          if (field_book$Rep[obs_index[ij[1]]] != field_book$Rep[obs_index[ij[2]]]) next
          
          cand <- best_trt
          cand[ij] <- cand[rev(ij)]
          sc <- score_dispersion(cand, is_check_obs, Ksub, pairs)
          
          if (sc < best_score) {
            best_score <- sc
            best_trt <- cand
          }
        }
        
        field_book$Treatment[obs_index] <- best_trt
        field_book$Family[obs_index] <- family_lookup[field_book$Treatment[obs_index]]
        field_book$Gcluster[obs_index] <- gcluster_lookup[field_book$Treatment[obs_index]]
        field_book$Check[obs_index] <- field_book$Treatment[obs_index] %in% check_treatments
      }
    })
  }
  
  # ============================================================
  # 10. LAYOUT MATRIX
  # ============================================================
  layout_matrix <- matrix(NA_character_, nrow = n_rows, ncol = n_cols)
  for (i in seq_len(nrow(field_book))) {
    layout_matrix[field_book$Row[i], field_book$Column[i]] <- field_book$Treatment[i]
  }
  
  # ============================================================
  # 11. EFFICIENCY EVALUATION
  # ============================================================
  efficiency <- NULL
  
  if (isTRUE(eval_efficiency) && (treatment_effect == "fixed" || prediction_type != "none")) {
    fb <- field_book[!is.na(field_book$Treatment), , drop = FALSE]
    nn <- nrow(fb)
    
    trt <- as.character(fb$Treatment)
    repf <- as.character(fb$Rep)
    ibf  <- as.character(fb$IBlock)
    rw   <- as.character(fb$Row)
    cl   <- as.character(fb$Column)
    is_check <- trt %in% check_treatments
    
    sigma_e2   <- varcomp$sigma_e2
    sigma_g2   <- varcomp$sigma_g2
    sigma_rep2 <- varcomp$sigma_rep2
    sigma_ib2  <- varcomp$sigma_ib2
    sigma_r2   <- varcomp$sigma_r2
    sigma_c2   <- varcomp$sigma_c2
    
    if (any(c(sigma_e2, sigma_rep2, sigma_ib2, sigma_r2, sigma_c2) <= 0)) {
      stop("sigma_e2, sigma_rep2, sigma_ib2, sigma_r2, sigma_c2 must be > 0.")
    }
    if (treatment_effect == "random" && sigma_g2 <= 0) {
      stop("sigma_g2 must be > 0 when treatment_effect='random'.")
    }
    
    if (residual_structure == "IID") {
      Q <- (1 / sigma_e2) * Matrix::Diagonal(nn, 1)
    } else {
      Qrow <- ar1_precision_sparse(n_rows, rho_row)
      Qcol <- if (residual_structure == "AR1") Matrix::Diagonal(n_cols, 1) else ar1_precision_sparse(n_cols, rho_col)
      Qgrid <- Matrix::kronecker(Qcol, Qrow) * (1 / sigma_e2)
      grid_index <- (as.integer(rw) - 1) * n_cols + as.integer(cl)
      Q <- Qgrid[grid_index, grid_index, drop = FALSE]
    }
    
    X <- Matrix::Matrix(rep(1, nn), ncol = 1, sparse = TRUE)
    colnames(X) <- "(Intercept)"
    
    Xrep <- make_sparse_incidence(repf)$M
    if (ncol(Xrep) > 1) {
      Xrep <- Xrep[, -1, drop = FALSE]
      colnames(Xrep) <- paste0("Rep_", seq_len(ncol(Xrep)))
      X <- cbind(X, Xrep)
    }
    
    if (isTRUE(check_as_fixed)) {
      chk_vec <- ifelse(is_check, trt, NA)
      Xchk <- make_sparse_incidence(chk_vec)$M
      if (ncol(Xchk) > 0) {
        colnames(Xchk) <- paste0("Check_", colnames(Xchk))
        X <- cbind(X, Xchk)
      }
    }
    
    if (treatment_effect == "fixed") {
      trt_fix_vec <- ifelse(!is_check, trt, NA)
      Xtreat <- make_sparse_incidence(trt_fix_vec)$M
      if (ncol(Xtreat) < 2) stop("Not enough fixed non-check treatments.")
      colnames(Xtreat) <- paste0("Line_", colnames(Xtreat))
      X <- cbind(X, Xtreat)
    }
    
    Zib <- make_sparse_incidence(ibf)$M
    Zr  <- make_sparse_incidence(rw)$M
    Zc  <- make_sparse_incidence(cl)$M
    colnames(Zib) <- paste0("IB_", colnames(Zib))
    colnames(Zr)  <- paste0("Row_", colnames(Zr))
    colnames(Zc)  <- paste0("Col_", colnames(Zc))
    
    Z_list <- list(Zib = Zib, Zr = Zr, Zc = Zc)
    Zg <- NULL
    
    if (treatment_effect == "random") {
      g_vec <- ifelse(!is_check, trt, NA)
      Zg <- make_sparse_incidence(g_vec)$M
      if (ncol(Zg) < 2) stop("Not enough random non-check treatments.")
      colnames(Zg) <- paste0("Line_", colnames(Zg))
      Z_list$Zg <- Zg
    }
    
    Z <- do.call(cbind, Z_list)
    
    XtQX <- Matrix::crossprod(X, Q %*% X)
    XtQZ <- Matrix::crossprod(X, Q %*% Z)
    ZtQX <- t(XtQZ)
    ZtQZ <- Matrix::crossprod(Z, Q %*% Z)
    
    Ginv <- Matrix::Diagonal(ncol(Z), 0)
    
    idx0 <- 0
    pIB <- ncol(Zib)
    pR  <- ncol(Zr)
    pC  <- ncol(Zc)
    
    if (pIB > 0) {
      ii <- (idx0 + 1):(idx0 + pIB)
      Ginv[ii, ii] <- Matrix::Diagonal(pIB, 1 / sigma_ib2)
      idx0 <- idx0 + pIB
    }
    if (pR > 0) {
      ii <- (idx0 + 1):(idx0 + pR)
      Ginv[ii, ii] <- Matrix::Diagonal(pR, 1 / sigma_r2)
      idx0 <- idx0 + pR
    }
    if (pC > 0) {
      ii <- (idx0 + 1):(idx0 + pC)
      Ginv[ii, ii] <- Matrix::Diagonal(pC, 1 / sigma_c2)
      idx0 <- idx0 + pC
    }
    
    trt_idx <- integer(0)
    
    if (treatment_effect == "random") {
      pG <- ncol(Zg)
      trt_idx <- (idx0 + 1):(idx0 + pG)
      
      if (prediction_type == "IID") {
        Ginv[trt_idx, trt_idx] <- Matrix::Diagonal(pG, 1 / sigma_g2)
      } else {
        if (is.null(K)) stop("K must be supplied for GBLUP/PBLUP.")
        if (is.null(rownames(K)) || is.null(colnames(K))) stop("K must have rownames and colnames.")
        
        line_trt_order <- sub("^Line_", "", colnames(Zg))
        if (is.null(line_id_map)) {
          line_ids <- setNames(line_trt_order, line_trt_order)
        } else {
          if (!is.data.frame(line_id_map) || !all(c("Treatment", "LineID") %in% names(line_id_map))) {
            stop("line_id_map must contain Treatment and LineID.")
          }
          line_ids <- setNames(line_id_map$LineID, line_id_map$Treatment)
        }
        
        ids <- unname(line_ids[line_trt_order])
        miss <- setdiff(ids, rownames(K))
        if (length(miss) > 0) stop("Some treatment IDs are missing in K.")
        
        Ksub2 <- K[ids, ids, drop = FALSE]
        Kinv_try <- try({
          Km <- Matrix::Matrix(Ksub2, sparse = FALSE)
          facK <- Matrix::Cholesky(Km, LDL = FALSE, Imult = 0)
          Matrix::solve(facK, Matrix::Diagonal(nrow(Km), 1))
        }, silent = TRUE)
        
        if (!inherits(Kinv_try, "try-error")) {
          Kinv <- Kinv_try
        } else {
          Kinv <- Matrix::Matrix(pinv_sym_dense(Ksub2), sparse = FALSE)
        }
        
        Ginv[trt_idx, trt_idx] <- (1 / sigma_g2) * Kinv
      }
    }
    
    Cmat <- rbind(cbind(XtQX, XtQZ), cbind(ZtQX, ZtQZ + Ginv))
    Cmat <- Matrix::Matrix(Cmat, sparse = TRUE)
    
    eff <- list(
      treatment_effect = treatment_effect,
      prediction_type = if (treatment_effect == "random") prediction_type else NA_character_,
      residual_structure = residual_structure,
      n_reps = n_reps,
      n_blocks_total = nrow(block_meta),
      order = order,
      serpentine = serpentine
    )
    
    if (treatment_effect == "fixed") {
      xnames <- colnames(X)
      trt_cols <- grep("^Line_", xnames)
      p <- length(trt_cols)
      
      if (p <= eff_full_max) {
        B <- Matrix::sparseMatrix(i = trt_cols, j = seq_along(trt_cols), x = 1, dims = c(nrow(Cmat), p))
        Xsol <- solve_C(Cmat, B)
        Vsub <- as.matrix(Xsol[trt_cols, , drop = FALSE])
        
        mean_var_diff <- pairwise_diff_mean_var(Vsub)
        A_opt <- 1 / mean_var_diff
        
        H <- diag(p) - matrix(1 / p, p, p)
        Vctr <- H %*% Vsub %*% H
        logdet <- safe_logdet_psd_dense(Vctr)
        D_opt <- if (is.finite(logdet)) exp(-logdet / (p - 1)) else NA_real_
        
        eff$mode <- "FIXED_TREATMENT_BLUE_CONTRAST"
        eff$A <- A_opt
        eff$D <- D_opt
        eff$mean_VarDiff <- mean_var_diff
        eff$n_trt <- p
      } else {
        tr_est <- trace_subinv_est(Cmat, trt_cols, m = eff_trace_samples, seed_local = 1)
        mean_var_coef <- tr_est / p
        eff$mode <- "FIXED_TREATMENT_BLUE_APPROX"
        eff$A <- 1 / mean_var_coef
        eff$D <- NA_real_
        eff$mean_VarCoef <- mean_var_coef
        eff$n_trt <- p
      }
      
    } else {
      p <- length(trt_idx)
      
      if (p <= eff_full_max) {
        B <- Matrix::sparseMatrix(i = trt_idx, j = seq_along(trt_idx), x = 1, dims = c(nrow(Cmat), p))
        Xsol <- solve_C(Cmat, B)
        PEVsub <- as.matrix(Xsol[trt_idx, , drop = FALSE])
        
        eff$mode <- paste0("RANDOM_TREATMENT_BLUP_", prediction_type)
        eff$A <- 1 / mean(diag(PEVsub))
        logdet <- safe_logdet_psd_dense(PEVsub)
        eff$D <- if (is.finite(logdet)) exp(-logdet / p) else NA_real_
        eff$mean_PEV <- mean(diag(PEVsub))
        eff$n_lines <- p
      } else {
        tr_est <- trace_subinv_est(Cmat, trt_idx, m = eff_trace_samples, seed_local = 1)
        mean_pev <- tr_est / p
        eff$mode <- paste0("RANDOM_TREATMENT_BLUP_", prediction_type, "_APPROX")
        eff$A <- 1 / mean_pev
        eff$D <- NA_real_
        eff$mean_PEV <- mean_pev
        eff$n_lines <- p
      }
    }
    
    efficiency <- eff
  }
  
  # ============================================================
  # 12. FINAL FIELD BOOK COLUMN ORDER
  # ============================================================
  field_book <- field_book[order(field_book$PlotStream), ]
  field_book <- field_book[, c(
    "Treatment", "Family", "Gcluster", "Check",
    "PlotStream", "Rep", "IBlock", "BlockInRep",
    "Row", "Column"
  )]
  
  # ============================================================
  # 13. RETURN
  # ============================================================
  design_info <- list(
    n_reps = n_reps,
    n_rows = n_rows,
    n_cols = n_cols,
    total_plots = total_plots,
    used_plots = total_used_plots,
    unused_plots = total_plots - total_used_plots,
    rep_sizes = rep_sizes,
    n_blocks_per_rep = n_blocks_per_rep,
    block_meta = block_meta,
    n_checks = n_checks,
    n_entries = n_entries,
    min_entry_slots_per_block = min_entry_slots_per_block,
    order = order,
    serpentine = serpentine,
    warn_and_correct = warn_and_correct,
    fix_rows = fix_rows
  )
  
  list(
    layout_matrix = layout_matrix,
    field_book = field_book,
    efficiency = efficiency,
    seed_used = seed_used,
    design_info = design_info
  )
}
