#' Example data for OptiDesign
#'
#' Synthetic example data shipped with `OptiDesign` for illustrating and
#' testing the main design constructors, relationship-matrix workflows,
#' optional dispersion optimisation, and efficiency evaluation across both
#' design families.
#'
#' @name OptiDesign_example_data
#' @docType data
#' @keywords datasets
#'
#' @format A named list with the following components:
#' \describe{
#'
#'   \item{`OptiDesign_lines`}{A data frame with columns `Treatment` and
#'   `Family` representing a synthetic pool of line identifiers and their
#'   family labels. Used as a starting point for constructing both
#'   `prep_famoptg()` and `alpha_rc_stream()` treatment vectors.}
#'
#'   \item{`OptiDesign_id_map`}{A data frame with columns `Treatment` and
#'   `LineID` used to map treatment labels to the rownames and colnames of
#'   relationship matrices when those names differ from the treatment
#'   identifiers. Relevant for `cluster_source %in% c("GRM", "A")` and for
#'   `prediction_type %in% c("GBLUP", "PBLUP")`.}
#'
#'   \item{`OptiDesign_GRM`}{A synthetic genomic relationship matrix with
#'   rownames and colnames corresponding to line IDs. Can be used for
#'   matrix-based grouping (`cluster_source = "GRM"`), dispersion optimisation
#'   (`dispersion_source = "GRM"`), or as a `K` matrix for GBLUP efficiency
#'   evaluation.}
#'
#'   \item{`OptiDesign_A`}{A synthetic pedigree-style relationship matrix with
#'   rownames and colnames corresponding to line IDs. Can be used for
#'   pedigree-based grouping (`cluster_source = "A"`), dispersion optimisation
#'   (`dispersion_source = "A"`), or as a `K` matrix for PBLUP efficiency
#'   evaluation.}
#'
#'   \item{`OptiDesign_K`}{A synthetic genomic kernel matrix with rownames and
#'   colnames corresponding to line IDs. Intended for workflows that supply a
#'   `K` matrix directly, such as `prediction_type = "GBLUP"` in
#'   `evaluate_famoptg_efficiency()` or `evaluate_alpha_efficiency()`, CDmean
#'   computation in `optimize_famoptg()` or `optimize_alpha_rc()`, and
#'   dispersion optimisation with `dispersion_source = "K"`.}
#'
#'   \item{`OptiDesign_famoptg_example`}{A named list containing the core
#'   treatment vectors and field dimensions needed to call `prep_famoptg()`.
#'   Includes check treatments, p-rep treatments with replication counts,
#'   unreplicated treatments, family labels for all treatment classes, number
#'   of blocks, and field dimensions (`n_rows`, `n_cols`). Combine with one of
#'   `OptiDesign_famoptg_args_family` or `OptiDesign_famoptg_args_grm` to form
#'   a complete argument list.}
#'
#'   \item{`OptiDesign_alpha_example`}{A named list containing the core
#'   treatment vectors and field dimensions needed to call
#'   `alpha_rc_stream()`. Includes check treatments, entry treatments,
#'   family labels, number of replicates, and field dimensions (`n_rows`,
#'   `n_cols`). Combine with one of `OptiDesign_alpha_args_family` or
#'   `OptiDesign_alpha_args_grm` to form a complete argument list.}
#'
#'   \item{`OptiDesign_famoptg_args_family`}{A named list of supplementary
#'   arguments for a family-based call to `prep_famoptg()`. Specifies
#'   `cluster_source = "Family"`, field traversal settings, block placement
#'   options, and does not include efficiency evaluation arguments (those are
#'   now passed to `evaluate_famoptg_efficiency()` separately).}
#'
#'   \item{`OptiDesign_famoptg_args_grm`}{A named list of supplementary
#'   arguments for a GRM-based call to `prep_famoptg()`. Specifies
#'   `cluster_source = "GRM"`, includes `GRM` and `id_map` references,
#'   and illustrates matrix-based grouping and optional dispersion settings.
#'   Efficiency evaluation arguments are passed to
#'   `evaluate_famoptg_efficiency()` separately.}
#'
#'   \item{`OptiDesign_alpha_args_family`}{A named list of supplementary
#'   arguments for a family-based call to `alpha_rc_stream()`. Specifies
#'   `cluster_source = "Family"`, field traversal settings, and block size
#'   constraints. Efficiency evaluation arguments are passed to
#'   `evaluate_alpha_efficiency()` separately.}
#'
#'   \item{`OptiDesign_alpha_args_grm`}{A named list of supplementary
#'   arguments for a GRM-based call to `alpha_rc_stream()`. Specifies
#'   `cluster_source = "GRM"`, includes `GRM` and `id_map` references,
#'   and illustrates matrix-based grouping and optional dispersion settings.
#'   Efficiency evaluation arguments are passed to
#'   `evaluate_alpha_efficiency()` separately.}
#'
#' }
#'
#' @details
#' The objects in this dataset are fully synthetic and are intended for:
#'
#' - package examples and vignettes,
#' - unit tests,
#' - workflow demonstrations,
#' - and user exploration of argument combinations.
#'
#' The dataset covers both design families in the package and is structured
#' so that `OptiDesign_famoptg_example` and `OptiDesign_alpha_example`
#' provide the treatment and field inputs, while the `*_args_*` lists provide
#' the algorithmic and model settings. This separation mirrors the
#' construct-then-evaluate architecture of the package.
#'
#' ## Illustrative workflows
#'
#' **Repeated-check block design (family-based):**
#' ```r
#' x <- OptiDesign_example_data
#'
#' design <- do.call(prep_famoptg,
#'   c(x$OptiDesign_famoptg_example, x$OptiDesign_famoptg_args_family)
#' )
#'
#' eff <- evaluate_famoptg_efficiency(
#'   field_book       = design$field_book,
#'   n_rows           = x$OptiDesign_famoptg_example$n_rows,
#'   n_cols           = x$OptiDesign_famoptg_example$n_cols,
#'   check_treatments = x$OptiDesign_famoptg_example$check_treatments,
#'   treatment_effect = "fixed"
#' )
#' ```
#'
#' **Alpha row-column stream design (GRM-based):**
#' ```r
#' x <- OptiDesign_example_data
#'
#' design <- do.call(alpha_rc_stream,
#'   c(x$OptiDesign_alpha_example, x$OptiDesign_alpha_args_grm)
#' )
#'
#' eff <- evaluate_alpha_efficiency(
#'   field_book       = design$field_book,
#'   n_rows           = x$OptiDesign_alpha_example$n_rows,
#'   n_cols           = x$OptiDesign_alpha_example$n_cols,
#'   check_treatments = x$OptiDesign_alpha_example$check_treatments,
#'   treatment_effect = "fixed"
#' )
#' ```
#'
#' **Optimised repeated-check block design:**
#' ```r
#' x <- OptiDesign_example_data
#'
#' opt <- do.call(optimize_famoptg,
#'   c(
#'     x$OptiDesign_famoptg_example,
#'     x$OptiDesign_famoptg_args_family,
#'     list(
#'       treatment_effect = "fixed",
#'       criterion        = "A",
#'       n_restarts       = 20
#'     )
#'   )
#' )
#' opt$optimization$best_score
#' ```
#'
#' **Optimised alpha row-column stream design:**
#' ```r
#' x <- OptiDesign_example_data
#'
#' opt <- do.call(optimize_alpha_rc,
#'   c(
#'     x$OptiDesign_alpha_example,
#'     x$OptiDesign_alpha_args_family,
#'     list(
#'       treatment_effect = "fixed",
#'       method           = "RS",
#'       criterion        = "A",
#'       n_restarts       = 20
#'     )
#'   )
#' )
#' opt$optimization$best_score
#' ```
#'
#' @source Generated internally by `data-raw/generate_example_data.R`.
NULL
