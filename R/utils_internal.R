# ==============================================================================
# utils_internal.R
# Internal utilities for OptiDesign
#
# Contents:
#   - Global variable declarations (suppresses R CMD CHECK NOTEs)
#   - .optidesign_dbg()        - conditional debug message helper
#   - .optidesign_z()          - z-score normalisation helper
#   - .validate_square_matrix() - square matrix validator
# ==============================================================================

# nocov start
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    
    # -- Field book columns ---------------------------------------------------
    # Columns present in field books returned by prep_famoptg() and
    # alpha_rc_stream(). Declared to suppress R CMD CHECK notes when these
    # names appear unquoted in NSE contexts (e.g. dplyr pipelines in examples).
    "Plot",
    "Row",
    "Col",
    "Column",
    "Block",
    "IBlock",
    "BlockInRep",
    "Rep",
    "Replicate",
    "Entry",
    "Treatment",
    "Family",
    "Gcluster",
    "Check",
    "IsCheck",
    "StreamPos",
    
    # -- Optimizer / efficiency internals -------------------------------------
    # Names used internally by optimize_famoptg() and optimize_alpha_rc()
    # that may appear in NSE expressions.
    "score",
    "best_score",
    "n_failed",
    
    # -- Relationship matrix helpers ------------------------------------------
    # Column-selection helper used in matrix subsetting utilities.
    "..keep_cols",
    
    # -- Dispersion optimisation ----------------------------------------------
    # Variables used internally during the swap-based local search in
    # prep_famoptg() and alpha_rc_stream().
    "SpatialGroup",
    "RelatednessScore",
    
    # -- pracma import --------------------------------------------------------
    # mod() is imported from pracma for serpentine traversal parity in
    # prep_famoptg(). Declared here to suppress NOTE about binding.
    "mod"
  ))
}
# nocov end


# ------------------------------------------------------------------------------
# .optidesign_dbg
# ------------------------------------------------------------------------------
#' Internal conditional debug message helper
#'
#' Emits a formatted message via [base::message()] when `debug = TRUE`.
#' Used internally to provide optional verbose output during design
#' construction and optimisation without cluttering normal output.
#'
#' @param debug Logical. If `TRUE`, the message is emitted. If `FALSE` or
#'   `NULL`, the function returns invisibly with no output.
#' @param ... Arguments passed to [base::sprintf()]. The first argument
#'   should be a format string; subsequent arguments are substituted into
#'   the format string in order.
#'
#' @return Invisibly returns `NULL`.
#' @keywords internal
.optidesign_dbg <- function(debug, ...) {
  if (isTRUE(debug)) message(sprintf(...))
}


# ------------------------------------------------------------------------------
# .optidesign_z
# ------------------------------------------------------------------------------
#' Internal z-score normalisation helper
#'
#' Computes the z-score (standard normal transformation) of a numeric vector,
#' handling edge cases where the standard deviation is zero or non-finite.
#' Used internally when combining or scaling multiple optimality criteria
#' prior to comparison.
#'
#' If the standard deviation is zero or non-finite (e.g. all values identical,
#' or all `NA`), returns a zero vector of the same length rather than
#' producing `NaN` or `Inf`.
#'
#' @param x Numeric vector. `NA` values are ignored when computing the mean
#'   and standard deviation but are returned as `NA` in the output.
#'
#' @return Numeric vector of the same length as `x` containing z-scores,
#'   or a zero vector if the standard deviation is zero or non-finite.
#' @keywords internal
.optidesign_z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}


# ------------------------------------------------------------------------------
# .validate_square_matrix
# ------------------------------------------------------------------------------
#' Internal square matrix validator
#'
#' Checks that an object is a base R matrix with the expected square
#' dimensions. Stops with an informative error if either condition is not met.
#' Used internally to validate relationship matrices (`GRM`, `A`, `K`) before
#' they are used for clustering, dispersion optimisation, or efficiency
#' evaluation.
#'
#' @param M Object to validate. Must be a base R `matrix`.
#' @param p Expected dimension. The matrix must have `nrow(M) == p` and
#'   `ncol(M) == p`.
#' @param nm Character scalar. Object name used in error messages to identify
#'   which argument failed validation. Default `"matrix"`.
#'
#' @return Invisibly returns `TRUE` if validation passes.
#' @keywords internal
.validate_square_matrix <- function(M, p, nm = "matrix") {
  if (!is.matrix(M))
    stop(sprintf("%s must be a matrix.", nm), call. = FALSE)
  if (nrow(M) != p || ncol(M) != p)
    stop(sprintf("%s must be a %d x %d matrix.", nm, p, p), call. = FALSE)
  invisible(TRUE)
}
