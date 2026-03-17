#' Example data for OptiDesign
#'
#' Synthetic example data shipped with `OptiDesign` for illustrating and testing
#' the main design constructors, relationship-matrix workflows, optional
#' dispersion optimization, and optional efficiency settings.
#'
#' @name OptiDesign_example_data
#' @docType data
#' @keywords datasets
#'
#' @format A named list with the following components:
#' \describe{
#'   \item{`OptiDesign_lines`}{A data frame with columns `Treatment` and `Family`
#'   representing a synthetic pool of line identifiers and their family labels.}
#'
#'   \item{`OptiDesign_id_map`}{A data frame with columns `Treatment` and `LineID`
#'   used to map treatment labels to row and column names of relationship matrices
#'   when those names differ.}
#'
#'   \item{`OptiDesign_GRM`}{A synthetic genomic relationship matrix with rownames
#'   and colnames corresponding to line IDs. It can be used for matrix-based
#'   grouping or genomic dispersion examples.}
#'
#'   \item{`OptiDesign_A`}{A synthetic pedigree-style relationship matrix with
#'   rownames and colnames corresponding to line IDs. It can be used for
#'   pedigree-based grouping or dispersion examples.}
#'
#'   \item{`OptiDesign_K`}{A synthetic relationship or kernel matrix with rownames
#'   and colnames corresponding to line IDs. It is intended for workflows that
#'   require a `K` matrix, such as BLUP-style prediction settings or dispersion
#'   optimization using `dispersion_source = "K"`.}
#'
#'   \item{`OptiDesign_famoptg_example`}{A named list containing the core treatment
#'   vectors and field dimensions needed to run `prep_famoptg()`. This includes
#'   checks, p-rep treatments, unreplicated treatments, and basic field dimensions.}
#'
#'   \item{`OptiDesign_alpha_example`}{A named list containing the core treatment
#'   vectors and field dimensions needed to run `alpha_rc_stream()`.
#'   This includes checks, entries, number of replicates, and fixed field dimensions.}
#'
#'   \item{`OptiDesign_famoptg_args_family`}{A named list of additional arguments
#'   for a family-based call to `prep_famoptg()`. It illustrates a simple design
#'   setup without efficiency evaluation.}
#'
#'   \item{`OptiDesign_famoptg_args_grm`}{A named list of additional arguments for
#'   a GRM-based call to `prep_famoptg()`. It illustrates matrix-based grouping and
#'   optional dispersion settings.}
#'
#'   \item{`OptiDesign_alpha_args_family`}{A named list of additional arguments
#'   for a family-based call to `alpha_rc_stream()`. It illustrates
#'   a stream-based row-column design without efficiency evaluation.}
#'
#'   \item{`OptiDesign_alpha_args_grm`}{A named list of additional arguments for
#'   a GRM-based call to `alpha_rc_stream()`. It illustrates
#'   matrix-based grouping and optional dispersion settings.}
#' }
#'
#' @details
#' The objects in this dataset are synthetic and are intended for:
#'
#' - package examples,
#' - tests,
#' - workflow demonstrations,
#' - and user exploration of argument combinations.
#'
#' The dataset is especially useful for understanding how the package handles:
#'
#' - repeated checks,
#' - p-rep and unreplicated materials,
#' - fixed-grid stream designs,
#' - family-based grouping,
#' - GRM- or pedigree-based grouping,
#' - and matrix-based dispersion settings.
#'
#' Load the dataset with:
#'
#' `data("OptiDesign_example_data", package = "OptiDesign")`
#'
#' Then extract individual components as needed, for example:
#'
#' `x <- OptiDesign_example_data`
#'
#' `x$OptiDesign_famoptg_example`
#'
#' `x$OptiDesign_GRM`
#'
#' `x$OptiDesign_alpha_args_family`
#'
#' @source Generated internally by `data-raw/generate_example_data.R`.
NULL