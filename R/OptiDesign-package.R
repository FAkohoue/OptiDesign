#' OptiDesign: Experimental Field Design Utilities for Optimized Layout Construction
#'
#' `OptiDesign` provides tools for constructing optimized experimental field
#' designs for plant breeding and related applications. The package focuses on
#' layouts where field structure, treatment grouping, and optional mixed-model
#' evaluation should be integrated into the design stage rather than treated as
#' separate downstream tasks.
#'
#' @details
#' The package currently includes two main design constructors:
#'
#' - `prep_famoptg()` for partially replicated (p-rep) designs with repeated
#'   checks, optional p-rep and unreplicated entries, optional grouping from
#'   family or relationship matrices, optional dispersion optimization, and
#'   optional efficiency evaluation.
#'
#' - `prep_alpha_checks_rc_stream()` for fixed-grid alpha row-column stream
#'   designs where the field is treated as a global traversal stream, replicates
#'   are contiguous stream segments, incomplete blocks may be unequal, checks are
#'   repeated in every block, and unused cells appear only at the end of the
#'   stream.
#'
#' The package is intended for users who need designs that are not only feasible
#' in the field, but also informed by:
#'
#' - family structure,
#' - pedigree or genomic relationship matrices,
#' - spacing of related materials,
#' - and mixed-model design efficiency considerations.
#'
#' @section Main objectives of the package:
#' `OptiDesign` is particularly useful when the user wants to:
#'
#' - construct designs with repeated checks and heterogeneous entry replication,
#' - avoid local clustering of related entries,
#' - use family, pedigree, or genomic structure to guide placement,
#' - preserve realistic field-book order,
#' - evaluate a design under fixed- or random-effect mixed models,
#' - compare alternative layout strategies before field implementation.
#'
#' @section Relationship matrices and grouping:
#' Several design routines in the package can use one of three grouping sources:
#'
#' - direct family labels,
#' - a genomic relationship matrix (`GRM`),
#' - a pedigree relationship matrix (`A`).
#'
#' In addition, some workflows may use a separate relationship matrix `K` for:
#'
#' - BLUP-style prediction efficiency calculations,
#' - or dispersion optimization.
#'
#' When treatment names do not match matrix row names directly, mapping tables
#' such as `id_map` or `line_id_map` can be used to reconcile identifiers.
#'
#' @section Included example dataset:
#' The package ships with the dataset `OptiDesign_example_data`, which contains
#' synthetic example objects for:
#'
#' - treatment and family lists,
#' - example relationship matrices,
#' - example argument lists for both main design functions.
#'
#' Load it with:
#'
#' `data("OptiDesign_example_data", package = "OptiDesign")`
#'
#' @section Imported functionality:
#' The package imports `mod()` from **pracma** for serpentine parity handling in
#' layouts where alternating row or column direction is required.
#'
#' Matrix computations used in efficiency evaluation and some optimization steps
#' rely on the **Matrix** package.
#'
#' @keywords internal
"_PACKAGE"