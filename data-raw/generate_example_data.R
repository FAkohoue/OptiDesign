# ==============================================================================
# generate_example_data.R
# Generates OptiDesign_example_data and saves it to data/
#
# Run from the package root with:
#   source("data-raw/generate_example_data.R")
#
# Changes from previous version:
#   - eval_efficiency removed from all *_args_* lists (now passed separately
#     to evaluate_famoptg_efficiency() or evaluate_alpha_efficiency())
#   - min_entry_slots_per_block / max_blocks_per_rep replaced with the correct
#     alpha_rc_stream() arguments: min_block_size / max_block_size
#   - treatment_effect / prediction_type / varcomp / residual args removed from
#     *_args_* lists (now passed to the evaluation functions directly)
#   - verbose = FALSE retained in alpha args lists
# ==============================================================================

set.seed(123)

# ==============================================================================
# 1. LINE IDENTIFIERS AND FAMILIES
# ==============================================================================
n_lines  <- 220
line_ids <- paste0("L", sprintf("%03d", seq_len(n_lines)))
family_ids <- paste0("F", sprintf("%02d", sample(1:24, n_lines, replace = TRUE)))

OptiDesign_lines <- data.frame(
  Treatment = line_ids,
  Family    = family_ids,
  stringsAsFactors = FALSE
)

# id_map: identity mapping (Treatment == LineID)
# Used when treatment labels already match relationship matrix rownames.
OptiDesign_id_map <- data.frame(
  Treatment = line_ids,
  LineID    = line_ids,
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. RELATIONSHIP MATRICES
# ==============================================================================
make_psd_matrix <- function(ids, n_features = 20, diag_scale = 1) {
  X <- matrix(rnorm(length(ids) * n_features),
              nrow = length(ids), ncol = n_features)
  K <- tcrossprod(scale(X, center = TRUE, scale = TRUE)) / n_features
  diag(K) <- diag(K) + diag_scale
  rownames(K) <- ids
  colnames(K) <- ids
  K
}

OptiDesign_GRM <- make_psd_matrix(line_ids, n_features = 30, diag_scale = 1.0)
OptiDesign_A   <- make_psd_matrix(line_ids, n_features = 15, diag_scale = 1.2)
OptiDesign_K   <- make_psd_matrix(line_ids, n_features = 25, diag_scale = 1.1)

# ==============================================================================
# 3. TREATMENT VECTORS FOR prep_famoptg()
# ==============================================================================
# 4 checks + 24 p-rep entries (each replicated twice) + 48 unreplicated entries
# Total plots = 8 blocks x 4 checks + 24 x 2 + 48 = 32 + 48 + 48 = 128
# Field = 16 x 8 = 128 plots — exact fit
famopt_check_ids <- line_ids[1:4]
famopt_prep_ids  <- line_ids[5:28]
famopt_unrep_ids <- line_ids[29:76]

OptiDesign_famoptg_example <- list(
  check_treatments        = famopt_check_ids,
  check_families          = OptiDesign_lines$Family[
    match(famopt_check_ids, OptiDesign_lines$Treatment)],
  p_rep_treatments        = famopt_prep_ids,
  p_rep_reps              = rep(2L, length(famopt_prep_ids)),
  p_rep_families          = OptiDesign_lines$Family[
    match(famopt_prep_ids, OptiDesign_lines$Treatment)],
  unreplicated_treatments = famopt_unrep_ids,
  unreplicated_families   = OptiDesign_lines$Family[
    match(famopt_unrep_ids, OptiDesign_lines$Treatment)],
  n_blocks = 8L,
  n_rows   = 16L,
  n_cols   = 8L
)

# ==============================================================================
# 4. TREATMENT VECTORS FOR alpha_rc_stream()
# ==============================================================================
# 3 checks + 60 entries, 2 reps
# Total plots = 2 x (60 + blocks_per_rep x 3)
# With min_block_size = 10, max_block_size = 12:
#   max entry slots = 12 - 3 = 9  -> min blocks = ceil(60/9) = 7
#   min entry slots = 10 - 3 = 7  -> max blocks = floor(60/7) = 8
# Used plots = 2 x (60 + 8 x 3) = 2 x 84 = 168; field = 12 x 14 = 168
alpha_check_ids <- line_ids[101:103]
alpha_entry_ids <- line_ids[104:163]

OptiDesign_alpha_example <- list(
  check_treatments = alpha_check_ids,
  check_families   = OptiDesign_lines$Family[
    match(alpha_check_ids, OptiDesign_lines$Treatment)],
  entry_treatments = alpha_entry_ids,
  entry_families   = OptiDesign_lines$Family[
    match(alpha_entry_ids, OptiDesign_lines$Treatment)],
  n_reps = 2L,
  n_rows = 12L,
  n_cols = 14L
)

# ==============================================================================
# 5. SUPPLEMENTARY ARGS FOR prep_famoptg() — FAMILY CLUSTERING
# ==============================================================================
# Contains only construction arguments. Evaluation arguments
# (treatment_effect, prediction_type, varcomp, residual_structure, etc.)
# are passed directly to evaluate_famoptg_efficiency() at call time.
OptiDesign_famoptg_args_family <- list(
  order              = "column",
  serpentine         = FALSE,
  seed               = 123L,
  attempts           = 1000L,
  warn_and_correct   = TRUE,
  fix_rows           = TRUE,
  cluster_source     = "Family",
  GRM                = NULL,
  A                  = NULL,
  id_map             = NULL,
  cluster_method     = "kmeans",
  cluster_seed       = 1L,
  cluster_attempts   = 25L,
  n_pcs_use          = Inf,
  check_placement    = "optimal",
  check_opt_attempts = 20L,
  use_dispersion     = FALSE,
  dispersion_source  = "K",
  dispersion_radius  = 1L,
  dispersion_iters   = 200L,
  dispersion_seed    = 123L,
  K                  = NULL,
  line_id_map        = NULL,
  verbose            = FALSE
)

# ==============================================================================
# 6. SUPPLEMENTARY ARGS FOR prep_famoptg() — GRM CLUSTERING + DISPERSION
# ==============================================================================
OptiDesign_famoptg_args_grm <- list(
  order              = "column",
  serpentine         = FALSE,
  seed               = 123L,
  attempts           = 1000L,
  warn_and_correct   = TRUE,
  fix_rows           = TRUE,
  cluster_source     = "GRM",
  GRM                = OptiDesign_GRM,
  A                  = NULL,
  id_map             = OptiDesign_id_map,
  cluster_method     = "kmeans",
  cluster_seed       = 1L,
  cluster_attempts   = 25L,
  n_pcs_use          = 5L,
  check_placement    = "random",
  check_opt_attempts = 20L,
  use_dispersion     = TRUE,
  dispersion_source  = "K",
  dispersion_radius  = 1L,
  dispersion_iters   = 100L,
  dispersion_seed    = 123L,
  K                  = OptiDesign_K,
  line_id_map        = OptiDesign_id_map,
  verbose            = FALSE
)

# ==============================================================================
# 7. SUPPLEMENTARY ARGS FOR alpha_rc_stream() — FAMILY CLUSTERING
# ==============================================================================
# min_block_size / max_block_size replace the stale min_entry_slots_per_block /
# max_blocks_per_rep arguments that were in the previous version.
# eval_efficiency removed — evaluation is now done via evaluate_alpha_efficiency().
OptiDesign_alpha_args_family <- list(
  order                  = "row",
  serpentine             = TRUE,
  seed                   = 123L,
  attempts               = 3000L,
  warn_and_correct       = TRUE,
  fix_rows               = TRUE,
  cluster_source         = "Family",
  GRM                    = NULL,
  A                      = NULL,
  id_map                 = NULL,
  cluster_method         = "kmeans",
  cluster_seed           = 1L,
  cluster_attempts       = 25L,
  n_pcs_use              = Inf,
  min_block_size         = 10L,
  max_block_size         = 12L,
  check_placement        = "random",
  check_position_pattern = "spread",
  use_dispersion         = FALSE,
  dispersion_source      = "K",
  dispersion_radius      = 1L,
  dispersion_iters       = 200L,
  dispersion_seed        = 123L,
  K                      = NULL,
  line_id_map            = NULL,
  verbose                = FALSE
)

# ==============================================================================
# 8. SUPPLEMENTARY ARGS FOR alpha_rc_stream() — GRM CLUSTERING + DISPERSION
# ==============================================================================
OptiDesign_alpha_args_grm <- list(
  order                  = "row",
  serpentine             = TRUE,
  seed                   = 123L,
  attempts               = 3000L,
  warn_and_correct       = TRUE,
  fix_rows               = TRUE,
  cluster_source         = "GRM",
  GRM                    = OptiDesign_GRM,
  A                      = NULL,
  id_map                 = OptiDesign_id_map,
  cluster_method         = "kmeans",
  cluster_seed           = 1L,
  cluster_attempts       = 25L,
  n_pcs_use              = 5L,
  min_block_size         = 10L,
  max_block_size         = 12L,
  check_placement        = "random",
  check_position_pattern = "spread",
  use_dispersion         = TRUE,
  dispersion_source      = "K",
  dispersion_radius      = 1L,
  dispersion_iters       = 100L,
  dispersion_seed        = 123L,
  K                      = OptiDesign_K,
  line_id_map            = OptiDesign_id_map,
  verbose                = FALSE
)

# ==============================================================================
# 9. ASSEMBLE AND SAVE
# ==============================================================================
OptiDesign_example_data <- list(
  OptiDesign_lines               = OptiDesign_lines,
  OptiDesign_id_map              = OptiDesign_id_map,
  OptiDesign_GRM                 = OptiDesign_GRM,
  OptiDesign_A                   = OptiDesign_A,
  OptiDesign_K                   = OptiDesign_K,
  OptiDesign_famoptg_example     = OptiDesign_famoptg_example,
  OptiDesign_alpha_example       = OptiDesign_alpha_example,
  OptiDesign_famoptg_args_family = OptiDesign_famoptg_args_family,
  OptiDesign_famoptg_args_grm    = OptiDesign_famoptg_args_grm,
  OptiDesign_alpha_args_family   = OptiDesign_alpha_args_family,
  OptiDesign_alpha_args_grm      = OptiDesign_alpha_args_grm
)

save(
  OptiDesign_example_data,
  file     = "data/OptiDesign_example_data.RData",
  compress = "bzip2",
  version  = 2
)

message("OptiDesign_example_data saved to data/OptiDesign_example_data.RData")