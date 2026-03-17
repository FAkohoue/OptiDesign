set.seed(123)

n_lines <- 220
line_ids <- paste0("L", sprintf("%03d", seq_len(n_lines)))
family_ids <- paste0("F", sprintf("%02d", sample(1:24, n_lines, replace = TRUE)))

OptiDesign_lines <- data.frame(
  Treatment = line_ids,
  Family = family_ids,
  stringsAsFactors = FALSE
)

OptiDesign_id_map <- data.frame(
  Treatment = line_ids,
  LineID = line_ids,
  stringsAsFactors = FALSE
)

make_psd_matrix <- function(ids, n_features = 20, diag_scale = 1) {
  X <- matrix(rnorm(length(ids) * n_features), nrow = length(ids), ncol = n_features)
  K <- tcrossprod(scale(X, center = TRUE, scale = TRUE)) / n_features
  diag(K) <- diag(K) + diag_scale
  rownames(K) <- ids
  colnames(K) <- ids
  K
}

OptiDesign_GRM <- make_psd_matrix(line_ids, n_features = 30, diag_scale = 1.0)
OptiDesign_A   <- make_psd_matrix(line_ids, n_features = 15, diag_scale = 1.2)
OptiDesign_K   <- make_psd_matrix(line_ids, n_features = 25, diag_scale = 1.1)

famopt_check_ids <- line_ids[1:4]
famopt_prep_ids  <- line_ids[5:28]
famopt_unrep_ids <- line_ids[29:76]

OptiDesign_famoptg_example <- list(
  check_treatments = famopt_check_ids,
  check_families = OptiDesign_lines$Family[match(famopt_check_ids, OptiDesign_lines$Treatment)],
  p_rep_treatments = famopt_prep_ids,
  p_rep_reps = rep(2L, length(famopt_prep_ids)),
  p_rep_families = OptiDesign_lines$Family[match(famopt_prep_ids, OptiDesign_lines$Treatment)],
  unreplicated_treatments = famopt_unrep_ids,
  unreplicated_families = OptiDesign_lines$Family[match(famopt_unrep_ids, OptiDesign_lines$Treatment)],
  n_blocks = 8L,
  n_rows = 16L,
  n_cols = 8L
)

alpha_check_ids <- line_ids[101:103]
alpha_entry_ids <- line_ids[104:163]

OptiDesign_alpha_example <- list(
  check_treatments = alpha_check_ids,
  check_families = OptiDesign_lines$Family[match(alpha_check_ids, OptiDesign_lines$Treatment)],
  entry_treatments = alpha_entry_ids,
  entry_families = OptiDesign_lines$Family[match(alpha_entry_ids, OptiDesign_lines$Treatment)],
  n_reps = 2L,
  n_rows = 12L,
  n_cols = 14L
)

OptiDesign_famoptg_args_family <- list(
  order = "column",
  serpentine = FALSE,
  seed = 123,
  attempts = 1000,
  warn_and_correct = TRUE,
  fix_rows = TRUE,
  cluster_source = "Family",
  GRM = NULL,
  A = NULL,
  id_map = NULL,
  cluster_method = "kmeans",
  cluster_seed = 1,
  cluster_attempts = 25,
  n_pcs_use = Inf,
  eval_efficiency = FALSE,
  treatment_effect = "random",
  prediction_type = "none",
  K = NULL,
  line_id_map = NULL,
  varcomp = list(
    sigma_e2 = 1,
    sigma_g2 = 1,
    sigma_b2 = 1,
    sigma_r2 = 1,
    sigma_c2 = 1
  ),
  check_as_fixed = TRUE,
  residual_structure = "IID",
  rho_row = 0,
  rho_col = 0,
  spatial_engine = "auto",
  dense_max_n = 5000,
  eff_trace_samples = 80,
  eff_full_max = 400,
  check_placement = "optimal",
  check_opt_attempts = 20,
  use_dispersion = FALSE,
  dispersion_source = "K",
  dispersion_radius = 1,
  dispersion_iters = 200,
  dispersion_seed = 123
)

OptiDesign_famoptg_args_grm <- list(
  order = "column",
  serpentine = FALSE,
  seed = 123,
  attempts = 1000,
  warn_and_correct = TRUE,
  fix_rows = TRUE,
  cluster_source = "GRM",
  GRM = OptiDesign_GRM,
  A = NULL,
  id_map = OptiDesign_id_map,
  cluster_method = "kmeans",
  cluster_seed = 1,
  cluster_attempts = 25,
  n_pcs_use = 5,
  eval_efficiency = FALSE,
  treatment_effect = "random",
  prediction_type = "none",
  K = OptiDesign_K,
  line_id_map = OptiDesign_id_map,
  varcomp = list(
    sigma_e2 = 1,
    sigma_g2 = 1,
    sigma_b2 = 1,
    sigma_r2 = 1,
    sigma_c2 = 1
  ),
  check_as_fixed = TRUE,
  residual_structure = "IID",
  rho_row = 0,
  rho_col = 0,
  spatial_engine = "auto",
  dense_max_n = 5000,
  eff_trace_samples = 40,
  eff_full_max = 400,
  check_placement = "random",
  check_opt_attempts = 20,
  use_dispersion = TRUE,
  dispersion_source = "K",
  dispersion_radius = 1,
  dispersion_iters = 100,
  dispersion_seed = 123
)

OptiDesign_alpha_args_family <- list(
  order = "row",
  serpentine = TRUE,
  seed = 123,
  attempts = 3000,
  warn_and_correct = TRUE,
  fix_rows = TRUE,
  cluster_source = "Family",
  GRM = NULL,
  A = NULL,
  id_map = NULL,
  cluster_method = "kmeans",
  cluster_seed = 1,
  cluster_attempts = 25,
  n_pcs_use = Inf,
  min_entry_slots_per_block = 10,
  max_blocks_per_rep = NULL,
  eval_efficiency = FALSE,
  treatment_effect = "random",
  prediction_type = "none",
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
  residual_structure = "IID",
  rho_row = 0,
  rho_col = 0,
  spatial_engine = "auto",
  dense_max_n = 5000,
  eff_trace_samples = 80,
  eff_full_max = 400,
  check_placement = "random",
  check_position_pattern = "corners_first",
  use_dispersion = FALSE,
  dispersion_source = "K",
  dispersion_radius = 1,
  dispersion_iters = 200,
  dispersion_seed = 123,
  verbose = FALSE
)

OptiDesign_alpha_args_grm <- list(
  order = "row",
  serpentine = TRUE,
  seed = 123,
  attempts = 3000,
  warn_and_correct = TRUE,
  fix_rows = TRUE,
  cluster_source = "GRM",
  GRM = OptiDesign_GRM,
  A = NULL,
  id_map = OptiDesign_id_map,
  cluster_method = "kmeans",
  cluster_seed = 1,
  cluster_attempts = 25,
  n_pcs_use = 5,
  min_entry_slots_per_block = 10,
  max_blocks_per_rep = NULL,
  eval_efficiency = FALSE,
  treatment_effect = "random",
  prediction_type = "none",
  K = OptiDesign_K,
  line_id_map = OptiDesign_id_map,
  varcomp = list(
    sigma_e2 = 1,
    sigma_g2 = 1,
    sigma_rep2 = 1,
    sigma_ib2 = 1,
    sigma_r2 = 1,
    sigma_c2 = 1
  ),
  check_as_fixed = TRUE,
  residual_structure = "IID",
  rho_row = 0,
  rho_col = 0,
  spatial_engine = "auto",
  dense_max_n = 5000,
  eff_trace_samples = 40,
  eff_full_max = 400,
  check_placement = "random",
  check_position_pattern = "corners_first",
  use_dispersion = TRUE,
  dispersion_source = "K",
  dispersion_radius = 1,
  dispersion_iters = 100,
  dispersion_seed = 123,
  verbose = FALSE
)

OptiDesign_example_data <- list(
  OptiDesign_lines = OptiDesign_lines,
  OptiDesign_id_map = OptiDesign_id_map,
  OptiDesign_GRM = OptiDesign_GRM,
  OptiDesign_A = OptiDesign_A,
  OptiDesign_K = OptiDesign_K,
  OptiDesign_famoptg_example = OptiDesign_famoptg_example,
  OptiDesign_alpha_example = OptiDesign_alpha_example,
  OptiDesign_famoptg_args_family = OptiDesign_famoptg_args_family,
  OptiDesign_famoptg_args_grm = OptiDesign_famoptg_args_grm,
  OptiDesign_alpha_args_family = OptiDesign_alpha_args_family,
  OptiDesign_alpha_args_grm = OptiDesign_alpha_args_grm
)

save(
  OptiDesign_example_data,
  file = "data/OptiDesign_example_data.RData",
  compress = "bzip2",
  version = 2
)