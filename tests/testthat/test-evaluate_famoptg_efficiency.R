# ── Shared fixture ────────────────────────────────────────────────────────────
make_famoptg_design <- function(seed = 1L) {
  prep_famoptg(
    check_treatments        = c("CHK1", "CHK2"),
    check_families          = c("CHECK", "CHECK"),
    p_rep_treatments        = paste0("P", 1:10),
    p_rep_reps              = rep(2L, 10L),
    p_rep_families          = rep(c("F1", "F2"), 5L),
    unreplicated_treatments = paste0("U", 1:20),
    unreplicated_families   = rep(c("F1", "F2", "F3", "F4"), 5L),
    n_blocks = 4L, n_rows = 6L,  n_cols = 8L,
    seed = seed
  )
}

base_eval_args <- function(design) {
  list(
    field_book       = design$field_book,
    n_rows           = 6L,
    n_cols           = 8L,
    check_treatments = c("CHK1", "CHK2"),
    treatment_effect = "fixed"
  )
}

# ── Return value structure ─────────────────────────────────────────────────────

test_that("evaluate_famoptg_efficiency returns a named list", {
  design <- make_famoptg_design()
  eff    <- do.call(evaluate_famoptg_efficiency, base_eval_args(design))

  expect_type(eff, "list")
  expect_true(all(c("model", "treatment_effect", "mode",
                    "A_criterion", "A_efficiency", "n_trt") %in% names(eff)))
})

test_that("fixed effects: A_criterion and D_criterion are positive numerics", {
  design <- make_famoptg_design()
  eff    <- do.call(evaluate_famoptg_efficiency, base_eval_args(design))

  expect_true(is.numeric(eff$A_criterion))
  expect_true(eff$A_criterion > 0)
  expect_true(is.numeric(eff$D_criterion))
  expect_true(eff$D_criterion > 0)
})

test_that("A_efficiency = 1 / A_criterion", {
  design <- make_famoptg_design()
  eff    <- do.call(evaluate_famoptg_efficiency, base_eval_args(design))
  expect_equal(eff$A_efficiency, 1 / eff$A_criterion, tolerance = 1e-10)
})

test_that("D_efficiency = 1 / D_criterion", {
  design <- make_famoptg_design()
  eff    <- do.call(evaluate_famoptg_efficiency, base_eval_args(design))
  expect_equal(eff$D_efficiency, 1 / eff$D_criterion, tolerance = 1e-10)
})

test_that("backward-compatible aliases A and D match efficiency forms", {
  design <- make_famoptg_design()
  eff    <- do.call(evaluate_famoptg_efficiency, base_eval_args(design))
  expect_equal(eff$A, eff$A_efficiency)
  expect_equal(eff$D, eff$D_efficiency)
})

# ── Residual structures ───────────────────────────────────────────────────────

test_that("IID residual structure runs without error", {
  design <- make_famoptg_design()
  args   <- base_eval_args(design)
  args$residual_structure <- "IID"
  expect_no_error(do.call(evaluate_famoptg_efficiency, args))
})

test_that("AR1 residual structure runs without error", {
  design <- make_famoptg_design()
  args   <- base_eval_args(design)
  args$residual_structure <- "AR1"
  args$rho_row            <- 0.3
  expect_no_error(do.call(evaluate_famoptg_efficiency, args))
})

test_that("AR1xAR1 residual structure runs without error", {
  design <- make_famoptg_design()
  args   <- base_eval_args(design)
  args$residual_structure <- "AR1xAR1"
  args$rho_row            <- 0.3
  args$rho_col            <- 0.2
  expect_no_error(do.call(evaluate_famoptg_efficiency, args))
})

test_that("spatial correlation changes criterion values", {
  design   <- make_famoptg_design()
  eff_iid  <- do.call(evaluate_famoptg_efficiency, base_eval_args(design))
  args_ar1 <- base_eval_args(design)
  args_ar1$residual_structure <- "AR1xAR1"
  args_ar1$rho_row            <- 0.5
  args_ar1$rho_col            <- 0.5
  eff_ar1  <- do.call(evaluate_famoptg_efficiency, args_ar1)
  expect_false(isTRUE(all.equal(eff_iid$A_criterion, eff_ar1$A_criterion)))
})

# ── Random effects and CDmean ─────────────────────────────────────────────────

test_that("random effects IID returns PEV_criterion and CDmean", {
  design <- make_famoptg_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  eff <- do.call(evaluate_famoptg_efficiency, args)

  expect_true(is.numeric(eff$PEV_criterion))
  expect_true(eff$PEV_criterion > 0)
  expect_true(is.numeric(eff$CDmean))
  expect_true(eff$CDmean >= 0 && eff$CDmean <= 1)
})

test_that("CDmean = 1 - mean_PEV / sigma_g2", {
  design <- make_famoptg_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  args$varcomp          <- list(sigma_e2 = 1, sigma_g2 = 2,
                                sigma_b2 = 1, sigma_r2 = 1, sigma_c2 = 1)
  eff <- do.call(evaluate_famoptg_efficiency, args)
  expect_equal(eff$CDmean, 1 - eff$mean_PEV / 2, tolerance = 1e-10)
})

test_that("random effects: D_criterion is NA", {
  design <- make_famoptg_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  eff <- do.call(evaluate_famoptg_efficiency, args)
  expect_true(is.na(eff$D_criterion))
})

test_that("CD_per_line has length equal to n_trt in full mode", {
  design <- make_famoptg_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  args$eff_full_max     <- 1000L
  eff <- do.call(evaluate_famoptg_efficiency, args)
  expect_length(eff$CD_per_line, eff$n_trt)
})

# ── Pre-flight guards ─────────────────────────────────────────────────────────

test_that("stops when treatment_effect = random and prediction_type = none", {
  design <- make_famoptg_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "none"
  expect_error(do.call(evaluate_famoptg_efficiency, args))
})

test_that("stops when GBLUP requested but K is NULL", {
  design <- make_famoptg_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "GBLUP"
  args$K                <- NULL
  expect_error(do.call(evaluate_famoptg_efficiency, args), "K must be provided")
})

test_that("stops when field_book is missing required columns", {
  design <- make_famoptg_design()
  args   <- base_eval_args(design)
  args$field_book <- args$field_book[, setdiff(names(args$field_book), "Block")]
  expect_error(do.call(evaluate_famoptg_efficiency, args), "missing columns")
})

test_that("stops when varcomp is missing sigma_b2", {
  design <- make_famoptg_design()
  args   <- base_eval_args(design)
  args$varcomp <- list(sigma_e2 = 1, sigma_g2 = 1,
                       sigma_r2 = 1, sigma_c2 = 1)   # sigma_b2 missing
  expect_error(do.call(evaluate_famoptg_efficiency, args))
})
