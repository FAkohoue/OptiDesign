# ── Shared fixture ────────────────────────────────────────────────────────────
make_alpha_design <- function(seed = 1L) {
  alpha_rc_stream(
    check_treatments = c("CHK1", "CHK2", "CHK3"),
    check_families   = c("CHECK", "CHECK", "CHECK"),
    entry_treatments = paste0("E", 1:30),
    entry_families   = rep(paste0("F", 1:5), 6L),
    n_reps           = 2L,
    n_rows           = 12L,
    n_cols           = 10L,
    min_block_size   = 8L,
    max_block_size   = 10L,
    seed             = seed
  )
}

base_eval_args <- function(design) {
  list(
    field_book       = design$field_book,
    n_rows           = 12L,
    n_cols           = 10L,
    check_treatments = c("CHK1", "CHK2", "CHK3"),
    treatment_effect = "fixed"
  )
}

# ── Return value structure ─────────────────────────────────────────────────────

test_that("evaluate_alpha_efficiency returns a named list", {
  design <- make_alpha_design()
  eff    <- do.call(evaluate_alpha_efficiency, base_eval_args(design))

  expect_type(eff, "list")
  expect_true(all(c("model", "treatment_effect", "mode",
                    "A_criterion", "A_efficiency", "n_trt") %in% names(eff)))
})

test_that("model string references Rep and IBlock structure", {
  design <- make_alpha_design()
  eff    <- do.call(evaluate_alpha_efficiency, base_eval_args(design))
  expect_true(grepl("IBlock", eff$model))
  expect_true(grepl("Rep",    eff$model))
})

test_that("fixed effects: A and D criteria are positive", {
  design <- make_alpha_design()
  eff    <- do.call(evaluate_alpha_efficiency, base_eval_args(design))
  expect_true(eff$A_criterion > 0)
  expect_true(eff$D_criterion > 0)
})

test_that("A_efficiency = 1 / A_criterion", {
  design <- make_alpha_design()
  eff    <- do.call(evaluate_alpha_efficiency, base_eval_args(design))
  expect_equal(eff$A_efficiency, 1 / eff$A_criterion, tolerance = 1e-10)
})

test_that("D_efficiency = 1 / D_criterion", {
  design <- make_alpha_design()
  eff    <- do.call(evaluate_alpha_efficiency, base_eval_args(design))
  expect_equal(eff$D_efficiency, 1 / eff$D_criterion, tolerance = 1e-10)
})

test_that("backward-compatible aliases A and D are correct", {
  design <- make_alpha_design()
  eff    <- do.call(evaluate_alpha_efficiency, base_eval_args(design))
  expect_equal(eff$A, eff$A_efficiency)
  expect_equal(eff$D, eff$D_efficiency)
})

# ── Residual structures ───────────────────────────────────────────────────────

test_that("IID residual structure runs without error", {
  design <- make_alpha_design()
  args   <- base_eval_args(design)
  args$residual_structure <- "IID"
  expect_no_error(do.call(evaluate_alpha_efficiency, args))
})

test_that("AR1 residual structure runs without error", {
  design <- make_alpha_design()
  args   <- base_eval_args(design)
  args$residual_structure <- "AR1"
  args$rho_row            <- 0.3
  expect_no_error(do.call(evaluate_alpha_efficiency, args))
})

test_that("AR1xAR1 residual structure runs without error", {
  design <- make_alpha_design()
  args   <- base_eval_args(design)
  args$residual_structure <- "AR1xAR1"
  args$rho_row            <- 0.3
  args$rho_col            <- 0.2
  expect_no_error(do.call(evaluate_alpha_efficiency, args))
})

test_that("spatial correlation changes criterion values vs IID", {
  design   <- make_alpha_design()
  eff_iid  <- do.call(evaluate_alpha_efficiency, base_eval_args(design))
  args_ar1 <- base_eval_args(design)
  args_ar1$residual_structure <- "AR1xAR1"
  args_ar1$rho_row            <- 0.5
  args_ar1$rho_col            <- 0.5
  eff_ar1  <- do.call(evaluate_alpha_efficiency, args_ar1)
  expect_false(isTRUE(all.equal(eff_iid$A_criterion, eff_ar1$A_criterion)))
})

# ── varcomp difference from famoptg ──────────────────────────────────────────

test_that("varcomp uses sigma_rep2 and sigma_ib2 (not sigma_b2)", {
  design <- make_alpha_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  # Correct varcomp for alpha designs
  args$varcomp <- list(sigma_e2 = 1, sigma_g2 = 1,
                       sigma_rep2 = 1, sigma_ib2 = 1,
                       sigma_r2 = 1, sigma_c2 = 1)
  expect_no_error(do.call(evaluate_alpha_efficiency, args))
})

# ── Random effects and CDmean ─────────────────────────────────────────────────

test_that("random IID returns PEV_criterion and CDmean", {
  design <- make_alpha_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  eff <- do.call(evaluate_alpha_efficiency, args)

  expect_true(is.numeric(eff$PEV_criterion))
  expect_true(eff$PEV_criterion > 0)
  expect_true(is.numeric(eff$CDmean))
  expect_true(eff$CDmean >= 0 && eff$CDmean <= 1)
})

test_that("CDmean = 1 - mean_PEV / sigma_g2", {
  design <- make_alpha_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  args$varcomp <- list(sigma_e2 = 1, sigma_g2 = 3,
                       sigma_rep2 = 1, sigma_ib2 = 1,
                       sigma_r2 = 1, sigma_c2 = 1)
  eff <- do.call(evaluate_alpha_efficiency, args)
  expect_equal(eff$CDmean, 1 - eff$mean_PEV / 3, tolerance = 1e-10)
})

test_that("random: D_criterion is NA", {
  design <- make_alpha_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  eff <- do.call(evaluate_alpha_efficiency, args)
  expect_true(is.na(eff$D_criterion))
})

# ── Pre-flight guards ─────────────────────────────────────────────────────────

test_that("stops when treatment_effect = random and prediction_type = none", {
  design <- make_alpha_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "none"
  expect_error(do.call(evaluate_alpha_efficiency, args))
})

test_that("stops when GBLUP requested but K is NULL", {
  design <- make_alpha_design()
  args   <- base_eval_args(design)
  args$treatment_effect <- "random"
  args$prediction_type  <- "GBLUP"
  args$K                <- NULL
  expect_error(do.call(evaluate_alpha_efficiency, args), "K must be provided")
})

test_that("stops when field_book is missing required columns", {
  design <- make_alpha_design()
  args   <- base_eval_args(design)
  args$field_book <- args$field_book[, setdiff(names(args$field_book), "IBlock")]
  expect_error(do.call(evaluate_alpha_efficiency, args), "missing columns")
})
