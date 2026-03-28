# ── Shared fixture ────────────────────────────────────────────────────────────
base_opt_args <- function() {
  list(
    check_treatments = c("CHK1", "CHK2", "CHK3"),
    check_families   = c("CHECK", "CHECK", "CHECK"),
    entry_treatments = paste0("E", 1:30),
    entry_families   = rep(paste0("F", 1:5), 6L),
    n_reps           = 2L,
    n_rows           = 12L,
    n_cols           = 10L,
    min_block_size   = 8L,
    max_block_size   = 10L,
    treatment_effect = "fixed",
    criterion        = "A",
    n_restarts       = 5L,
    verbose_opt      = FALSE
  )
}

# ── Return value structure ─────────────────────────────────────────────────────

test_that("optimize_alpha_rc returns correct list structure", {
  out <- do.call(optimize_alpha_rc, base_opt_args())

  expect_type(out, "list")
  expect_true(is.matrix(out$layout_matrix))
  expect_true(is.data.frame(out$field_book))
  expect_true(is.list(out$efficiency))
  expect_true(is.list(out$optimization))
})

test_that("optimization metadata has required fields", {
  out <- do.call(optimize_alpha_rc, base_opt_args())
  opt <- out$optimization

  expect_true(opt$method %in% c("RS", "SA", "GA"))
  expect_equal(opt$criterion, "A")
  expect_true(is.numeric(opt$best_score))
  expect_true(is.numeric(opt$master_seed))
  expect_true(is.integer(opt$n_failed) || is.numeric(opt$n_failed))
})

test_that("score_history length equals n_restarts for RS", {
  args           <- base_opt_args()
  args$method    <- "RS"
  args$n_restarts <- 7L
  out <- do.call(optimize_alpha_rc, args)
  expect_length(out$optimization$score_history, 7L)
})

# ── RS method ────────────────────────────────────────────────────────────────

test_that("RS: best_score equals minimum of score_history", {
  args        <- base_opt_args()
  args$method <- "RS"
  out    <- do.call(optimize_alpha_rc, args)
  scores <- out$optimization$score_history
  valid  <- scores[!is.na(scores)]
  if (length(valid) > 0)
    expect_equal(out$optimization$best_score, min(valid), tolerance = 1e-10)
})

# ── SA method ────────────────────────────────────────────────────────────────

test_that("SA method runs without error", {
  args              <- base_opt_args()
  args$method       <- "SA"
  args$n_restarts   <- 2L
  args$sa_max_iter  <- 10L
  args$sa_temp_start <- 1.0
  args$sa_temp_end  <- 0.01
  expect_no_error(do.call(optimize_alpha_rc, args))
})

test_that("SA: optimization list contains sa_params", {
  args              <- base_opt_args()
  args$method       <- "SA"
  args$n_restarts   <- 2L
  args$sa_max_iter  <- 5L
  args$sa_temp_start <- 1.0
  args$sa_temp_end  <- 0.01
  out <- do.call(optimize_alpha_rc, args)
  expect_true(is.list(out$optimization$sa_params))
  expect_null(out$optimization$ga_params)
})

test_that("SA: score_history is a list of length n_restarts", {
  args              <- base_opt_args()
  args$method       <- "SA"
  args$n_restarts   <- 3L
  args$sa_max_iter  <- 5L
  args$sa_temp_start <- 1.0
  args$sa_temp_end  <- 0.01
  out <- do.call(optimize_alpha_rc, args)
  expect_true(is.list(out$optimization$score_history))
  expect_length(out$optimization$score_history, 3L)
})

# ── GA method ────────────────────────────────────────────────────────────────

test_that("GA method runs without error", {
  args                  <- base_opt_args()
  args$method           <- "GA"
  args$ga_pop_size      <- 4L
  args$ga_n_generations <- 3L
  expect_no_error(do.call(optimize_alpha_rc, args))
})

test_that("GA: optimization list contains ga_params", {
  args                  <- base_opt_args()
  args$method           <- "GA"
  args$ga_pop_size      <- 4L
  args$ga_n_generations <- 3L
  out <- do.call(optimize_alpha_rc, args)
  expect_true(is.list(out$optimization$ga_params))
  expect_null(out$optimization$sa_params)
})

test_that("GA: score_history has length equal to n_generations", {
  args                  <- base_opt_args()
  args$method           <- "GA"
  args$ga_pop_size      <- 4L
  args$ga_n_generations <- 5L
  out <- do.call(optimize_alpha_rc, args)
  expect_length(out$optimization$score_history, 5L)
})

# ── Criterion directions ──────────────────────────────────────────────────────

test_that("criterion D returns finite best_score", {
  args           <- base_opt_args()
  args$criterion <- "D"
  out <- do.call(optimize_alpha_rc, args)
  expect_true(is.finite(out$optimization$best_score))
})

test_that("criterion both returns finite best_score", {
  args           <- base_opt_args()
  args$criterion <- "both"
  out <- do.call(optimize_alpha_rc, args)
  expect_true(is.finite(out$optimization$best_score))
})

test_that("CDmean requires random treatment effect", {
  args                  <- base_opt_args()
  args$criterion        <- "CDmean"
  args$treatment_effect <- "fixed"
  expect_error(do.call(optimize_alpha_rc, args), "CDmean")
})

test_that("CDmean: best_score is positive and score_history values are positive", {
  args                  <- base_opt_args()
  args$criterion        <- "CDmean"
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  out    <- do.call(optimize_alpha_rc, args)
  scores <- out$optimization$score_history
  valid  <- scores[!is.na(scores)]

  expect_true(out$optimization$best_score > 0)
  if (length(valid) > 0) expect_true(all(valid > 0))
})

# ── Constraint preservation ───────────────────────────────────────────────────

test_that("optimised design: each entry appears once per rep", {
  out <- do.call(optimize_alpha_rc, base_opt_args())
  fb  <- out$field_book[!out$field_book$Check & !is.na(out$field_book$Treatment), ]

  entry_rep_counts <- table(fb$Rep, fb$Treatment)
  expect_true(all(entry_rep_counts <= 1L))
})

test_that("optimised design: checks appear in every block", {
  args_list <- base_opt_args()
  out  <- do.call(optimize_alpha_rc, args_list)
  fb   <- out$field_book[!is.na(out$field_book$Treatment), ]

  n_blocks <- max(fb$IBlock, na.rm = TRUE)
  for (chk in args_list$check_treatments) {
    blocks_with_chk <- unique(fb$IBlock[fb$Treatment == chk])
    expect_setequal(blocks_with_chk, seq_len(n_blocks))
  }
})

# ── Pre-flight guards ─────────────────────────────────────────────────────────

test_that("stops when random + prediction_type = none", {
  args                  <- base_opt_args()
  args$treatment_effect <- "random"
  args$prediction_type  <- "none"
  expect_error(do.call(optimize_alpha_rc, args))
})

test_that("stops when n_restarts < 1", {
  args            <- base_opt_args()
  args$n_restarts <- 0L
  expect_error(do.call(optimize_alpha_rc, args))
})

test_that("SA stops when sa_temp_end >= sa_temp_start", {
  args               <- base_opt_args()
  args$method        <- "SA"
  args$sa_temp_start <- 0.5
  args$sa_temp_end   <- 1.0
  expect_error(do.call(optimize_alpha_rc, args))
})

test_that("GA stops when ga_pop_size < 4", {
  args             <- base_opt_args()
  args$method      <- "GA"
  args$ga_pop_size <- 3L
  expect_error(do.call(optimize_alpha_rc, args))
})
