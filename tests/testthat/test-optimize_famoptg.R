# ── Shared fixture ────────────────────────────────────────────────────────────
base_opt_args <- function() {
  list(
    check_treatments        = c("CHK1", "CHK2"),
    check_families          = c("CHECK", "CHECK"),
    p_rep_treatments        = paste0("P", 1:10),
    p_rep_reps              = rep(2L, 10L),
    p_rep_families          = rep(c("F1", "F2"), 5L),
    unreplicated_treatments = paste0("U", 1:20),
    unreplicated_families   = rep(c("F1", "F2", "F3", "F4"), 5L),
    n_blocks         = 4L,
    n_rows           = 6L,
    n_cols           = 8L,
    treatment_effect = "fixed",
    criterion        = "A",
    n_restarts       = 5L,    # small for speed
    verbose_opt      = FALSE
  )
}

# ── Return value structure ─────────────────────────────────────────────────────

test_that("optimize_famoptg returns correct list structure", {
  out <- do.call(optimize_famoptg, base_opt_args())

  expect_type(out, "list")
  expect_true(is.matrix(out$layout_matrix))
  expect_true(is.data.frame(out$field_book))
  expect_true(is.list(out$efficiency))
  expect_true(is.list(out$optimization))
})

test_that("optimization metadata has required fields", {
  out <- do.call(optimize_famoptg, base_opt_args())
  opt <- out$optimization

  expect_equal(opt$method, "RS")
  expect_equal(opt$criterion, "A")
  expect_length(opt$score_history, 5L)
  expect_true(is.numeric(opt$best_score))
  expect_true(is.integer(opt$n_restarts) || is.numeric(opt$n_restarts))
  expect_true(is.integer(opt$n_failed)   || is.numeric(opt$n_failed))
  expect_true(is.numeric(opt$master_seed))
})

test_that("score_history has length equal to n_restarts", {
  args <- base_opt_args()
  args$n_restarts <- 8L
  out  <- do.call(optimize_famoptg, args)
  expect_length(out$optimization$score_history, 8L)
})

# ── Criterion directions ──────────────────────────────────────────────────────

test_that("best_score for criterion A is the minimum of score_history", {
  out    <- do.call(optimize_famoptg, base_opt_args())
  scores <- out$optimization$score_history
  valid  <- scores[!is.na(scores)]
  if (length(valid) > 0)
    expect_equal(out$optimization$best_score, min(valid), tolerance = 1e-10)
})

test_that("criterion D runs and returns finite best_score", {
  args           <- base_opt_args()
  args$criterion <- "D"
  out <- do.call(optimize_famoptg, args)
  expect_true(is.finite(out$optimization$best_score))
})

test_that("criterion CDmean requires random treatment effect", {
  args                   <- base_opt_args()
  args$criterion         <- "CDmean"
  args$treatment_effect  <- "fixed"
  expect_error(do.call(optimize_famoptg, args), "CDmean")
})

test_that("CDmean criterion: best_score is positive and <= 1", {
  args                  <- base_opt_args()
  args$criterion        <- "CDmean"
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  out <- do.call(optimize_famoptg, args)
  expect_true(out$optimization$best_score > 0)
  expect_true(out$optimization$best_score <= 1)
})

test_that("CDmean score_history values are positive", {
  args                  <- base_opt_args()
  args$criterion        <- "CDmean"
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  out    <- do.call(optimize_famoptg, args)
  scores <- out$optimization$score_history
  valid  <- scores[!is.na(scores)]
  expect_true(all(valid > 0))
})

# ── Constraint preservation ───────────────────────────────────────────────────

test_that("optimised design satisfies p-rep constraint", {
  out <- do.call(optimize_famoptg, base_opt_args())
  fb  <- out$field_book
  args <- base_opt_args()

  for (trt in args$p_rep_treatments) {
    trt_blocks <- fb$Block[fb$Treatment == trt]
    expect_equal(length(trt_blocks), length(unique(trt_blocks)),
                 label = paste0("p-rep constraint for ", trt))
  }
})

test_that("checks appear in every block of optimised design", {
  args_list <- base_opt_args()
  out  <- do.call(optimize_famoptg, args_list)
  fb   <- out$field_book

  for (chk in args_list$check_treatments) {
    blocks_with_chk <- unique(fb$Block[fb$Treatment == chk])
    expect_setequal(blocks_with_chk, seq_len(args_list$n_blocks))
  }
})

test_that("unreplicated treatments appear exactly once in optimised design", {
  args_list <- base_opt_args()
  out  <- do.call(optimize_famoptg, args_list)
  fb   <- out$field_book

  for (trt in args_list$unreplicated_treatments) {
    expect_equal(sum(fb$Treatment == trt), 1L)
  }
})

# ── Pre-flight guards ─────────────────────────────────────────────────────────

test_that("stops when random + prediction_type = none", {
  args                  <- base_opt_args()
  args$treatment_effect <- "random"
  args$prediction_type  <- "none"
  expect_error(do.call(optimize_famoptg, args))
})

test_that("stops when n_restarts < 1", {
  args             <- base_opt_args()
  args$n_restarts  <- 0L
  expect_error(do.call(optimize_famoptg, args))
})
