# ── Shared fixtures ───────────────────────────────────────────────────────────
make_famoptg_args <- function() {
  list(
    check_treatments        = c("CHK1", "CHK2"),
    check_families          = c("CHECK", "CHECK"),
    p_rep_treatments        = paste0("P", 1:10),
    p_rep_reps              = rep(2L, 10L),
    p_rep_families          = rep(c("F1", "F2"), 5L),
    unreplicated_treatments = paste0("U", 1:20),
    unreplicated_families   = rep(c("F1", "F2", "F3", "F4"), 5L),
    n_blocks                = 4L,
    n_rows                  = 6L,
    n_cols                  = 8L     # 4*2 + 10*2 + 20 = 48 = 6*8
  )
}

# ── Return value structure ─────────────────────────────────────────────────────

test_that("prep_famoptg returns correct list structure", {
  out <- do.call(prep_famoptg, make_famoptg_args())

  expect_type(out, "list")
  expect_named(out, c("layout_matrix", "field_book", "seed_used"),
               ignore.order = TRUE)
  expect_false("efficiency" %in% names(out))  # efficiency removed from constructor
})

test_that("layout_matrix has correct dimensions", {
  args <- make_famoptg_args()
  out  <- do.call(prep_famoptg, args)
  expect_equal(dim(out$layout_matrix), c(args$n_rows, args$n_cols))
})

test_that("field_book has required columns", {
  out <- do.call(prep_famoptg, make_famoptg_args())
  required_cols <- c("Treatment", "Family", "Gcluster", "Block", "Plot",
                     "Row", "Column")
  expect_true(all(required_cols %in% names(out$field_book)))
})

test_that("assigned plots match non-NA layout cells", {
  out <- do.call(prep_famoptg, make_famoptg_args())
  expect_equal(nrow(out$field_book), sum(!is.na(c(out$layout_matrix))))
})

test_that("seed_used is a single integer", {
  out <- do.call(prep_famoptg, make_famoptg_args())
  expect_length(out$seed_used, 1L)
  expect_true(is.numeric(out$seed_used))
})

test_that("design is reproducible with fixed seed", {
  args       <- make_famoptg_args()
  args$seed  <- 42L
  out1 <- do.call(prep_famoptg, args)
  out2 <- do.call(prep_famoptg, args)
  expect_identical(out1$field_book$Treatment, out2$field_book$Treatment)
})

# ── Block structure constraints ───────────────────────────────────────────────

test_that("checks appear in every block", {
  args <- make_famoptg_args()
  out  <- do.call(prep_famoptg, args)
  fb   <- out$field_book

  for (chk in args$check_treatments) {
    blocks_with_chk <- unique(fb$Block[fb$Treatment == chk])
    expect_setequal(blocks_with_chk, seq_len(args$n_blocks))
  }
})

test_that("p-rep treatments appear in exactly p_rep_reps[i] blocks", {
  args <- make_famoptg_args()
  out  <- do.call(prep_famoptg, args)
  fb   <- out$field_book

  for (i in seq_along(args$p_rep_treatments)) {
    trt    <- args$p_rep_treatments[i]
    actual <- sum(fb$Treatment == trt)
    expect_equal(actual, args$p_rep_reps[i],
                 label = paste0("replication count for ", trt))
  }
})

test_that("p-rep treatments never appear twice in the same block", {
  args <- make_famoptg_args()
  out  <- do.call(prep_famoptg, args)
  fb   <- out$field_book

  for (trt in args$p_rep_treatments) {
    trt_blocks <- fb$Block[fb$Treatment == trt]
    expect_equal(length(trt_blocks), length(unique(trt_blocks)),
                 label = paste0("no within-block duplicate for ", trt))
  }
})

test_that("unreplicated treatments appear exactly once", {
  args <- make_famoptg_args()
  out  <- do.call(prep_famoptg, args)
  fb   <- out$field_book

  for (trt in args$unreplicated_treatments) {
    expect_equal(sum(fb$Treatment == trt), 1L,
                 label = paste0("exactly one occurrence of ", trt))
  }
})

# ── Three design classes ──────────────────────────────────────────────────────

test_that("augmented design: no p-rep treatments", {
  args <- make_famoptg_args()
  args$p_rep_treatments        <- character(0)
  args$p_rep_reps               <- integer(0)
  args$p_rep_families           <- character(0)
  # Adjust field size: 4 blocks * 2 checks + 20 unrep = 28 plots = 4*7
  args$n_rows <- 4L; args$n_cols <- 7L

  out <- do.call(prep_famoptg, args)
  expect_true(is.matrix(out$layout_matrix))
  expect_equal(dim(out$layout_matrix), c(4L, 7L))
})

test_that("RCBD-type design: all non-check entries equally replicated", {
  # 2 checks, 10 entries each replicated 4 times, 4 blocks
  # Plots = 4*2 + 10*4 = 8 + 40 = 48 = 6*8
  args <- list(
    check_treatments        = c("CHK1", "CHK2"),
    check_families          = c("CHECK", "CHECK"),
    p_rep_treatments        = paste0("E", 1:10),
    p_rep_reps              = rep(4L, 10L),
    p_rep_families          = rep(c("F1", "F2"), 5L),
    unreplicated_treatments = character(0),
    unreplicated_families   = character(0),
    n_blocks = 4L, n_rows = 6L, n_cols = 8L
  )
  out <- do.call(prep_famoptg, args)
  fb  <- out$field_book

  for (trt in args$p_rep_treatments) {
    expect_equal(sum(fb$Treatment == trt), 4L)
    trt_blocks <- fb$Block[fb$Treatment == trt]
    expect_equal(length(unique(trt_blocks)), 4L)
  }
})

# ── Input validation ──────────────────────────────────────────────────────────

test_that("stops on duplicate check_treatments", {
  args <- make_famoptg_args()
  args$check_treatments <- c("CHK1", "CHK1")
  args$check_families   <- c("CHECK", "CHECK")
  expect_error(do.call(prep_famoptg, args), "Duplicate")
})

test_that("stops when treatment appears in two classes", {
  args <- make_famoptg_args()
  args$p_rep_treatments[1] <- args$check_treatments[1]
  expect_error(do.call(prep_famoptg, args))
})

test_that("stops when p_rep_reps exceeds n_blocks", {
  args <- make_famoptg_args()
  args$p_rep_reps[1] <- args$n_blocks + 1L
  expect_error(do.call(prep_famoptg, args))
})

test_that("stops when field size mismatch and warn_and_correct = FALSE", {
  args                   <- make_famoptg_args()
  args$n_cols            <- args$n_cols + 1L   # deliberate mismatch
  args$warn_and_correct  <- FALSE
  expect_error(do.call(prep_famoptg, args))
})

test_that("adjusts dimensions when warn_and_correct = TRUE", {
  args                  <- make_famoptg_args()
  args$n_cols           <- args$n_cols + 1L
  args$warn_and_correct <- TRUE
  expect_warning(out <- do.call(prep_famoptg, args))
  expect_true(is.matrix(out$layout_matrix))
})

# ── Grouping ──────────────────────────────────────────────────────────────────

test_that("family clustering leaves Gcluster as NA", {
  args <- make_famoptg_args()
  out  <- do.call(prep_famoptg, args)
  expect_true(all(is.na(out$field_book$Gcluster)))
})

test_that("check plots always have Gcluster = NA regardless of cluster_source", {
  args <- make_famoptg_args()
  out  <- do.call(prep_famoptg, args)
  check_rows <- out$field_book[
    out$field_book$Treatment %in% args$check_treatments, ]
  expect_true(all(is.na(check_rows$Gcluster)))
})

# ── check_placement modes ─────────────────────────────────────────────────────

test_that("check_placement = 'systematic' runs without error", {
  args                <- make_famoptg_args()
  args$check_placement <- "systematic"
  expect_no_error(do.call(prep_famoptg, args))
})

test_that("check_placement = 'optimal' runs without error", {
  args                    <- make_famoptg_args()
  args$check_placement     <- "optimal"
  args$check_opt_attempts  <- 5L
  expect_no_error(do.call(prep_famoptg, args))
})
