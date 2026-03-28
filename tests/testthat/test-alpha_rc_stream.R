# ── Shared fixture ────────────────────────────────────────────────────────────
make_alpha_args <- function() {
  list(
    check_treatments = c("CHK1", "CHK2", "CHK3"),
    check_families   = c("CHECK", "CHECK", "CHECK"),
    entry_treatments = paste0("E", 1:30),
    entry_families   = rep(paste0("F", 1:5), 6L),
    n_reps           = 2L,
    n_rows           = 12L,
    n_cols           = 10L,
    min_block_size   = 8L,
    max_block_size   = 10L
  )
}

# ── Return value structure ─────────────────────────────────────────────────────

test_that("alpha_rc_stream returns correct list structure", {
  out <- do.call(alpha_rc_stream, make_alpha_args())

  expect_type(out, "list")
  expect_true("layout_matrix"  %in% names(out))
  expect_true("field_book"     %in% names(out))
  expect_true("design_info"    %in% names(out))
  expect_true("seed_used"      %in% names(out))
  expect_false("efficiency"    %in% names(out))  # evaluation is separate
})

test_that("layout_matrix has correct dimensions", {
  args <- make_alpha_args()
  out  <- do.call(alpha_rc_stream, args)
  expect_equal(dim(out$layout_matrix), c(args$n_rows, args$n_cols))
})

test_that("field_book has required columns", {
  out <- do.call(alpha_rc_stream, make_alpha_args())
  required <- c("Plot", "Row", "Column", "Rep", "IBlock",
                "BlockInRep", "Treatment", "Family", "Gcluster", "Check")
  expect_true(all(required %in% names(out$field_book)))
})

test_that("used plots match non-NA layout cells", {
  out <- do.call(alpha_rc_stream, make_alpha_args())
  expect_equal(
    sum(!is.na(c(out$layout_matrix))),
    out$design_info$total_used_plots
  )
})

test_that("seed_used is a single numeric", {
  out <- do.call(alpha_rc_stream, make_alpha_args())
  expect_length(out$seed_used, 1L)
})

test_that("design is reproducible with fixed seed", {
  args      <- make_alpha_args()
  args$seed <- 99L
  out1 <- do.call(alpha_rc_stream, args)
  out2 <- do.call(alpha_rc_stream, args)
  expect_identical(out1$field_book$Treatment, out2$field_book$Treatment)
})

# ── Replicate and block structure ─────────────────────────────────────────────

test_that("n_reps is recorded correctly in design_info", {
  args <- make_alpha_args()
  out  <- do.call(alpha_rc_stream, args)
  expect_equal(out$design_info$n_reps, args$n_reps)
})

test_that("each entry appears exactly once per replicate", {
  args <- make_alpha_args()
  out  <- do.call(alpha_rc_stream, args)
  fb   <- out$field_book[!out$field_book$Check & !is.na(out$field_book$Treatment), ]

  entry_rep_counts <- table(fb$Rep, fb$Treatment)
  expect_true(all(entry_rep_counts <= 1L))

  # Every entry should appear in every rep exactly once
  expect_true(all(entry_rep_counts == 1L))
})

test_that("checks appear in every incomplete block", {
  out <- do.call(alpha_rc_stream, make_alpha_args())
  fb  <- out$field_book[!is.na(out$field_book$Treatment), ]

  n_blocks <- max(fb$IBlock, na.rm = TRUE)
  for (chk in make_alpha_args()$check_treatments) {
    blocks_with_chk <- unique(fb$IBlock[fb$Treatment == chk])
    expect_setequal(blocks_with_chk, seq_len(n_blocks))
  }
})

test_that("design_info block sizes are within min/max constraints", {
  args <- make_alpha_args()
  out  <- do.call(alpha_rc_stream, args)

  block_sizes <- out$design_info$block_plan$BlockSize
  expect_true(all(block_sizes >= args$min_block_size))
  expect_true(all(block_sizes <= args$max_block_size))
})

# ── Field traversal ───────────────────────────────────────────────────────────

test_that("column order traversal works", {
  args       <- make_alpha_args()
  args$order <- "column"
  expect_no_error(do.call(alpha_rc_stream, args))
})

test_that("row order traversal works", {
  args       <- make_alpha_args()
  args$order <- "row"
  expect_no_error(do.call(alpha_rc_stream, args))
})

test_that("serpentine traversal works", {
  args            <- make_alpha_args()
  args$serpentine <- TRUE
  expect_no_error(do.call(alpha_rc_stream, args))
})

# ── Grouping ──────────────────────────────────────────────────────────────────

test_that("family clustering: Gcluster is NA for all plots", {
  args                <- make_alpha_args()
  args$cluster_source <- "Family"
  out <- do.call(alpha_rc_stream, args)
  used <- out$field_book[!is.na(out$field_book$Treatment), ]
  expect_true(all(is.na(used$Gcluster)))
})

test_that("check plots always have Gcluster = NA", {
  out       <- do.call(alpha_rc_stream, make_alpha_args())
  chk_rows  <- out$field_book[out$field_book$Check & !is.na(out$field_book$Treatment), ]
  expect_true(all(is.na(chk_rows$Gcluster)))
})

# ── Input validation ──────────────────────────────────────────────────────────

test_that("stops on duplicate check_treatments", {
  args <- make_alpha_args()
  args$check_treatments <- c("CHK1", "CHK1", "CHK3")
  args$check_families   <- c("CHECK", "CHECK", "CHECK")
  expect_error(do.call(alpha_rc_stream, args), "Duplicate")
})

test_that("stops when check and entry overlap", {
  args <- make_alpha_args()
  args$entry_treatments[1] <- args$check_treatments[1]
  expect_error(do.call(alpha_rc_stream, args))
})

test_that("stops when min_block_size > max_block_size", {
  args                <- make_alpha_args()
  args$min_block_size <- 10L
  args$max_block_size <- 8L
  expect_error(do.call(alpha_rc_stream, args))
})
