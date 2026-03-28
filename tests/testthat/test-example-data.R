# Tests using the shipped example dataset.
# These test the full do.call workflow for both design families.

data("OptiDesign_example_data", package = "OptiDesign")
x <- OptiDesign_example_data

# ── prep_famoptg ──────────────────────────────────────────────────────────────

test_that("prep_famoptg works with shipped family-based example", {
  out <- do.call(
    prep_famoptg,
    c(x$OptiDesign_famoptg_example, x$OptiDesign_famoptg_args_family)
  )

  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_true(is.data.frame(out$field_book))
  expect_length(out$seed_used, 1L)

  # No efficiency slot in refactored constructor
  expect_false("efficiency" %in% names(out))

  # All assigned plots appear in layout matrix
  expect_equal(
    nrow(out$field_book),
    sum(!is.na(c(out$layout_matrix)))
  )

  # Expected field book columns
  expect_true(all(c("Treatment", "Family", "Gcluster", "Block",
                    "Plot", "Row", "Column") %in% names(out$field_book)))
})

test_that("prep_famoptg works with shipped GRM-based example", {
  out <- do.call(
    prep_famoptg,
    c(x$OptiDesign_famoptg_example, x$OptiDesign_famoptg_args_grm)
  )

  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_true(is.data.frame(out$field_book))
  expect_length(out$seed_used, 1L)
  expect_false("efficiency" %in% names(out))
  expect_equal(
    nrow(out$field_book),
    sum(!is.na(c(out$layout_matrix)))
  )

  # GRM clustering should populate Gcluster column for non-checks
  non_check_rows <- out$field_book[
    !out$field_book$Treatment %in% x$OptiDesign_famoptg_example$check_treatments, ]
  expect_true(any(!is.na(non_check_rows$Gcluster)))
})

# ── alpha_rc_stream ───────────────────────────────────────────────────────────

test_that("alpha_rc_stream works with shipped family-based example", {
  out <- do.call(
    alpha_rc_stream,
    c(x$OptiDesign_alpha_example, x$OptiDesign_alpha_args_family)
  )

  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_true(is.data.frame(out$field_book))
  expect_true(is.list(out$design_info))
  expect_false("efficiency" %in% names(out))
  expect_length(out$seed_used, 1L)

  expect_equal(out$design_info$n_reps, 2L)
  expect_equal(
    sum(!is.na(c(out$layout_matrix))),
    out$design_info$total_used_plots   # corrected from used_plots
  )

  # Expected field book columns
  expect_true(all(c("Plot", "Row", "Column", "Rep", "IBlock",
                    "BlockInRep", "Treatment", "Family",
                    "Gcluster", "Check") %in% names(out$field_book)))
})

test_that("alpha_rc_stream works with shipped GRM-based example", {
  out <- do.call(
    alpha_rc_stream,
    c(x$OptiDesign_alpha_example, x$OptiDesign_alpha_args_grm)
  )

  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_true(is.data.frame(out$field_book))
  expect_true(is.list(out$design_info))
  expect_false("efficiency" %in% names(out))
  expect_equal(out$design_info$n_reps, 2L)
  expect_equal(
    sum(!is.na(c(out$layout_matrix))),
    out$design_info$total_used_plots
  )

  # GRM clustering populates Gcluster for non-check entries
  non_check_rows <- out$field_book[!out$field_book$Check, ]
  expect_true(any(!is.na(non_check_rows$Gcluster)))
})
