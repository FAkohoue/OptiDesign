test_that("prep_famoptg works with shipped family-based example", {
  data("OptiDesign_example_data", package = "OptiDesign")
  x <- OptiDesign_example_data
  
  out <- do.call(
    prep_famoptg,
    c(
      x$OptiDesign_famoptg_example,
      x$OptiDesign_famoptg_args_family
    )
  )
  
  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_true(is.data.frame(out$field_book))
  expect_true(is.null(out$efficiency))
  expect_length(out$seed_used, 1)
  expect_equal(
    nrow(out$field_book),
    sum(!is.na(c(out$layout_matrix)))
  )
})

test_that("prep_famoptg works with shipped GRM example", {
  data("OptiDesign_example_data", package = "OptiDesign")
  x <- OptiDesign_example_data
  
  out <- do.call(
    prep_famoptg,
    c(
      x$OptiDesign_famoptg_example,
      x$OptiDesign_famoptg_args_grm
    )
  )
  
  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_true(is.data.frame(out$field_book))
  expect_true(is.null(out$efficiency))
  expect_length(out$seed_used, 1)
  expect_equal(
    nrow(out$field_book),
    sum(!is.na(c(out$layout_matrix)))
  )
})

test_that("prep_alpha_checks_rc_stream works with shipped family-based example", {
  data("OptiDesign_example_data", package = "OptiDesign")
  x <- OptiDesign_example_data
  
  out <- do.call(
    prep_alpha_checks_rc_stream,
    c(
      x$OptiDesign_alpha_example,
      x$OptiDesign_alpha_args_family
    )
  )
  
  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_true(is.data.frame(out$field_book))
  expect_true(is.list(out$design_info))
  expect_true(is.null(out$efficiency))
  expect_equal(out$design_info$n_reps, 2)
  expect_equal(
    sum(!is.na(c(out$layout_matrix))),
    out$design_info$used_plots
  )
})

test_that("prep_alpha_checks_rc_stream works with shipped GRM example", {
  data("OptiDesign_example_data", package = "OptiDesign")
  x <- OptiDesign_example_data
  
  out <- do.call(
    prep_alpha_checks_rc_stream,
    c(
      x$OptiDesign_alpha_example,
      x$OptiDesign_alpha_args_grm
    )
  )
  
  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_true(is.data.frame(out$field_book))
  expect_true(is.list(out$design_info))
  expect_true(is.null(out$efficiency))
  expect_equal(out$design_info$n_reps, 2)
  expect_equal(
    sum(!is.na(c(out$layout_matrix))),
    out$design_info$used_plots
  )
})