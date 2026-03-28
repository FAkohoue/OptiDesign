test_that("package loads and all exported functions exist", {
  expect_true(isNamespaceLoaded("OptiDesign"))

  exported <- getNamespaceExports("OptiDesign")
  expect_true("prep_famoptg"                %in% exported)
  expect_true("evaluate_famoptg_efficiency" %in% exported)
  expect_true("optimize_famoptg"            %in% exported)
  expect_true("alpha_rc_stream"             %in% exported)
  expect_true("evaluate_alpha_efficiency"   %in% exported)
  expect_true("optimize_alpha_rc"           %in% exported)
})

test_that("internal helpers exist and are not exported", {
  # Helpers should be accessible via ::: but not via ::
  expect_true(is.function(OptiDesign:::.make_sparse_incidence))
  expect_true(is.function(OptiDesign:::.ar1_precision_sparse))
  expect_true(is.function(OptiDesign:::.solve_C))
  expect_true(is.function(OptiDesign:::.pinv_sym_dense))
  expect_true(is.function(OptiDesign:::.safe_logdet_psd_dense))
  expect_true(is.function(OptiDesign:::.pairwise_diff_mean_var))
  expect_true(is.function(OptiDesign:::.trace_subinv_est))
  expect_true(is.function(OptiDesign:::.with_local_seed))
  expect_true(is.function(OptiDesign:::.build_neighbor_pairs))
  expect_true(is.function(OptiDesign:::.score_dispersion))

  # Confirm they are not in the public namespace
  expect_false(".make_sparse_incidence" %in% getNamespaceExports("OptiDesign"))
  expect_false(".ar1_precision_sparse"  %in% getNamespaceExports("OptiDesign"))
})

test_that("example dataset loads and has expected structure", {
  data("OptiDesign_example_data", package = "OptiDesign")
  expect_true(exists("OptiDesign_example_data"))
  expect_true(is.list(OptiDesign_example_data))

  expected_names <- c(
    "OptiDesign_lines",
    "OptiDesign_id_map",
    "OptiDesign_GRM",
    "OptiDesign_A",
    "OptiDesign_K",
    "OptiDesign_famoptg_example",
    "OptiDesign_alpha_example",
    "OptiDesign_famoptg_args_family",
    "OptiDesign_famoptg_args_grm",
    "OptiDesign_alpha_args_family",
    "OptiDesign_alpha_args_grm"
  )
  expect_true(all(expected_names %in% names(OptiDesign_example_data)))
})
