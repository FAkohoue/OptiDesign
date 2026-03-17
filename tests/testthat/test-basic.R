test_that("package example data loads correctly", {
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