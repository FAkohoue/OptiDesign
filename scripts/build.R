# Build helper script for OptiDesign
# Example usage in R:

source("data-raw/generate_example_data.R")
devtools::document()
devtools::check()
devtools::build()


devtools::document()
devtools::install()
devtools::check()
pkgdown::build_site()

devtools::build_vignettes()
pkgdown::build_site()
