# Build helper script for OptiDesign
# Example usage in R:
#unlink("pkgdown", recursive = TRUE)

source("data-raw/generate_example_data.R")
devtools::document()
devtools::test()

devtools::check()

devtools::install()

unloadNamespace("OptiDesign")
pkgdown::build_site()

devtools::build()

devtools::build_vignettes()
pkgdown::build_site()

install.packages("rsvg")
rsvg::rsvg_png(
  "man/figures/logo.svg",
  "man/figures/logo.png",
  width  = 1360,
  height = 560
)

# Then build favicons from the PNG instead
pkgdown::build_favicons(overwrite = TRUE)

options(timeout = 3000)  # 5 minutes
pkgdown::build_favicons(overwrite = TRUE)

# Regenerate favicons from the correct logo
pkgdown::build_favicons(overwrite = TRUE)
pkgdown::build_site()
pkgdown::build_site(override = list(template = list(favicon = FALSE)))

library(magick)

# Load your logo
img <- image_read("man/figures/logo.png")

# Create favicon sizes
sizes <- c(16, 32, 48, 64, 180, 192, 512)

dir.create("pkgdown/assets", recursive = TRUE, showWarnings = FALSE)

for (s in sizes) {
  resized <- image_resize(img, paste0(s, "x", s))
  image_write(resized, paste0("pkgdown/assets/favicon-", s, ".png"))
}

# Create favicon.ico (multi-size)
ico <- image_resize(img, "64x64")
image_write(ico, "pkgdown/assets/favicon.ico")

writeLines(
  c(
    '<link rel="icon" type="image/svg+xml" href="logo.svg">',
    '<link rel="icon" type="image/png" sizes="32x32" href="favicon-32.png">',
    '<link rel="apple-touch-icon" href="favicon-180.png">'
  ),
  "pkgdown/assets/favicon.html"
)


# Create the correct folder
dir.create("pkgdown/favicon", showWarnings = FALSE)

# Copy all favicon files from assets/ to favicon/
file.copy(
  from      = list.files("pkgdown/assets", full.names = TRUE),
  to        = "pkgdown/favicon/",
  overwrite = TRUE
)

list.files("pkgdown/favicon")
