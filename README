# OptiDesign

`OptiDesign` provides field-design construction utilities for breeding and agronomic experiments.

## Exported functions

- `prep_famoptg()`: partially replicated (p-rep) design construction with optional grouping, dispersion optimization, and efficiency diagnostics.
- `prep_alpha_checks_rc_stream()`: fixed-grid alpha row-column design with repeated checks in every incomplete block, optional grouping, dispersion optimization, and efficiency diagnostics.

## Installation

After creating the GitHub repository under `FAkohoue/OptiDesign`, install with:

```r
# install.packages("remotes")
remotes::install_github("FAkohoue/OptiDesign")
```

## Important note about `mod()`

The function `prep_famoptg()` was kept exactly as provided and still calls `mod()` inside the serpentine position builder. As requested, `mod()` is **not** defined in this package skeleton.

That means:

- the package can be built and versioned as provided,
- but calling `prep_famoptg()` with code paths that reach `mod()` requires `mod()` to exist in the package namespace or to be added later.

A reminder note is also stored in `inst/NOTICE_MOD.txt`.
