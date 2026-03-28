# Internal z-score normalisation helper

Computes the z-score (standard normal transformation) of a numeric
vector, handling edge cases where the standard deviation is zero or
non-finite. Used internally when combining or scaling multiple
optimality criteria prior to comparison.

## Usage

``` r
.optidesign_z(x)
```

## Arguments

- x:

  Numeric vector. `NA` values are ignored when computing the mean and
  standard deviation but are returned as `NA` in the output.

## Value

Numeric vector of the same length as `x` containing z-scores, or a zero
vector if the standard deviation is zero or non-finite.

## Details

If the standard deviation is zero or non-finite (e.g. all values
identical, or all `NA`), returns a zero vector of the same length rather
than producing `NaN` or `Inf`.
