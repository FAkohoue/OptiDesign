# Internal square matrix validator

Checks that an object is a base R matrix with the expected square
dimensions. Stops with an informative error if either condition is not
met. Used internally to validate relationship matrices (`GRM`, `A`, `K`)
before they are used for clustering, dispersion optimisation, or
efficiency evaluation.

## Usage

``` r
.validate_square_matrix(M, p, nm = "matrix")
```

## Arguments

- M:

  Object to validate. Must be a base R `matrix`.

- p:

  Expected dimension. The matrix must have `nrow(M) == p` and
  `ncol(M) == p`.

- nm:

  Character scalar. Object name used in error messages to identify which
  argument failed validation. Default `"matrix"`.

## Value

Invisibly returns `TRUE` if validation passes.
