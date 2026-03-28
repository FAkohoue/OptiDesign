# Internal conditional debug message helper

Emits a formatted message via
[`base::message()`](https://rdrr.io/r/base/message.html) when
`debug = TRUE`. Used internally to provide optional verbose output
during design construction and optimisation without cluttering normal
output.

## Usage

``` r
.optidesign_dbg(debug, ...)
```

## Arguments

- debug:

  Logical. If `TRUE`, the message is emitted. If `FALSE` or `NULL`, the
  function returns invisibly with no output.

- ...:

  Arguments passed to
  [`base::sprintf()`](https://rdrr.io/r/base/sprintf.html). The first
  argument should be a format string; subsequent arguments are
  substituted into the format string in order.

## Value

Invisibly returns `NULL`.
