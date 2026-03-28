---
name: Bug report
about: Report a reproducible problem in OptiDesign
title: "[BUG] "
labels: bug
assignees: FAkohoue
---

## What happened?

<!-- Describe the problem clearly. What did the function do that it should not have done,
     or what did it fail to do? Include any error messages or warnings in full. -->

## What did you expect?

<!-- Describe what you expected to happen instead. -->

## Which function is affected?

<!-- Tick the relevant function(s). -->

- [ ] `prep_famoptg()`
- [ ] `evaluate_famoptg_efficiency()`
- [ ] `optimize_famoptg()`
- [ ] `alpha_rc_stream()`
- [ ] `evaluate_alpha_efficiency()`
- [ ] `optimize_alpha_rc()`
- [ ] Other / unsure

## Minimal reproducible example

<!--
Please provide the smallest self-contained code that reproduces the problem.
If your data is sensitive, replace it with synthetic data of the same structure.
The example below uses the shipped dataset as a starting point.
-->

```r
library(OptiDesign)

# Using the shipped example data as a starting point:
data("OptiDesign_example_data", package = "OptiDesign")
x <- OptiDesign_example_data

# Replace with the code that triggers the bug:

```

## Observed output or error message

<!-- Paste the full output, error message, or warning here. -->

```
# paste output here
```

## Session info

```r
sessionInfo()
```

<!-- Paste the output of sessionInfo() here: -->

```
# paste sessionInfo() output here
```

## Additional context

<!-- Any other information that might be relevant: field dimensions,
     number of treatments, variance component values, optimisation method used,
     operating system, or R version constraints. -->