# OptiDesign: Optimized Experimental Field Design for Plant Breeding

``` r
library(OptiDesign)
```

## Introduction

Experimental design is a foundational component of plant breeding and
agronomic research. The ability to accurately estimate genetic effects,
compare treatments, and predict breeding values depends critically on
how field trials are constructed.

Classical designs — randomised complete block designs (RCBD),
alpha-lattice designs, and augmented designs — provide well-understood
statistical properties but make assumptions that modern breeding
programs increasingly cannot meet: regular block structures, balanced
replication across all entries, and no prior knowledge of genetic
relationships. Contemporary trials involve large numbers of candidates,
limited field capacity, spatial heterogeneity, and rich genomic
information that classical frameworks treat as external to the design
process.

`OptiDesign` addresses these challenges through a unified framework that
integrates:

- flexible field layout construction for two design families
- genetic structure from family labels, pedigree matrices, or genomic
  relationship matrices
- optional spatial dispersion optimisation to reduce clustering of
  related entries
- mixed-model efficiency evaluation under A, D, and CDmean optimality
  criteria
- **criterion-driven design search** that returns the statistically best
  design across many randomisations

The package follows a **single-responsibility architecture**:
construction, evaluation, and optimisation are separated into distinct
functions that can be called independently or chained, making each step
transparent and reproducible.

------------------------------------------------------------------------

## Package Architecture

`OptiDesign` provides six exported functions organised into two design
families:

| Function                                                                                                          | Role                        | Design family           |
|-------------------------------------------------------------------------------------------------------------------|-----------------------------|-------------------------|
| [`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)                               | Construction                | Repeated-check block    |
| [`evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiDesign/reference/evaluate_famoptg_efficiency.md) | Evaluation                  | Repeated-check block    |
| [`optimize_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/optimize_famoptg.md)                       | Optimisation (RS)           | Repeated-check block    |
| [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md)                         | Construction                | Alpha row-column stream |
| [`evaluate_alpha_efficiency()`](https://FAkohoue.github.io/OptiDesign/reference/evaluate_alpha_efficiency.md)     | Evaluation                  | Alpha row-column stream |
| [`optimize_alpha_rc()`](https://FAkohoue.github.io/OptiDesign/reference/optimize_alpha_rc.md)                     | Optimisation (RS / SA / GA) | Alpha row-column stream |

The two families differ in their blocking structure, replication model,
and the depth of their optimisation capabilities:

| Feature          | `prep_famoptg` family       | `alpha_rc_stream` family       |
|------------------|-----------------------------|--------------------------------|
| Blocking         | Flat blocks                 | Replicates → incomplete blocks |
| Replication      | Flexible per-entry          | Uniform across entries         |
| Design types     | Augmented, p-rep, RCBD-type | Alpha-lattice                  |
| Block variance   | `sigma_b2`                  | `sigma_rep2` + `sigma_ib2`     |
| Optimisation     | RS only                     | RS, SA, GA                     |
| P-rep constraint | Enforced by construction    | Not applicable                 |

------------------------------------------------------------------------

## Statistical Framework

### The mixed model

Both design families are evaluated under the same general mixed model:

$$y = X\beta + Zu + e$$

| Symbol  | Description                                                 |
|---------|-------------------------------------------------------------|
| $y$     | Vector of observed phenotypes                               |
| $X$     | Fixed effects design matrix                                 |
| $\beta$ | Fixed effects (intercept, checks, entry effects when fixed) |
| $Z$     | Incidence matrix linking random effects to plots            |
| $u$     | Random effects (blocks, rows, columns, entries when random) |
| $e$     | Residual vector                                             |

Random effects: $u \sim N(0,\, G)$ and $e \sim N(0,\, R)$.

For the **repeated-check block** family:
$$G^{- 1} = \text{blockdiag}\left( \sigma_{b}^{- 2}I,\;\sigma_{r}^{- 2}I,\;\sigma_{c}^{- 2}I,\;\sigma_{g}^{- 2}K^{- 1} \right)$$

For the **alpha row-column stream** family:
$$G^{- 1} = \text{blockdiag}\left( \sigma_{\text{rep}}^{- 2}I,\;\sigma_{\text{ib}}^{- 2}I,\;\sigma_{r}^{- 2}I,\;\sigma_{c}^{- 2}I,\;\sigma_{g}^{- 2}K^{- 1} \right)$$

### Mixed model coefficient matrix

Efficiency criteria are derived from the mixed model coefficient matrix:

$$C = \begin{pmatrix}
{X^{\top}QX} & {X^{\top}QZ} \\
{Z^{\top}QX} & {Z^{\top}QZ + G^{- 1}}
\end{pmatrix}$$

where $Q = R^{- 1}$ is the residual precision matrix.

### Residual structures

Three residual structures are supported. For an AR1 process of length
$n$ with autocorrelation $\rho$, the precision matrix $Q_{\text{AR1}}$
is tridiagonal with interior diagonal entries
$\left( 1 + \rho^{2} \right)/\left( 1 - \rho^{2} \right)$, edge diagonal
entries $1/\left( 1 - \rho^{2} \right)$, and off-diagonal entries
$- \rho/\left( 1 - \rho^{2} \right)$.

| Structure | Formula                                                                                                                        | Parameters           |
|-----------|--------------------------------------------------------------------------------------------------------------------------------|----------------------|
| IID       | $R = \sigma_{e}^{2}I$                                                                                                          | `sigma_e2`           |
| AR1       | $R^{- 1} = \sigma_{e}^{- 2}\left( Q_{\text{AR1}}\left( \rho_{r} \right) \otimes I_{c} \right)$                                 | `rho_row`            |
| AR1×AR1   | $R^{- 1} = \sigma_{e}^{- 2}\left( Q_{\text{AR1}}\left( \rho_{c} \right) \otimes Q_{\text{AR1}}\left( \rho_{r} \right) \right)$ | `rho_row`, `rho_col` |

### Optimality criteria

**A-criterion** (lower is better): minimises the mean pairwise contrast
variance under fixed treatment effects, or the mean prediction error
variance (PEV) under random treatment effects.

$$A_{\text{criterion}} = \frac{2}{p(p - 1)}\sum\limits_{i < j}\text{Var}\left( {\widehat{\tau}}_{i} - {\widehat{\tau}}_{j} \right)$$

**D-criterion** (lower is better): minimises the geometric mean of the
contrast covariance eigenvalues (fixed effects only).

$$D_{\text{criterion}} = \exp\!\left( \frac{\log\det(HVH)}{p - 1} \right)$$

where $H = I_{p} - p^{- 1}J_{p}$ is the centering matrix and $V$ is the
treatment variance-covariance submatrix of $C^{- 1}$.

**CDmean** (higher is better): the mean coefficient of determination for
genomic breeding value (GEBV) prediction (Rincent et al. 2012). Measures
the proportion of genetic variance explained by prediction on average
across lines.

$$\text{CDmean} = 1 - \frac{\text{mean PEV}}{\sigma_{g}^{2}}$$

CDmean ranges from 0 (no information) to 1 (perfect prediction). It is
the primary criterion for optimising training population designs in
genomic selection.

### Large-design approximation

When the number of treatments exceeds `eff_full_max` (default 400),
exact inversion of the $C$ submatrix is replaced by the Hutchinson
stochastic trace estimator (Hutchinson 1990), which approximates
$\text{trace}\left( C^{- 1}\left\lbrack \text{idx},\text{idx} \right\rbrack \right)$
using $m$ Rademacher random vectors:

$$\text{trace}\left( C_{\text{idx}}^{- 1} \right) \approx \frac{1}{m}\sum\limits_{k = 1}^{m}z_{k}^{\top}C^{- 1}z_{k},\quad z_{k} \sim \text{Rademacher}$$

The result carries the `_APPROX` mode suffix and `D_criterion` is `NA`.

------------------------------------------------------------------------

## Design Family 1: Repeated-Check Block Designs

### Construction with `prep_famoptg()`

[`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)
builds a repeated-check block design in which checks appear in every
block and non-check entries are allocated across blocks according to
their specified replication levels.

**Three design classes share the same function:**

**Augmented design** — all test entries unreplicated, checks repeated in
every block. Suitable for very large entry sets where resources only
allow a single observation per candidate.

**P-rep design** — some entries replicated across multiple distinct
blocks, others unreplicated. The most general case: a mixture of
candidate entries with different priority levels.

**RCBD-type design** — all non-check entries equally replicated. When
the replication number equals `n_blocks`, every entry appears in every
block — the closest repeated-check analogue of a classical RCBD.

**The p-rep constraint** — the core structural rule: no replicated
treatment ever appears twice in the same block. Enforced by construction
at every call, not by post-hoc checking.

Total required plots:
$$\text{total} = n_{\text{blocks}} \times n_{\text{checks}} + \sum\limits_{i = 1}^{v_{p}}r_{i} + v_{u}$$

where $r_{i}$ is the replication count of p-rep entry $i$ and $v_{u}$ is
the number of unreplicated entries.

### Evaluation with `evaluate_famoptg_efficiency()`

Takes the `field_book` returned by
[`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)
and computes A, D, and CDmean criteria. Fully decoupled from
construction — the same field book can be evaluated multiple times under
different model assumptions.

The random effect model for this family uses `sigma_b2` for the flat
block structure (no replicate or incomplete-block nesting).

### Optimisation with `optimize_famoptg()`

Runs `n_restarts` independent calls to
[`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)
with different seeds and returns the design with the best criterion
value. **Random Restart (RS) is the only method** because the p-rep
constraint is enforced by construction at every call — permutation-based
methods (SA, GA) would require block-aware swap logic to preserve it.

------------------------------------------------------------------------

## Design Family 2: Alpha Row-Column Stream Designs

### Construction with `alpha_rc_stream()`

[`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md)
builds a fixed-grid alpha row-column design using a stream-based layout.
The field is converted to a one-dimensional planting stream, partitioned
into $n_{\text{reps}}$ contiguous replicate segments, and each segment
is divided into incomplete blocks. Checks appear in every incomplete
block; each entry appears exactly once per replicate. Unused cells
appear only at the end of the stream.

Block-size constraints are expressed in total block size (checks +
entries) via `min_block_size` and `max_block_size`. The number of
incomplete blocks per replicate $b$ must satisfy:

$$\left\lceil \frac{v}{\text{max\_block\_size} - c} \right\rceil \leq b \leq \left\lfloor \frac{v}{\text{min\_block\_size} - c} \right\rfloor$$

where $v$ is the number of entries and $c$ is the number of checks.

### Evaluation with `evaluate_alpha_efficiency()`

Computes A, D, and CDmean criteria for a design produced by
[`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md).
The model contains Rep + IBlock(Rep) + Row + Column random effects,
using `sigma_rep2` and `sigma_ib2` — distinct from
[`evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiDesign/reference/evaluate_famoptg_efficiency.md)
which uses `sigma_b2`.

### Optimisation with `optimize_alpha_rc()`

Wraps
[`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md)
and
[`evaluate_alpha_efficiency()`](https://FAkohoue.github.io/OptiDesign/reference/evaluate_alpha_efficiency.md)
in an optimisation loop with three search strategies:

**RS (Random Restart)** — generate `n_restarts` independent designs,
return the best. Guaranteed validity, simple, easily parallelisable.

**SA (Simulated Annealing)** — iterative entry permutation swaps with
Metropolis acceptance:
$$P\left( \text{accept worse} \right) = \exp\!\left( - \frac{\Delta}{T_{k}} \right)$$
where $T_{k}$ cools from `sa_temp_start` to `sa_temp_end`. Better at
escaping local optima than RS. Invalid swap proposals are treated as
neutral events and do not affect the acceptance rate.

**GA (Genetic Algorithm)** — population of entry permutations evolved
via Order Crossover (OX1), random swap mutation, tournament selection,
and elitism. Most powerful for global search.

------------------------------------------------------------------------

## Integrity Checking

Both optimisers implement a four-point integrity checking strategy that
guarantees the returned design is structurally valid:

1.  **Post-construction** — every candidate validated immediately after
    [`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)
    or
    [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md)
    returns, before scoring.
2.  **Pre-storage** — candidate re-checked before updating the running
    best.
3.  **Pre-return** — stored best re-checked one final time before
    returning to the user.
4.  **Emergency fallback** — if no valid design is found after all
    iterations, up to 10 fresh random designs are attempted before
    stopping with an informative error.

For the `prep_famoptg` family, five structural constraints are verified:

- No non-check entry appears more than once in a single block
- Each p-rep treatment appears in exactly `p_rep_reps[i]` blocks
- Each unreplicated treatment appears exactly once
- All checks appear in every block
- No p-rep treatment occupies the same block twice (the core p-rep
  guarantee)

------------------------------------------------------------------------

## Example Dataset

The package ships with a built-in example dataset for both design
families:

``` r
data("OptiDesign_example_data", package = "OptiDesign")
x <- OptiDesign_example_data
names(x)
#>  [1] "OptiDesign_lines"               "OptiDesign_id_map"             
#>  [3] "OptiDesign_GRM"                 "OptiDesign_A"                  
#>  [5] "OptiDesign_K"                   "OptiDesign_famoptg_example"    
#>  [7] "OptiDesign_alpha_example"       "OptiDesign_famoptg_args_family"
#>  [9] "OptiDesign_famoptg_args_grm"    "OptiDesign_alpha_args_family"  
#> [11] "OptiDesign_alpha_args_grm"
```

The dataset contains treatment vectors, field dimensions, relationship
matrices, and ready-to-use argument lists structured for
[`do.call()`](https://rdrr.io/r/base/do.call.html) workflows.

------------------------------------------------------------------------

## Workflow 1: Repeated-Check Block Design (Family-Based)

### Step 1 — Construct

``` r
design_fam <- do.call(
  prep_famoptg,
  c(x$OptiDesign_famoptg_example, x$OptiDesign_famoptg_args_family)
)

dim(design_fam$layout_matrix)
#> [1] 16  8
head(design_fam$field_book)
#>   Treatment Family Gcluster Block Plot Row Column
#> 1      L027    F21     <NA>     1    1   1      1
#> 2      L016    F03     <NA>     1    2   2      1
#> 3      L002    F19     <NA>     1    3   3      1
#> 4      L006    F18     <NA>     1    4   4      1
#> 5      L008    F11     <NA>     1    5   5      1
#> 6      L001    F15     <NA>     1    6   6      1
```

The field book contains one row per assigned plot with treatment
identity, family label, genomic cluster (if applicable), block, plot
number, row, and column. There is no `efficiency` slot — evaluation is a
separate step.

### Step 2 — Evaluate

``` r
eff_fam <- evaluate_famoptg_efficiency(
  field_book         = design_fam$field_book,
  n_rows             = x$OptiDesign_famoptg_example$n_rows,
  n_cols             = x$OptiDesign_famoptg_example$n_cols,
  check_treatments   = x$OptiDesign_famoptg_example$check_treatments,
  treatment_effect   = "fixed",
  residual_structure = "IID"
)

cat("A-criterion (lower is better):", round(eff_fam$A_criterion, 4), "\n")
#> A-criterion (lower is better): 2.4774
cat("D-criterion (lower is better):", round(eff_fam$D_criterion, 4), "\n")
#> D-criterion (lower is better): 0.9643
cat("A-efficiency (higher is better):", round(eff_fam$A_efficiency, 4), "\n")
#> A-efficiency (higher is better): 0.4036
cat("Mode:", eff_fam$mode, "\n")
#> Mode: FIXED_TREATMENT_BLUE_CONTRAST
```

The same field book can be re-evaluated under a spatial model without
rebuilding the design:

``` r
eff_fam_ar1 <- evaluate_famoptg_efficiency(
  field_book         = design_fam$field_book,
  n_rows             = x$OptiDesign_famoptg_example$n_rows,
  n_cols             = x$OptiDesign_famoptg_example$n_cols,
  check_treatments   = x$OptiDesign_famoptg_example$check_treatments,
  treatment_effect   = "fixed",
  residual_structure = "AR1xAR1",
  rho_row            = 0.3,
  rho_col            = 0.2
)

cat("A-criterion under AR1xAR1:", round(eff_fam_ar1$A_criterion, 4), "\n")
#> A-criterion under AR1xAR1: 2.1032
```

### Step 3 — Optimise (optional)

``` r
opt_fam <- optimize_famoptg(
  # Construction arguments
  check_treatments        = x$OptiDesign_famoptg_example$check_treatments,
  check_families          = x$OptiDesign_famoptg_example$check_families,
  p_rep_treatments        = x$OptiDesign_famoptg_example$p_rep_treatments,
  p_rep_reps              = x$OptiDesign_famoptg_example$p_rep_reps,
  p_rep_families          = x$OptiDesign_famoptg_example$p_rep_families,
  unreplicated_treatments = x$OptiDesign_famoptg_example$unreplicated_treatments,
  unreplicated_families   = x$OptiDesign_famoptg_example$unreplicated_families,
  n_blocks                = x$OptiDesign_famoptg_example$n_blocks,
  n_rows                  = x$OptiDesign_famoptg_example$n_rows,
  n_cols                  = x$OptiDesign_famoptg_example$n_cols,
  # Evaluation arguments
  treatment_effect   = "fixed",
  residual_structure = "IID",
  # Optimiser arguments
  criterion   = "A",
  n_restarts  = 20,
  verbose_opt = FALSE
)

cat("Best A-criterion:", round(opt_fam$optimization$best_score, 4), "\n")
cat("Valid restarts:", opt_fam$optimization$n_restarts -
                        opt_fam$optimization$n_failed, "/",
                        opt_fam$optimization$n_restarts, "\n")
```

------------------------------------------------------------------------

## Workflow 2: Alpha Row-Column Stream Design (Family-Based)

### Step 1 — Construct

``` r
design_alpha <- do.call(
  alpha_rc_stream,
  c(x$OptiDesign_alpha_example, x$OptiDesign_alpha_args_family)
)

dim(design_alpha$layout_matrix)
#> [1] 12 14
cat("Blocks per rep:", design_alpha$design_info$n_blocks_per_rep, "\n")
#> Blocks per rep: 8
cat("Total used plots:", design_alpha$design_info$total_used_plots, "\n")
#> Total used plots: 168
cat("Trailing NA plots:", design_alpha$design_info$trailing_na_plots, "\n")
#> Trailing NA plots: 0
head(design_alpha$field_book)
#>   Plot Row Column Rep IBlock BlockInRep Treatment Family Gcluster Check
#> 1    1   1      1   1      1          1      L146    F08     <NA> FALSE
#> 2    2   1      2   1      1          1      L118    F14     <NA> FALSE
#> 3    3   1      3   1      1          1      L134    F09     <NA> FALSE
#> 4    4   1      4   1      1          1      L153    F20     <NA> FALSE
#> 5    5   1      5   1      1          1      L145    F01     <NA> FALSE
#> 6    6   1      6   1      1          1      L117    F06     <NA> FALSE
```

### Step 2 — Evaluate

``` r
eff_alpha <- evaluate_alpha_efficiency(
  field_book         = design_alpha$field_book,
  n_rows             = x$OptiDesign_alpha_example$n_rows,
  n_cols             = x$OptiDesign_alpha_example$n_cols,
  check_treatments   = x$OptiDesign_alpha_example$check_treatments,
  treatment_effect   = "fixed",
  residual_structure = "IID"
)

cat("A-criterion (lower is better):", round(eff_alpha$A_criterion, 4), "\n")
#> A-criterion (lower is better): 1.3017
cat("D-criterion (lower is better):", round(eff_alpha$D_criterion, 4), "\n")
#> D-criterion (lower is better): 0.5794
cat("Number of treatments evaluated:", eff_alpha$n_trt, "\n")
#> Number of treatments evaluated: 63
```

### Step 3 — Optimise with Random Restart

``` r
opt_rs <- optimize_alpha_rc(
  check_treatments   = x$OptiDesign_alpha_example$check_treatments,
  check_families     = x$OptiDesign_alpha_example$check_families,
  entry_treatments   = x$OptiDesign_alpha_example$entry_treatments,
  entry_families     = x$OptiDesign_alpha_example$entry_families,
  n_reps             = x$OptiDesign_alpha_example$n_reps,
  n_rows             = x$OptiDesign_alpha_example$n_rows,
  n_cols             = x$OptiDesign_alpha_example$n_cols,
  min_block_size     = 10L,
  max_block_size     = 12L,
  treatment_effect   = "fixed",
  residual_structure = "IID",
  method             = "RS",
  criterion          = "A",
  n_restarts         = 20,
  verbose_opt        = FALSE
)

cat("Best A-criterion (RS):", round(opt_rs$optimization$best_score, 4), "\n")
plot(opt_rs$optimization$score_history, type = "b",
     xlab = "Restart", ylab = "A-criterion",
     main = "RS: A-criterion across restarts")
```

### Step 4 — Optimise with Simulated Annealing

``` r
opt_sa <- optimize_alpha_rc(
  check_treatments   = x$OptiDesign_alpha_example$check_treatments,
  check_families     = x$OptiDesign_alpha_example$check_families,
  entry_treatments   = x$OptiDesign_alpha_example$entry_treatments,
  entry_families     = x$OptiDesign_alpha_example$entry_families,
  n_reps             = x$OptiDesign_alpha_example$n_reps,
  n_rows             = x$OptiDesign_alpha_example$n_rows,
  n_cols             = x$OptiDesign_alpha_example$n_cols,
  min_block_size     = 10L,
  max_block_size     = 12L,
  treatment_effect   = "fixed",
  residual_structure = "IID",
  method             = "SA",
  criterion          = "A",
  n_restarts         = 3,
  sa_max_iter        = 200,
  sa_temp_start      = 1.0,
  sa_temp_end        = 0.001,
  sa_cooling         = "exponential",
  sa_swap_scope      = "global",
  verbose_opt        = FALSE
)

cat("Best A-criterion (SA):", round(opt_sa$optimization$best_score, 4), "\n")
```

------------------------------------------------------------------------

## Workflow 3: GRM-Based Design with Dispersion

When a genomic relationship matrix is available it replaces family
labels as the grouping source, and can also drive the optional
dispersion optimisation:

``` r
design_grm <- do.call(
  alpha_rc_stream,
  c(x$OptiDesign_alpha_example, x$OptiDesign_alpha_args_grm)
)

# GRM clustering populates the Gcluster column for non-check entries
non_check <- design_grm$field_book[!design_grm$field_book$Check, ]
cat("Unique genomic clusters:", length(unique(non_check$Gcluster[
  !is.na(non_check$Gcluster)])), "\n")
#> Unique genomic clusters: 23
# alpha_rc_stream field books use IBlock and Rep, not Block
head(non_check[, c("Treatment", "Family", "Gcluster", "Rep",
                    "IBlock", "Row", "Column")])
#>   Treatment Family Gcluster Rep IBlock Row Column
#> 1      L145    F01      G14   1      1   1      1
#> 2      L146    F08       G7   1      1   1      2
#> 3      L154    F08       G1   1      1   1      3
#> 7      L134    F09       G2   1      1   1      7
#> 8      L104    F22      G20   1      1   1      8
#> 9      L117    F06      G22   1      1   1      9
```

------------------------------------------------------------------------

## Workflow 4: CDmean Optimisation for Genomic Selection

CDmean is the primary criterion when the objective is maximising GEBV
prediction reliability rather than contrast precision. It requires
random treatment effects and a genomic prediction model.

``` r
# CDmean optimisation with GBLUP requires a K matrix.
# Here we use the example K from the shipped dataset.
opt_cdmean <- optimize_alpha_rc(
  check_treatments   = x$OptiDesign_alpha_example$check_treatments,
  check_families     = x$OptiDesign_alpha_example$check_families,
  entry_treatments   = x$OptiDesign_alpha_example$entry_treatments,
  entry_families     = x$OptiDesign_alpha_example$entry_families,
  n_reps             = x$OptiDesign_alpha_example$n_reps,
  n_rows             = x$OptiDesign_alpha_example$n_rows,
  n_cols             = x$OptiDesign_alpha_example$n_cols,
  min_block_size     = 10L,
  max_block_size     = 12L,
  # Genomic prediction model
  treatment_effect   = "random",
  prediction_type    = "GBLUP",
  K                  = x$OptiDesign_K,
  line_id_map        = x$OptiDesign_id_map,
  varcomp = list(
    sigma_g2   = 0.4,
    sigma_e2   = 0.6,
    sigma_rep2 = 0.1,
    sigma_ib2  = 0.05,
    sigma_r2   = 0.02,
    sigma_c2   = 0.02
  ),
  # CDmean criterion
  method      = "RS",
  criterion   = "CDmean",
  n_restarts  = 20,
  verbose_opt = FALSE
)

cat("Best CDmean:", round(opt_cdmean$optimization$best_score, 4),
    "(higher is better)\n")
cat("CDmean from efficiency slot:",
    round(opt_cdmean$efficiency$CDmean, 4), "\n")
```

CDmean and A-criterion can point in different directions because they
measure different objectives — contrast precision versus prediction
reliability. Choose based on your trial objective.

------------------------------------------------------------------------

## Comparing Optimality Criteria

Both criterion types can be computed on the same design to understand
the trade-off:

``` r
# Evaluate the same design under both fixed and random models
eff_fixed <- evaluate_alpha_efficiency(
  field_book         = design_alpha$field_book,
  n_rows             = x$OptiDesign_alpha_example$n_rows,
  n_cols             = x$OptiDesign_alpha_example$n_cols,
  check_treatments   = x$OptiDesign_alpha_example$check_treatments,
  treatment_effect   = "fixed",
  residual_structure = "IID"
)

eff_random <- evaluate_alpha_efficiency(
  field_book         = design_alpha$field_book,
  n_rows             = x$OptiDesign_alpha_example$n_rows,
  n_cols             = x$OptiDesign_alpha_example$n_cols,
  check_treatments   = x$OptiDesign_alpha_example$check_treatments,
  treatment_effect   = "random",
  prediction_type    = "IID",
  residual_structure = "IID"
)

results <- data.frame(
  Criterion  = c("A-criterion", "D-criterion",
                 "A-efficiency", "D-efficiency",
                 "Mean PEV", "CDmean"),
  Value      = c(
    round(eff_fixed$A_criterion,  4),
    round(eff_fixed$D_criterion,  4),
    round(eff_fixed$A_efficiency, 4),
    round(eff_fixed$D_efficiency, 4),
    round(eff_random$mean_PEV,    4),
    round(eff_random$CDmean,      4)
  ),
  Direction  = c("lower=better", "lower=better",
                 "higher=better", "higher=better",
                 "lower=better", "higher=better"),
  Model      = c("fixed", "fixed", "fixed", "fixed", "random", "random")
)
knitr::kable(results, caption = "Efficiency criteria for the example alpha design")
```

| Criterion    |  Value | Direction     | Model  |
|:-------------|-------:|:--------------|:-------|
| A-criterion  | 1.3017 | lower=better  | fixed  |
| D-criterion  | 0.5794 | lower=better  | fixed  |
| A-efficiency | 0.7682 | higher=better | fixed  |
| D-efficiency | 1.7260 | higher=better | fixed  |
| Mean PEV     | 0.3883 | lower=better  | random |
| CDmean       | 0.6117 | higher=better | random |

Efficiency criteria for the example alpha design

------------------------------------------------------------------------

## Grouping and Dispersion Options

### When to use each grouping strategy

| Strategy      | `cluster_source` | Use when                                                                           |
|---------------|------------------|------------------------------------------------------------------------------------|
| Family labels | `"Family"`       | Family structure is meaningful and interpretable; no relationship matrix available |
| Genomic (GRM) | `"GRM"`          | Genomic data is available; relatedness precision matters more than family labels   |
| Pedigree (A)  | `"A"`            | Only pedigree information is available                                             |

### Dispersion optimisation

The dispersion step minimises the total genomic relatedness among
neighbouring plots:

$$S = \sum\limits_{{(i,j)} \in \mathcal{N}}K_{ij}$$

where $\mathcal{N}$ is the set of plot pairs within Chebyshev distance
`dispersion_radius`. The local swap search accepts a swap if and only if
it reduces $S$. Key parameters:

| Parameter           | Effect                                                               |
|---------------------|----------------------------------------------------------------------|
| `dispersion_radius` | Neighbourhood size: 1 = 8-connected, 2 = 24-connected                |
| `dispersion_iters`  | Number of swap proposals; more iterations = lower $S$ at linear cost |
| `dispersion_source` | Which matrix to score against: `"K"`, `"GRM"`, or `"A"`              |

------------------------------------------------------------------------

## Practical Guidelines

### Choosing between design families

| Situation                                                | Recommended function                                                                            |
|----------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| Many entries, not all need replication                   | [`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)             |
| Some entries need priority replication                   | [`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md) (p-rep)     |
| All entries equally replicated, checks needed everywhere | [`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md) (RCBD-type) |
| Fixed field dimensions, operational field-book order     | [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md)       |
| Alpha-lattice structure needed                           | [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md)       |

### Choosing an optimality criterion

| Objective                            | Criterion   | Function argument      |
|--------------------------------------|-------------|------------------------|
| Maximise contrast precision          | A-criterion | `criterion = "A"`      |
| Minimise overall estimation volume   | D-criterion | `criterion = "D"`      |
| Balance both                         | Combined    | `criterion = "both"`   |
| Maximise GEBV prediction reliability | CDmean      | `criterion = "CDmean"` |

CDmean requires `treatment_effect = "random"` and
`prediction_type %in% c("IID", "GBLUP", "PBLUP")`. For genomic selection
training population optimisation, GBLUP with a real kinship matrix $K$
is strongly recommended.

### Choosing an optimisation method (`optimize_alpha_rc` only)

| Method | Best for                                     | Cost   |
|--------|----------------------------------------------|--------|
| RS     | Quick exploration, guaranteed valid designs  | Low    |
| SA     | Moderate search depth, escaping local optima | Medium |
| GA     | Thorough global search                       | High   |

For production use, SA or GA with 5–10 restarts and 500–1000 iterations/
generations typically provides a good balance of quality and computation
time.

### Variance component specification

Variance components have a large effect on efficiency values but a small
effect on the *ranking* of designs. For design comparison purposes,
equal variance components (`sigma_* = 1`) are a reasonable default. For
absolute criterion values that are meaningful on the scale of real data,
use heritability-consistent components:

``` r
# Example: h^2 = 0.5, moderate spatial correlation
varcomp <- list(
  sigma_g2   = 0.5,
  sigma_e2   = 0.5,
  sigma_rep2 = 0.1,
  sigma_ib2  = 0.05,
  sigma_r2   = 0.02,
  sigma_c2   = 0.02
)
```

------------------------------------------------------------------------

## Session Information

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] OptiDesign_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] Matrix_1.7-4      xfun_0.57         lattice_0.22-9    cachem_1.1.0     
#>  [9] knitr_1.51        htmltools_0.5.9   rmarkdown_2.31    lifecycle_1.0.5  
#> [13] cli_3.6.5         grid_4.5.3        sass_0.4.10       pkgdown_2.2.0    
#> [17] textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2 compiler_4.5.3   
#> [21] tools_4.5.3       ragg_1.5.2        pracma_2.4.6      evaluate_1.0.5   
#> [25] bslib_0.10.0      yaml_2.3.12       jsonlite_2.0.0    rlang_1.1.7      
#> [29] fs_2.0.1
```

------------------------------------------------------------------------

## References

Rincent, R., Laloë, D., Nicolas, S., Altmann, T., Brunel, D., Revilla,
P., …, & Moreau, L. (2012). Maximizing the reliability of genomic
selection by optimizing the calibration set of reference individuals:
comparison of methods in two diverse groups of maize inbreds (*Zea mays*
L.). *Genetics*, 192(2), 715–728.
<https://doi.org/10.1534/genetics.112.141473>

Jones, B., Allen-Moyer, K., & Goos, P. (2021). A-optimal versus
D-optimal design of screening experiments. *Journal of Quality
Technology*, 53(4), 369–382.
<https://doi.org/10.1080/00224065.2020.1757391>

Hutchinson, M. F. (1990). A stochastic estimator of the trace of the
influence matrix for Laplacian smoothing splines. *Communications in
Statistics — Simulation and Computation*, 19(2), 433–450.
<https://doi.org/10.1080/03610919008812866>

Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983). Optimization by
simulated annealing. *Science*, 220(4598), 671–680.
<https://doi.org/10.1126/science.220.4598.671>

Holland, J. H. (1992). *Adaptation in Natural and Artificial Systems*.
MIT Press.
