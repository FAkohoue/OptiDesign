# Create a fixed-grid alpha row-column design with unequal incomplete blocks and checks in every block

`alpha_rc_stream()` is intended for practical field situations where:

- the overall field size is fixed in advance,

- replicate boundaries are determined by field-book order rather than by
  rigid rectangular subfields,

- replicate 2 begins exactly where replicate 1 ends in the traversal
  stream,

- incomplete blocks may differ in size,

- checks must be repeated in every incomplete block,

- all entries must appear exactly once in each replicate,

- any leftover field cells should remain unused and appear at the end of
  the stream.

The function:

1.  Builds the full field stream using `order` and `serpentine`.

2.  Determines how many incomplete blocks can be supported in each
    replicate.

3.  Divides each replicate into incomplete blocks.

4.  Allocates checks and entries subject to the block plan.

5.  Arranges entries heuristically to reduce clustering of similar
    materials.

6.  Leaves any surplus cells as trailing `NA`.

7.  Optionally performs dispersion optimization.

8.  Optionally computes mixed-model efficiency diagnostics.

## Usage

``` r
alpha_rc_stream(
  check_treatments,
  check_families,
  entry_treatments,
  entry_families,
  n_reps,
  n_rows,
  n_cols,
  order = "column",
  serpentine = FALSE,
  seed = NULL,
  attempts = 5000,
  warn_and_correct = TRUE,
  fix_rows = TRUE,
  cluster_source = c("Family", "GRM", "A"),
  GRM = NULL,
  A = NULL,
  id_map = NULL,
  cluster_method = c("kmeans", "hclust"),
  cluster_seed = 1,
  cluster_attempts = 25,
  n_pcs_use = Inf,
  min_entry_slots_per_block = 8,
  max_blocks_per_rep = NULL,
  eval_efficiency = FALSE,
  treatment_effect = c("random", "fixed"),
  prediction_type = c("none", "IID", "GBLUP", "PBLUP"),
  K = NULL,
  line_id_map = NULL,
  varcomp = list(sigma_e2 = 1, sigma_g2 = 1, sigma_rep2 = 1, sigma_ib2 = 1, sigma_r2 = 1,
    sigma_c2 = 1),
  check_as_fixed = TRUE,
  residual_structure = c("IID", "AR1", "AR1xAR1"),
  rho_row = 0,
  rho_col = 0,
  spatial_engine = c("auto", "sparse", "dense"),
  dense_max_n = 5000,
  eff_trace_samples = 80,
  eff_full_max = 400,
  check_placement = c("systematic", "random"),
  check_position_pattern = c("spread", "corners_first"),
  use_dispersion = FALSE,
  dispersion_source = c("K", "A", "GRM"),
  dispersion_radius = 1,
  dispersion_iters = 2000,
  dispersion_seed = 1,
  verbose = TRUE
)
```

## Arguments

- check_treatments:

  Character vector of check IDs.

  Every check is placed exactly once in every incomplete block.

  This argument is always active and is central to the design because
  repeated checks define the benchmark structure across blocks and
  replicates.

  Use this when the trial requires standard controls or benchmark
  entries in every incomplete block.

  Example use case: a trial with 3 commercial checks that must appear in
  all incomplete blocks.

- check_families:

  Character vector of the same length as `check_treatments`.

  Family labels for checks.

  These are used directly when `cluster_source = "Family"` and are also
  stored in the output field book under all settings.

  This argument must align exactly with `check_treatments`.

  Use it when checks belong to known families or control groups, or when
  the user wants family labels preserved in the output.

- entry_treatments:

  Character vector of distinct entries.

  Each entry appears exactly once **per replicate**.

  This is the main set of test entries being evaluated.

  Use this for breeding lines, hybrids, progeny, or other candidate
  materials that should be represented once in each replicate.

  Entries must not overlap with `check_treatments`.

- entry_families:

  Character vector of the same length as `entry_treatments`.

  Family labels for entries.

  These are used directly when `cluster_source = "Family"`. When
  `cluster_source %in% c("GRM", "A")`, they still help determine the
  target number of clusters among entries.

  Use this when family structure matters for spacing, grouping, or
  interpretation.

- n_reps:

  Integer giving the number of replicates.

  Each replicate is defined as a contiguous segment of the global
  stream.

  Use more replicates when stronger replication is needed for entry
  precision, provided the field can still support repeated checks in
  every block.

  This argument affects:

  - total number of used plots,

  - the length of each replicate segment,

  - the total number of times entries appear.

- n_rows:

  Integer giving the number of field rows.

  The field size is fixed in this function and is never altered.

  Use this to reflect the true row dimension of the physical field.

  It directly affects:

  - total field capacity,

  - the global traversal stream,

  - spatial coordinates in the returned `layout_matrix` and
    `field_book`.

- n_cols:

  Integer giving the number of field columns.

  The field size is fixed in this function and is never altered.

  Use this to reflect the true column dimension of the physical field.

  Together with `n_rows`, it determines field capacity and stream
  structure.

- order:

  Character specifying the global traversal order.

  Allowed values:

  - `"column"`: column-major traversal,

  - `"row"`: row-major traversal.

  This determines:

  - the order of the global field stream,

  - where replicate boundaries occur,

  - where incomplete block boundaries occur,

  - where trailing `NA` cells are placed.

  Use `"row"` when operational movement follows rows. Use `"column"`
  when operational movement follows columns.

  This argument interacts directly with `serpentine`.

- serpentine:

  Logical indicating whether alternate rows or columns reverse direction
  during stream generation.

  If `TRUE`, traversal alternates direction:

  - by row when `order = "row"`,

  - by column when `order = "column"`.

  Use this when field-book order should mimic serpentine movement in the
  field.

  This changes the stream order, which in turn changes replicate and
  incomplete block boundaries.

  Depends on:

  - `order`

- seed:

  Optional integer seed for reproducibility.

  If `NULL`, a seed is generated internally and returned as `seed_used`.

  Use this when the layout should be reproducible.

  It affects:

  - entry allocation,

  - within-block arrangement,

  - random check placement,

  - optional dispersion optimization (unless a separate dispersion seed
    is supplied).

- attempts:

  Integer giving the number of swap proposals used when improving entry
  allocation across blocks within each replicate.

  Larger values increase search effort but do not change the design
  algorithm itself.

  Use larger values when:

  - family structure is strongly imbalanced,

  - the user wants better separation of similar materials,

  - the number of entries per replicate is large.

- warn_and_correct:

  Logical retained for interface continuity.

  In this function, the field dimensions are fixed and are **not**
  altered.

  This argument is included so the interface stays similar to related
  design functions, but it does not drive field resizing here.

- fix_rows:

  Logical retained for interface continuity.

  In this function, the field dimensions are fixed and are **not**
  altered.

  Included mainly for consistency with related interfaces.

- cluster_source:

  Character specifying the grouping source used during allocation and
  arrangement.

  Allowed values:

  - `"Family"`

  - `"GRM"`

  - `"A"`

  Use `"Family"` when user-supplied family labels are adequate and
  interpretable.

  Use `"GRM"` when genomic structure should define similarity more
  precisely.

  Use `"A"` when pedigree relationships are available and should define
  grouping.

  This argument determines whether `GRM`, `A`, `id_map`,
  `cluster_method`, `cluster_seed`, `cluster_attempts`, and `n_pcs_use`
  become active.

- GRM:

  Optional genomic relationship matrix.

  Required when:

  - `cluster_source = "GRM"`, or

  - `use_dispersion = TRUE` and `dispersion_source = "GRM"`.

  Ignored otherwise.

  Use this when genomic relatedness should influence grouping or
  dispersion.

  Matrix row and column names must match entry IDs, or be reachable
  through `id_map`.

- A:

  Optional pedigree relationship matrix.

  Required when:

  - `cluster_source = "A"`, or

  - `use_dispersion = TRUE` and `dispersion_source = "A"`.

  Ignored otherwise.

  Use this when pedigree relatedness should influence grouping or
  dispersion.

- id_map:

  Optional `data.frame` with columns `Treatment` and `LineID`.

  Used only when `cluster_source %in% c("GRM", "A")` and treatment IDs
  do not already match matrix row/column names.

  Use this when field-book treatment labels are different from matrix
  IDs.

  Example: treatment labels may be breeder-friendly names while `GRM` or
  `A` uses internal IDs.

- cluster_method:

  Character specifying the clustering method applied after PCA in
  matrix-based grouping.

  Allowed values:

  - `"kmeans"`

  - `"hclust"`

  Use `"kmeans"` when reproducible compact clustering is preferred.

  Use `"hclust"` when a hierarchical approach is preferred.

  Active only when `cluster_source %in% c("GRM", "A")`.

- cluster_seed:

  Integer seed for k-means initialization.

  Active only when:

  - `cluster_source %in% c("GRM", "A")`, and

  - `cluster_method = "kmeans"`.

  Use this to make matrix-based grouping reproducible.

- cluster_attempts:

  Integer number of random starts for k-means clustering.

  Active only when:

  - `cluster_source %in% c("GRM", "A")`, and

  - `cluster_method = "kmeans"`.

  Use larger values when clustering stability is important.

- n_pcs_use:

  Integer or `Inf` giving the number of principal components used for
  matrix-based clustering.

  Ignored when `cluster_source = "Family"`.

  Use smaller values when only broad structure should drive grouping.
  Use `Inf` when the function should use as many informative components
  as possible.

- min_entry_slots_per_block:

  Integer giving the minimum number of entry slots allowed in any
  incomplete block.

  This argument constrains the automatically chosen number of blocks per
  replicate.

  Use smaller values when more blocks per replicate are acceptable. Use
  larger values when blocks should not become too small once checks are
  inserted.

  This is particularly important when checks are numerous, because
  repeated checks consume capacity in every incomplete block.

- max_blocks_per_rep:

  Optional integer giving an upper bound on the number of incomplete
  blocks per replicate.

  If `NULL`, the number of blocks per replicate is derived from capacity
  constraints.

  Use this when the user wants to prevent the design from creating too
  many incomplete blocks, even if more could fit mathematically.

  Useful when operational simplicity or analysis preference favors fewer
  blocks.

- eval_efficiency:

  Logical indicating whether efficiency diagnostics should be computed.

  Use `TRUE` when the user wants to compare design quality under a mixed
  model.

  Use `FALSE` when only the layout itself is needed.

  If `FALSE`, all efficiency-related arguments are ignored.

- treatment_effect:

  Character indicating how entries are treated in efficiency evaluation.

  Allowed values:

  - `"random"`

  - `"fixed"`

  Active only when `eval_efficiency = TRUE`.

  Use `"fixed"` when interest lies in contrast precision among entry
  means.

  Use `"random"` when interest lies in BLUP-style prediction quality.

- prediction_type:

  Character controlling the random-effect efficiency model.

  Allowed values:

  - `"none"`

  - `"IID"`

  - `"GBLUP"`

  - `"PBLUP"`

  Active only when:

  - `eval_efficiency = TRUE`, and

  - `treatment_effect = "random"`.

  Use `"none"` to skip random-effect efficiency.

  Use `"IID"` when entry random effects are assumed independent.

  Use `"GBLUP"` or `"PBLUP"` when prediction should use a relationship
  matrix `K`.

- K:

  Optional relationship matrix used in two contexts:

  1.  random-effect efficiency when:

      - `eval_efficiency = TRUE`,

      - `treatment_effect = "random"`,

      - `prediction_type %in% c("GBLUP", "PBLUP")`;

  2.  dispersion optimization when:

      - `use_dispersion = TRUE`,

      - `dispersion_source = "K"`.

  Ignored otherwise.

  Use this when prediction or dispersion should be based on an explicit
  relationship matrix.

- line_id_map:

  Optional `data.frame` with columns `Treatment` and `LineID`.

  Used only when `K` is active and treatment labels differ from
  `rownames(K)`.

  This serves the same role for `K` that `id_map` serves for `GRM` or
  `A`.

- varcomp:

  Named list of variance components used only in efficiency evaluation.

  Must contain:

  - `sigma_e2`

  - `sigma_g2`

  - `sigma_rep2`

  - `sigma_ib2`

  - `sigma_r2`

  - `sigma_c2`

  Use realistic values when efficiency should reflect a plausible trial
  model.

  Example: use larger `sigma_ib2` when incomplete blocks are expected to
  differ strongly.

- check_as_fixed:

  Logical indicating whether checks are included as fixed indicators
  during efficiency evaluation.

  Active only when `eval_efficiency = TRUE`.

  Use `TRUE` when checks should be explicitly modeled as benchmark fixed
  effects.

- residual_structure:

  Character specifying the residual correlation model.

  Allowed values:

  - `"IID"`

  - `"AR1"`

  - `"AR1xAR1"`

  Active only when `eval_efficiency = TRUE`.

  Use `"IID"` when no spatial residual correlation is assumed.

  Use `"AR1"` when row-wise spatial correlation is expected.

  Use `"AR1xAR1"` when both row and column spatial correlation are
  expected.

- rho_row:

  Numeric AR1 row parameter.

  Active only when:

  - `eval_efficiency = TRUE`, and

  - `residual_structure %in% c("AR1", "AR1xAR1")`.

  Use values near 0 for weak row correlation and larger absolute values
  for stronger correlation.

- rho_col:

  Numeric AR1 column parameter.

  Active only when:

  - `eval_efficiency = TRUE`, and

  - `residual_structure = "AR1xAR1"`.

  Use this when column-wise correlation is expected in addition to
  row-wise correlation.

- spatial_engine:

  Character retained for interface compatibility.

  Allowed values:

  - `"auto"`

  - `"sparse"`

  - `"dense"`

  In this function the efficiency code primarily uses sparse Matrix
  operations, but the argument is retained so the interface remains
  aligned with related functions.

- dense_max_n:

  Integer retained for interface compatibility.

  This argument is included for consistency with related functions.

- eff_trace_samples:

  Integer number of Hutchinson trace samples used when approximate
  efficiency is needed because the target treatment dimension exceeds
  `eff_full_max`.

  Larger values improve approximation stability but increase runtime.

- eff_full_max:

  Integer maximum target dimension for exact efficiency extraction.

  Above this threshold, the function switches to an approximate
  trace-based method.

  Use larger values for more exact computation when memory and runtime
  allow.

- check_placement:

  Character specifying how check positions are chosen within blocks.

  Allowed values:

  - `"systematic"`

  - `"random"`

  Use `"systematic"` when checks should be approximately evenly spread
  in each block.

  Use `"random"` when check positions may vary freely.

- check_position_pattern:

  Argument retained for interface compatibility.

  It is not used by the current stream-based placement logic.

- use_dispersion:

  Logical; if `TRUE`, apply post-hoc dispersion optimization among
  non-check treatments.

  Use this when genetically or pedigree-similar entries should be less
  clustered within replicates.

  If `FALSE`, all dispersion-related arguments are ignored.

- dispersion_source:

  Character selecting which matrix is used for dispersion scoring.

  Allowed values:

  - `"K"`

  - `"A"`

  - `"GRM"`

  Active only when `use_dispersion = TRUE`.

  Use `"K"` when a dedicated relationship matrix is preferred. Use `"A"`
  for pedigree-based dispersion. Use `"GRM"` for genomic-based
  dispersion.

- dispersion_radius:

  Integer neighborhood radius for dispersion scoring.

  Two plots are considered neighbors when:

  `max(|dr|, |dc|) <= dispersion_radius`

  Use `1` for immediate neighbors only. Use larger values when broader
  local dispersion matters.

- dispersion_iters:

  Integer number of swap proposals used in the dispersion search.

  Larger values increase optimization effort and runtime.

- dispersion_seed:

  Optional integer seed used in dispersion search.

  Use this when dispersion optimization should be reproducible
  independently of the initial design seed.

  If left at its default, the function uses the supplied value in the
  interface.

- verbose:

  Logical; if `TRUE`, print the derived replicate and block structure.

  Use this for diagnostics, debugging, or understanding how the function
  chose replicate sizes and incomplete block sizes.

## Value

A list with:

- layout_matrix:

  A character matrix of size `n_rows × n_cols`. Unused cells are `NA`.

- field_book:

  A plot-level data frame with one row per field cell, including unused
  cells.

- efficiency:

  `NULL` or a list of efficiency summaries computed from non-`NA` plots.

- seed_used:

  The actual random seed used internally.

- design_info:

  A list summarizing replicate sizes, incomplete block sizes, used and
  unused cells, and key settings.

## Details

Construct a fixed-field row-column design in which the full
`n_rows × n_cols` grid is first converted into a single ordered stream
of positions, then split into contiguous replicate segments, then
subdivided into incomplete blocks of possibly unequal size. Within each
replicate, every entry appears exactly once, every incomplete block
contains all checks, and any unused field cells are left as trailing
`NA` positions at the end of the global stream.

Let:

- `E` = number of entries,

- `C` = number of checks,

- `R` = number of replicates,

- `b` = number of incomplete blocks per replicate.

Then the number of **used plots per replicate** is:

`E + b * C`

and the **total number of used plots** is:

`R * (E + b * C)`

This total must fit within the fixed field capacity:

`n_rows * n_cols`

The function chooses a block plan that respects:

- the number of entries,

- the number of checks,

- the minimum entry capacity per block,

- any user-supplied cap on blocks per replicate,

- the total fixed field capacity.

In efficiency evaluation, the function uses only the observed non-`NA`
plots. Depending on `treatment_effect`, the returned metrics represent
either:

- precision of fixed treatment contrasts, or

- average prediction error variance of random treatment effects.

The arguments `warn_and_correct`, `fix_rows`, `spatial_engine`,
`dense_max_n`, and `check_position_pattern` are retained for interface
continuity, even where they do not materially alter the stream-based
design logic.

## Conceptual design

The key feature of this function is that it is **stream-based** rather
than **rectangle-based**.

Instead of cutting the field into rectangular replicate areas, the
function:

- lists all grid cells in a single global order,

- cuts that order into `n_reps` contiguous replicate segments,

- then cuts each replicate segment into incomplete blocks.

This is useful when:

- planting follows a continuous field-book order,

- operational field movement is row-wise or column-wise,

- replicate boundaries are administrative or logistical rather than
  geometric,

- the user wants all unused cells to appear only at the tail of the
  field stream.

For each replicate:

- each entry is used exactly once,

- each incomplete block receives every check,

- the remaining block capacity is filled with entries,

- allocation is then improved heuristically to reduce local grouping
  conflicts.

## What this function is most useful for

This function is particularly useful when:

- the field geometry is fixed and should never be auto-resized,

- the design requires repeated checks in every incomplete block,

- the user wants a row-column design without forcing rectangular
  replicate areas,

- practical planting or harvesting follows a stream order,

- family, pedigree, or genomic structure should influence entry
  placement,

- the user wants to evaluate the resulting design using a mixed-model
  framework.

## Dependency guide

Many arguments are active only under particular modes.

**Grouping source**

- If `cluster_source = "Family"`:

  - grouping comes directly from `check_families` and `entry_families`,

  - `GRM`, `A`, `id_map`, `cluster_method`, `cluster_seed`,
    `cluster_attempts`, and `n_pcs_use` are ignored.

- If `cluster_source = "GRM"`:

  - `GRM` is required,

  - `id_map` is only needed if entry names differ from `rownames(GRM)`,

  - grouping labels are derived by PCA followed by clustering.

- If `cluster_source = "A"`:

  - `A` is required,

  - `id_map` is only needed if entry names differ from `rownames(A)`,

  - grouping labels are derived by PCA followed by clustering.

**Derived block structure**

- `min_entry_slots_per_block` controls the minimum number of non-check
  entry slots permitted in an incomplete block.

- `max_blocks_per_rep` optionally limits the number of blocks per
  replicate.

- If `max_blocks_per_rep = NULL`, block count is derived from capacity
  and `min_entry_slots_per_block`.

**Efficiency evaluation**

- If `eval_efficiency = FALSE`, all efficiency arguments are ignored.

- If `eval_efficiency = TRUE` and `treatment_effect = "fixed"`:

  - fixed-treatment contrast precision is computed,

  - `prediction_type`, `K`, and `line_id_map` are ignored.

- If `eval_efficiency = TRUE` and `treatment_effect = "random"`:

  - `prediction_type` becomes active.

- If `prediction_type = "none"`:

  - random-effect efficiency is skipped.

- If `prediction_type = "IID"`:

  - entry random effects are treated as independent,

  - `K` and `line_id_map` are ignored.

- If `prediction_type %in% c("GBLUP", "PBLUP")`:

  - `K` is required,

  - `line_id_map` may be required if treatment IDs differ from
    `rownames(K)`.

**Residual structure**

- If `residual_structure = "IID"`:

  - `rho_row` and `rho_col` are ignored.

- If `residual_structure = "AR1"`:

  - only `rho_row` is used.

- If `residual_structure = "AR1xAR1"`:

  - both `rho_row` and `rho_col` are used.

**Dispersion optimization**

- If `use_dispersion = FALSE`, all dispersion arguments are ignored.

- If `use_dispersion = TRUE`, `dispersion_source` selects which matrix
  is used.

- If `dispersion_source = "K"`:

  - `K` is required,

  - `line_id_map` may be needed if treatment IDs differ from matrix row
    names.

- If `dispersion_source = "A"`:

  - `A` is required.

- If `dispersion_source = "GRM"`:

  - `GRM` is required.

**Check placement**

- `check_placement = "systematic"` spreads checks approximately evenly
  within each block stream.

- `check_placement = "random"` randomizes check positions.

- `check_position_pattern` is retained for interface continuity but is
  not used by the current stream-based placement logic.

## Examples

``` r
data("OptiDesign_example_data", package = "OptiDesign")
x <- OptiDesign_example_data

## Family-based row-column stream design
out_alpha_family <- do.call(
  alpha_rc_stream,
  c(
    x$OptiDesign_alpha_example,
    x$OptiDesign_alpha_args_family
  )
)

dim(out_alpha_family$layout_matrix)
#> [1] 12 14
head(out_alpha_family$field_book)
#>   Treatment Family Gcluster Check PlotStream Rep IBlock BlockInRep Row Column
#> 1      L140    F24     <NA> FALSE          1   1      1          1   1      1
#> 2      L103    F14     <NA>  TRUE          2   1      1          1   1      2
#> 3      L102    F17     <NA>  TRUE          3   1      1          1   1      3
#> 4      L117    F06     <NA> FALSE          4   1      1          1   1      4
#> 5      L101    F14     <NA>  TRUE          5   1      1          1   1      5
#> 6      L133    F22     <NA> FALSE          6   1      1          1   1      6
out_alpha_family$design_info$n_blocks_per_rep
#> [1] 6

if (FALSE) { # \dontrun{
## GRM-based row-column stream design with dispersion
out_alpha_grm <- do.call(
  alpha_rc_stream,
  c(
    x$OptiDesign_alpha_example,
    x$OptiDesign_alpha_args_grm
  )
)

dim(out_alpha_grm$layout_matrix)
head(out_alpha_grm$field_book)
out_alpha_grm$efficiency
} # }
```
