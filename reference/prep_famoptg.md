# Create a repeated-check block design with flexible replication

`prep_famoptg()` is a general repeated-check block design constructor in
which **check treatments are included in every block**, while non-check
treatments are allocated across blocks according to user-specified
replication levels.

Depending on the supplied treatment structure, the same function can
generate:

- an **augmented design**, when only checks are repeated and all other
  entries are unreplicated;

- a **partially replicated (p-rep) design**, when some non-check entries
  are replicated and others are unreplicated;

- an **RCBD-type repeated-check design**, when all non-check entries are
  given the same replication number greater than 1, especially when that
  replication equals `n_blocks`.

The basic allocation rules are:

- **Checks** are included in **every block**.

- **Replicated non-check treatments** are assigned to a subset of blocks
  and may appear more than once overall, but **at most once within any
  single block**.

- **Unreplicated treatments** appear exactly once in the full design.

This means, for example, that if a treatment has replication 2 and
`n_blocks >= 2`, its two replicates are placed in **two different
blocks** rather than repeated within the same block.

The function returns:

- a **layout matrix** representing the field grid, and

- a **field book** giving one row per assigned plot with treatment
  identity, grouping metadata, block membership, and field coordinates.

Optionally, the function can also:

- derive treatment groups from a genomic relationship matrix (`GRM`) or
  a pedigree relationship matrix (`A`);

- optimize the spatial dispersion of checks;

- apply a swap-based local search to reduce local relatedness among
  nearby non-check plots;

- compute mixed-model design efficiency metrics for fixed or random
  treatment effects.

## Usage

``` r
prep_famoptg(
  check_treatments,
  check_families,
  p_rep_treatments,
  p_rep_reps,
  p_rep_families,
  unreplicated_treatments,
  unreplicated_families,
  n_blocks,
  n_rows,
  n_cols,
  order = "column",
  serpentine = FALSE,
  seed = NULL,
  attempts = 1000,
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
  eval_efficiency = FALSE,
  treatment_effect = c("random", "fixed"),
  prediction_type = c("none", "IID", "GBLUP", "PBLUP"),
  K = NULL,
  line_id_map = NULL,
  varcomp = list(sigma_e2 = 1, sigma_g2 = 1, sigma_b2 = 1, sigma_r2 = 1, sigma_c2 = 1),
  check_as_fixed = TRUE,
  residual_structure = c("IID", "AR1", "AR1xAR1"),
  rho_row = 0,
  rho_col = 0,
  spatial_engine = c("auto", "sparse", "dense"),
  dense_max_n = 5000,
  eff_trace_samples = 80,
  eff_full_max = 400,
  check_placement = c("random", "systematic", "optimal"),
  check_opt_attempts = 200,
  use_dispersion = FALSE,
  dispersion_source = c("K", "A", "GRM"),
  dispersion_radius = 1,
  dispersion_iters = 2000,
  dispersion_seed = 1
)
```

## Arguments

- check_treatments:

  Character vector of check treatment identifiers.

  These treatments are included in **every block**, so they define the
  repeated reference structure of the design.

  Use this argument when you have standard checks, control varieties,
  benchmark cultivars, or fixed reference entries that must be present
  everywhere to improve block comparability.

  This argument is always required.

  Example use case: a breeding trial with 3 control cultivars included
  in all blocks.

- check_families:

  Character vector of the same length as `check_treatments`.

  Provides the family or group label for each check treatment.

  Even when `cluster_source != "Family"`, checks still need family
  labels because checks are always assigned a group label in the output
  and remain part of the grouping metadata used internally.

  This argument always depends on `check_treatments` and must align
  exactly in order and length.

  Example: if the three checks all belong to the same control class,
  `check_families` might be `c("CHECK", "CHECK", "CHECK")`.

- p_rep_treatments:

  Character vector of replicated non-check treatment IDs.

  These are the entries that appear more than once overall, but at most
  once in any given block.

  Use this argument when some or all non-check treatments should be
  replicated.

  Depending on how it is used:

  - if it is empty and only `unreplicated_treatments` are supplied, the
    function behaves as an augmented repeated-check constructor;

  - if it contains only some non-check entries, the function behaves as
    a p-rep repeated-check constructor;

  - if it contains all non-check entries with common replication greater
    than 1, the function behaves as an RCBD-type repeated-check
    constructor.

  Can be `NULL` or `character(0)` if no replicated non-check treatments
  are present.

  Depends on:

  - `p_rep_reps`

  - `p_rep_families`

  Example use case: elite lines replicated twice while many other lines
  remain unreplicated.

- p_rep_reps:

  Integer vector giving the total number of replicates for each entry in
  `p_rep_treatments`.

  This controls how many times each replicated non-check treatment
  appears in the full design.

  Must:

  - have the same length as `p_rep_treatments`;

  - align element-wise with `p_rep_treatments`;

  - satisfy `p_rep_reps[i] <= n_blocks`, because a treatment cannot
    occur twice in the same block.

  A treatment with replication 2 and `n_blocks >= 2` will be placed in
  two different blocks.

  Use a common vector such as `rep(2, length(p_rep_treatments))` when
  all replicated treatments should appear twice.

  Use unequal values when some entries deserve more replication than
  others.

- p_rep_families:

  Character vector of the same length as `p_rep_treatments`.

  Family or group labels for replicated non-check treatments.

  These labels are required whenever replicated non-check treatments
  exist. Even when `cluster_source != "Family"`, they remain useful
  because the function uses non-check family structure to determine the
  target number of clusters in matrix-based grouping modes.

  Depends on:

  - `p_rep_treatments`

  Example use case: replicated candidate lines belong to multiple
  families and the user wants adjacency control to reflect that
  structure.

- unreplicated_treatments:

  Character vector of treatments that appear exactly once.

  These are single-plot entries used when broad screening is needed and
  full replication is not possible.

  Use `NULL` or `character(0)` if no unreplicated treatments are
  present.

  In combination with repeated checks and no replicated non-checks,
  these entries define an augmented repeated-check design.

  Depends on:

  - `unreplicated_families`

  Example use case: a large early-generation nursery where each new line
  is observed once.

- unreplicated_families:

  Character vector of the same length as `unreplicated_treatments`.

  Family or group labels for unreplicated treatments.

  Required whenever `unreplicated_treatments` is provided.

  These labels are especially useful when:

  - `cluster_source = "Family"`;

  - matrix-based grouping is used and the target number of clusters
    should reflect the family structure among non-checks.

- n_blocks:

  Integer giving the number of experimental blocks.

  This determines how many times checks are repeated, and it also
  constrains the maximum allowable replication count for any replicated
  non-check treatment.

  Use more blocks when:

  - many checks must be accommodated;

  - finer local control is desired;

  - field heterogeneity suggests stronger blocking.

  This argument affects:

  - total plot count;

  - feasibility of replicated-treatment assignment;

  - distribution of unreplicated treatments.

  When the same replicated treatment should appear in every block, its
  replication must equal `n_blocks`.

- n_rows:

  Integer giving the number of rows in the final field grid.

  Controls the physical field geometry returned in `layout_matrix`.

  Use this to match the intended number of field rows.

  If `warn_and_correct = TRUE`, this value may be retained while
  `n_cols` is adjusted when `fix_rows = TRUE`.

- n_cols:

  Integer giving the number of columns in the final field grid.

  Controls the physical field geometry returned in `layout_matrix`.

  Use this together with `n_rows` to represent the intended field
  dimensions.

  If `warn_and_correct = TRUE`, this value may be retained while
  `n_rows` is adjusted when `fix_rows = FALSE`.

- order:

  Character specifying how the field grid is filled.

  Allowed values:

  - `"row"`: fill row by row;

  - `"column"`: fill column by column.

  Use `"row"` when field-book order or planting order follows rows. Use
  `"column"` when field-book order follows columns.

  This argument interacts directly with `serpentine`.

- serpentine:

  Logical indicating whether alternate rows or columns should reverse
  direction during grid filling.

  If `TRUE`:

  - with `order = "row"`, alternate rows are reversed;

  - with `order = "column"`, alternate columns are reversed.

  Use this when planting, labeling, or harvesting follows a serpentine
  movement pattern.

  This affects only how the final ordered treatments are mapped into
  space; it does not change block composition.

  Depends on:

  - `order`

- seed:

  Optional integer seed controlling stochastic components of the design.

  Use this when reproducibility is important.

  If `NULL`, a seed is generated internally and returned as `seed_used`.

  This seed affects:

  - replicated-treatment block assignment;

  - shuffling within blocks;

  - candidate generation for optimal check placement;

  - optional dispersion optimization when `dispersion_seed` is `NULL`.

- attempts:

  Integer giving the maximum number of shuffle attempts per block.

  Used to reduce adjacency of same-group entries in the 1D block
  ordering.

  Larger values may improve adjacency avoidance but increase runtime.

  Use larger values when:

  - family structure is highly imbalanced;

  - many similar entries occur in the same block;

  - stricter adjacency avoidance is desired.

- warn_and_correct:

  Logical controlling what happens when the supplied field size does not
  match the required number of plots.

  If `FALSE`, the function stops with an error.

  If `TRUE`, the function adjusts one dimension upward so the field can
  contain all required plots.

  Use `FALSE` when the field dimensions are fixed and must not be
  changed. Use `TRUE` when a near-feasible layout should be repaired
  automatically.

  Depends on:

  - `fix_rows`

- fix_rows:

  Logical used only when `warn_and_correct = TRUE`.

  If `TRUE`, `n_rows` is kept fixed and `n_cols` is adjusted. If
  `FALSE`, `n_cols` is kept fixed and `n_rows` is adjusted.

  Use this according to which field dimension is physically fixed in
  practice.

  Example: if the field has a fixed number of rows but flexible row
  length, use `fix_rows = TRUE`.

- cluster_source:

  Character specifying the source of grouping used for non-check
  adjacency control.

  Allowed values:

  - `"Family"`

  - `"GRM"`

  - `"A"`

  Use `"Family"` when family labels are trusted and easy to interpret.

  Use `"GRM"` when genomic relatedness should define grouping more
  accurately.

  Use `"A"` when pedigree structure is available and pedigree-based
  grouping is desired.

  This argument determines which of the following become active:

  - `GRM`

  - `A`

  - `id_map`

  - `cluster_method`

  - `cluster_seed`

  - `cluster_attempts`

  - `n_pcs_use`

- GRM:

  Optional square genomic relationship matrix with rownames and
  colnames.

  Required when:

  - `cluster_source = "GRM"`, or

  - `use_dispersion = TRUE` and `dispersion_source = "GRM"`.

  Ignored otherwise.

  Use this when genomic similarity should drive grouping or dispersion
  scoring.

  Row and column names must correspond to treatment IDs or to IDs
  reachable via `id_map`.

- A:

  Optional square pedigree relationship matrix with rownames and
  colnames.

  Required when:

  - `cluster_source = "A"`, or

  - `use_dispersion = TRUE` and `dispersion_source = "A"`.

  Ignored otherwise.

  Use this when pedigree structure should drive grouping or dispersion
  scoring.

- id_map:

  Optional data frame with columns `Treatment` and `LineID`.

  Used only when `cluster_source %in% c("GRM", "A")` and treatment
  labels in the design do not exactly match the rownames of the selected
  relationship matrix.

  Use this when:

  - treatment names in the design are breeder-friendly labels;

  - but the matrix uses internal IDs, numeric codes, or standardized
    line names.

  Ignored when `cluster_source = "Family"`.

- cluster_method:

  Clustering method used after PCA in matrix-based grouping.

  Allowed values:

  - `"kmeans"`

  - `"hclust"`

  Use `"kmeans"` when compact clusters are acceptable and reproducible
  initialization is desired.

  Use `"hclust"` when hierarchical grouping is preferred.

  Active only when `cluster_source %in% c("GRM", "A")`.

- cluster_seed:

  Integer seed used when `cluster_method = "kmeans"`.

  Ignored otherwise.

  Use this when k-means clustering must be reproducible.

- cluster_attempts:

  Integer number of random starts passed to
  [`kmeans()`](https://rdrr.io/r/stats/kmeans.html).

  Active only when:

  - `cluster_source %in% c("GRM", "A")`, and

  - `cluster_method = "kmeans"`.

  Larger values may improve clustering stability at higher runtime cost.

- n_pcs_use:

  Integer or `Inf` giving the number of principal components used for
  PCA-based clustering in matrix-based grouping.

  Use a smaller value when:

  - only the leading structure should define grouping;

  - smaller PCs are considered mostly noise.

  Use `Inf` to let the function use as many informative PCs as
  available.

  Ignored when `cluster_source = "Family"`.

- eval_efficiency:

  Logical; if `TRUE`, compute design efficiency metrics.

  Use `TRUE` when the design should be evaluated under a mixed-model
  framework.

  Use `FALSE` when only the layout itself is needed.

  This argument activates:

  - `treatment_effect`

  - `prediction_type`

  - `K`

  - `line_id_map`

  - `varcomp`

  - `check_as_fixed`

  - `residual_structure`

  - `rho_row`

  - `rho_col`

  - `spatial_engine`

  - `dense_max_n`

  - `eff_trace_samples`

  - `eff_full_max`

- treatment_effect:

  Character indicating whether non-check treatments are treated as
  `"fixed"` or `"random"` for efficiency evaluation.

  Active only when `eval_efficiency = TRUE`.

  Use `"fixed"` when interest is in the precision of estimated treatment
  contrasts.

  Use `"random"` when interest is in prediction quality of treatment
  effects.

- prediction_type:

  Character controlling the random-effect prediction model.

  Allowed values:

  - `"none"`

  - `"IID"`

  - `"GBLUP"`

  - `"PBLUP"`

  Active only when:

  - `eval_efficiency = TRUE`, and

  - `treatment_effect = "random"`.

  Use `"none"` to skip random-effect efficiency.

  Use `"IID"` when all non-check effects are assumed independent.

  Use `"GBLUP"` or `"PBLUP"` when prediction should use a relationship
  matrix `K`.

- K:

  Optional square relationship matrix.

  Used in two contexts:

  1.  for efficiency evaluation when random effects are predicted using
      `"GBLUP"` or `"PBLUP"`;

  2.  for dispersion optimization when `dispersion_source = "K"`.

  Ignored otherwise.

  Use this when prediction or dispersion should be based on an explicit
  relationship matrix.

- line_id_map:

  Optional data frame with columns `Treatment` and `LineID`.

  Used only when `K` is active and treatment names differ from
  `rownames(K)`.

  This is analogous to `id_map`, but specifically for the matrix `K`.

  Use this when the design uses one set of treatment labels while `K`
  uses another.

- varcomp:

  Named list of variance components used for efficiency evaluation.

  Must contain:

  - `sigma_e2`

  - `sigma_g2`

  - `sigma_b2`

  - `sigma_r2`

  - `sigma_c2`

  Active only when `eval_efficiency = TRUE`.

  Use realistic values when the efficiency metric should reflect a
  plausible data-generating model.

  Example use case: a user may set larger `sigma_e2` when expecting
  noisy phenotypes, or larger `sigma_b2` when block-to-block
  heterogeneity is substantial.

- check_as_fixed:

  Logical indicating whether checks are included as fixed indicator
  columns during efficiency evaluation.

  Active only when `eval_efficiency = TRUE`.

  Use `TRUE` when checks should be explicitly modeled as fixed benchmark
  effects.

- residual_structure:

  Character specifying the residual correlation model.

  Allowed values:

  - `"IID"`

  - `"AR1"`

  - `"AR1xAR1"`

  Use `"IID"` when no spatial residual correlation is assumed.

  Use `"AR1"` when correlation is expected mainly along rows.

  Use `"AR1xAR1"` when correlation is expected along both rows and
  columns.

  Active only when `eval_efficiency = TRUE`.

- rho_row:

  Numeric AR1 row correlation parameter.

  Active only when:

  - `eval_efficiency = TRUE`, and

  - `residual_structure %in% c("AR1", "AR1xAR1")`.

  Use values near 0 for weak row correlation and values closer to 1 in
  magnitude for stronger row correlation.

- rho_col:

  Numeric AR1 column correlation parameter.

  Active only when:

  - `eval_efficiency = TRUE`, and

  - `residual_structure = "AR1xAR1"`.

  Use this when column-wise spatial correlation is part of the assumed
  residual model.

- spatial_engine:

  Character controlling the computational implementation of spatial
  residual calculations.

  Allowed values:

  - `"auto"`

  - `"dense"`

  - `"sparse"`

  Active only when:

  - `eval_efficiency = TRUE`, and

  - `residual_structure != "IID"`.

  Use `"auto"` in most cases.

  Use `"dense"` for smaller problems when dense linear algebra is
  acceptable.

  Use `"sparse"` for larger problems where sparse precision matrices are
  preferred.

- dense_max_n:

  Integer threshold used when `spatial_engine = "auto"`.

  If the number of observed plots is at most this threshold, dense
  computation is used; otherwise sparse computation is used.

- eff_trace_samples:

  Integer giving the number of Hutchinson trace samples used when
  approximate efficiency evaluation is needed.

  Active only when:

  - `eval_efficiency = TRUE`, and

  - the target dimension exceeds `eff_full_max`.

  Larger values improve stability of the approximation but increase
  runtime.

- eff_full_max:

  Integer giving the maximum target dimension for exact inverse
  sub-block extraction in efficiency evaluation.

  Active only when `eval_efficiency = TRUE`.

  Use larger values when exact evaluation is desired and computation is
  affordable. Use smaller values to switch to approximation earlier and
  reduce runtime.

- check_placement:

  Character specifying how checks are placed within blocks.

  Allowed values:

  - `"random"`

  - `"systematic"`

  - `"optimal"`

  Use `"random"` when no explicit spacing pattern is required.

  Use `"systematic"` when checks should be spread approximately evenly
  within each block.

  Use `"optimal"` when spatial dispersion of checks in the final 2D grid
  is a priority.

  If `"optimal"` is selected, `check_opt_attempts` becomes active.

- check_opt_attempts:

  Integer number of candidate layouts evaluated when
  `check_placement = "optimal"`.

  Larger values may improve the final spatial spread of checks but
  increase runtime.

  Ignored unless `check_placement = "optimal"`.

- use_dispersion:

  Logical; if `TRUE`, apply post-layout swap optimization to reduce
  local relatedness among nearby non-check treatments.

  Use this when genetically or pedigree-similar entries should be less
  clustered in the final field.

  If `FALSE`, all dispersion-related arguments are ignored.

- dispersion_source:

  Character specifying which matrix is used during dispersion
  optimization.

  Allowed values:

  - `"K"`

  - `"A"`

  - `"GRM"`

  Active only when `use_dispersion = TRUE`.

  Use `"K"` when a dedicated relationship matrix for dispersion is
  available. Use `"A"` for pedigree-based spacing. Use `"GRM"` for
  genomic-based spacing.

- dispersion_radius:

  Integer neighborhood radius used in dispersion scoring.

  Neighboring plots are defined using Chebyshev distance, meaning plots
  are neighbors when `max(|dr|, |dc|) <= dispersion_radius`.

  Use `1` for immediate neighbors only. Use larger values when broader
  local neighborhoods matter.

  Active only when `use_dispersion = TRUE`.

- dispersion_iters:

  Integer number of swap proposals used during local search.

  Larger values allow a more extensive search for reduced local
  relatedness but increase runtime.

  Active only when `use_dispersion = TRUE`.

- dispersion_seed:

  Optional integer seed used for dispersion optimization.

  Active only when `use_dispersion = TRUE`.

  If `NULL`, the function uses `seed_used`.

  Use this when the user wants reproducibility of the dispersion step
  independently from the initial design seed.

## Value

A list with:

- layout_matrix:

  A character matrix of treatment IDs with dimensions `n_rows × n_cols`.
  Unused cells are `NA`.

- field_book:

  A data frame with columns `Treatment`, `Family`, `Gcluster`, `Block`,
  `Plot`, `Row`, and `Column`. Each row corresponds to an assigned plot.

- efficiency:

  `NULL` or a list of efficiency metrics and notes, depending on
  `eval_efficiency`.

- seed_used:

  The actual random seed used internally.

## Details

Construct an augmented, partially replicated, or RCBD-type
repeated-check field design with optional treatment grouping, optional
post-layout dispersion optimization, and optional mixed-model efficiency
diagnostics.

The function uses `mod()` imported from the **pracma** package for
serpentine traversal logic.

Field dimensions do not have to match the exact number of required plots
if `warn_and_correct = TRUE`. In that case, the function expands one
dimension so the field can contain all required plots. Extra cells
remain `NA` in the returned `layout_matrix`, but the `field_book`
includes only assigned plots.

The function separates:

- **construction logic** (how treatments are placed);

- **grouping logic** (how similarity is defined);

- **dispersion logic** (how nearby related entries are discouraged);

- **efficiency logic** (how the final design is evaluated).

This separation is important because a user may, for example:

- use family labels for adjacency control;

- use a genomic matrix for dispersion optimization;

- and choose not to compute efficiency;

or alternatively:

- use GRM-based clustering;

- skip dispersion optimization;

- and compute GBLUP-based efficiency using a separate matrix `K`.

## Design classes represented by the same function

**1. Augmented repeated-check design**

This is obtained when:

- `p_rep_treatments = NULL` or `character(0)`,

- `p_rep_reps = NULL` or `integer(0)`,

- `p_rep_families = NULL` or `character(0)`,

- non-check entries are supplied only through `unreplicated_treatments`.

In this case:

- checks are repeated in every block;

- all test entries appear once;

- the design behaves as an augmented repeated-check layout.

This setting is especially useful for early-stage screening when many
entries must be observed but only checks can be repeated systematically.

**2. Partially replicated repeated-check design**

This is obtained when:

- some non-check treatments are supplied in `p_rep_treatments`;

- those treatments have replication counts greater than 1;

- other entries may remain unreplicated.

This is the classical use case of the function: a mixture of repeated
checks, replicated candidate entries, and single-plot candidate entries.

**3. RCBD-type repeated-check design**

This is obtained when:

- all non-check treatments are supplied through `p_rep_treatments`;

- they all receive the same replication number greater than 1.

If that common replication equals `n_blocks`, then every non-check
treatment appears once in every block, which is the closest
repeated-check analogue of a classical RCBD under this framework.

If the common replication is less than `n_blocks`, the design remains
balanced and repeated, but it is not a strict classical RCBD because not
every treatment appears in every block.

## Conceptual workflow

The function proceeds in several stages:

1.  **Validate inputs** and reconcile field size with the required
    number of plots.

2.  **Prepare grouping information** from family labels or from
    clustering on `GRM` / `A`.

3.  **Assign replicated non-check entries to blocks** subject to the
    rule that no treatment appears twice in the same block.

4.  **Distribute unreplicated entries** across blocks.

5.  **Insert checks into each block** using the chosen placement
    strategy.

6.  **Shuffle treatment order within blocks** to reduce adjacency of
    same-group entries in the 1D block ordering.

7.  **Map the ordered treatments to the field grid** according to
    `order` and `serpentine`.

8.  Optionally apply **genetic dispersion optimization**.

9.  Optionally compute **design efficiency metrics**.

## What this function is most useful for

This function is especially useful for:

- early- to intermediate-stage breeding trials where many entries must
  be screened and not all can be equally replicated;

- augmented repeated-check designs with a large number of single-plot
  test entries;

- p-rep layouts mixing replicated and unreplicated entries;

- balanced repeated-check block designs in which all candidate entries
  are replicated;

- experiments where checks must appear in every block;

- layouts where related lines should not be clustered too closely;

- applications where family, pedigree, or genomic relationships should
  influence physical layout construction;

- cases where the user wants to compare alternative design choices using
  model-based efficiency metrics.

## Grouping logic and why it matters

Several parts of the function rely on the idea of a **group**:

- During block construction, groups are used to reduce adjacency of
  similar entries in the 1D ordering within each block.

- During optional dispersion optimization, relationship matrices are
  used to discourage close spatial placement of related non-check
  entries.

Group labels may come from:

- user-supplied families (`cluster_source = "Family"`);

- clusters derived from a genomic relationship matrix
  (`cluster_source = "GRM"`);

- clusters derived from a pedigree relationship matrix
  (`cluster_source = "A"`).

Family-based grouping is the simplest option and is usually appropriate
when family labels are meaningful and easy to interpret.

Matrix-based grouping is more useful when:

- family labels are too coarse;

- genomic or pedigree structure is more informative than nominal family
  labels;

- the objective is to spread genetically similar materials more
  effectively.

## Dependency guide

Many arguments are active only in certain modes.

**Treatment structure**

- If `p_rep_treatments` is empty (or `NULL`) and
  `unreplicated_treatments` is non-empty, the function behaves as an
  **augmented repeated-check design** constructor.

- If both `p_rep_treatments` and `unreplicated_treatments` are present,
  the function behaves as a **p-rep repeated-check design** constructor.

- If all non-check treatments are supplied in `p_rep_treatments` with a
  common replication number greater than 1, the function behaves as an
  **RCBD-type repeated-check design** constructor.

**Grouping mode**

- If `cluster_source = "Family"`:

  - `check_families`, `p_rep_families`, and `unreplicated_families` are
    used.

  - `GRM`, `A`, `id_map`, `cluster_method`, `cluster_seed`,
    `cluster_attempts`, and `n_pcs_use` are ignored.

- If `cluster_source = "GRM"`:

  - `GRM` is required.

  - `id_map` is needed only if treatment IDs do not match
    `rownames(GRM)`.

  - `cluster_method`, `cluster_seed`, `cluster_attempts`, and
    `n_pcs_use` are used.

  - `A` is ignored.

- If `cluster_source = "A"`:

  - `A` is required.

  - `id_map` is needed only if treatment IDs do not match `rownames(A)`.

  - `cluster_method`, `cluster_seed`, `cluster_attempts`, and
    `n_pcs_use` are used.

  - `GRM` is ignored.

**Efficiency evaluation**

- If `eval_efficiency = FALSE`, all efficiency-related arguments are
  ignored.

- If `eval_efficiency = TRUE` and `treatment_effect = "fixed"`:

  - fixed-effect design efficiency is computed;

  - `prediction_type`, `K`, and `line_id_map` are ignored.

- If `eval_efficiency = TRUE` and `treatment_effect = "random"`:

  - `prediction_type` becomes active.

- If `prediction_type = "none"`:

  - no random-effect prediction efficiency is computed.

- If `prediction_type = "IID"`:

  - non-check treatments are treated as independent random effects;

  - `K` and `line_id_map` are ignored.

- If `prediction_type %in% c("GBLUP", "PBLUP")`:

  - `K` is required;

  - `line_id_map` may be needed if treatment IDs do not match
    `rownames(K)`.

**Residual structure**

- If `residual_structure = "IID"`:

  - residuals are assumed independent;

  - `rho_row` and `rho_col` are ignored.

- If `residual_structure = "AR1"`:

  - row-wise correlation is used;

  - `rho_row` is active;

  - `rho_col` is ignored.

- If `residual_structure = "AR1xAR1"`:

  - row and column correlations are both used;

  - `rho_row` and `rho_col` are both active.

**Check placement**

- If `check_placement = "optimal"`, `check_opt_attempts` is used.

- Otherwise, `check_opt_attempts` is ignored.

**Dispersion optimization**

- If `use_dispersion = FALSE`, all dispersion arguments are ignored.

- If `use_dispersion = TRUE` and `dispersion_source = "K"`:

  - `K` is required;

  - `line_id_map` may be needed if treatment IDs differ from
    `rownames(K)`.

- If `use_dispersion = TRUE` and `dispersion_source = "A"`:

  - `A` is required.

- If `use_dispersion = TRUE` and `dispersion_source = "GRM"`:

  - `GRM` is required.

## Examples

``` r
data("OptiDesign_example_data", package = "OptiDesign")
x <- OptiDesign_example_data

## ---------------------------------------------------------
## Example 1: Family-based repeated-check design
## This example uses the shipped family-based arguments.
## Depending on the supplied treatment lists, the same function
## can represent augmented, p-rep, or RCBD-type repeated-check
## layouts.
## ---------------------------------------------------------
out_family <- do.call(
  prep_famoptg,
  c(
    x$OptiDesign_famoptg_example,
    x$OptiDesign_famoptg_args_family
  )
)

dim(out_family$layout_matrix)
#> [1] 16  8
head(out_family$field_book)
#>        Treatment Family Gcluster Block Plot Row Column
#> plot_1      L076    F03     <NA>     1    1   1      1
#> plot_2      L061    F16     <NA>     1    2   2      1
#> plot_3      L028    F12     <NA>     1    3   3      1
#> plot_4      L030    F10     <NA>     1    4   4      1
#> plot_5      L027    F21     <NA>     1    5   5      1
#> plot_6      L051    F09     <NA>     1    6   6      1
out_family$seed_used
#> [1] 123

if (FALSE) { # \dontrun{
## ---------------------------------------------------------
## Example 2: GRM-based grouping with dispersion optimization
## ---------------------------------------------------------
out_grm <- do.call(
  prep_famoptg,
  c(
    x$OptiDesign_famoptg_example,
    x$OptiDesign_famoptg_args_grm
  )
)

dim(out_grm$layout_matrix)
head(out_grm$field_book)
out_grm$efficiency

## ---------------------------------------------------------
## Example 3: Augmented repeated-check use pattern
## In practice, this is obtained by leaving p-rep arguments empty
## and supplying candidate entries only as unreplicated entries.
## ---------------------------------------------------------
aug_args <- x$OptiDesign_famoptg_example
aug_args$p_rep_treatments <- character(0)
aug_args$p_rep_reps <- integer(0)
aug_args$p_rep_families <- character(0)

out_aug <- do.call(
  prep_famoptg,
  c(aug_args, x$OptiDesign_famoptg_args_family)
)

dim(out_aug$layout_matrix)
head(out_aug$field_book)

## ---------------------------------------------------------
## Example 4: RCBD-type repeated-check use pattern
## Here all non-check entries are treated as replicated entries.
## If their common replication equals n_blocks, each entry appears
## once in every block.
## ---------------------------------------------------------
rcbd_args <- x$OptiDesign_famoptg_example
all_noncheck <- c(rcbd_args$p_rep_treatments, rcbd_args$unreplicated_treatments)
all_noncheck_fam <- c(rcbd_args$p_rep_families, rcbd_args$unreplicated_families)

rcbd_args$p_rep_treatments <- all_noncheck
rcbd_args$p_rep_reps <- rep(2L, length(all_noncheck))
rcbd_args$p_rep_families <- all_noncheck_fam
rcbd_args$unreplicated_treatments <- character(0)
rcbd_args$unreplicated_families <- character(0)

out_rcbd_type <- do.call(
  prep_famoptg,
  c(rcbd_args, x$OptiDesign_famoptg_args_family)
)

dim(out_rcbd_type$layout_matrix)
head(out_rcbd_type$field_book)
} # }
```
