# OptiDesign: Advanced Experimental Design for Breeding Trials

## Introduction

Experimental design is a fundamental component of plant breeding and
agronomic research. The ability to accurately estimate genetic effects,
compare treatments, and predict field performance depends critically on
how trials are constructed.

Traditional designs such as randomized complete block designs (RCBD) and
alpha-lattice designs typically assume regular block structures,
balanced replication, and independence between spatial and genetic
effects. Modern breeding programs face increasingly complex constraints
that these classical frameworks do not address well: large numbers of
genotypes, limited field capacity, heterogeneous environments, and the
availability of genomic and pedigree information.

`OptiDesign` addresses these challenges by integrating flexible field
layout construction, genetic structure (family, pedigree, genomic),
spatial-aware arrangement, and optional mixed-model efficiency
evaluation into a unified design framework.

------------------------------------------------------------------------

## Conceptual Framework

### The mixed model

Let $y$ be the vector of observed phenotypes, $g$ the genetic effects,
and $e$ the residuals. The general mixed model underlying design
evaluation is:

$$y = X\beta + Zg + e$$

where:

| Symbol  | Description                                                   |
|---------|---------------------------------------------------------------|
| $X$     | Design matrix for fixed effects                               |
| $\beta$ | Vector of fixed effects (e.g., block, check)                  |
| $Z$     | Incidence matrix linking genotypes to plots                   |
| $g$     | Genetic effects, $g \sim N\left( 0,\,\sigma_{g}^{2}K \right)$ |
| $K$     | Genetic relationship matrix (GRM or A)                        |
| $e$     | Residuals, $e \sim N(0,\, R)$                                 |

The design problem is to choose $X$ and $Z$ — through field layout
construction — to optimise the precision of $\widehat{\beta}$ (BLUEs)
and the prediction accuracy of $\widehat{g}$ (BLUPs).

### Field as a spatial process

Field plots are spatially structured. Residuals are rarely independent
owing to row and column effects, soil gradients, and management
variation. Three residual structures are supported:

- **IID** — independent and identically distributed residuals (classical
  assumption)
- **AR1** — first-order autoregressive structure along one direction
- **AR1×AR1** — separable two-dimensional autoregressive structure (rows
  and columns)

### Genetic structure in design

Entries can be grouped based on family labels, the pedigree relationship
matrix $A$, or the genomic relationship matrix (GRM). Incorporating
genetic structure at the design stage allows the layout to reduce
clustering of related genotypes, improve diversity spread across blocks,
and reduce confounding between spatial and genetic signals.

------------------------------------------------------------------------

## Design Strategies

### Block-based designs with repeated checks — `prep_famoptg()`

[`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)
constructs a block-based design in which checks are repeated in every
block, while non-check entries are allocated across blocks according to
their specified replication levels. Depending on the replication pattern
supplied, the same function generates:

- an **augmented design** — unreplicated test entries with checks in
  every block
- a **partially replicated (p-rep) design** — some entries replicated,
  others unreplicated
- an **RCBD-type design** — all non-check entries replicated equally,
  with `n_blocks` equal to the replication number

The key design rule is that **a treatment may appear at most once per
block**, regardless of its total replication. This ensures replicates
are distributed across blocks rather than clustered.

**Use
[`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)
when:**

- checks must appear in every block
- many entries must be screened with limited resources
- the design may be augmented, p-rep, or RCBD-type depending on
  replication
- early- to intermediate-stage breeding evaluation is the objective

### Stream-based alpha designs — `alpha_rc_stream()`

[`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md)
constructs alpha row–column designs on a fixed grid using a stream-based
layout. The field is converted into a single ordered stream of
positions, split into contiguous replicate segments, and further divided
into incomplete blocks.

A key property of this approach is that **unused cells are placed only
at the end of the field stream** — not scattered across the grid — which
is critical for practical field implementation.

**Use
[`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md)
when:**

- field dimensions are fixed and cannot change
- planting follows field-book (serpentine) order
- replicates are defined operationally as contiguous segments
- incomplete blocks may be unequal in size
- classical rectangular alpha designs are not feasible

------------------------------------------------------------------------

## Example Dataset

The package ships with a built-in example dataset containing inputs for
both main functions:

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

------------------------------------------------------------------------

## Example 1 — Family-based block design with repeated checks

This example calls
[`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)
with family-based grouping. Depending on the replication pattern in the
example inputs, the function acts as an augmented, p-rep, or RCBD-type
constructor:

``` r
out_fam <- do.call(
  OptiDesign::prep_famoptg,
  c(
    x$OptiDesign_famoptg_example,
    x$OptiDesign_famoptg_args_family
  )
)

dim(out_fam$layout_matrix)
#> [1] 16  8
head(out_fam$field_book)
#>        Treatment Family Gcluster Block Plot Row Column
#> plot_1      L076    F03     <NA>     1    1   1      1
#> plot_2      L061    F16     <NA>     1    2   2      1
#> plot_3      L028    F12     <NA>     1    3   3      1
#> plot_4      L030    F10     <NA>     1    4   4      1
#> plot_5      L027    F21     <NA>     1    5   5      1
#> plot_6      L051    F09     <NA>     1    6   6      1
```

The result contains:

- `layout_matrix` — a spatial grid of plot assignments
  (`n_rows × n_cols`)
- `field_book` — a row-by-row record of plot, block, treatment, and
  grouping information suitable for field implementation

In this design: checks are present in every block; non-check entries are
allocated across blocks according to their replication levels; no
non-check entry appears twice in the same block; and grouping is based
on family labels.

------------------------------------------------------------------------

## Example 2 — Stream-based alpha row–column design

``` r
out_alpha <- do.call(
  OptiDesign::alpha_rc_stream,
  c(
    x$OptiDesign_alpha_example,
    x$OptiDesign_alpha_args_family
  )
)

dim(out_alpha$layout_matrix)
#> [1] 12 14
head(out_alpha$field_book)
#>   Treatment Family Gcluster Check PlotStream Rep IBlock BlockInRep Row Column
#> 1      L140    F24     <NA> FALSE          1   1      1          1   1      1
#> 2      L103    F14     <NA>  TRUE          2   1      1          1   1      2
#> 3      L102    F17     <NA>  TRUE          3   1      1          1   1      3
#> 4      L117    F06     <NA> FALSE          4   1      1          1   1      4
#> 5      L101    F14     <NA>  TRUE          5   1      1          1   1      5
#> 6      L133    F22     <NA> FALSE          6   1      1          1   1      6
```

In this design: the field is traversed as a single ordered stream;
replicate segments are contiguous within that stream; checks are placed
in every incomplete block; and any unused cells appear only at the end
of the stream.

------------------------------------------------------------------------

## Example 3 — GRM-based stream design

When a genomic relationship matrix is available it can replace family
labels as the grouping source, allowing the layout to reflect genomic
similarity more precisely:

``` r
out_grm <- do.call(
  OptiDesign::alpha_rc_stream,
  c(
    x$OptiDesign_alpha_example,
    x$OptiDesign_alpha_args_grm
  )
)

head(out_grm$field_book)
#>   Treatment Family Gcluster Check PlotStream Rep IBlock BlockInRep Row Column
#> 1      L140    F24      G21 FALSE          1   1      1          1   1      1
#> 2      L103    F14     <NA>  TRUE          2   1      1          1   1      2
#> 3      L102    F17     <NA>  TRUE          3   1      1          1   1      3
#> 4      L117    F06      G22 FALSE          4   1      1          1   1      4
#> 5      L101    F14     <NA>  TRUE          5   1      1          1   1      5
#> 6      L133    F22      G19 FALSE          6   1      1          1   1      6
```

------------------------------------------------------------------------

## Grouping Options

Both functions support three grouping strategies:

| Strategy       | Source                    | Best for                            |
|----------------|---------------------------|-------------------------------------|
| Family-based   | User-defined labels       | Simple, interpretable grouping      |
| GRM-based      | Genomic similarity matrix | Captures real genomic relationships |
| Pedigree-based | $A$ matrix                | When genomic data is unavailable    |

------------------------------------------------------------------------

## Dispersion Optimization

An optional step that improves the spatial distribution of genetically
related entries. If $d(i,j)$ denotes the similarity between plots $i$
and $j$, the optimizer seeks to reduce:

$$\sum\limits_{{(i,j)} \in \mathcal{N}}d(i,j)$$

where $\mathcal{N}$ is the set of neighbouring plot pairs. This reduces
local clustering of related materials, improves robustness to spatial
confounding, and promotes better distribution of genetic diversity
across the field.

Use dispersion optimization when genomic or pedigree structure is
important and strong local spatial correlation is expected.

------------------------------------------------------------------------

## Efficiency Evaluation

An optional diagnostic step that evaluates a candidate design under a
mixed-model framework before field implementation.

**Fixed-effect precision** — the variance of BLUE estimates:

$${Var}\left( \widehat{\beta} \right) = \left( X^{\top}V^{- 1}X \right)^{- 1}$$

Lower diagonal values indicate higher precision for the corresponding
fixed effects.

**Random-effect prediction** — the prediction error variance (PEV):

$${PEV} = {Var}\left( \widehat{g} - g \right)$$

Lower PEV values indicate higher prediction accuracy for genetic
effects.

Supported models include IID, AR1, and AR1×AR1 spatial structures, and
prediction modes include BLUP, GBLUP, and PBLUP. Use efficiency
evaluation when comparing alternative designs or studying the impact of
block structure and replication on estimability.

------------------------------------------------------------------------

## Practical Guidelines

| Criterion             | [`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md) | [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md) |
|-----------------------|-------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------|
| Checks in every block | Yes — by design                                                                     | Yes — by design                                                                           |
| Field dimensions      | Flexible                                                                            | Fixed grid required                                                                       |
| Replicate structure   | Block-based                                                                         | Contiguous stream segments                                                                |
| Unequal block sizes   | Not applicable                                                                      | Supported                                                                                 |
| Unused cell placement | Depends on layout                                                                   | End of stream only                                                                        |
| Typical use           | Early/mid breeding stage                                                            | Operational field-book trials                                                             |

**Grouping choice:** use family labels for simple interpretability, GRM
when genomic data is available and relatedness precision matters, and
the $A$ matrix when only pedigree information is available.

------------------------------------------------------------------------

## Summary

`OptiDesign` provides two complementary design strategies — block-based
repeated-check layouts via
[`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)
and stream-based alpha row–column layouts via
[`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md)
— together with optional dispersion optimization and mixed-model
efficiency evaluation. Both functions support family, genomic, and
pedigree-based grouping and are designed to be practically implementable
in field trials with real-world constraints.

------------------------------------------------------------------------

## Session Information

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
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
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] Matrix_1.7-4      xfun_0.57         lattice_0.22-9    OptiDesign_0.1.0 
#>  [9] cachem_1.1.0      knitr_1.51        htmltools_0.5.9   rmarkdown_2.30   
#> [13] lifecycle_1.0.5   cli_3.6.5         grid_4.5.3        sass_0.4.10      
#> [17] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [21] compiler_4.5.3    tools_4.5.3       ragg_1.5.2        pracma_2.4.6     
#> [25] evaluate_1.0.5    bslib_0.10.0      yaml_2.3.12       jsonlite_2.0.0   
#> [29] rlang_1.1.7       fs_2.0.0
```
