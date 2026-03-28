# OptiDesign 0.1.0

## Initial release

### New design constructors

- Added `prep_famoptg()`: constructs repeated-check block designs with
  flexible replication, supporting augmented, partially replicated (p-rep),
  and RCBD-type repeated-check layouts. The p-rep constraint — no replicated
  treatment appears twice in the same block — is enforced by construction.

- Added `alpha_rc_stream()`: constructs fixed-grid alpha row-column stream
  designs with repeated checks. The field is treated as a global traversal
  stream partitioned into contiguous replicate segments and incomplete blocks.
  Block-size constraints are expressed as total block size (checks + entries)
  via `min_block_size` and `max_block_size`.

### New evaluation functions

- Added `evaluate_famoptg_efficiency()`: standalone mixed-model efficiency
  evaluator for designs produced by `prep_famoptg()`. Model contains Block +
  Row + Column random effects (`sigma_b2`). Computes A-criterion, D-criterion,
  and CDmean. Fully decoupled from construction.

- Added `evaluate_alpha_efficiency()`: standalone mixed-model efficiency
  evaluator for designs produced by `alpha_rc_stream()`. Model contains
  Rep + IBlock(Rep) + Row + Column random effects (`sigma_rep2`, `sigma_ib2`).
  Computes A-criterion, D-criterion, and CDmean. Fully decoupled from
  construction.

### New optimisation functions

- Added `optimize_famoptg()`: Random Restart (RS) optimiser for
  `prep_famoptg()` designs. Searches across `n_restarts` independent
  randomisations and returns the design with the best criterion value.
  Supports A, D, both, and CDmean criteria. RS is used exclusively because
  the p-rep constraint is preserved by construction at every call.

- Added `optimize_alpha_rc()`: criterion-driven optimiser for
  `alpha_rc_stream()` designs with three search strategies:
  - **RS** (Random Restart): independent randomisations, return best.
  - **SA** (Simulated Annealing): iterative entry-permutation swaps with
    temperature-governed acceptance, repeated across multiple restarts.
  - **GA** (Genetic Algorithm): population of entry permutations evolved
    via Order Crossover (OX1), swap mutation, and elitism.
  Supports A, D, both, and CDmean criteria. All methods preserve structural
  constraints by construction.

### Optimality criteria

Both evaluation and optimisation functions support:

- **A-criterion**: mean pairwise contrast variance (fixed effects) or mean
  prediction error variance / PEV (random effects). Lower is better.
- **D-criterion**: geometric mean of contrast covariance eigenvalues (fixed
  effects only). Lower is better.
- **CDmean**: mean coefficient of determination for genomic breeding value
  prediction, `CDmean = 1 - mean_PEV / sigma_g2` (Rincent et al. 2012).
  Higher is better. Requires random treatment effects and a genomic prediction
  model (`prediction_type = "GBLUP"` or `"PBLUP"`).

### Integrity checking

Both optimisers implement a four-point integrity checking strategy:
post-construction validation, re-check before storing as best, final check
before return, and emergency fallback if no valid design is found across all
iterations.

### Internal architecture

- Added `alpha_rc_helpers.R`: internal helper functions shared across all six
  exported functions — sparse incidence matrix construction, AR1 precision
  matrices, the mixed model solver, the Hutchinson stochastic trace estimator
  (Hutchinson 1990) for scalable criterion approximation in large designs,
  Chebyshev neighbourhood enumeration, and the dispersion scoring function.

- Construction, evaluation, and optimisation follow a
  **single-responsibility architecture**: each step is a separate function
  that can be called independently or chained.

### Supported residual structures

All evaluation and optimisation functions support:

- `"IID"`: independent residuals.
- `"AR1"`: row-only AR1 autocorrelation.
- `"AR1xAR1"`: separable row × column AR1 autocorrelation.

### Grouping and dispersion

Both design families support:

- Family-based adjacency scoring (`cluster_source = "Family"`).
- GRM-based clustering (`cluster_source = "GRM"`).
- Pedigree-based clustering (`cluster_source = "A"`).
- Optional post-layout genetic dispersion optimisation (`use_dispersion = TRUE`)
  using swap-based local search to reduce spatial clustering of related entries.

### Package infrastructure

- Added package documentation (`OptiDesign-package.R`) covering both design
  families, typical workflows, function overview table, and key references.
- Added `OptiDesign_example_data`: synthetic dataset with treatment vectors,
  field dimensions, relationship matrices (GRM, A, K), and ready-to-use
  argument lists for both design families.
- Added pkgdown site configuration (`_pkgdown.yml`) with reference sections
  organised by design family.
- Added GitHub Actions workflows for R CMD CHECK and pkgdown site deployment.
- Added unit tests (`testthat >= 3.0.0`).