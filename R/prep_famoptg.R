#' Create a repeated-check block design with flexible replication
#'
#' Construct an augmented, partially replicated, or RCBD-type repeated-check
#' field design with optional treatment grouping, optional post-layout
#' dispersion optimization, and optional mixed-model efficiency diagnostics.
#'
#' @description
#' `prep_famoptg()` is a general repeated-check block design constructor in
#' which **check treatments are included in every block**, while non-check
#' treatments are allocated across blocks according to user-specified
#' replication levels.
#'
#' Depending on the supplied treatment structure, the same function can generate:
#'
#' - an **augmented design**, when only checks are repeated and all other entries
#'   are unreplicated;
#' - a **partially replicated (p-rep) design**, when some non-check entries are
#'   replicated and others are unreplicated;
#' - an **RCBD-type repeated-check design**, when all non-check entries are given
#'   the same replication number greater than 1, especially when that replication
#'   equals `n_blocks`.
#'
#' The basic allocation rules are:
#'
#' - **Checks** are included in **every block**.
#' - **Replicated non-check treatments** are assigned to a subset of blocks and
#'   may appear more than once overall, but **at most once within any single block**.
#' - **Unreplicated treatments** appear exactly once in the full design.
#'
#' This means, for example, that if a treatment has replication 2 and
#' `n_blocks >= 2`, its two replicates are placed in **two different blocks**
#' rather than repeated within the same block.
#'
#' The function returns:
#'
#' - a **layout matrix** representing the field grid, and
#' - a **field book** giving one row per assigned plot with treatment identity,
#'   grouping metadata, block membership, and field coordinates.
#'
#' Optionally, the function can also:
#'
#' - derive treatment groups from a genomic relationship matrix (`GRM`) or a
#'   pedigree relationship matrix (`A`);
#' - optimize the spatial dispersion of checks;
#' - apply a swap-based local search to reduce local relatedness among nearby
#'   non-check plots;
#' - compute mixed-model design efficiency metrics for fixed or random treatment
#'   effects.
#'
#' @section Design classes represented by the same function:
#'
#' **1. Augmented repeated-check design**
#'
#' This is obtained when:
#'
#' - `p_rep_treatments = NULL` or `character(0)`,
#' - `p_rep_reps = NULL` or `integer(0)`,
#' - `p_rep_families = NULL` or `character(0)`,
#' - non-check entries are supplied only through `unreplicated_treatments`.
#'
#' In this case:
#'
#' - checks are repeated in every block;
#' - all test entries appear once;
#' - the design behaves as an augmented repeated-check layout.
#'
#' This setting is especially useful for early-stage screening when many entries
#' must be observed but only checks can be repeated systematically.
#'
#' **2. Partially replicated repeated-check design**
#'
#' This is obtained when:
#'
#' - some non-check treatments are supplied in `p_rep_treatments`;
#' - those treatments have replication counts greater than 1;
#' - other entries may remain unreplicated.
#'
#' This is the classical use case of the function: a mixture of repeated checks,
#' replicated candidate entries, and single-plot candidate entries.
#'
#' **3. RCBD-type repeated-check design**
#'
#' This is obtained when:
#'
#' - all non-check treatments are supplied through `p_rep_treatments`;
#' - they all receive the same replication number greater than 1.
#'
#' If that common replication equals `n_blocks`, then every non-check treatment
#' appears once in every block, which is the closest repeated-check analogue of
#' a classical RCBD under this framework.
#'
#' If the common replication is less than `n_blocks`, the design remains balanced
#' and repeated, but it is not a strict classical RCBD because not every
#' treatment appears in every block.
#'
#' @section Conceptual workflow:
#' The function proceeds in several stages:
#'
#' 1. **Validate inputs** and reconcile field size with the required number of plots.
#' 2. **Prepare grouping information** from family labels or from clustering on
#'    `GRM` / `A`.
#' 3. **Assign replicated non-check entries to blocks** subject to the rule that
#'    no treatment appears twice in the same block.
#' 4. **Distribute unreplicated entries** across blocks.
#' 5. **Insert checks into each block** using the chosen placement strategy.
#' 6. **Shuffle treatment order within blocks** to reduce adjacency of same-group
#'    entries in the 1D block ordering.
#' 7. **Map the ordered treatments to the field grid** according to `order` and
#'    `serpentine`.
#' 8. Optionally apply **genetic dispersion optimization**.
#' 9. Optionally compute **design efficiency metrics**.
#'
#' @section What this function is most useful for:
#' This function is especially useful for:
#'
#' - early- to intermediate-stage breeding trials where many entries must be
#'   screened and not all can be equally replicated;
#' - augmented repeated-check designs with a large number of single-plot test entries;
#' - p-rep layouts mixing replicated and unreplicated entries;
#' - balanced repeated-check block designs in which all candidate entries are replicated;
#' - experiments where checks must appear in every block;
#' - layouts where related lines should not be clustered too closely;
#' - applications where family, pedigree, or genomic relationships should influence
#'   physical layout construction;
#' - cases where the user wants to compare alternative design choices using
#'   model-based efficiency metrics.
#'
#' @section Grouping logic and why it matters:
#' Several parts of the function rely on the idea of a **group**:
#'
#' - During block construction, groups are used to reduce adjacency of similar
#'   entries in the 1D ordering within each block.
#' - During optional dispersion optimization, relationship matrices are used to
#'   discourage close spatial placement of related non-check entries.
#'
#' Group labels may come from:
#'
#' - user-supplied families (`cluster_source = "Family"`);
#' - clusters derived from a genomic relationship matrix (`cluster_source = "GRM"`);
#' - clusters derived from a pedigree relationship matrix (`cluster_source = "A"`).
#'
#' Family-based grouping is the simplest option and is usually appropriate when
#' family labels are meaningful and easy to interpret.
#'
#' Matrix-based grouping is more useful when:
#'
#' - family labels are too coarse;
#' - genomic or pedigree structure is more informative than nominal family labels;
#' - the objective is to spread genetically similar materials more effectively.
#'
#' @section Dependency guide:
#' Many arguments are active only in certain modes.
#'
#' **Treatment structure**
#'
#' - If `p_rep_treatments` is empty (or `NULL`) and
#'   `unreplicated_treatments` is non-empty, the function behaves as an
#'   **augmented repeated-check design** constructor.
#'
#' - If both `p_rep_treatments` and `unreplicated_treatments` are present,
#'   the function behaves as a **p-rep repeated-check design** constructor.
#'
#' - If all non-check treatments are supplied in `p_rep_treatments` with a common
#'   replication number greater than 1, the function behaves as an
#'   **RCBD-type repeated-check design** constructor.
#'
#' **Grouping mode**
#'
#' - If `cluster_source = "Family"`:
#'   - `check_families`, `p_rep_families`, and `unreplicated_families` are used.
#'   - `GRM`, `A`, `id_map`, `cluster_method`, `cluster_seed`,
#'     `cluster_attempts`, and `n_pcs_use` are ignored.
#'
#' - If `cluster_source = "GRM"`:
#'   - `GRM` is required.
#'   - `id_map` is needed only if treatment IDs do not match `rownames(GRM)`.
#'   - `cluster_method`, `cluster_seed`, `cluster_attempts`, and `n_pcs_use`
#'     are used.
#'   - `A` is ignored.
#'
#' - If `cluster_source = "A"`:
#'   - `A` is required.
#'   - `id_map` is needed only if treatment IDs do not match `rownames(A)`.
#'   - `cluster_method`, `cluster_seed`, `cluster_attempts`, and `n_pcs_use`
#'     are used.
#'   - `GRM` is ignored.
#'
#' **Efficiency evaluation**
#'
#' - If `eval_efficiency = FALSE`, all efficiency-related arguments are ignored.
#'
#' - If `eval_efficiency = TRUE` and `treatment_effect = "fixed"`:
#'   - fixed-effect design efficiency is computed;
#'   - `prediction_type`, `K`, and `line_id_map` are ignored.
#'
#' - If `eval_efficiency = TRUE` and `treatment_effect = "random"`:
#'   - `prediction_type` becomes active.
#'
#' - If `prediction_type = "none"`:
#'   - no random-effect prediction efficiency is computed.
#'
#' - If `prediction_type = "IID"`:
#'   - non-check treatments are treated as independent random effects;
#'   - `K` and `line_id_map` are ignored.
#'
#' - If `prediction_type %in% c("GBLUP", "PBLUP")`:
#'   - `K` is required;
#'   - `line_id_map` may be needed if treatment IDs do not match `rownames(K)`.
#'
#' **Residual structure**
#'
#' - If `residual_structure = "IID"`:
#'   - residuals are assumed independent;
#'   - `rho_row` and `rho_col` are ignored.
#'
#' - If `residual_structure = "AR1"`:
#'   - row-wise correlation is used;
#'   - `rho_row` is active;
#'   - `rho_col` is ignored.
#'
#' - If `residual_structure = "AR1xAR1"`:
#'   - row and column correlations are both used;
#'   - `rho_row` and `rho_col` are both active.
#'
#' **Check placement**
#'
#' - If `check_placement = "optimal"`, `check_opt_attempts` is used.
#' - Otherwise, `check_opt_attempts` is ignored.
#'
#' **Dispersion optimization**
#'
#' - If `use_dispersion = FALSE`, all dispersion arguments are ignored.
#'
#' - If `use_dispersion = TRUE` and `dispersion_source = "K"`:
#'   - `K` is required;
#'   - `line_id_map` may be needed if treatment IDs differ from `rownames(K)`.
#'
#' - If `use_dispersion = TRUE` and `dispersion_source = "A"`:
#'   - `A` is required.
#'
#' - If `use_dispersion = TRUE` and `dispersion_source = "GRM"`:
#'   - `GRM` is required.
#'
#' @details
#' The function uses `mod()` imported from the **pracma** package for serpentine
#' traversal logic.
#'
#' Field dimensions do not have to match the exact number of required plots if
#' `warn_and_correct = TRUE`. In that case, the function expands one dimension
#' so the field can contain all required plots. Extra cells remain `NA` in the
#' returned `layout_matrix`, but the `field_book` includes only assigned plots.
#'
#' The function separates:
#'
#' - **construction logic** (how treatments are placed);
#' - **grouping logic** (how similarity is defined);
#' - **dispersion logic** (how nearby related entries are discouraged);
#' - **efficiency logic** (how the final design is evaluated).
#'
#' This separation is important because a user may, for example:
#'
#' - use family labels for adjacency control;
#' - use a genomic matrix for dispersion optimization;
#' - and choose not to compute efficiency;
#'
#' or alternatively:
#'
#' - use GRM-based clustering;
#' - skip dispersion optimization;
#' - and compute GBLUP-based efficiency using a separate matrix `K`.
#'
#' @param check_treatments Character vector of check treatment identifiers.
#'
#' These treatments are included in **every block**, so they define the repeated
#' reference structure of the design.
#'
#' Use this argument when you have standard checks, control varieties, benchmark
#' cultivars, or fixed reference entries that must be present everywhere to
#' improve block comparability.
#'
#' This argument is always required.
#'
#' Example use case:
#' a breeding trial with 3 control cultivars included in all blocks.
#'
#' @param check_families Character vector of the same length as `check_treatments`.
#'
#' Provides the family or group label for each check treatment.
#'
#' Even when `cluster_source != "Family"`, checks still need family labels
#' because checks are always assigned a group label in the output and remain part
#' of the grouping metadata used internally.
#'
#' This argument always depends on `check_treatments` and must align exactly in
#' order and length.
#'
#' Example:
#' if the three checks all belong to the same control class,
#' `check_families` might be `c("CHECK", "CHECK", "CHECK")`.
#'
#' @param p_rep_treatments Character vector of replicated non-check treatment IDs.
#'
#' These are the entries that appear more than once overall, but at most once in
#' any given block.
#'
#' Use this argument when some or all non-check treatments should be replicated.
#'
#' Depending on how it is used:
#'
#' - if it is empty and only `unreplicated_treatments` are supplied, the function
#'   behaves as an augmented repeated-check constructor;
#' - if it contains only some non-check entries, the function behaves as a p-rep
#'   repeated-check constructor;
#' - if it contains all non-check entries with common replication greater than 1,
#'   the function behaves as an RCBD-type repeated-check constructor.
#'
#' Can be `NULL` or `character(0)` if no replicated non-check treatments are present.
#'
#' Depends on:
#' - `p_rep_reps`
#' - `p_rep_families`
#'
#' Example use case:
#' elite lines replicated twice while many other lines remain unreplicated.
#'
#' @param p_rep_reps Integer vector giving the total number of replicates for each
#' entry in `p_rep_treatments`.
#'
#' This controls how many times each replicated non-check treatment appears in
#' the full design.
#'
#' Must:
#' - have the same length as `p_rep_treatments`;
#' - align element-wise with `p_rep_treatments`;
#' - satisfy `p_rep_reps[i] <= n_blocks`, because a treatment cannot occur twice
#'   in the same block.
#'
#' A treatment with replication 2 and `n_blocks >= 2` will be placed in two
#' different blocks.
#'
#' Use a common vector such as `rep(2, length(p_rep_treatments))` when all
#' replicated treatments should appear twice.
#'
#' Use unequal values when some entries deserve more replication than others.
#'
#' @param p_rep_families Character vector of the same length as `p_rep_treatments`.
#'
#' Family or group labels for replicated non-check treatments.
#'
#' These labels are required whenever replicated non-check treatments exist.
#' Even when `cluster_source != "Family"`, they remain useful because the
#' function uses non-check family structure to determine the target number of
#' clusters in matrix-based grouping modes.
#'
#' Depends on:
#' - `p_rep_treatments`
#'
#' Example use case:
#' replicated candidate lines belong to multiple families and the user wants
#' adjacency control to reflect that structure.
#'
#' @param unreplicated_treatments Character vector of treatments that appear exactly once.
#'
#' These are single-plot entries used when broad screening is needed and full
#' replication is not possible.
#'
#' Use `NULL` or `character(0)` if no unreplicated treatments are present.
#'
#' In combination with repeated checks and no replicated non-checks, these
#' entries define an augmented repeated-check design.
#'
#' Depends on:
#' - `unreplicated_families`
#'
#' Example use case:
#' a large early-generation nursery where each new line is observed once.
#'
#' @param unreplicated_families Character vector of the same length as
#' `unreplicated_treatments`.
#'
#' Family or group labels for unreplicated treatments.
#'
#' Required whenever `unreplicated_treatments` is provided.
#'
#' These labels are especially useful when:
#' - `cluster_source = "Family"`;
#' - matrix-based grouping is used and the target number of clusters should
#'   reflect the family structure among non-checks.
#'
#' @param n_blocks Integer giving the number of experimental blocks.
#'
#' This determines how many times checks are repeated, and it also constrains the
#' maximum allowable replication count for any replicated non-check treatment.
#'
#' Use more blocks when:
#' - many checks must be accommodated;
#' - finer local control is desired;
#' - field heterogeneity suggests stronger blocking.
#'
#' This argument affects:
#' - total plot count;
#' - feasibility of replicated-treatment assignment;
#' - distribution of unreplicated treatments.
#'
#' When the same replicated treatment should appear in every block, its
#' replication must equal `n_blocks`.
#'
#' @param n_rows Integer giving the number of rows in the final field grid.
#'
#' Controls the physical field geometry returned in `layout_matrix`.
#'
#' Use this to match the intended number of field rows.
#'
#' If `warn_and_correct = TRUE`, this value may be retained while `n_cols` is
#' adjusted when `fix_rows = TRUE`.
#'
#' @param n_cols Integer giving the number of columns in the final field grid.
#'
#' Controls the physical field geometry returned in `layout_matrix`.
#'
#' Use this together with `n_rows` to represent the intended field dimensions.
#'
#' If `warn_and_correct = TRUE`, this value may be retained while `n_rows` is
#' adjusted when `fix_rows = FALSE`.
#'
#' @param order Character specifying how the field grid is filled.
#'
#' Allowed values:
#' - `"row"`: fill row by row;
#' - `"column"`: fill column by column.
#'
#' Use `"row"` when field-book order or planting order follows rows.
#' Use `"column"` when field-book order follows columns.
#'
#' This argument interacts directly with `serpentine`.
#'
#' @param serpentine Logical indicating whether alternate rows or columns should
#' reverse direction during grid filling.
#'
#' If `TRUE`:
#' - with `order = "row"`, alternate rows are reversed;
#' - with `order = "column"`, alternate columns are reversed.
#'
#' Use this when planting, labeling, or harvesting follows a serpentine movement pattern.
#'
#' This affects only how the final ordered treatments are mapped into space; it
#' does not change block composition.
#'
#' Depends on:
#' - `order`
#'
#' @param seed Optional integer seed controlling stochastic components of the design.
#'
#' Use this when reproducibility is important.
#'
#' If `NULL`, a seed is generated internally and returned as `seed_used`.
#'
#' This seed affects:
#' - replicated-treatment block assignment;
#' - shuffling within blocks;
#' - candidate generation for optimal check placement;
#' - optional dispersion optimization when `dispersion_seed` is `NULL`.
#'
#' @param attempts Integer giving the maximum number of shuffle attempts per block.
#'
#' Used to reduce adjacency of same-group entries in the 1D block ordering.
#'
#' Larger values may improve adjacency avoidance but increase runtime.
#'
#' Use larger values when:
#' - family structure is highly imbalanced;
#' - many similar entries occur in the same block;
#' - stricter adjacency avoidance is desired.
#'
#' @param warn_and_correct Logical controlling what happens when the supplied field
#' size does not match the required number of plots.
#'
#' If `FALSE`, the function stops with an error.
#'
#' If `TRUE`, the function adjusts one dimension upward so the field can contain
#' all required plots.
#'
#' Use `FALSE` when the field dimensions are fixed and must not be changed.
#' Use `TRUE` when a near-feasible layout should be repaired automatically.
#'
#' Depends on:
#' - `fix_rows`
#'
#' @param fix_rows Logical used only when `warn_and_correct = TRUE`.
#'
#' If `TRUE`, `n_rows` is kept fixed and `n_cols` is adjusted.
#' If `FALSE`, `n_cols` is kept fixed and `n_rows` is adjusted.
#'
#' Use this according to which field dimension is physically fixed in practice.
#'
#' Example:
#' if the field has a fixed number of rows but flexible row length,
#' use `fix_rows = TRUE`.
#'
#' @param cluster_source Character specifying the source of grouping used for
#' non-check adjacency control.
#'
#' Allowed values:
#' - `"Family"`
#' - `"GRM"`
#' - `"A"`
#'
#' Use `"Family"` when family labels are trusted and easy to interpret.
#'
#' Use `"GRM"` when genomic relatedness should define grouping more accurately.
#'
#' Use `"A"` when pedigree structure is available and pedigree-based grouping is desired.
#'
#' This argument determines which of the following become active:
#' - `GRM`
#' - `A`
#' - `id_map`
#' - `cluster_method`
#' - `cluster_seed`
#' - `cluster_attempts`
#' - `n_pcs_use`
#'
#' @param GRM Optional square genomic relationship matrix with rownames and colnames.
#'
#' Required when:
#' - `cluster_source = "GRM"`, or
#' - `use_dispersion = TRUE` and `dispersion_source = "GRM"`.
#'
#' Ignored otherwise.
#'
#' Use this when genomic similarity should drive grouping or dispersion scoring.
#'
#' Row and column names must correspond to treatment IDs or to IDs reachable via `id_map`.
#'
#' @param A Optional square pedigree relationship matrix with rownames and colnames.
#'
#' Required when:
#' - `cluster_source = "A"`, or
#' - `use_dispersion = TRUE` and `dispersion_source = "A"`.
#'
#' Ignored otherwise.
#'
#' Use this when pedigree structure should drive grouping or dispersion scoring.
#'
#' @param id_map Optional data frame with columns `Treatment` and `LineID`.
#'
#' Used only when `cluster_source %in% c("GRM", "A")` and treatment labels in
#' the design do not exactly match the rownames of the selected relationship matrix.
#'
#' Use this when:
#' - treatment names in the design are breeder-friendly labels;
#' - but the matrix uses internal IDs, numeric codes, or standardized line names.
#'
#' Ignored when `cluster_source = "Family"`.
#'
#' @param cluster_method Clustering method used after PCA in matrix-based grouping.
#'
#' Allowed values:
#' - `"kmeans"`
#' - `"hclust"`
#'
#' Use `"kmeans"` when compact clusters are acceptable and reproducible
#' initialization is desired.
#'
#' Use `"hclust"` when hierarchical grouping is preferred.
#'
#' Active only when `cluster_source %in% c("GRM", "A")`.
#'
#' @param cluster_seed Integer seed used when `cluster_method = "kmeans"`.
#'
#' Ignored otherwise.
#'
#' Use this when k-means clustering must be reproducible.
#'
#' @param cluster_attempts Integer number of random starts passed to `kmeans()`.
#'
#' Active only when:
#' - `cluster_source %in% c("GRM", "A")`, and
#' - `cluster_method = "kmeans"`.
#'
#' Larger values may improve clustering stability at higher runtime cost.
#'
#' @param n_pcs_use Integer or `Inf` giving the number of principal components used
#' for PCA-based clustering in matrix-based grouping.
#'
#' Use a smaller value when:
#' - only the leading structure should define grouping;
#' - smaller PCs are considered mostly noise.
#'
#' Use `Inf` to let the function use as many informative PCs as available.
#'
#' Ignored when `cluster_source = "Family"`.
#'
#' @param eval_efficiency Logical; if `TRUE`, compute design efficiency metrics.
#'
#' Use `TRUE` when the design should be evaluated under a mixed-model framework.
#'
#' Use `FALSE` when only the layout itself is needed.
#'
#' This argument activates:
#' - `treatment_effect`
#' - `prediction_type`
#' - `K`
#' - `line_id_map`
#' - `varcomp`
#' - `check_as_fixed`
#' - `residual_structure`
#' - `rho_row`
#' - `rho_col`
#' - `spatial_engine`
#' - `dense_max_n`
#' - `eff_trace_samples`
#' - `eff_full_max`
#'
#' @param treatment_effect Character indicating whether non-check treatments are
#' treated as `"fixed"` or `"random"` for efficiency evaluation.
#'
#' Active only when `eval_efficiency = TRUE`.
#'
#' Use `"fixed"` when interest is in the precision of estimated treatment contrasts.
#'
#' Use `"random"` when interest is in prediction quality of treatment effects.
#'
#' @param prediction_type Character controlling the random-effect prediction model.
#'
#' Allowed values:
#' - `"none"`
#' - `"IID"`
#' - `"GBLUP"`
#' - `"PBLUP"`
#'
#' Active only when:
#' - `eval_efficiency = TRUE`, and
#' - `treatment_effect = "random"`.
#'
#' Use `"none"` to skip random-effect efficiency.
#'
#' Use `"IID"` when all non-check effects are assumed independent.
#'
#' Use `"GBLUP"` or `"PBLUP"` when prediction should use a relationship matrix `K`.
#'
#' @param K Optional square relationship matrix.
#'
#' Used in two contexts:
#'
#' 1. for efficiency evaluation when random effects are predicted using
#'    `"GBLUP"` or `"PBLUP"`;
#' 2. for dispersion optimization when `dispersion_source = "K"`.
#'
#' Ignored otherwise.
#'
#' Use this when prediction or dispersion should be based on an explicit
#' relationship matrix.
#'
#' @param line_id_map Optional data frame with columns `Treatment` and `LineID`.
#'
#' Used only when `K` is active and treatment names differ from `rownames(K)`.
#'
#' This is analogous to `id_map`, but specifically for the matrix `K`.
#'
#' Use this when the design uses one set of treatment labels while `K` uses another.
#'
#' @param varcomp Named list of variance components used for efficiency evaluation.
#'
#' Must contain:
#' - `sigma_e2`
#' - `sigma_g2`
#' - `sigma_b2`
#' - `sigma_r2`
#' - `sigma_c2`
#'
#' Active only when `eval_efficiency = TRUE`.
#'
#' Use realistic values when the efficiency metric should reflect a plausible
#' data-generating model.
#'
#' Example use case:
#' a user may set larger `sigma_e2` when expecting noisy phenotypes, or larger
#' `sigma_b2` when block-to-block heterogeneity is substantial.
#'
#' @param check_as_fixed Logical indicating whether checks are included as fixed
#' indicator columns during efficiency evaluation.
#'
#' Active only when `eval_efficiency = TRUE`.
#'
#' Use `TRUE` when checks should be explicitly modeled as fixed benchmark effects.
#'
#' @param residual_structure Character specifying the residual correlation model.
#'
#' Allowed values:
#' - `"IID"`
#' - `"AR1"`
#' - `"AR1xAR1"`
#'
#' Use `"IID"` when no spatial residual correlation is assumed.
#'
#' Use `"AR1"` when correlation is expected mainly along rows.
#'
#' Use `"AR1xAR1"` when correlation is expected along both rows and columns.
#'
#' Active only when `eval_efficiency = TRUE`.
#'
#' @param rho_row Numeric AR1 row correlation parameter.
#'
#' Active only when:
#' - `eval_efficiency = TRUE`, and
#' - `residual_structure %in% c("AR1", "AR1xAR1")`.
#'
#' Use values near 0 for weak row correlation and values closer to 1 in magnitude
#' for stronger row correlation.
#'
#' @param rho_col Numeric AR1 column correlation parameter.
#'
#' Active only when:
#' - `eval_efficiency = TRUE`, and
#' - `residual_structure = "AR1xAR1"`.
#'
#' Use this when column-wise spatial correlation is part of the assumed residual model.
#'
#' @param spatial_engine Character controlling the computational implementation of
#' spatial residual calculations.
#'
#' Allowed values:
#' - `"auto"`
#' - `"dense"`
#' - `"sparse"`
#'
#' Active only when:
#' - `eval_efficiency = TRUE`, and
#' - `residual_structure != "IID"`.
#'
#' Use `"auto"` in most cases.
#'
#' Use `"dense"` for smaller problems when dense linear algebra is acceptable.
#'
#' Use `"sparse"` for larger problems where sparse precision matrices are preferred.
#'
#' @param dense_max_n Integer threshold used when `spatial_engine = "auto"`.
#'
#' If the number of observed plots is at most this threshold, dense computation
#' is used; otherwise sparse computation is used.
#'
#' @param eff_trace_samples Integer giving the number of Hutchinson trace samples
#' used when approximate efficiency evaluation is needed.
#'
#' Active only when:
#' - `eval_efficiency = TRUE`, and
#' - the target dimension exceeds `eff_full_max`.
#'
#' Larger values improve stability of the approximation but increase runtime.
#'
#' @param eff_full_max Integer giving the maximum target dimension for exact inverse
#' sub-block extraction in efficiency evaluation.
#'
#' Active only when `eval_efficiency = TRUE`.
#'
#' Use larger values when exact evaluation is desired and computation is affordable.
#' Use smaller values to switch to approximation earlier and reduce runtime.
#'
#' @param check_placement Character specifying how checks are placed within blocks.
#'
#' Allowed values:
#' - `"random"`
#' - `"systematic"`
#' - `"optimal"`
#'
#' Use `"random"` when no explicit spacing pattern is required.
#'
#' Use `"systematic"` when checks should be spread approximately evenly within each block.
#'
#' Use `"optimal"` when spatial dispersion of checks in the final 2D grid is a priority.
#'
#' If `"optimal"` is selected, `check_opt_attempts` becomes active.
#'
#' @param check_opt_attempts Integer number of candidate layouts evaluated when
#' `check_placement = "optimal"`.
#'
#' Larger values may improve the final spatial spread of checks but increase runtime.
#'
#' Ignored unless `check_placement = "optimal"`.
#'
#' @param use_dispersion Logical; if `TRUE`, apply post-layout swap optimization to
#' reduce local relatedness among nearby non-check treatments.
#'
#' Use this when genetically or pedigree-similar entries should be less clustered
#' in the final field.
#'
#' If `FALSE`, all dispersion-related arguments are ignored.
#'
#' @param dispersion_source Character specifying which matrix is used during
#' dispersion optimization.
#'
#' Allowed values:
#' - `"K"`
#' - `"A"`
#' - `"GRM"`
#'
#' Active only when `use_dispersion = TRUE`.
#'
#' Use `"K"` when a dedicated relationship matrix for dispersion is available.
#' Use `"A"` for pedigree-based spacing.
#' Use `"GRM"` for genomic-based spacing.
#'
#' @param dispersion_radius Integer neighborhood radius used in dispersion scoring.
#'
#' Neighboring plots are defined using Chebyshev distance, meaning plots are neighbors
#' when `max(|dr|, |dc|) <= dispersion_radius`.
#'
#' Use `1` for immediate neighbors only.
#' Use larger values when broader local neighborhoods matter.
#'
#' Active only when `use_dispersion = TRUE`.
#'
#' @param dispersion_iters Integer number of swap proposals used during local search.
#'
#' Larger values allow a more extensive search for reduced local relatedness but
#' increase runtime.
#'
#' Active only when `use_dispersion = TRUE`.
#'
#' @param dispersion_seed Optional integer seed used for dispersion optimization.
#'
#' Active only when `use_dispersion = TRUE`.
#'
#' If `NULL`, the function uses `seed_used`.
#'
#' Use this when the user wants reproducibility of the dispersion step independently
#' from the initial design seed.
#'
#' @return A list with:
#' \describe{
#'   \item{layout_matrix}{A character matrix of treatment IDs with dimensions
#'   `n_rows × n_cols`. Unused cells are `NA`.}
#'   \item{field_book}{A data frame with columns `Treatment`, `Family`, `Gcluster`,
#'   `Block`, `Plot`, `Row`, and `Column`. Each row corresponds to an assigned plot.}
#'   \item{efficiency}{`NULL` or a list of efficiency metrics and notes,
#'   depending on `eval_efficiency`.}
#'   \item{seed_used}{The actual random seed used internally.}
#' }
#'
#' @examples
#' data("OptiDesign_example_data", package = "OptiDesign")
#' x <- OptiDesign_example_data
#'
#' ## ---------------------------------------------------------
#' ## Example 1: Family-based repeated-check design
#' ## This example uses the shipped family-based arguments.
#' ## Depending on the supplied treatment lists, the same function
#' ## can represent augmented, p-rep, or RCBD-type repeated-check
#' ## layouts.
#' ## ---------------------------------------------------------
#' out_family <- do.call(
#'   prep_famoptg,
#'   c(
#'     x$OptiDesign_famoptg_example,
#'     x$OptiDesign_famoptg_args_family
#'   )
#' )
#'
#' dim(out_family$layout_matrix)
#' head(out_family$field_book)
#' out_family$seed_used
#'
#' \dontrun{
#' ## ---------------------------------------------------------
#' ## Example 2: GRM-based grouping with dispersion optimization
#' ## ---------------------------------------------------------
#' out_grm <- do.call(
#'   prep_famoptg,
#'   c(
#'     x$OptiDesign_famoptg_example,
#'     x$OptiDesign_famoptg_args_grm
#'   )
#' )
#'
#' dim(out_grm$layout_matrix)
#' head(out_grm$field_book)
#' out_grm$efficiency
#'
#' ## ---------------------------------------------------------
#' ## Example 3: Augmented repeated-check use pattern
#' ## In practice, this is obtained by leaving p-rep arguments empty
#' ## and supplying candidate entries only as unreplicated entries.
#' ## ---------------------------------------------------------
#' aug_args <- x$OptiDesign_famoptg_example
#' aug_args$p_rep_treatments <- character(0)
#' aug_args$p_rep_reps <- integer(0)
#' aug_args$p_rep_families <- character(0)
#'
#' out_aug <- do.call(
#'   prep_famoptg,
#'   c(aug_args, x$OptiDesign_famoptg_args_family)
#' )
#'
#' dim(out_aug$layout_matrix)
#' head(out_aug$field_book)
#'
#' ## ---------------------------------------------------------
#' ## Example 4: RCBD-type repeated-check use pattern
#' ## Here all non-check entries are treated as replicated entries.
#' ## If their common replication equals n_blocks, each entry appears
#' ## once in every block.
#' ## ---------------------------------------------------------
#' rcbd_args <- x$OptiDesign_famoptg_example
#' all_noncheck <- c(rcbd_args$p_rep_treatments, rcbd_args$unreplicated_treatments)
#' all_noncheck_fam <- c(rcbd_args$p_rep_families, rcbd_args$unreplicated_families)
#'
#' rcbd_args$p_rep_treatments <- all_noncheck
#' rcbd_args$p_rep_reps <- rep(2L, length(all_noncheck))
#' rcbd_args$p_rep_families <- all_noncheck_fam
#' rcbd_args$unreplicated_treatments <- character(0)
#' rcbd_args$unreplicated_families <- character(0)
#'
#' out_rcbd_type <- do.call(
#'   prep_famoptg,
#'   c(rcbd_args, x$OptiDesign_famoptg_args_family)
#' )
#'
#' dim(out_rcbd_type$layout_matrix)
#' head(out_rcbd_type$field_book)
#' }
#'
#' @importFrom stats runif setNames
#' @importFrom utils head
#' @importFrom pracma mod
#' @export
#' 
prep_famoptg <- function(
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
    varcomp = list(
      sigma_e2 = 1,
      sigma_g2 = 1,
      sigma_b2 = 1,
      sigma_r2 = 1,
      sigma_c2 = 1
    ),
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
) {
  
  # ============================================================
  # 0. RNG POLICY
  # ============================================================
  # Reproducible if seed provided; otherwise random but record the seed used.
  seed_used <- seed
  if (is.null(seed_used)) {
    seed_used <- sample.int(.Machine$integer.max, 1)
  }
  set.seed(seed_used)
  
  # Utility: set a seed locally without hijacking the global RNG stream
  .with_local_seed <- function(seed_local, expr) {
    old_seed <- NULL
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    
    set.seed(seed_local)
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }
    }, add = TRUE)
    
    force(expr)
  }
  
  # ============================================================
  # 1. INPUT VALIDATION AND ARGUMENT NORMALIZATION
  # ============================================================
  cluster_source <- match.arg(cluster_source)
  cluster_method <- match.arg(cluster_method)
  treatment_effect <- match.arg(treatment_effect)
  prediction_type <- match.arg(prediction_type)
  residual_structure <- match.arg(residual_structure)
  spatial_engine <- match.arg(spatial_engine)
  check_placement <- match.arg(check_placement)
  dispersion_source <- match.arg(dispersion_source)
  
  if (length(check_treatments) != length(check_families)) {
    stop("Length of check_families must match length of check_treatments.")
  }
  if (length(p_rep_treatments) != length(p_rep_reps) ||
      length(p_rep_treatments) != length(p_rep_families)) {
    stop("Lengths of p_rep_treatments, p_rep_reps, and p_rep_families must all match.")
  }
  if (length(unreplicated_treatments) != length(unreplicated_families)) {
    stop("Length of unreplicated_families must match length of unreplicated_treatments.")
  }
  if (any(p_rep_reps > n_blocks)) {
    stop("Each p-rep treatment replication count must not exceed n_blocks.")
  }
  if (n_blocks < 1) {
    stop("n_blocks must be at least 1.")
  }
  if (!order %in% c("column", "row")) {
    stop("Invalid 'order'. Use 'row' or 'column'.")
  }
  if (!(
    is.numeric(n_pcs_use) && length(n_pcs_use) == 1 &&
    (is.finite(n_pcs_use) || is.infinite(n_pcs_use)) && n_pcs_use > 0
  )) {
    stop("n_pcs_use must be a single positive number or Inf.")
  }
  if (!is.list(varcomp) ||
      !all(c("sigma_e2", "sigma_g2", "sigma_b2", "sigma_r2", "sigma_c2") %in% names(varcomp))) {
    stop("varcomp must be a named list with sigma_e2, sigma_g2, sigma_b2, sigma_r2, sigma_c2.")
  }
  if (residual_structure != "IID") {
    if (!is.numeric(rho_row) || length(rho_row) != 1 || abs(rho_row) >= 1) {
      stop("rho_row must satisfy |rho_row| < 1.")
    }
    if (!is.numeric(rho_col) || length(rho_col) != 1 || abs(rho_col) >= 1) {
      stop("rho_col must satisfy |rho_col| < 1.")
    }
  }
  if (dispersion_radius < 1) stop("dispersion_radius must be >= 1.")
  if (dispersion_iters < 0) stop("dispersion_iters must be >= 0.")
  if (check_opt_attempts < 1) stop("check_opt_attempts must be >= 1.")
  
  # ============================================================
  # 2. FIELD DIMENSION ACCOUNTING AND CORRECTION
  # ============================================================
  total_checks <- n_blocks * length(check_treatments)
  total_prep <- sum(p_rep_reps)
  total_unrep <- length(unreplicated_treatments)
  total_required <- total_checks + total_prep + total_unrep
  
  field_size <- n_rows * n_cols
  if (field_size != total_required) {
    if (warn_and_correct) {
      warning(
        paste0(
          "Field size (", n_rows, " x ", n_cols, " = ", field_size,
          ") does not match required (", total_required, "). Adjusting dimensions."
        )
      )
      if (fix_rows) {
        n_cols <- ceiling(total_required / n_rows)
      } else {
        n_rows <- ceiling(total_required / n_cols)
      }
      field_size <- n_rows * n_cols
    } else {
      stop(
        paste0(
          "Provided field size (", field_size,
          ") does not match required (", total_required,
          "). Adjust n_rows/n_cols or enable warn_and_correct."
        )
      )
    }
  }
  
  # ============================================================
  # 3. CLUSTER PREPARATION (FAMILY OR GRM/A-BASED)
  # ============================================================
  family_lookup <- setNames(
    c(check_families, p_rep_families, unreplicated_families),
    c(check_treatments, p_rep_treatments, unreplicated_treatments)
  )
  
  checks_trt <- check_treatments
  noncheck_trt <- c(p_rep_treatments, unreplicated_treatments)
  
  gcluster_lookup <- setNames(
    rep(NA_character_, length(c(checks_trt, noncheck_trt))),
    c(checks_trt, noncheck_trt)
  )
  
  if (cluster_source %in% c("GRM", "A")) {
    Kc <- if (cluster_source == "GRM") GRM else A
    if (is.null(Kc)) stop(paste0("cluster_source='", cluster_source, "' selected but matrix is NULL."))
    if (is.null(rownames(Kc)) || is.null(colnames(Kc))) stop("GRM/A must have rownames and colnames.")
    
    noncheck_fams <- c(p_rep_families, unreplicated_families)
    k_clusters <- length(unique(noncheck_fams))
    if (k_clusters < 2) stop("Non-check treatments have <2 unique families; clustering is not meaningful.")
    
    if (is.null(id_map)) {
      line_ids <- setNames(noncheck_trt, noncheck_trt)
    } else {
      if (!is.data.frame(id_map) || !all(c("Treatment", "LineID") %in% names(id_map))) {
        stop("id_map must be a data.frame with columns: Treatment, LineID")
      }
      line_ids <- setNames(id_map$LineID, id_map$Treatment)
    }
    
    noncheck_line_ids <- unname(line_ids[noncheck_trt])
    missing <- setdiff(noncheck_line_ids, rownames(Kc))
    if (length(missing) > 0) stop("Some non-check LineIDs are not found in GRM/A rownames.")
    
    Ksub <- Kc[noncheck_line_ids, noncheck_line_ids, drop = FALSE]
    eg <- eigen(Ksub, symmetric = TRUE)
    pos <- which(eg$values > 1e-10)
    if (length(pos) < 2) stop("GRM/A has too few positive eigenvalues for PCA clustering.")
    
    max_possible_pcs <- min(length(pos), nrow(Ksub) - 1)
    if (is.infinite(n_pcs_use)) {
      n_pcs <- max_possible_pcs
    } else {
      n_pcs_req <- as.integer(n_pcs_use)
      n_pcs <- min(n_pcs_req, max_possible_pcs)
      if (n_pcs_req > max_possible_pcs) {
        warning(
          paste0(
            "Requested n_pcs_use=", n_pcs_use,
            " but only ", max_possible_pcs, " PCs are available; using ", n_pcs, "."
          )
        )
      }
    }
    if (n_pcs < 2) stop("After bounding, fewer than 2 PCs are available for clustering.")
    
    pcs <- eg$vectors[, pos[seq_len(n_pcs)], drop = FALSE]
    pcs <- sweep(pcs, 2, sqrt(eg$values[pos[seq_len(n_pcs)]]), `*`)
    
    if (cluster_method == "kmeans") {
      clust <- .with_local_seed(cluster_seed, {
        stats::kmeans(pcs, centers = k_clusters, nstart = cluster_attempts)$cluster
      })
    } else {
      d <- stats::dist(pcs)
      hc <- stats::hclust(d, method = "ward.D2")
      clust <- stats::cutree(hc, k = k_clusters)
    }
    
    prefix <- if (cluster_source == "GRM") "G" else "A"
    gcluster_lookup[noncheck_trt] <- paste0(prefix, clust)
  }
  
  get_adj_group <- function(trt) {
    if (cluster_source == "Family") {
      family_lookup[trt]
    } else if (trt %in% check_treatments) {
      family_lookup[trt]
    } else {
      gcluster_lookup[trt]
    }
  }
  
  place_checks_in_block <- function(checks, others, mode = c("random", "systematic")) {
    mode <- match.arg(mode)
    total_len <- length(checks) + length(others)
    
    if (mode == "random") {
      return(sample(c(checks, others)))
    }
    
    pos <- round(seq(1, total_len, length.out = length(checks) + 2))[2:(length(checks) + 1)]
    pos <- pmin(pmax(pos, 1), total_len)
    
    out <- rep(NA_character_, total_len)
    out[pos] <- checks
    out[is.na(out)] <- sample(others)
    out
  }
  
  # ============================================================
  # 4. BLOCK CONSTRUCTION
  # ============================================================
  p_rep_assignments <- vector("list", length(p_rep_treatments))
  names(p_rep_assignments) <- p_rep_treatments
  
  available_blocks <- vector("list", length(p_rep_treatments))
  for (i in seq_along(p_rep_treatments)) {
    available_blocks[[i]] <- sample(seq_len(n_blocks))
  }
  
  for (i in seq_along(p_rep_treatments)) {
    trt_i <- p_rep_treatments[i]
    reps_i <- p_rep_reps[i]
    if (reps_i > length(available_blocks[[i]])) {
      stop(paste0("Not enough unique blocks available for p-rep treatment '", trt_i, "'."))
    }
    assigned <- sample(available_blocks[[i]], reps_i)
    available_blocks[[i]] <- setdiff(available_blocks[[i]], assigned)
    p_rep_assignments[[trt_i]] <- assigned
  }
  
  unrep_treatments_shuffled <- sample(unreplicated_treatments)
  block_unrep_list <- split(
    unrep_treatments_shuffled,
    rep(seq_len(n_blocks), length.out = length(unrep_treatments_shuffled))
  )
  
  build_blocks_once <- function(seed_offset = 0) {
    if (!is.null(seed_used)) {
      set.seed(seed_used + seed_offset)
    }
    
    blocks <- vector("list", n_blocks)
    
    for (b in seq_len(n_blocks)) {
      others <- character(0)
      
      for (trt in names(p_rep_assignments)) {
        if (is.element(b, p_rep_assignments[[trt]])) {
          others <- c(others, trt)
        }
      }
      if (b <= length(block_unrep_list)) {
        others <- c(others, block_unrep_list[[b]])
      }
      
      if (check_placement %in% c("random", "systematic")) {
        block_treatments <- place_checks_in_block(check_treatments, others, mode = check_placement)
      } else {
        block_treatments <- place_checks_in_block(check_treatments, others, mode = "random")
      }
      
      block_family <- vapply(block_treatments, function(trt) family_lookup[trt], character(1))
      block_gcluster <- vapply(block_treatments, function(trt) gcluster_lookup[trt], character(1))
      block_adj <- vapply(block_treatments, function(trt) get_adj_group(trt), character(1))
      
      valid_order <- FALSE
      attempt <- 1
      trt_sh <- block_treatments
      fam_sh <- block_family
      gcl_sh <- block_gcluster
      
      while (!valid_order && attempt <= attempts) {
        ord <- sample(seq_along(block_treatments))
        trt_sh <- block_treatments[ord]
        fam_sh <- block_family[ord]
        gcl_sh <- block_gcluster[ord]
        adj_sh <- block_adj[ord]
        
        if (!any(adj_sh[-1] == adj_sh[-length(adj_sh)])) {
          valid_order <- TRUE
        }
        attempt <- attempt + 1
      }
      
      if (!valid_order) {
        warning(
          paste0(
            "Could not fully avoid adjacent groups in block ", b,
            " after ", attempts, " attempts."
          )
        )
      }
      
      blocks[[b]] <- data.frame(
        Treatment = trt_sh,
        Family = fam_sh,
        Gcluster = gcl_sh,
        Block = b,
        stringsAsFactors = FALSE
      )
    }
    
    blocks
  }
  
  build_positions <- function(n_rows, n_cols, order, serpentine) {
    positions <- vector("list", n_rows * n_cols)
    kpos <- 1
    
    if (order == "row") {
      for (r in seq_len(n_rows)) {
        cols <- seq_len(n_cols)
        if (serpentine && mod(r, 2) == 0) cols <- rev(cols)
        for (c in cols) {
          positions[[kpos]] <- c(Row = r, Column = c)
          kpos <- kpos + 1
        }
      }
    } else {
      for (c in seq_len(n_cols)) {
        rows <- seq_len(n_rows)
        if (serpentine && mod(c, 2) == 0) rows <- rev(rows)
        for (r in rows) {
          positions[[kpos]] <- c(Row = r, Column = c)
          kpos <- kpos + 1
        }
      }
    }
    
    pos_mat <- as.data.frame(do.call(rbind, positions))
    pos_mat$Row <- as.integer(pos_mat$Row)
    pos_mat$Column <- as.integer(pos_mat$Column)
    pos_mat
  }
  
  # ============================================================
  # 5. CHECK PLACEMENT STRATEGY (INCLUDING OPTIMIZATION)
  # ============================================================
  if (check_placement != "optimal") {
    blocks <- build_blocks_once(seed_offset = 0)
  } else {
    best_score <- -Inf
    best_blocks <- NULL
    
    for (k in seq_len(check_opt_attempts)) {
      cand_blocks <- build_blocks_once(seed_offset = k)
      fb_cand <- do.call(rbind, cand_blocks)
      rownames(fb_cand) <- paste0("plot_", seq_len(nrow(fb_cand)))
      
      n_assigned_c <- nrow(fb_cand)
      pos_mat <- build_positions(n_rows, n_cols, order, serpentine)
      fb_cand$Row <- pos_mat$Row[seq_len(n_assigned_c)]
      fb_cand$Column <- pos_mat$Column[seq_len(n_assigned_c)]
      
      is_chk <- fb_cand$Treatment %in% check_treatments
      chk_rc <- unique(cbind(fb_cand$Row[is_chk], fb_cand$Column[is_chk]))
      
      if (nrow(chk_rc) <= 1) {
        score <- 0
      } else {
        d <- as.matrix(stats::dist(chk_rc))
        diag(d) <- Inf
        score <- mean(apply(d, 1, min))
      }
      
      if (score > best_score) {
        best_score <- score
        best_blocks <- cand_blocks
      }
    }
    blocks <- best_blocks
  }
  
  final_data <- do.call(rbind, blocks)
  n_assigned <- nrow(final_data)
  rownames(final_data) <- paste0("plot_", seq_len(n_assigned))
  
  pos_mat <- build_positions(n_rows, n_cols, order, serpentine)
  final_data$Plot <- seq_len(n_assigned)
  final_data$Row <- pos_mat$Row[seq_len(n_assigned)]
  final_data$Column <- pos_mat$Column[seq_len(n_assigned)]
  
  layout_matrix <- matrix(NA, nrow = n_rows, ncol = n_cols)
  for (i in seq_len(n_assigned)) {
    layout_matrix[final_data$Row[i], final_data$Column[i]] <- final_data$Treatment[i]
  }
  
  # ============================================================
  # 6. OPTIONAL GENETIC DISPERSION LOCAL SEARCH
  # ============================================================
  build_neighbor_pairs <- function(row, col, radius = 1) {
    n <- length(row)
    if (n < 2) return(matrix(integer(0), ncol = 2, dimnames = list(NULL, c("i", "j"))))
    
    pairs <- vector("list", 0)
    k <- 1
    for (i in seq_len(n - 1)) {
      dr <- abs(row[i] - row[(i + 1):n])
      dc <- abs(col[i] - col[(i + 1):n])
      ok <- pmax(dr, dc) <= radius
      if (any(ok)) {
        jj <- (i + 1):n
        jj <- jj[ok]
        pairs[[k]] <- cbind(i = rep(i, length(jj)), j = jj)
        k <- k + 1
      }
    }
    if (length(pairs) == 0) return(matrix(integer(0), ncol = 2, dimnames = list(NULL, c("i", "j"))))
    out <- do.call(rbind, pairs)
    colnames(out) <- c("i", "j")
    out
  }
  
  score_dispersion <- function(trt_vec, is_check, Ksub, pairs) {
    if (nrow(pairs) == 0) return(0)
    
    line_levels <- colnames(Ksub)
    line_idx <- rep(NA_integer_, length(trt_vec))
    line_idx[!is_check] <- match(trt_vec[!is_check], line_levels)
    
    ii <- pairs[, "i"]; jj <- pairs[, "j"]
    li <- line_idx[ii]; lj <- line_idx[jj]
    ok <- !is.na(li) & !is.na(lj)
    if (!any(ok)) return(0)
    
    sum(Ksub[cbind(li[ok], lj[ok])])
  }
  
  apply_genetic_dispersion <- function(
    fb,
    check_treatments,
    Kdisp,
    line_id_map,
    radius,
    iters,
    seed_local
  ) {
    .with_local_seed(seed_local, {
      trt <- as.character(fb$Treatment)
      is_check <- trt %in% check_treatments
      movable <- which(!is_check)
      
      if (length(movable) < 2 || iters <= 0) return(fb)
      
      non_trt <- unique(trt[!is_check])
      
      if (is.null(line_id_map)) {
        line_ids <- setNames(non_trt, non_trt)
      } else {
        if (!is.data.frame(line_id_map) || !all(c("Treatment", "LineID") %in% names(line_id_map))) {
          stop("line_id_map must be a data.frame with columns: Treatment, LineID")
        }
        line_ids <- setNames(line_id_map$LineID, line_id_map$Treatment)
      }
      
      ids <- unname(line_ids[non_trt])
      miss <- setdiff(ids, rownames(Kdisp))
      if (length(miss) > 0) stop("Some non-check LineIDs are missing in dispersion matrix rownames/colnames.")
      
      Ksub <- Kdisp[ids, ids, drop = FALSE]
      colnames(Ksub) <- non_trt
      rownames(Ksub) <- non_trt
      
      pairs <- build_neighbor_pairs(fb$Row, fb$Column, radius = radius)
      best_trt <- trt
      best_score <- score_dispersion(best_trt, is_check, Ksub, pairs)
      
      for (it in seq_len(iters)) {
        ij <- sample(movable, 2, replace = FALSE)
        cand <- best_trt
        cand[ij] <- cand[rev(ij)]
        sc <- score_dispersion(cand, is_check, Ksub, pairs)
        if (sc < best_score) {
          best_score <- sc
          best_trt <- cand
        }
      }
      
      fb$Treatment <- best_trt
      fb
    })
  }
  
  if (isTRUE(use_dispersion)) {
    Kdisp <- switch(
      dispersion_source,
      "K" = K,
      "A" = A,
      "GRM" = GRM
    )
    
    if (is.null(Kdisp)) stop("use_dispersion=TRUE but selected dispersion matrix is NULL.")
    if (is.null(rownames(Kdisp)) || is.null(colnames(Kdisp))) stop("Dispersion matrix must have rownames and colnames.")
    
    seed_disp_used <- dispersion_seed
    if (is.null(seed_disp_used)) seed_disp_used <- seed_used
    
    final_data <- apply_genetic_dispersion(
      fb = final_data,
      check_treatments = check_treatments,
      Kdisp = Kdisp,
      line_id_map = line_id_map,
      radius = dispersion_radius,
      iters = dispersion_iters,
      seed_local = seed_disp_used
    )
    
    layout_matrix[,] <- NA
    for (i in seq_len(n_assigned)) {
      layout_matrix[final_data$Row[i], final_data$Column[i]] <- final_data$Treatment[i]
    }
  }
  
  # ============================================================
  # 7. OPTIONAL EFFICIENCY EVALUATION
  # ============================================================
  efficiency <- NULL
  
  if (isTRUE(eval_efficiency) && (treatment_effect == "fixed" || prediction_type != "none")) {
    
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' is required for efficiency evaluation.")
    }
    
    ar1_precision_sparse <- function(nn, rho) {
      if (nn <= 0) stop("n must be >= 1")
      if (nn == 1) return(Matrix::Diagonal(1, 1))
      if (abs(rho) >= 1) stop("AR1 requires |rho| < 1")
      
      a <- 1 / (1 - rho^2)
      d <- rep((1 + rho^2) * a, nn)
      d[1] <- a
      d[nn] <- a
      o <- rep(-rho * a, nn - 1)
      
      Matrix::sparseMatrix(
        i = c(seq_len(nn), seq_len(nn - 1), 2:nn),
        j = c(seq_len(nn), 2:nn, seq_len(nn - 1)),
        x = c(d, o, o),
        dims = c(nn, nn)
      )
    }
    
    make_sparse_incidence <- function(levels_vec) {
      nn <- length(levels_vec)
      lv <- unique(levels_vec[!is.na(levels_vec)])
      if (length(lv) == 0) {
        return(list(M = Matrix::Matrix(0, nrow = nn, ncol = 0, sparse = TRUE), levels = character(0)))
      }
      j <- match(levels_vec, lv)
      ok <- !is.na(j)
      M <- Matrix::sparseMatrix(i = which(ok), j = j[ok], x = 1, dims = c(nn, length(lv)))
      colnames(M) <- lv
      list(M = M, levels = lv)
    }
    
    pinv_sym_dense <- function(A, tol = 1e-10) {
      eg <- eigen(A, symmetric = TRUE)
      vals <- eg$values
      vecs <- eg$vectors
      keep <- vals > tol
      if (!any(keep)) return(matrix(0, nrow(A), ncol(A)))
      vecs[, keep, drop = FALSE] %*% diag(1 / vals[keep]) %*% t(vecs[, keep, drop = FALSE])
    }
    
    solve_C <- function(Cmat, B) {
      out <- try({
        fac <- Matrix::Cholesky(Cmat, LDL = TRUE, Imult = 0)
        Matrix::solve(fac, B)
      }, silent = TRUE)
      if (inherits(out, "try-error")) out <- Matrix::solve(Cmat, B)
      out
    }
    
    trace_subinv_est <- function(Cmat, idx, m = 80, seed_local = 1) {
      .with_local_seed(seed_local, {
        p <- length(idx)
        if (p == 0) return(NA_real_)
        nn <- nrow(Cmat)
        acc <- 0
        for (k in seq_len(m)) {
          z <- sample(c(-1, 1), p, replace = TRUE)
          u <- Matrix::sparseVector(i = idx, x = z, length = nn)
          x <- solve_C(Cmat, u)
          acc <- acc + as.numeric(Matrix::crossprod(u, x))
        }
        acc / m
      })
    }
    
    safe_logdet_psd_dense <- function(A, tol = 1e-10) {
      eg <- eigen(A, symmetric = TRUE)
      lam <- eg$values
      lam_pos <- lam[lam > tol]
      if (length(lam_pos) == 0) return(-Inf)
      sum(log(lam_pos))
    }
    
    ar1_cov_dense <- function(nn, rho) {
      idx <- seq_len(nn)
      outer(idx, idx, function(i, j) rho^abs(i - j))
    }
    
    pairwise_diff_mean_var <- function(V) {
      p <- nrow(V)
      if (p < 2) return(NA_real_)
      one <- rep(1, p)
      num <- p * sum(diag(V)) - as.numeric(t(one) %*% V %*% one)
      (2 * num) / (p * (p - 1))
    }
    
    fb <- final_data
    nn <- nrow(fb)
    
    trt <- as.character(fb$Treatment)
    blk <- as.character(fb$Block)
    rw  <- as.character(fb$Row)
    cl  <- as.character(fb$Column)
    
    is_check <- trt %in% check_treatments
    
    sigma_e2 <- varcomp$sigma_e2
    sigma_g2 <- varcomp$sigma_g2
    sigma_b2 <- varcomp$sigma_b2
    sigma_r2 <- varcomp$sigma_r2
    sigma_c2 <- varcomp$sigma_c2
    
    if (any(c(sigma_e2, sigma_b2, sigma_r2, sigma_c2) <= 0)) {
      stop("Variance components sigma_e2, sigma_b2, sigma_r2, sigma_c2 must be > 0.")
    }
    if (treatment_effect == "random" && sigma_g2 <= 0) {
      stop("sigma_g2 must be > 0 when treatment_effect='random'.")
    }
    
    has_holes <- n_assigned < (n_rows * n_cols)
    
    spatial_engine_use <- spatial_engine
    if (spatial_engine_use == "auto") {
      spatial_engine_use <- if (nn <= dense_max_n) "dense" else "sparse"
    }
    
    eff_notes <- character(0)
    residual_use <- residual_structure
    if (residual_structure != "IID" && has_holes) {
      eff_notes <- c(eff_notes, "Spatial residual requested and layout has holes; AR1 submatrix on observed plots was used.")
    }
    
    if (residual_use == "IID") {
      Q <- (1 / sigma_e2) * Matrix::Diagonal(nn, 1)
      
    } else if (spatial_engine_use == "dense") {
      
      if (residual_use == "AR1") {
        Rrow <- ar1_cov_dense(n_rows, rho_row)
        Rcol <- diag(n_cols)
      } else if (residual_use == "AR1xAR1") {
        Rrow <- ar1_cov_dense(n_rows, rho_row)
        Rcol <- ar1_cov_dense(n_cols, rho_col)
      } else stop("Unsupported residual_structure.")
      
      Rgrid <- sigma_e2 * kronecker(Rcol, Rrow)
      grid_index <- (as.integer(rw) - 1) * n_cols + as.integer(cl)
      Robs <- Rgrid[grid_index, grid_index, drop = FALSE]
      Q <- Matrix::Matrix(solve(Robs), sparse = FALSE)
      
    } else {
      
      Qrow <- ar1_precision_sparse(n_rows, rho_row)
      if (residual_use == "AR1") {
        Qcol <- Matrix::Diagonal(n_cols, 1)
      } else if (residual_use == "AR1xAR1") {
        Qcol <- ar1_precision_sparse(n_cols, rho_col)
      } else stop("Unsupported residual_structure.")
      
      Qgrid <- Matrix::kronecker(Qcol, Qrow) * (1 / sigma_e2)
      grid_index <- (as.integer(rw) - 1) * n_cols + as.integer(cl)
      Q <- Qgrid[grid_index, grid_index, drop = FALSE]
    }
    
    X_int <- Matrix::Matrix(rep(1, nn), ncol = 1, sparse = TRUE)
    colnames(X_int) <- "(Intercept)"
    X <- X_int
    
    if (isTRUE(check_as_fixed)) {
      chk_vec <- ifelse(is_check, trt, NA)
      chk_inc <- make_sparse_incidence(chk_vec)
      if (ncol(chk_inc$M) > 0) {
        keep <- intersect(check_treatments, colnames(chk_inc$M))
        chkM <- chk_inc$M[, keep, drop = FALSE]
        chkM <- chkM[, check_treatments[check_treatments %in% keep], drop = FALSE]
        colnames(chkM) <- paste0("Check_", colnames(chkM))
        X <- cbind(X, chkM)
      }
    }
    
    if (treatment_effect == "fixed") {
      trt_fix_vec <- ifelse(!is_check, trt, NA)
      trt_inc <- make_sparse_incidence(trt_fix_vec)
      if (ncol(trt_inc$M) < 2) stop("Not enough fixed non-check treatments to compute efficiency.")
      colnames(trt_inc$M) <- paste0("Line_", colnames(trt_inc$M))
      X <- cbind(X, trt_inc$M)
    }
    
    Zb <- make_sparse_incidence(blk)$M
    Zr <- make_sparse_incidence(rw)$M
    Zc <- make_sparse_incidence(cl)$M
    colnames(Zb) <- paste0("Block_", colnames(Zb))
    colnames(Zr) <- paste0("Row_",   colnames(Zr))
    colnames(Zc) <- paste0("Col_",   colnames(Zc))
    
    Z_list <- list(Zb = Zb, Zr = Zr, Zc = Zc)
    
    Zg <- NULL
    if (treatment_effect == "random") {
      g_vec <- ifelse(!is_check, trt, NA)
      g_inc <- make_sparse_incidence(g_vec)
      Zg <- g_inc$M
      if (ncol(Zg) < 2) stop("Not enough random non-check treatments to compute efficiency.")
      colnames(Zg) <- paste0("Line_", colnames(Zg))
      Z_list$Zg <- Zg
    }
    
    Z <- do.call(cbind, Z_list)
    
    XtQX <- Matrix::crossprod(X, Q %*% X)
    XtQZ <- Matrix::crossprod(X, Q %*% Z)
    ZtQX <- t(XtQZ)
    ZtQZ <- Matrix::crossprod(Z, Q %*% Z)
    
    Ginv <- Matrix::Diagonal(ncol(Z), 0)
    
    idx0 <- 0
    pB <- ncol(Zb); pR <- ncol(Zr); pC <- ncol(Zc)
    
    if (pB > 0) {
      ii <- (idx0 + 1):(idx0 + pB)
      Ginv[ii, ii] <- Matrix::Diagonal(pB, 1 / sigma_b2)
      idx0 <- idx0 + pB
    }
    if (pR > 0) {
      ii <- (idx0 + 1):(idx0 + pR)
      Ginv[ii, ii] <- Matrix::Diagonal(pR, 1 / sigma_r2)
      idx0 <- idx0 + pR
    }
    if (pC > 0) {
      ii <- (idx0 + 1):(idx0 + pC)
      Ginv[ii, ii] <- Matrix::Diagonal(pC, 1 / sigma_c2)
      idx0 <- idx0 + pC
    }
    
    trt_idx <- integer(0)
    
    if (treatment_effect == "random") {
      pG <- ncol(Zg)
      trt_idx <- (idx0 + 1):(idx0 + pG)
      
      if (prediction_type == "IID") {
        Ginv[trt_idx, trt_idx] <- Matrix::Diagonal(pG, 1 / sigma_g2)
      } else {
        if (is.null(K)) stop("K must be provided for prediction_type='GBLUP'/'PBLUP'.")
        if (is.null(rownames(K)) || is.null(colnames(K))) stop("K must have rownames and colnames.")
        
        line_trt_order <- sub("^Line_", "", colnames(Zg))
        
        if (is.null(line_id_map)) {
          line_ids <- setNames(line_trt_order, line_trt_order)
        } else {
          if (!is.data.frame(line_id_map) || !all(c("Treatment", "LineID") %in% names(line_id_map))) {
            stop("line_id_map must be a data.frame with columns: Treatment, LineID")
          }
          line_ids <- setNames(line_id_map$LineID, line_id_map$Treatment)
        }
        
        ids <- unname(line_ids[line_trt_order])
        miss <- setdiff(ids, rownames(K))
        if (length(miss) > 0) stop("Some line IDs are missing in K rownames/colnames.")
        
        Ksub2 <- K[ids, ids, drop = FALSE]
        
        Kinv_try <- try({
          Km <- Matrix::Matrix(Ksub2, sparse = FALSE)
          facK <- Matrix::Cholesky(Km, LDL = FALSE, Imult = 0)
          Matrix::solve(facK, Matrix::Diagonal(nrow(Km), 1))
        }, silent = TRUE)
        
        if (!inherits(Kinv_try, "try-error")) {
          Kinv <- Kinv_try
        } else {
          if (nrow(Ksub2) > 2500) warning("K is large and not Cholesky-factorable; dense generalized inverse may be slow.")
          Kinv <- Matrix::Matrix(pinv_sym_dense(Ksub2), sparse = FALSE)
        }
        
        Ginv[trt_idx, trt_idx] <- (1 / sigma_g2) * Kinv
      }
    }
    
    Cmat <- rbind(cbind(XtQX, XtQZ), cbind(ZtQX, ZtQZ + Ginv))
    if (spatial_engine_use == "dense") {
      Cmat <- as.matrix(Cmat)
    } else {
      Cmat <- Matrix::Matrix(Cmat, sparse = TRUE)
    }
    
    eff <- list(
      model = "mu + Treatment + Block + Row + Column + e",
      treatment_effect = treatment_effect,
      prediction_type = if (treatment_effect == "random") prediction_type else NA_character_,
      residual_structure_requested = residual_structure,
      residual_structure_used = residual_use,
      spatial_engine_used = spatial_engine_use,
      notes = eff_notes
    )
    
    if (treatment_effect == "fixed") {
      
      xnames <- colnames(X)
      trt_cols <- grep("^Line_", xnames)
      if (length(trt_cols) < 2) stop("Not enough fixed non-check treatments to compute efficiency.")
      p <- length(trt_cols)
      
      if (p <= eff_full_max) {
        B <- Matrix::sparseMatrix(i = trt_cols, j = seq_along(trt_cols), x = 1, dims = c(nrow(Cmat), p))
        Xsol <- solve_C(Cmat, B)
        Vsub <- as.matrix(Xsol[trt_cols, , drop = FALSE])
        
        mean_var_diff <- pairwise_diff_mean_var(Vsub)
        A_opt <- 1 / mean_var_diff
        
        H <- diag(p) - matrix(1 / p, p, p)
        Vctr <- H %*% Vsub %*% H
        logdet <- safe_logdet_psd_dense(Vctr)
        D_opt <- if (is.finite(logdet)) exp(-logdet / (p - 1)) else NA_real_
        
        eff$mode <- "FIXED_TREATMENT_BLUE_CONTRAST"
        eff$A <- A_opt
        eff$D <- D_opt
        eff$mean_VarDiff <- mean_var_diff
        eff$n_trt <- p
        
      } else {
        tr_est <- trace_subinv_est(Cmat, trt_cols, m = eff_trace_samples, seed_local = 1)
        mean_var_coef <- tr_est / p
        A_opt <- 1 / mean_var_coef
        
        eff$mode <- "FIXED_TREATMENT_BLUE_APPROX"
        eff$A <- A_opt
        eff$D <- NA_real_
        eff$mean_VarCoef <- mean_var_coef
        eff$mean_Var <- mean_var_coef
        eff$n_trt <- p
        eff$notes <- c(
          eff$notes,
          paste0(
            "Fixed-treatment target dimension (", p,
            ") > eff_full_max (", eff_full_max,
            "); A via stochastic trace estimator; D = NA."
          )
        )
      }
      
    } else {
      
      p <- length(trt_idx)
      if (p < 2) stop("Not enough random non-check treatments to compute efficiency.")
      
      if (p <= eff_full_max) {
        B <- Matrix::sparseMatrix(i = trt_idx, j = seq_along(trt_idx), x = 1, dims = c(nrow(Cmat), p))
        Xsol <- solve_C(Cmat, B)
        PEVsub <- as.matrix(Xsol[trt_idx, , drop = FALSE])
        
        A_pred <- 1 / mean(diag(PEVsub))
        logdet <- safe_logdet_psd_dense(PEVsub)
        D_pred <- if (is.finite(logdet)) exp(-logdet / p) else NA_real_
        
        eff$mode <- paste0("RANDOM_TREATMENT_BLUP_", prediction_type)
        eff$A <- A_pred
        eff$D <- D_pred
        eff$mean_PEV <- mean(diag(PEVsub))
        eff$n_lines <- p
        
      } else {
        tr_est <- trace_subinv_est(Cmat, trt_idx, m = eff_trace_samples, seed_local = 1)
        mean_pev <- tr_est / p
        A_pred <- 1 / mean_pev
        
        eff$mode <- paste0("RANDOM_TREATMENT_BLUP_", prediction_type, "_APPROX")
        eff$A <- A_pred
        eff$D <- NA_real_
        eff$mean_PEV <- mean_pev
        eff$n_lines <- p
        eff$notes <- c(
          eff$notes,
          paste0(
            "Random-treatment target dimension (", p,
            ") > eff_full_max (", eff_full_max,
            "); A via stochastic trace estimator; D = NA."
          )
        )
      }
    }
    
    efficiency <- eff
  }
  
  # ============================================================
  # 8. RETURN OBJECT
  # ============================================================
  list(
    layout_matrix = layout_matrix,
    field_book = final_data,
    efficiency = efficiency,
    seed_used = seed_used
  )
}
