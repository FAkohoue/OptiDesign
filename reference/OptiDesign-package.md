# OptiDesign: Experimental Field Design Utilities for Optimized Layout Construction

`OptiDesign` provides tools for constructing optimized experimental
field designs for plant breeding and related applications. The package
focuses on layouts where field structure, treatment grouping, and
optional mixed-model evaluation are integrated into the design stage
rather than treated as separate downstream tasks.

## Details

The package currently includes two main design constructors:

- [`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)
  for repeated-check block designs with flexible replication, supporting
  augmented designs, partially replicated (p-rep) designs, and RCBD-type
  repeated-check designs. This function allows optional grouping from
  family labels or relationship matrices, optional dispersion
  optimization, and optional mixed-model efficiency evaluation.

- [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md)
  for fixed-grid alpha row-column stream designs where the field is
  treated as a global traversal stream, replicates are contiguous stream
  segments, incomplete blocks may be unequal, checks are included in
  every block, and unused cells appear only at the end of the stream.

The package is intended for users who need designs that are not only
feasible in the field, but also informed by:

- family structure,

- pedigree or genomic relationship matrices,

- spatial distribution of related materials,

- and mixed-model design efficiency considerations.

## Main objectives of the package

`OptiDesign` is particularly useful when the user wants to:

- construct designs with repeated checks and heterogeneous entry
  replication,

- generate augmented, p-rep, or balanced repeated-check layouts,

- avoid local clustering of related entries,

- use family, pedigree, or genomic structure to guide placement,

- preserve realistic field-book order and layout constraints,

- evaluate a design under fixed- or random-effect mixed models,

- compare alternative layout strategies before field implementation.

## Relationship matrices and grouping

Several design routines in the package can use one of three grouping
sources:

- direct family labels,

- a genomic relationship matrix (`GRM`),

- a pedigree relationship matrix (`A`).

In addition, some workflows may use a separate relationship matrix `K`
for:

- BLUP-style prediction efficiency calculations,

- or dispersion optimization.

When treatment names do not match matrix row names directly, mapping
tables such as `id_map` or `line_id_map` can be used to reconcile
identifiers.

## Included example dataset

The package ships with the dataset `OptiDesign_example_data`, which
contains synthetic example objects for:

- treatment and family lists,

- example relationship matrices,

- example argument lists for both main design functions.

Load it with:

`data("OptiDesign_example_data", package = "OptiDesign")`

## Imported functionality

The package imports `mod()` from **pracma** for serpentine parity
handling in layouts where alternating row or column direction is
required.

Matrix computations used in efficiency evaluation and some optimization
steps rely on the **Matrix** package.

## See also

Useful links:

- <https://github.com/FAkohoue/OptiDesign>

- <https://FAkohoue.github.io/OptiDesign/>

- Report bugs at <https://github.com/FAkohoue/OptiDesign/issues>

## Author

**Maintainer**: Félicien Akohoue <akohoue.f@gmail.com>
([ORCID](https://orcid.org/0000-0002-2160-0182))
