# Package index

## Alpha row-column stream designs

Functions for constructing, evaluating, and optimising fixed-grid alpha
row-column designs with repeated checks. Optimisation supports Random
Restart, Simulated Annealing, and a Genetic Algorithm.

- [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiDesign/reference/alpha_rc_stream.md)
  : Construct a stream-based repeated-check alpha row-column design
- [`evaluate_alpha_efficiency()`](https://FAkohoue.github.io/OptiDesign/reference/evaluate_alpha_efficiency.md)
  : Evaluate the statistical efficiency of an alpha-lattice design
- [`optimize_alpha_rc()`](https://FAkohoue.github.io/OptiDesign/reference/optimize_alpha_rc.md)
  : Search for a criterion-optimal alpha-lattice design

## Repeated-check block designs

Functions for constructing, evaluating, and optimising augmented,
partially replicated (p-rep), and RCBD-type repeated-check block
designs. The p-rep constraint (no replicated entry in the same block
twice) is enforced by construction.

- [`evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiDesign/reference/evaluate_famoptg_efficiency.md)
  : Evaluate the statistical efficiency of a repeated-check block design
- [`optimize_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/optimize_famoptg.md)
  : Search for a criterion-optimal repeated-check block design
- [`prep_famoptg()`](https://FAkohoue.github.io/OptiDesign/reference/prep_famoptg.md)
  : Construct a repeated-check block design with flexible replication

## Datasets

Synthetic example data for illustrating and testing both design
families, including treatment vectors, field dimensions, relationship
matrices, and ready-to-use argument lists.

- [`OptiDesign_example_data`](https://FAkohoue.github.io/OptiDesign/reference/OptiDesign_example_data.md)
  : Example data for OptiDesign

## Package

Package-level documentation covering the full function overview, typical
workflows for both design families, and key references.

- [`OptiDesign`](https://FAkohoue.github.io/OptiDesign/reference/OptiDesign-package.md)
  [`OptiDesign-package`](https://FAkohoue.github.io/OptiDesign/reference/OptiDesign-package.md)
  : OptiDesign: Optimized Experimental Field Design for Plant Breeding
