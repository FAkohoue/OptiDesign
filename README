# OptiDesign

OptiDesign provides advanced tools for constructing optimized experimental field designs for plant breeding and agronomic experiments.

The package integrates:

- field layout construction  
- genetic structure (family, pedigree, genomic relationships)  
- optional spatial dispersion optimization  
- optional mixed-model efficiency evaluation  

into a single, flexible workflow.

---

## Why OptiDesign?

In many breeding programs, experimental design is treated as a purely logistical step. However, design choices strongly affect:

- precision of estimates  
- prediction accuracy  
- ability to separate genetic effects  
- robustness to field heterogeneity  

OptiDesign allows users to **explicitly control these aspects at the design stage**, rather than correcting them later during analysis.

---

## Core Design Concepts

OptiDesign is built around three key ideas:

### 1. Field realism

Designs reflect how trials are actually implemented:

- fixed grid layouts  
- field-book ordering  
- serpentine movement  
- contiguous replicate structure  

---

### 2. Genetic awareness

Entries can be arranged based on:

- family labels  
- pedigree relationships (A matrix)  
- genomic relationships (GRM)  

This allows:

- reduction of local relatedness  
- improved sampling of genetic diversity  
- better estimation of genetic effects  

---

### 3. Integrated evaluation

Designs can be evaluated directly using mixed-model principles:

- fixed-effect precision (BLUEs)  
- random-effect prediction (BLUP / GBLUP / PBLUP)  
- spatial residual structures (IID, AR1, AR1×AR1)  

---

## Main Functions

---

### `prep_famoptg()`

#### Purpose

Construct **repeated-check block designs with flexible replication**, including:

- **augmented designs**  
- **partially replicated (p-rep) designs**  
- **RCBD-type repeated-check designs**  

where:

- checks are included in every block  
- some entries may be replicated  
- others may be unreplicated  

---

#### When to use

Use this function when:

- you have **many entries** but limited resources  
- you want **checks in every block**  
- you want to **balance replication and coverage**  
- you are working in **early- or intermediate-stage breeding trials**  
- you want a **flexible design framework covering augmented, p-rep, or balanced layouts**  

---

#### Key capabilities

- Flexible replication per entry  
- Block-wise allocation of treatments (no treatment repeats within a block)  
- Supports:
  - augmented designs  
  - p-rep designs  
  - balanced repeated-check (RCBD-type) designs  
- Optional grouping using:
  - family labels  
  - genomic relationship matrix (GRM)  
  - pedigree matrix (A)  
- Optional **dispersion optimization** to reduce clustering of related entries  
- Optional **efficiency diagnostics** under mixed models  

---

#### Key design rule

For replicated entries:

> A treatment can appear multiple times overall, but **at most once per block**.

This ensures proper distribution of replicates across blocks.

---

#### Typical workflow

1. Define:
   - checks  
   - replicated entries (optional)  
   - unreplicated entries (optional)  
2. Specify:
   - number of blocks  
   - field dimensions  
3. Choose grouping source (family, GRM, or A)  
4. Optionally:
   - optimize dispersion  
   - compute efficiency  

---

---

### `alpha_rc_stream()`

#### Purpose

Construct **alpha row–column designs on a fixed grid**, using a **stream-based layout**.

---

#### Key idea: Stream-based design

Instead of dividing the field into rectangles, the field is:

1. Converted into a **single ordered stream of positions**
2. Split into **replicates as contiguous segments**
3. Split further into **incomplete blocks**

---

#### When to use

Use this function when:

- field dimensions are **fixed and cannot change**  
- planting follows **field-book order**  
- replicates are defined operationally (not geometrically)  
- you need **checks in every incomplete block**  
- classical rectangular alpha designs are not practical  

---

#### Key capabilities

- Fixed grid (`n_rows × n_cols`)  
- Replicates as contiguous field segments  
- Unequal incomplete block sizes  
- Checks in every block  
- One observation per entry per replicate  
- Optional:
  - grouping (family, GRM, A)  
  - dispersion optimization  
  - efficiency diagnostics  

---

#### Important consequence

Unused cells:

- are not scattered  
- are placed **only at the end of the field stream**  

This is critical for practical field implementation.

---

## Grouping Options

Both main functions support three grouping strategies:

### Family-based

- simplest and most interpretable  
- uses user-defined family labels  

### GRM-based (genomic)

- uses genomic similarity  
- better captures real genetic relationships  

### Pedigree-based (A matrix)

- uses pedigree structure  
- useful when genomic data is unavailable  

---

## Dispersion Optimization

Optional step to improve spatial distribution of entries.

### Purpose

Reduce the probability that:

- closely related genotypes  
- or highly similar entries  

are placed near each other in the field.

---

### When to use

- genomic or pedigree structure is important  
- strong local correlation is expected  
- you want to avoid confounding spatial and genetic effects  

---

## Efficiency Evaluation

Optional but powerful feature.

### What it does

Evaluates the design using a mixed-model framework:

- fixed-effect precision (BLUEs)  
- random-effect prediction (BLUP / GBLUP / PBLUP)  
- spatial models (AR1, AR1×AR1)  

---

### When to use

- comparing alternative designs  
- optimizing design before field implementation  
- studying impact of block structure and replication  

---

## Installation

**Build Vignettes**

```r
# install.packages("remotes")
remotes::install_github("FAkohoue/OptiDesign", build_vignettes = TRUE,
  dependencies = TRUE
)
```

**Without Vignettes** 

```r
# install.packages("remotes")
remotes::install_github("FAkohoue/OptiDesign", build_vignettes = FALSE,
  dependencies = TRUE
)
```
---

## Documentation

Full documentation and tutorials are available at: https://FAkohoue.github.io/OptiDesign/

---

# Citation

Akohoue, F. (2026). OptiDesign: Experimental Field Design Utilities for Optimized Layout Construction.
https://github.com/FAkohoue/OptiDesign


# Contributing

Issues and suggestions are welcome: https://github.com/FAkohoue/OptiDesign/issues