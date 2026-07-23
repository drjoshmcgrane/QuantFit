# QuantFit

Latent Structure Model Selection Framework for R

## Automated TI&D lattice

`select_model_ll()` now defaults to the successive Torres Irribarra–Diakow
lattice. It tests both UN→MON/IIO→DM paths edge by edge, then follows the
published DM→LCR→RM sequence. The DM-to-LCR scale comparison uses the fixed
minimum support-point count required by the paper's LCR/Rasch equivalence
argument (`ceiling((J + 1) / 2)` for dichotomous items), so a low selected
class count cannot block true Rasch data. A calibrated, class-profiled
RM-versus-LCR comparison then resolves continuity/grain size. Decisions use
corrected bootstrap p-values without a heuristic override by default. When a
class-count range is supplied, the UN-BIC ordinal enumeration is repeated in
every applicable null replicate so class selection is part of the calibration.
If finite-mixture optimisation leaves a general model below its nested child,
the general model receives one EM refinement initialized at the child solution;
the identical safeguard is used in bootstrap replicates, and its use is exposed
in the returned diagnostics.
The significance levels are per edge, not familywise over the entire path; a
true deep model can be stopped by a chance rejection at any preceding edge.
This is one reason exact recovery of Rasch data need not be perfect even when
every individual test is correctly calibrated.

Full calibration is computationally intensive because every lattice edge
refits multiple finite-mixture models in every bootstrap replicate. Runtime
grows sharply with test length and the class-count grid; long tests can take
hours even with parallel bootstrap workers. Use a small `B` only for code
shakedowns, and use `B = 99` or more for a reported 5% decision.

## Overview

**QuantFit** implements the model selection framework from Torres Irribarra & Diakow's paper "Categorization, Ordering and Quantification: Selecting a Latent Variable Model by Comparing Latent Structures."

The package enables researchers to compare six latent structure models to determine whether data supports classificatory, ordinal, or quantitative interpretations of a latent variable.

## The Six Models

| Model | Code | Structure | Key Constraints |
|-------|------|-----------|-----------------|
| Unconstrained Latent Class | `UN` | Qualitative | None (baseline) |
| Class Monotonicity | `MON` | Ordinal | β_ic ≤ β_ic' for c < c' |
| Invariant Item Ordering | `IIO` | Ordinal | β_ic ≤ β_i'c for i < i' |
| Double Monotonicity | `DM` | Ordinal | Both MON + IIO |
| Latent Class Rasch | `LCR` | Quantitative | θ_c - δ_i parameterization |
| Rasch Model | `RM` | Quantitative | Continuous θ |

## Installation

```r
# Install from GitHub (once available)
# devtools::install_github("quantfit/QuantFit")

# Or install locally
devtools::install("path/to/QuantFit")
```
## Quick Start

```r
library(QuantFit)

# Fit individual models
un_fit <- fit_un(data, n_classes = 3)
mon_fit <- fit_mon(data, n_classes = 3)
lcr_fit <- fit_lcr(data, n_classes = 3)

# Compare all models
comparison <- compare_models(data, n_classes = 3)
print(comparison)

# Successive comparison strategy
result <- successive_comparison(data, n_classes = 3)
summary(result)

# Visualize item response profiles
plot_irfs(mon_fit)
```

## Key Functions

### Model Fitting
- `fit_un()` - Unconstrained Latent Class Analysis
- `fit_mon()` - Class Monotonicity model
- `fit_iio()` - Invariant Item Ordering model
- `fit_dm()` - Double Monotonicity model
- `fit_lcr()` - Latent Class Rasch model
- `fit_rm()` - Rasch/partial-credit model (package MML engine)

### Model Comparison
- `compare_models()` - Compare multiple models using AIC/BIC
- `successive_comparison()` - Stepwise comparison following paper's strategy

### Visualization
- `plot_irfs()` - Item response function plots
- `plot_comparison()` - Model comparison plots

## Dependencies

**Required:**
- alabama (constrained optimization)
- nloptr (alternative optimizer)
- Matrix (sparse matrices)
- numDeriv (numerical derivatives)

**Suggested:**
- testthat (testing)
- mirt (optional cross-validation)
- poLCA (validation)
- eRm (validation)
- ggplot2 (enhanced plots)

## References

Torres Irribarra, D., & Diakow, R. Categorization, Ordering and Quantification: Selecting a Latent Variable Model by Comparing Latent Structures.

## License

MIT License
