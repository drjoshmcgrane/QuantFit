# QuantFit

Latent Structure Model Selection Framework for R

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
- `fit_rm()` - Rasch Model (mirt wrapper)

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
- mirt (Rasch model)
- Matrix (sparse matrices)
- numDeriv (numerical derivatives)

**Suggested:**
- testthat (testing)
- poLCA (validation)
- eRm (validation)
- ggplot2 (enhanced plots)

## References

Torres Irribarra, D., & Diakow, R. Categorization, Ordering and Quantification: Selecting a Latent Variable Model by Comparing Latent Structures.

## License

MIT License
