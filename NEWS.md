# QuantFit 0.2.0

First release of QuantFit as a unified package for checking quantitative
structure in item response data. It brings together three previously separate
lines of work into one framework with a common, optimised C++ core.

## Three approaches in one package

* **Conjoint cancellation checks** — Bayesian tests of the single, double, and
  triple cancellation axioms of additive conjoint measurement, extending the
  `ConjointChecks` package of Domingue (2014). Functions: `ConjointChecks()`,
  `SingleCancel()`, `DoubleCancel()`, `TripleCancel()`, `HiConjointChecks()`
  (hierarchical 4x4 test), and `PrepareChecks()` (build the score-by-item
  count matrices from raw responses).

* **Karabatsos omnibus test** — the Bayesian test of the additive conjoint
  measurement axioms using synthetic likelihood of Karabatsos (2018), ported
  to R/C++ and validated against the original MATLAB implementation
  (`ACMtest.m`) to a cell-level correlation > 0.99. Function: `KaraChecks()`.

* **Latent structure model selection** — the Torres Irribarra & Diakow
  framework comparing six latent structure models (UN, MON, IIO, DM, LCR, RM)
  to decide whether data support a classificatory, ordinal, or quantitative
  interpretation. Fit with `fit_un()`, `fit_mon()`, `fit_iio()`, `fit_dm()`,
  `fit_lcr()`, `fit_rm()`; compare and select with `compare_models()`,
  `successive_comparison()`, `select_model_constraint()`, and
  `select_model_ll()`.

## New in this release

* The triangulated verdict function is now `quant_fit()` - the package's
  namesake flagship function. `assess_quantitative()` is kept as a deprecated
  alias.

* `kara_bootstrap_null()` - applies the same Rasch bootstrap-null calibration
  to Karabatsos's global KL statistic that `cc_bootstrap_null()` applies to the
  cancellation-check violation rate, giving a per-dataset percentile p-value in
  place of a fixed KL > 0.01 cutoff. With this, all three routes of
  `assess_quantitative()` are calibrated the same way - a per-dataset
  parametric bootstrap with an accept/reject percentile - so their evidence is
  statistically consistent and comparable. The final LC step (RM vs LCR) is
  likewise now bootstrap-calibrated (the BIC difference against a Rasch null)
  rather than a raw BIC comparison, removing the last non-bootstrap criterion.

* `cc_bootstrap_null()` - calibrates the [ConjointChecks()] violation rate
  against a null distribution simulated from the Rasch model fitted to the
  data, following Student & Read (2025). The raw violation rate is not
  interpretable on its own (even interval-scalable Rasch data violate at a
  rate that depends on sample size, test length, and item parameters), so the
  observed rate's percentile within the Rasch null is treated as a p-value -
  interval scaling is rejected when it exceeds, say, the 95th percentile.
  Because observed and null data share the same pipeline, the null
  self-calibrates the baseline, so no fixed cutoff is needed. Includes `print`
  and `plot` methods and an under-power warning below ~1000 examinees.

* `assess_quantitative()` - a single triangulated judgement on whether data
  support a quantitative (interval / additive) interpretation, combining all
  three of the package's routes, each with its appropriate calibration:
  latent-structure model selection ([select_model_ll()], bootstrap
  chi-bar-squared) on the raw data; the cancellation checks calibrated against
  a Rasch bootstrap null ([cc_bootstrap_null()], per Student & Read 2025); and
  the Karabatsos synthetic-likelihood test ([KaraChecks()], his KL > 0.01
  criterion) on an ability-banded matrix. The Kara banding is deliberate - the
  KL test reads additive Rasch data as non-additive on raw sum-score groups,
  whereas ability bands with a real ability metric are well calibrated; the CC
  route needs no banding because its bootstrap null self-calibrates the
  pipeline. The verdict names the supporting routes and reports each route's
  statistics. It states its limits: the CC route is under-powered below ~1000
  examinees, and the (unidimensional) Kara banding can miss multidimensional
  departures.


* `select_model_ll()` and `ll_equivalence_test()` — a statistically calibrated
  model selection procedure. Because UN and the ordinal models (MON/IIO/DM)
  have identical parameter counts, information criteria cannot distinguish
  them; these functions instead use a parametric bootstrap of the
  likelihood-ratio statistic, whose asymptotic null is chi-bar-squared, to
  test each constrained model against its less-restricted parent. Selection
  proceeds down the hierarchy UN -> MON/IIO/DM -> LCR -> RM, using BIC only at
  the final step where parameter counts genuinely differ. A `method` argument
  offers two ways to test the ordinal layer: `"joint"` (default) tests the
  doubly-monotone model against UN directly, while `"lattice"` tests each
  constraint edge separately; a paired simulation study found them
  statistically indistinguishable in recovery, so the cheaper joint procedure
  is the default.

* `compute_se()` — observed-information (Hessian, with delta-method
  back-transform) and nonparametric bootstrap standard errors for all six
  models. RM standard errors are taken from `mirt`; parameters sitting on an
  active monotonicity constraint return `NA` with a warning, since Hessian
  standard errors are not valid on the boundary.

* `select_n_classes()` — chooses the number of latent classes by information
  criterion, fitting the unconstrained model across a range (default `C = 1`
  to `6`, where `C = 1` is the single-class independence baseline) and
  returning an enumeration table plus the BIC- or AIC-preferred count. This
  is the recommended first stage: decide how many classes the data support,
  then compare structure at that count. `select_model_ll()` now also accepts
  a vector for `n_classes` (e.g. `1:6`) and performs this selection
  internally, returning the enumeration as `n_classes_table`.

* Parallelisation. `ll_equivalence_test()`, `select_model_ll()`, and
  `select_n_classes()` gain an `mc.cores` argument that spreads the bootstrap
  refits (or the class-count fits) across cores with `parallel::mclapply()`
  on non-Windows platforms. Because every replicate and every fit is seeded
  independently, parallel runs are bit-identical to serial ones; a full
  `select_model_ll()` call runs about 5x faster on six cores.

* An optimised RcppArmadillo core for the latent-class EM (E-step, exact
  weighted-isotonic and Dykstra constrained M-steps), roughly 5x faster than
  the reference R implementation, which is retained and selectable via
  `use_cpp = FALSE` for validation.

## Statistical corrections

Several errors in the earlier latent-class estimation code were found and
fixed; together they had caused systematic over-selection of the Rasch model
and near-zero recovery of the ordinal models:

* The pool-adjacent-violators routine double-counted weights and did not solve
  the intended isotonic regression. It is rewritten as a correct
  block-merging weighted PAVA.

* The constrained M-steps now perform the exact expected-count-weighted
  isotonic regression that maximises the EM objective (per-item weighted PAVA
  for MON, per-class PAVA for IIO, Dykstra alternating projections for DM),
  replacing an unweighted projection and an optimiser path that could silently
  stall on active constraints.

* Parameter counts corrected: the Rasch model uses `J + 1` (item intercepts
  plus the latent variance, from `mirt`), and the latent-class Rasch model
  uses `2C + J - 2`.

* `g_squared()` now scales expected values by sample size; `lr_test()` returns
  `NA` with an explanatory note for equal-parameter-count (chi-bar-squared)
  comparisons rather than a spurious p-value.

## Notes

* Package relicensed GPL (>= 2) to match the absorbed `ConjointChecks` code.
* The bundled `rasch1000` dataset is a pre-computed `ConjointChecks` result on
  simulated Rasch data, useful for exploring the `checks` class and its
  `plot()` / `summary()` methods.
