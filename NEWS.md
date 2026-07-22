# QuantFit 0.3.2

* `select_model_ll()` now defaults to the complete adjacent-edge
  Torres Irribarra–Diakow lattice. MON and IIO are tested from UN, and DM must
  be supported on both parent paths before the published DM→LCR→RM sequence is
  entered.
* The LCR-vs-DM scale bridge is fit at the fixed minimum support-point count
  required by the paper's LCR/Rasch equivalence result (`ceiling((J + 1) / 2)`
  for dichotomous items). Grain is profiled separately for the calibrated
  RM-vs-LCR comparison, so a coarse selected LCR cannot veto true Rasch data.
  Observed and bootstrap fits use equal multi-start effort by default.
* When an ordinal class-count range is supplied, UN-BIC enumeration is now
  repeated inside every ordinal-edge null replicate. This calibrates the
  automated two-stage statistic rather than treating a data-selected class
  count as fixed.
* The former severity override is now optional and off by default. Public
  results expose LCR-vs-DM and RM-vs-LCR evidence separately.
* LCR multi-start fitting now uses genuinely dispersed mixture starts after the
  stable score start, fixing an extreme-class local optimum that materially
  depressed the LCR likelihood on correctly generated data.

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

* **Polytomous data support across all three routes.** Items scored with
  ordered integer categories `0, 1, ..., m` are now handled natively, with no
  loss of information. The latent-structure models (`fit_un()`, `fit_mon()`,
  `fit_iio()`, `fit_dm()`, `fit_lcr()`, `fit_rm()`) fit polytomous data directly
  via a multinomial EM engine: each item-by-class combination carries a full
  category-probability vector, class monotonicity becomes stochastic ordering of
  the category distributions across classes, invariant item ordering becomes
  non-crossing expected item-score curves, and the quantitative models use the
  partial credit model (latent-class PCM for `LCR`; `mirt`'s native PCM for
  `RM`). `select_model_ll()` and `quant_fit()` dispatch on the data automatically.
  For the conjoint cancellation routes, `prepare_polytomous()` recodes each item
  into its adjacent-category ("partial credit") dichotomisations and conditions
  on the total score before `ConjointChecks()` / `KaraChecks()`; the CC and Kara
  bootstrap nulls simulate from a fitted partial credit model. All six fits, the
  bootstrap tests, and the triangulated verdict work on dichotomous and
  polytomous data alike, choosing the representation appropriate to each route.
  The polytomous engine has a compiled RcppArmadillo core (multinomial E-step,
  expected counts, partial-credit probabilities); the ordered-model M-steps run
  NLopt's SLSQP directly from C++ (via `nloptrAPI.h`), with class monotonicity
  solved per item and an exact shortcut for the coupled IIO / DM cases, so
  constrained fits on constraint-consistent data are orders of magnitude faster.

* **The Rasch / partial-credit model is now estimated by the package's own
  engine, and the `mirt` dependency has been dropped** (moved to Suggests, used
  only to cross-check in the test suite). `fit_rm()` fits by marginal maximum
  likelihood over Gauss-Hermite quadrature, reusing the same EM engine as the
  latent-class models (the quadrature nodes play the role of classes); it
  handles dichotomous (Rasch) and polytomous (partial credit) data and
  reproduces `mirt`'s log-likelihood, item parameters, and latent variance to
  within rounding. `rm_scores()` (EAP / MAP / ML / WLE), `rm_itemfit()`,
  `rm_personfit()`, `rm_item_info()`, and the Rasch standard errors are all
  reimplemented against this engine.

* The triangulated verdict function is now `quant_fit()` - the package's
  namesake flagship function. `assess_quantitative()` is kept as a deprecated
  alias.

* `kara_bootstrap_null()` - applies the same Rasch bootstrap-null calibration
  to Karabatsos's global KL statistic that `cc_bootstrap_null()` applies to the
  cancellation-check violation rate, giving a per-dataset percentile p-value in
  place of a fixed KL > 0.01 cutoff. With this, all three routes of
  `quant_fit()` are calibrated the same way - a per-dataset
  parametric bootstrap with an accept/reject percentile - so their evidence is
  statistically consistent and comparable. The final LC step (RM vs LCR) is
  likewise now bootstrap-calibrated (the BIC difference against a Rasch null)
  rather than a raw BIC comparison, removing the last non-bootstrap criterion.

* **Missing data (MAR) is now supported across all three routes.** The
  latent-structure models and the Rasch/partial-credit model use a masked
  likelihood (each person contributes only their observed cells; binary data
  with missing responses runs as the m = 1 case of the masked multinomial
  engine; masked `fit_rm()` reproduces `mirt`'s native missing-data
  log-likelihood). The conjoint routes handle missingness by observation
  weighting - every cell of the conditioned matrix counts only the
  respondents who answered that item, exactly the same weighting that
  handles the structural out-of-play cells of the adjacent-category
  polytomous recoding, with the conditioning score taken over observed
  responses. Bootstrap null replicates inherit the observed missingness
  pattern by rank-matched assignment (the replicate person with the r-th
  smallest simulated total receives the mask of the r-th smallest observed
  total), which preserves ability-missingness dependence under the marginal
  redraw, so pipeline effects cancel in the calibration. Validated:
  verdicts concordant across complete, 15% MCAR, and MAR-by-ability versions
  of the same additive and non-additive datasets. ML/WLE person scores
  require complete data (score sufficiency) and say so; EAP/MAP work under
  missingness. Nonignorable (MNAR) missingness and heavily structured
  designs (booklets/matrix sampling) remain outside scope: with such
  designs, group-by-observed-total conditioning mixes very different
  patterns and results should be read with care.

* `cc_bootstrap_hierarchy()` - runs the calibrated cancellation checks in
  their logical order (single, then double, then triple), stopping at the
  first rejection. The cancellation axioms are hierarchical - double
  cancellation is only a distinct requirement once the single-cancellation
  orderings hold, and triple presupposes double - so the sequential procedure
  adds *attribution*: it reports the level at which additivity breaks, and
  deeper levels are not run (they would be uninformative given a shallower
  failure). `quant_fit()`'s CC route now uses this and reports `attribution`.

* Bootstrap nulls can now draw abilities from the *empirical* latent
  distribution (`latent = "empirical"`, the new default for
  `cc_bootstrap_null()`, `kara_bootstrap_null()`, and
  `cc_bootstrap_hierarchy()`): a Bock-Aitkin empirical-histogram refit
  re-estimates the latent distribution jointly with the item parameters, so
  null replicates reproduce the observed ability distribution (bimodal,
  skewed, censored samples) and its score-group structure. Additive conjoint
  structure is distribution-free, so the latent shape is a nuisance parameter
  here; making the null match it ensures rejections are attributable to
  non-additivity rather than population shape. (The RM-vs-LCR step in
  `select_model_ll()` deliberately keeps the normal density: there the shape
  of the latent distribution - continuous versus discrete - is the hypothesis
  being tested.) A type-I study on strongly bimodal additive data found the
  CC route well calibrated under both settings (sum-score conditioning is
  protected by Rasch sufficiency); the empirical default matters most for the
  ability-banded Kara route, which has no such protection.

* `select_model_ll()` rejections now pass an estimated-power / effect-size
  check before a constrained model is demoted: the LR distribution is also
  simulated under the fitted general model, and when it does not separate
  from the constrained null the data cannot distinguish the two models - the
  estimated departure is negligible - so parsimony keeps the constrained
  model. This eliminates the dominant recovery-audit error mode (true models
  demoted by unlucky upper-tail draws of their own null) while leaving
  genuine violations, where the distributions separate decisively, untouched.

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

* `quant_fit()` - a single triangulated judgement on whether data
  support a quantitative (interval / additive) interpretation, combining all
  three of the package's routes, each calibrated the same way: latent-structure
  model selection ([select_model_ll()], bootstrap chi-bar-squared) on the raw
  data; the cancellation checks calibrated against a Rasch bootstrap null
  ([cc_bootstrap_null()], per Student & Read 2025); and the Karabatsos
  synthetic-likelihood test calibrated against a Rasch bootstrap null
  ([kara_bootstrap_null()], the global KL statistic's percentile within the
  null rather than a fixed KL > 0.01 cutoff) on an ability-banded matrix. The
  Kara banding is deliberate - the KL test reads additive Rasch data as
  non-additive on raw sum-score groups, whereas ability bands with a real
  ability metric are well calibrated; the CC route needs no banding because its
  bootstrap null self-calibrates the pipeline. All three bootstrap nulls are
  marginal parametric bootstraps: each replicate redraws abilities
  theta ~ N(0, sigma^2) from the fitted Rasch latent variance rather than
  reusing fixed person estimates, and the percentile p-value carries the
  (1 + #{null >= observed}) / (B + 1) continuity correction. The verdict names
  the supporting routes and reports each route's statistics. It states its
  limits: the CC route is under-powered below ~1000 examinees, and the
  (unidimensional) Kara banding can miss multidimensional departures.


* `select_model_ll()` and `ll_equivalence_test()` — a statistically calibrated
  model selection procedure. Because UN and the ordinal models (MON/IIO/DM)
  have identical parameter counts, information criteria cannot distinguish
  them; these functions instead use a parametric bootstrap of the
  likelihood-ratio statistic, whose asymptotic null is chi-bar-squared, to
  test each constrained model against its less-restricted parent. The ordinal
  layer identifies the most restrictive supported model among UN/MON/IIO/DM,
  and a single quantitative gate then tests the parametric latent-class Rasch
  model directly against the unconstrained model (LCR vs UN): a one-step gate
  rather than a sequential DM-then-LCR path, so genuinely quantitative data are
  not lost to the ordinal layer by the compounded false-rejection rate of two
  tests. This gate uses a separate `alpha_quant` (default 0.05, the
  conventional level; true quantitative models are protected from chance
  rejections by the estimated-power check below, so the conventional level
  does not carry the compounded false-demotion cost of a plain sequential
  procedure) - the
  quantitative model is only demoted on strong evidence. A `method` argument
  offers two ways to test the ordinal layer: `"joint"` (default) tests the
  doubly-monotone model against UN directly, while `"lattice"` tests each
  constraint edge separately. Six-model recovery audits (`simulate_responses()`,
  N = 1500, J = 8, K = 30 per model, B = 99, with the estimated-power check)
  quantify the `alpha_quant` trade-off: at the default 0.05, false
  quantitative claims (ordinal data selected as LCR/RM) occur in 1.1% of
  ordinal datasets and quantitative-scale recovery is 88%; at 0.01, false
  quantitative claims rise to 4.4% while quantitative recovery reaches 97%.
  Nominal recovery is 100% and overall scale-type accuracy ~96% under both.
  The conventional default is the conservative choice about quantitativeness
  claims; set `alpha_quant = 0.01` when maximising quantitative recovery
  matters more than guarding against false quantitative verdicts.

* `simulate_responses()` - generate dichotomous or polytomous data from any of
  the six models (UN, MON, IIO, DM, LCR, RM) for validation and power studies;
  reproduces the generators of Torres Irribarra & Diakow for dichotomous data
  and generalises them to ordered categories. `inst/validation/selection_audit.R`
  uses it to reproduce the recovery study above.

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
