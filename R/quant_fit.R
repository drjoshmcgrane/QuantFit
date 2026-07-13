#' Assess whether item response data support a quantitative interpretation
#'
#' Combines the three routes QuantFit offers to the question *do these data
#' support a quantitative (interval / additive) interpretation?* into a single
#' triangulated verdict:
#'
#' \enumerate{
#'   \item \strong{LC} - latent-structure model selection
#'     ([select_model_ll()]) on the raw person-by-item data: does a
#'     quantitative model (LCR or RM) win over the classificatory and ordinal
#'     alternatives?
#'   \item \strong{CC} - Bayesian cancellation checks ([ConjointChecks()]),
#'     double (and optionally triple), calibrated against a Rasch bootstrap
#'     null via [cc_bootstrap_null()] (Student & Read, 2025): is the observed
#'     violation rate higher than expected under interval scaling?
#'   \item \strong{Kara} - the Karabatsos (2018) synthetic-likelihood global
#'     KL statistic ([KaraChecks()]) on an ability-banded matrix, calibrated
#'     against its own Rasch bootstrap null ([kara_bootstrap_null()]): is the
#'     KL departure from additivity higher than expected under interval
#'     scaling?
#' }
#'
#' @details
#' All three routes are calibrated the same way - a per-dataset parametric
#' bootstrap, with the statistic located in a null distribution and interval
#' scaling / the constrained model rejected above the `cc_cutoff` percentile.
#' \strong{LC} bootstraps the likelihood-ratio statistic against its
#' chi-bar-squared null (simulating from the fitted constrained model);
#' \strong{CC} bootstraps the [ConjointChecks()] violation rate against a Rasch
#' null (Student & Read, 2025); \strong{Kara} bootstraps the Karabatsos global
#' KL against a Rasch null. This makes the three routes statistically
#' consistent and their p-values comparable, and removes every fixed threshold.
#'
#' \strong{CC} runs on raw sum-score groups: because observed and null data
#' share the pipeline, the null self-calibrates the baseline, so no ability
#' banding is needed. \strong{Kara} does need banding - the KL test on raw
#' sum-score groups reads genuinely additive (Rasch) data as non-additive (a
#' sum-score-to-ability nonlinearity) - so persons are grouped into `n_bands`
#' ability bands by their Rasch ability estimate, with each band's mean ability
#' as the row metric. Its bootstrap null passes each simulated dataset through
#' the same banding. The verdict still reports Karabatsos's descriptive
#' summaries (global KL and the per-cell KL distribution) alongside the
#' bootstrap percentile.
#'
#' The only heuristic that remains is the ordinal / classificatory /
#' quantitative \emph{labelling} of the LC-selected model; every route's
#' accept/reject decision is a bootstrap percentile.
#'
#' Two limitations are worth stating. First, the axiom procedures are
#' under-powered below roughly 1000 examinees (Student & Read, 2025): at small
#' \eqn{N} a non-rejection may reflect low power rather than interval
#' scalability, and a note is printed. Second, Kara's ability banding is
#' unidimensional, so it can wash out \emph{multidimensional} departures from
#' additivity - read a low Kara statistic as within-scale additivity, not a
#' guarantee of unidimensionality.
#'
#' A quantitative reading is best supported when the model route selects
#' LCR/RM \emph{and} both axiom routes are additivity-consistent.
#'
#' @param data Binary response matrix (persons x items).
#' @param n_classes Class-count range for the LC route (default `1:6`).
#' @param n_bands Number of ability bands for the axiom routes (default 6;
#'   more bands re-introduce the sparse-extreme problem, fewer lose
#'   resolution).
#' @param cc_n_mat Submatrices sampled per [ConjointChecks()] run (default 50).
#' @param triple Also calibrate the triple-cancellation check against its own
#'   Rasch bootstrap null (default TRUE).
#' @param cc_B Number of Rasch-simulated null datasets for the CC route
#'   (default 100, passed to [cc_bootstrap_null()]).
#' @param cc_cutoff Null percentile above which the CC \emph{and} Kara routes
#'   reject interval scaling (default 0.95); it is passed as the `cutoff` to
#'   both [cc_bootstrap_null()] and [kara_bootstrap_null()].
#' @param kara_S,kara_N_synth Iterations and synthetic datasets for each
#'   [KaraChecks()] run (defaults 10000, 100), used identically for the
#'   observed statistic and every bootstrap replicate.
#' @param kara_B Number of Rasch-simulated null datasets for the Kara route
#'   (default 50; each is a full [KaraChecks()] run, so this is the most
#'   expensive step - raise it for a final analysis).
#' @param B Bootstrap replicates for the LC route (default 99).
#' @param mc.cores Cores for the parallelisable steps (default 1).
#' @param seed Optional integer seed.
#' @param verbose Print progress (default TRUE).
#' @param ... Passed to [select_model_ll()].
#'
#' @return An object of class `quantverdict`: a list with `verdict`
#'   (character), `support` (routes supporting a quantitative reading),
#'   `n_available` (routes that returned a result), per-route sub-lists `lc`,
#'   `cc`, `kara`, `n_bands`, `N` (sample size), and `thresholds`.
#'
#' @references
#' Karabatsos, G. (2018). On Bayesian testing of additive conjoint measurement
#' axioms using synthetic likelihood. \emph{Psychometrika}, 83(2), 321-332.
#'
#' Domingue, B. (2014). Evaluating the equal-interval hypothesis with test
#' score scales. \emph{Psychometrika}, 79(1), 1-19.
#'
#' Student, S. R., & Read, W. S. (2025). Applying Bayesian checks of
#' cancellation axioms for interval scaling in limited samples.
#' \emph{Behavior Research Methods}, 57, 305.
#'
#' @seealso [select_model_ll()], [cc_bootstrap_null()], [ConjointChecks()],
#'   [KaraChecks()].
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' theta <- rnorm(1500); beta <- seq(-2, 2, length.out = 12)
#' dat <- matrix(rbinom(1500 * 12, 1, plogis(outer(theta, beta, "-"))), 1500, 12)
#' v <- quant_fit(dat, mc.cores = 4, seed = 1)
#' v
#' }
#' @export
quant_fit <- function(data, n_classes = 1:6, n_bands = 6L,
                                cc_n_mat = 50, triple = TRUE, cc_B = 100,
                                cc_cutoff = 0.95,
                                kara_S = 10000, kara_N_synth = 100,
                                kara_B = 50, B = 99,
                                mc.cores = 1L, seed = NULL, verbose = TRUE,
                                ...) {

  data <- validate_data_any(data)

  # -- LC route (raw person x item data) ----------------------------------
  if (verbose) cat("[LC]   latent-structure model selection...\n")
  lc <- tryCatch({
    sel <- select_model_ll(data, n_classes = n_classes, B = B,
                           mc.cores = mc.cores, seed = seed, verbose = FALSE,
                           ...)
    list(available = TRUE, selected = sel$selected,
         interpretation = sel$interpretation, n_classes = sel$n_classes,
         scale = if (sel$selected %in% c("LCR", "RM")) "quantitative"
                 else if (sel$selected == "UN") "classificatory" else "ordinal",
         supports_quant = sel$selected %in% c("LCR", "RM"))
  }, error = function(e) list(available = FALSE, msg = conditionMessage(e),
                              supports_quant = NA))

  # -- CC route: bootstrapped null (Student & Read 2025) on the raw data ---
  # The Rasch-simulated null self-calibrates the sum-score pipeline, so this
  # runs on raw sum-score groups (not the bands) and needs no fixed threshold.
  # Levels run SEQUENTIALLY in their logical order (single -> double -> triple),
  # stopping at the first rejection: deeper cancellation conditions are only
  # distinct requirements once the shallower ones hold, so the first failing
  # level is where additivity breaks (`attribution`).
  if (verbose) cat("[CC]   bootstrapped cancellation checks (Student & Read)...\n")
  cc <- tryCatch({
    hier <- cc_bootstrap_hierarchy(
      data, levels = if (triple) c("single", "double", "triple")
                     else c("single", "double"),
      n.mat = cc_n_mat, B = cc_B, cutoff = cc_cutoff, mc.cores = mc.cores,
      seed = seed, verbose = FALSE)
    res <- list(available = TRUE, hierarchy = hier,
                attribution = hier$attribution)
    sng <- hier$levels$single
    if (!is.null(sng)) {
      res$single_rate <- sng$observed; res$single_percentile <- sng$percentile
      res$single_p <- sng$p_value; res$single_reject <- sng$reject
    }
    dbl <- hier$levels$double
    if (!is.null(dbl)) {
      res$double_rate <- dbl$observed
      res$double_null_mean <- mean(dbl$null)
      res$double_percentile <- dbl$percentile; res$double_p <- dbl$p_value
      res$double_reject <- dbl$reject
    }
    tri <- hier$levels$triple
    if (!is.null(tri)) {
      res$triple_rate <- tri$observed; res$triple_percentile <- tri$percentile
      res$triple_p <- tri$p_value; res$triple_reject <- tri$reject
    }
    # CC supports additivity when no tested level rejects
    res$supports_quant <- hier$supports_quant
    res
  }, error = function(e) list(available = FALSE, msg = conditionMessage(e),
                              supports_quant = NA))

  # -- Kara route: bootstrapped null for the global KL --------------------
  # Same per-dataset parametric-bootstrap logic as CC, applied to Karabatsos's
  # global KL on the ability-banded matrix (see kara_bootstrap_null()).
  if (verbose) cat("[Kara] bootstrapped Karabatsos KL test (banded)...\n")
  kara <- tryCatch({
    kn <- kara_bootstrap_null(data, n_bands = n_bands, B = kara_B,
                              cutoff = cc_cutoff, S = kara_S,
                              N_synth = kara_N_synth, mc.cores = mc.cores,
                              seed = if (!is.null(seed)) seed + 2L else NULL,
                              verbose = FALSE)
    list(available = TRUE, global_KL = kn$observed,
         null_mean = mean(kn$null), percentile = kn$percentile,
         p_value = kn$p_value, reject = kn$reject,
         kl_median = kn$kl_median, kl_q3 = kn$kl_q3, kl_p90 = kn$kl_p90,
         kl_max = kn$kl_max, supports_quant = !kn$reject)
  }, error = function(e) list(available = FALSE, msg = conditionMessage(e),
                              supports_quant = NA))

  # -- synthesise ---------------------------------------------------------
  votes <- c(lc$supports_quant, cc$supports_quant, kara$supports_quant)
  support <- sum(votes, na.rm = TRUE)
  n_avail <- sum(!is.na(votes))
  axiom_ok <- sum(c(cc$supports_quant, kara$supports_quant), na.rm = TRUE)
  axiom_avail <- sum(!is.na(c(cc$supports_quant, kara$supports_quant)))
  lc_scale <- if (isTRUE(lc$available)) lc$scale else NA_character_

  verdict <- if (n_avail == 0) {
    "INCONCLUSIVE - no method returned a usable result."
  } else if (identical(lc_scale, "classificatory")) {
    "CLASSIFICATORY - model selection favours unstructured latent classes; a quantitative interpretation is not supported."
  } else if (isTRUE(lc$supports_quant) && axiom_avail == 0) {
    "QUANTITATIVE MODEL, AXIOM CHECKS UNAVAILABLE - a quantitative model is selected, but neither conjoint-axiom route could be computed, so the additivity evidence is missing (not failed)."
  } else if (isTRUE(lc$supports_quant) && axiom_ok == axiom_avail) {
    "QUANTITATIVE (well supported) - a quantitative model is selected and the conjoint axioms are satisfied; the routes converge."
  } else if (isTRUE(lc$supports_quant) && axiom_ok > 0) {
    "QUANTITATIVE (qualified) - a quantitative model is selected, but one axiom check shows departures from additivity; interpret with caution."
  } else if (isTRUE(lc$supports_quant)) {
    "QUANTITATIVE MODEL, WEAK AXIOM SUPPORT - a quantitative model fits best, but the conjoint axioms are not satisfied; the interval interpretation is not tenable on the axiom evidence."
  } else if (identical(lc_scale, "ordinal") && axiom_avail == 0) {
    "ORDINAL - model selection favours an ordinal structure (the conjoint-axiom routes could not be computed)."
  } else if (identical(lc_scale, "ordinal") && axiom_ok == axiom_avail) {
    "BORDERLINE - model selection favours an ordinal structure, yet the conjoint axioms hold; possibly quantitative with weak class separation."
  } else if (identical(lc_scale, "ordinal")) {
    "ORDINAL - model selection favours an ordinal structure and the axiom checks do not support additivity."
  } else {
    "INCONCLUSIVE - the methods do not point clearly in one direction."
  }

  structure(list(verdict = verdict, support = support, n_available = n_avail,
                 lc = lc, cc = cc, kara = kara, n_bands = n_bands, N = nrow(data),
                 thresholds = c(cc_cutoff = cc_cutoff)),
            class = "quantverdict")
}

#' Print method for quantverdict objects
#'
#' @param x A quantverdict object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#' @export
print.quantverdict <- function(x, ...) {
  cat("\nQuantitative structure - triangulated judgement\n")
  cat("===============================================\n")
  cat(strwrap(x$verdict, width = 68), sep = "\n")
  cat(sprintf("\n\nRoutes supporting a quantitative reading: %d of %d\n\n",
              x$support, x$n_available))

  yn <- function(f) if (is.na(f)) "  ?  " else if (f) " yes " else "  no "

  cat("[LC]   Latent-structure model selection (raw data)\n")
  if (isTRUE(x$lc$available)) {
    cat(sprintf("       selected %s (%d classes) - %s   [%s]\n",
                x$lc$selected, x$lc$n_classes, x$lc$scale, yn(x$lc$supports_quant)))
  } else cat("       unavailable:", x$lc$msg, "\n")

  cat("[CC]   Cancellation checks vs Rasch bootstrap null (Student & Read)\n")
  if (isTRUE(x$cc$available)) {
    lvls <- character(0)
    if (!is.null(x$cc$single_percentile))
      lvls <- c(lvls, sprintf("single %.0f%%ile (p = %.3f)",
                              100 * x$cc$single_percentile, x$cc$single_p))
    if (!is.null(x$cc$double_percentile))
      lvls <- c(lvls, sprintf("double %.0f%%ile (p = %.3f)",
                              100 * x$cc$double_percentile, x$cc$double_p))
    if (!is.null(x$cc$triple_percentile))
      lvls <- c(lvls, sprintf("triple %.0f%%ile (p = %.3f)",
                              100 * x$cc$triple_percentile, x$cc$triple_p))
    cat(sprintf("       %s   [%s]\n", paste(lvls, collapse = "; "),
                yn(x$cc$supports_quant)))
    if (!identical(x$cc$attribution, "none"))
      cat(sprintf("       additivity fails at the %s-cancellation level\n",
                  x$cc$attribution))
  } else cat("       unavailable:", x$cc$msg, "\n")

  cat("[Kara] Karabatsos KL vs Rasch bootstrap null (banded)\n")
  if (isTRUE(x$kara$available)) {
    cat(sprintf("       global KL %.3f, %.0f%%ile of null (p = %.3f)   [%s]\n",
                x$kara$global_KL, 100 * x$kara$percentile, x$kara$p_value,
                yn(x$kara$supports_quant)))
    cat(sprintf("       per-cell KL median %.3f, Q3 %.3f, 90%% %.3f, max %.3f\n",
                x$kara$kl_median, x$kara$kl_q3, x$kara$kl_p90, x$kara$kl_max))
  } else cat("       unavailable:", x$kara$msg, "\n")

  cat(sprintf("\nAll three routes are bootstrap-calibrated (reject > %.0f%%ile of the null):",
              100 * x$thresholds["cc_cutoff"]),
      "\nLC by a chi-bar-squared LR bootstrap, CC and Kara by Rasch bootstrap nulls.")
  if (!is.null(x$N) && x$N < 1000) {
    cat(sprintf("\nNote: N = %d; the axiom procedures are under-powered below ~1000 (Student & Read 2025).", x$N))
  }
  cat("\nAxiom banding (Kara) is unidimensional, so it can miss multidimensional",
      "\ndepartures - read a low violation rate as within-scale additivity.\n")
  invisible(x)
}

#' @rdname quant_fit
#' @details `assess_quantitative()` is a deprecated alias for `quant_fit()`.
#' @export
assess_quantitative <- function(...) {
  .Deprecated("quant_fit")
  quant_fit(...)
}
