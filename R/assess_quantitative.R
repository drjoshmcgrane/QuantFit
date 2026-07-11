#' Overall judgement on quantitative structure from three methods
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
#'   \item \strong{Kara} - the Karabatsos (2018) synthetic-likelihood test
#'     ([KaraChecks()]) on an ability-banded matrix: do the additive axioms
#'     hold, as measured by Kullback-Leibler departures from additivity?
#' }
#'
#' @details
#' The three routes use the calibration appropriate to each. \strong{CC} runs
#' [cc_bootstrap_null()]: the observed [ConjointChecks()] violation rate is
#' located in a null distribution of rates simulated from the Rasch model
#' fitted to the data, and interval scaling is rejected when the observed rate
#' exceeds the `cc_cutoff` percentile of that null (Student & Read, 2025).
#' Because observed and null data pass through the same sum-score pipeline, the
#' null self-calibrates the baseline violation rate - no ability banding is
#' needed for CC, and there is no fixed violation threshold.
#'
#' \strong{Kara} does need help: applying the Karabatsos KL test to raw
#' sum-score groups reads genuinely additive (Rasch) data as non-additive (the
#' sum-score-to-ability nonlinearity injects spurious violations), so persons
#' are grouped into `n_bands` ability bands by their Rasch ability estimate
#' (from [fit_rm()]) and the KL test runs on that band-by-item matrix with each
#' band's mean ability as the row metric - which is well calibrated (simulated
#' Rasch data pass; data with dispersed item slopes are flagged).
#'
#' The Kara route is genuinely inferential in its own right: Karabatsos's
#' method performs approximate Bayesian inference (synthetic-likelihood
#' importance sampling) and rejects additivity when the per-cell
#' Kullback-Leibler divergence exceeds 0.01, a criterion he validated by
#' simulation. The verdict therefore reports his inferential quantities - the
#' global KL and the per-cell KL distribution - not just a flag. (Karabatsos
#' considered and rejected standardized-residual and credible-interval
#' statistics as, respectively, too liberal and too conservative, so the KL is
#' used here as he intended.)
#'
#' Calibration status of each route: \strong{LC} by bootstrap (a
#' chi-bar-squared p-value); \strong{CC} by the Rasch bootstrap null (a
#' percentile p-value, per Student & Read); \strong{Kara} by Karabatsos's
#' simulation-validated KL > 0.01 criterion. The only remaining heuristic is
#' the aggregation of Kara's per-cell rule into a single verdict, so the
#' verdict reports each route's raw statistics for inspection.
#'
#' Two limitations are worth stating. First, the CC procedure is under-powered
#' below roughly 1000 examinees (Student & Read, 2025): at small \eqn{N} a
#' non-rejection may reflect low power rather than interval scalability, and a
#' note is printed. Second, Kara's ability banding is unidimensional, so it can
#' wash out \emph{multidimensional} departures from additivity - read a low Kara
#' violation rate as within-scale additivity, not a guarantee of
#' unidimensionality.
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
#' @param cc_cutoff Null percentile above which the CC route rejects interval
#'   scaling (default 0.95).
#' @param kara_S,kara_N_synth Iterations and synthetic datasets for
#'   [KaraChecks()] (defaults 20000, 100).
#' @param B Bootstrap replicates for the LC route (default 99).
#' @param kara_max_viol Maximum proportion of Kara cells that may exceed
#'   Karabatsos's KL > 0.01 rejection criterion for the data to be judged
#'   additivity-consistent (default 0.15). The KL > 0.01 cutoff is
#'   Karabatsos's own, established by his simulation study (Rasch KL converges
#'   to 0, non-additive 2PL KL exceeds 0.01); the aggregation into a single
#'   proportion-of-cells tolerance is this package's, since his stated rule
#'   rejects additivity when one or more cells exceed the cutoff. The verdict
#'   also reports the global KL and the per-cell KL distribution
#'   (median, Q3, 90th percentile, max), which are his headline inferential
#'   summaries.
#' @param mc.cores Cores for the parallelisable steps (default 1).
#' @param seed Optional integer seed.
#' @param verbose Print progress (default TRUE).
#' @param ... Passed to [select_model_ll()].
#'
#' @return An object of class `quantverdict`: a list with `verdict`
#'   (character), `support` (0-3 routes supporting a quantitative reading),
#'   per-route sub-lists `lc`, `cc`, `kara`, the `n_bands` used, and the
#'   `thresholds`.
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
#' v <- assess_quantitative(dat, mc.cores = 4, seed = 1)
#' v
#' }
#' @export
assess_quantitative <- function(data, n_classes = 1:6, n_bands = 6L,
                                cc_n_mat = 50, triple = TRUE, cc_B = 100,
                                cc_cutoff = 0.95,
                                kara_S = 20000, kara_N_synth = 100, B = 99,
                                kara_max_viol = 0.15,
                                mc.cores = 1L, seed = NULL, verbose = TRUE,
                                ...) {

  data <- validate_data(data)

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

  # -- ability banding for the axiom routes -------------------------------
  band <- tryCatch({
    th <- rm_scores(suppressWarnings(fit_rm(data, verbose = FALSE)))$theta
    br <- stats::quantile(th, seq(0, 1, length.out = n_bands + 1))
    br[1] <- -Inf; br[length(br)] <- Inf
    grp <- cut(th, br, labels = FALSE)
    J <- ncol(data)
    Nb <- nb <- matrix(0, n_bands, J); ab <- numeric(n_bands)
    for (k in seq_len(n_bands)) {
      idx <- grp == k
      Nb[k, ] <- sum(idx); nb[k, ] <- colSums(data[idx, , drop = FALSE])
      ab[k] <- mean(th[idx])
    }
    ord <- order(colSums(data))                 # items by difficulty
    list(N = Nb[, ord], n = nb[, ord], ability = ab)
  }, error = function(e) NULL)

  # -- CC route: bootstrapped null (Student & Read 2025) on the raw data ---
  # The Rasch-simulated null self-calibrates the sum-score pipeline, so this
  # runs on raw sum-score groups (not the bands) and needs no fixed threshold.
  if (verbose) cat("[CC]   bootstrapped cancellation checks (Student & Read)...\n")
  cc <- tryCatch({
    dbl <- cc_bootstrap_null(data, check = "double", n.mat = cc_n_mat,
                             B = cc_B, cutoff = cc_cutoff, mc.cores = mc.cores,
                             seed = seed, verbose = FALSE)
    res <- list(available = TRUE, double_rate = dbl$observed,
                double_null_mean = mean(dbl$null),
                double_percentile = dbl$percentile, double_p = dbl$p_value,
                double_reject = dbl$reject)
    if (triple) {
      tri <- cc_bootstrap_null(data, check = "triple", n.mat = cc_n_mat,
                               B = cc_B, cutoff = cc_cutoff, mc.cores = mc.cores,
                               seed = if (!is.null(seed)) seed + 1L else NULL,
                               verbose = FALSE)
      res$triple_rate <- tri$observed; res$triple_percentile <- tri$percentile
      res$triple_p <- tri$p_value; res$triple_reject <- tri$reject
    }
    # CC supports additivity when interval scaling is NOT rejected
    res$supports_quant <- !isTRUE(res$double_reject) &&
      (is.null(res$triple_reject) || !isTRUE(res$triple_reject))
    res
  }, error = function(e) list(available = FALSE, msg = conditionMessage(e),
                              supports_quant = NA))

  # -- Kara route on the same banded matrix -------------------------------
  if (verbose) cat("[Kara] synthetic-likelihood additivity test on bands...\n")
  kara <- if (is.null(band)) {
    list(available = FALSE, msg = "ability banding failed", supports_quant = NA)
  } else tryCatch({
    nr <- nrow(band$N); nc <- ncol(band$N)
    kc <- KaraChecks(as.vector(band$N), as.vector(band$n), S = kara_S,
                     N_synth = kara_N_synth, mc.cores = mc.cores, verbose = FALSE,
                     testscore = rep(band$ability, nc), item = rep(1:nc, each = nr))
    # Karabatsos's inferential summaries: the global KL and the distribution
    # of per-cell KL, with his KL > 0.01 rejection criterion.
    kl <- as.vector(kc$KL)
    qs <- stats::quantile(kl, c(0.5, 0.75, 0.90, 1), names = FALSE)
    prop <- kc$n_violations / (nr * nc)
    list(available = TRUE, global_KL = kc$global_KL, n_violations = kc$n_violations,
         n_cells = nr * nc, prop_violations = prop,
         kl_median = qs[1], kl_q3 = qs[2], kl_p90 = qs[3], kl_max = qs[4],
         supports_quant = prop <= kara_max_viol)
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
  } else if (isTRUE(lc$supports_quant) && axiom_avail > 0 && axiom_ok == axiom_avail) {
    "QUANTITATIVE (well supported) - a quantitative model is selected and the conjoint axioms are satisfied; the three routes converge."
  } else if (isTRUE(lc$supports_quant) && axiom_ok > 0) {
    "QUANTITATIVE (qualified) - a quantitative model is selected, but one axiom check shows departures from additivity; interpret with caution."
  } else if (isTRUE(lc$supports_quant)) {
    "QUANTITATIVE MODEL, WEAK AXIOM SUPPORT - a quantitative model fits best, but the conjoint axioms are not satisfied; the interval interpretation is not tenable on the axiom evidence."
  } else if (identical(lc_scale, "ordinal") && axiom_avail > 0 && axiom_ok == axiom_avail) {
    "BORDERLINE - model selection favours an ordinal structure, yet the conjoint axioms hold; possibly quantitative with weak class separation."
  } else if (identical(lc_scale, "ordinal")) {
    "ORDINAL - model selection favours an ordinal structure and the axiom checks do not support additivity."
  } else {
    "INCONCLUSIVE - the methods do not point clearly in one direction."
  }

  structure(list(verdict = verdict, support = support, n_available = n_avail,
                 lc = lc, cc = cc, kara = kara, n_bands = n_bands, N = nrow(data),
                 thresholds = c(kara_max_viol = kara_max_viol,
                                cc_cutoff = cc_cutoff)),
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
    cat(sprintf("       double: rate %.3f, %.0f%%ile of null (p = %.3f)%s   [%s]\n",
                x$cc$double_rate, 100 * x$cc$double_percentile, x$cc$double_p,
                if (!is.null(x$cc$triple_percentile))
                  sprintf("; triple %.0f%%ile (p = %.3f)",
                          100 * x$cc$triple_percentile, x$cc$triple_p) else "",
                yn(x$cc$supports_quant)))
  } else cat("       unavailable:", x$cc$msg, "\n")

  cat("[Kara] Karabatsos synthetic-likelihood additivity test (banded)\n")
  if (isTRUE(x$kara$available)) {
    cat(sprintf("       global KL %.2f; per-cell KL median %.3f, Q3 %.3f, 90%% %.3f, max %.3f\n",
                x$kara$global_KL, x$kara$kl_median, x$kara$kl_q3,
                x$kara$kl_p90, x$kara$kl_max))
    cat(sprintf("       %d/%d cells exceed KL > 0.01 (%.0f%%)   [%s]\n",
                x$kara$n_violations, x$kara$n_cells,
                100 * x$kara$prop_violations, yn(x$kara$supports_quant)))
  } else cat("       unavailable:", x$kara$msg, "\n")

  cat("\nLC calibrated by bootstrap; CC by a Rasch bootstrap null",
      sprintf("(reject > %.0f%%ile);", 100 * x$thresholds["cc_cutoff"]),
      "Kara by Karabatsos's KL > 0.01 criterion.")
  if (!is.null(x$N) && x$N < 1000) {
    cat(sprintf("\nNote: N = %d; the CC procedure is under-powered below ~1000 (Student & Read 2025).", x$N))
  }
  cat("\nAxiom banding (Kara) is unidimensional, so it can miss multidimensional",
      "\ndepartures - read a low violation rate as within-scale additivity.\n")
  invisible(x)
}
