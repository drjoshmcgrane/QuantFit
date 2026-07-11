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
#'     double then triple, on an ability-banded score matrix: are the
#'     cancellation axioms of additive conjoint measurement satisfied?
#'   \item \strong{Kara} - the Karabatsos (2018) synthetic-likelihood test
#'     ([KaraChecks()]) on the same banded matrix: do the additive axioms
#'     hold, as measured by Kullback-Leibler departures from additivity?
#' }
#'
#' @details
#' Persons are grouped into `n_bands` ability bands by their Rasch ability
#' estimate (from [fit_rm()]); the CC and Kara routes both operate on the
#' resulting band-by-item matrix, with each band's mean ability as the row
#' metric. This banding is deliberate: applying these axiom checks to raw
#' sum-score groups makes them read genuinely additive (Rasch) data as
#' non-additive - the sum-score-to-ability nonlinearity and sparse extreme
#' groups inject spurious violations - whereas a small number of ability bands
#' with a real ability metric is well calibrated (simulated Rasch data pass;
#' data with dispersed item slopes are flagged).
#'
#' Two limitations are worth stating. First, banding on a unidimensional
#' ability estimate can wash out \emph{multidimensional} departures from
#' additivity, so a low Kara/CC violation rate is evidence about within-scale
#' additivity, not a guarantee of unidimensionality. Second, only the LC route
#' is calibrated in the strict sense (its bootstrap gives a chi-bar-squared
#' p-value); the CC and Kara routes use the fixed thresholds below, which are
#' heuristic - the verdict reports each route's raw statistic so the judgement
#' can be inspected, not just the flag.
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
#' @param triple Also run the triple-cancellation check (4x4), conditional on
#'   double (default TRUE). Requires `n_bands >= 4`.
#' @param kara_S,kara_N_synth Iterations and synthetic datasets for
#'   [KaraChecks()] (defaults 20000, 100).
#' @param B Bootstrap replicates for the LC route (default 99).
#' @param kara_max_viol Maximum proportion of Kara cells that may violate
#'   (KL > 0.01) to judge the data additivity-consistent (default 0.15).
#' @param cc_max_viol Maximum ConjointChecks weighted violation rate to judge
#'   the data additivity-consistent (default 0.10).
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
#' @seealso [select_model_ll()], [ConjointChecks()], [KaraChecks()].
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
                                cc_n_mat = 50, triple = TRUE,
                                kara_S = 20000, kara_N_synth = 100, B = 99,
                                kara_max_viol = 0.15, cc_max_viol = 0.10,
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

  # -- CC route (double, then triple) on the banded matrix ----------------
  if (verbose) cat("[CC]   cancellation checks on", n_bands, "ability bands...\n")
  cc <- if (is.null(band)) {
    list(available = FALSE, msg = "ability banding failed", supports_quant = NA)
  } else tryCatch({
    dbl <- ConjointChecks(band$N, band$n, n.mat = cc_n_mat, check = "double",
                          mc.cores = mc.cores)
    dbl_v <- dbl@means$weighted
    res <- list(available = TRUE, double_violation = dbl_v,
                double_ok = dbl_v <= cc_max_viol)
    if (triple && n_bands >= 4) {
      tri <- ConjointChecks(band$N, band$n, n.mat = cc_n_mat, check = "triple",
                            mc.cores = mc.cores)
      res$triple_violation <- tri@means$weighted
      res$triple_ok <- tri@means$weighted <= cc_max_viol
    }
    # CC supports additivity if double holds and (triple holds where run)
    res$supports_quant <- res$double_ok &&
      (is.null(res$triple_ok) || isTRUE(res$triple_ok))
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
    prop <- kc$n_violations / (nr * nc)
    list(available = TRUE, global_KL = kc$global_KL, n_violations = kc$n_violations,
         n_cells = nr * nc, prop_violations = prop,
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
                 lc = lc, cc = cc, kara = kara, n_bands = n_bands,
                 thresholds = c(kara_max_viol = kara_max_viol,
                                cc_max_viol = cc_max_viol)),
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

  cat(sprintf("[CC]   Cancellation checks (%d ability bands)\n", x$n_bands))
  if (isTRUE(x$cc$available)) {
    cat(sprintf("       double %.1f%%%s   [%s]\n",
                100 * x$cc$double_violation,
                if (!is.null(x$cc$triple_violation))
                  sprintf(", triple %.1f%%", 100 * x$cc$triple_violation) else "",
                yn(x$cc$supports_quant)))
  } else cat("       unavailable:", x$cc$msg, "\n")

  cat("[Kara] Synthetic-likelihood additivity test (banded)\n")
  if (isTRUE(x$kara$available)) {
    cat(sprintf("       global KL %.2f; %d/%d cells violate (%.0f%%)   [%s]\n",
                x$kara$global_KL, x$kara$n_violations, x$kara$n_cells,
                100 * x$kara$prop_violations, yn(x$kara$supports_quant)))
  } else cat("       unavailable:", x$kara$msg, "\n")

  cat("\nLC is bootstrap-calibrated; axiom thresholds are heuristic (Kara <=",
      sprintf("%.0f%% cells,", 100 * x$thresholds["kara_max_viol"]),
      "CC <=", sprintf("%.0f%%).", 100 * x$thresholds["cc_max_viol"]),
      "\nBanding is unidimensional, so the axiom routes can miss multidimensional\ndepartures - read a low violation rate as within-scale additivity.\n")
  invisible(x)
}
