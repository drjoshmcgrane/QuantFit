#' Bootstrapped null distribution for conjoint-cancellation violation rates
#'
#' Calibrates the [ConjointChecks()] weighted violation rate against a null
#' distribution simulated from the Rasch model fitted to the data, following
#' Student & Read (2025). The violation rate on its own is not interpretable -
#' even genuinely interval-scalable (Rasch) data produce violations at a rate
#' that depends on sample size, test length, and item parameters - so the
#' observed rate is located within a null distribution of rates from
#' Rasch-simulated data. Its percentile is treated as a p-value: an observed
#' rate above (say) the 95th percentile of the null is evidence that the data
#' are *not* compatible with an interval scale.
#'
#' @details
#' The procedure follows Student & Read (2025):
#' \enumerate{
#'   \item Compute the observed weighted mean proportion of cancellation-axiom
#'     violations with [ConjointChecks()].
#'   \item Fit the Rasch model to the data (via [fit_rm()]), giving the item
#'     difficulties and the latent variance.
#'   \item Simulate `B` datasets of the same size under the fitted *marginal*
#'     Rasch model - abilities are redrawn from N(0, sigma^2) afresh for each
#'     replicate - and compute the violation rate for each: the null
#'     distribution of rates *expected under interval scaling*.
#'   \item Locate the observed rate in that null distribution; its percentile
#'     is the (one-sided) evidence against interval scaling.
#' }
#' This uses a marginal parametric bootstrap (redrawing abilities) rather than
#' Student & Read's fixed plug-in of person-parameter estimates; the marginal
#' version includes the sampling variability of the latent distribution and is
#' the more defensible parametric bootstrap of the Rasch model. With the
#' default `latent = "empirical"`, abilities are redrawn from the Bock-Aitkin
#' empirical-histogram estimate of the latent distribution, which recovers the
#' distribution-free character of Student & Read's plug-in (their fixed person
#' estimates carry the observed ability distribution implicitly) while also
#' realising the error correction they propose as a future direction: their
#' plugged-in theta estimates carry measurement error that inflates the
#' person-distribution variance, and the joint EM re-estimation of the latent
#' weights with the item parameters deconvolves that error rather than
#' shrinking it ad hoc.
#' Because the observed and simulated datasets pass through the identical
#' sum-score [PrepareChecks()] + [ConjointChecks()] pipeline, any baseline
#' violation rate the pipeline induces on additive data (including the
#' sum-score-to-ability nonlinearity) is present in both and cancels - the
#' null is self-calibrating, so no ability banding is needed here.
#'
#' \strong{Power depends on sample size.} Student & Read find the procedure
#' well powered at \eqn{N = 1000} but not at \eqn{N = 250}; a warning is
#' issued for small samples. They recommend confirming, by a companion
#' simulation at your own \eqn{N} and parameter distribution, that the
#' procedure can detect the violations you care about before trusting a
#' non-rejection.
#'
#' @param data Binary response matrix (persons x items).
#' @param check Cancellation axiom passed to [ConjointChecks()]: `"double"`
#'   (default, as studied by Student & Read), `"single"`, or `"triple"`.
#' @param n.mat Number of submatrices sampled per [ConjointChecks()] run
#'   (default 50).
#' @param B Number of Rasch-simulated null datasets (default 100).
#' @param cutoff Percentile of the null above which interval scaling is
#'   rejected (default 0.95).
#' @param alpha Significance level: rejection is `p_value <= alpha` (the
#'   corrected Monte Carlo p, (1 + count of null draws at or above the observed)/(B+1)). Note B bounds the
#'   smallest attainable p at 1/(B+1): B >= 19 is needed for alpha = 0.05 to
#'   be attainable at all, and B >= 99 is recommended for decisions.
#' @param cutoff Retained for display only (null percentile reference); the
#'   decision uses `alpha`, not the percentile.
#' @param ss.lower Minimum sum-score-group size passed to [PrepareChecks()].
#' @param propagate_item_error Draw fresh item difficulties per replicate by
#'   refitting a provisional Rasch draw (full parametric bootstrap of the item
#'   parameters), so the null reflects item-estimation error rather than
#'   treating the plug-in estimates as exact. Costs one extra [fit_rm()] per
#'   replicate. Dichotomous data only (ignored, with a warning, for
#'   polytomous). Default `FALSE` (Student & Read's plug-in convention).
#' @param person_order Person ordering for score grouping when responses are
#'   missing (see [PrepareChecks()]): `"complete"` (default) uses complete
#'   cases only; `"facility"`/`"adjusted"` keep all respondents at the cost of
#'   an extra-ordinal commensuration assumption. Identical on complete data.
#' @param latent How person abilities are drawn in the null replicates.
#'   `"empirical"` (default) samples from the latent distribution *estimated
#'   from the data* (posterior mass at the quadrature nodes, the Bock-Aitkin
#'   empirical histogram), so the null reproduces the observed ability
#'   distribution - bimodal, skewed, or censored samples included - and its
#'   sum-score group structure; any observed-vs-null excess is then
#'   attributable to non-additivity rather than population shape (additive
#'   conjoint structure is itself distribution-free, so the latent shape is a
#'   nuisance parameter here). `"normal"` draws theta ~ N(0, sigma^2) with the
#'   fitted latent variance.
#' @param mc.cores Cores for the bootstrap (default 1). The `B` null datasets
#'   are processed in parallel with [parallel::mclapply()] on non-Windows
#'   platforms; each is seeded independently, so results are reproducible.
#' @param seed Optional integer seed.
#' @param verbose Print progress (default TRUE).
#'
#' @return An object of class `ccnull`: a list with `observed` (the observed
#'   weighted violation rate), `null` (vector of null rates), `percentile`
#'   (of `observed` within `null`, in \[0, 1]), `p_value` (upper-tail, with the
#'   `(1 + #{null >= observed}) / (B + 1)` continuity correction so it is never
#'   exactly 0), `reject` (logical, at `cutoff`), and the settings.
#'
#' @references
#' Student, S. R., & Read, W. S. (2025). Applying Bayesian checks of
#' cancellation axioms for interval scaling in limited samples.
#' \emph{Behavior Research Methods}, 57, 305.
#' \doi{10.3758/s13428-025-02844-7}
#'
#' Domingue, B. (2014). Evaluating the equal-interval hypothesis with test
#' score scales. \emph{Psychometrika}, 79(1), 1-19.
#'
#' @seealso [ConjointChecks()], [quant_fit()].
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' # Rasch data -> observed rate should sit mid-null (not rejected)
#' theta <- rnorm(1000); beta <- seq(-2, 2, length.out = 20)
#' dat <- matrix(rbinom(1000 * 20, 1, plogis(outer(theta, beta, "-"))), 1000, 20)
#' res <- cc_bootstrap_null(dat, B = 100, mc.cores = 4, seed = 1)
#' res
#' plot(res)
#' }
#' @section Relation to Student & Read (2025):
#' The design follows Student & Read: a per-dataset Rasch parametric
#' bootstrap of the [ConjointChecks()] weighted mean violation rate, observed
#' and null sharing an identical pipeline, judged by null percentile at a 95%
#' cutoff (B defaults to their 100). Two deliberate departures. First, they
#' hold each person's estimated ability fixed across replicates and note that
#' "estimated ... parameters ... inevitably contain random measurement error
#' that inflates the variance of both distributions", suggesting shrinkage as
#' a remedy; here abilities are REDRAWN each replicate from the
#' Bock-Aitkin empirical-histogram latent density (`latent = "empirical"`), a
#' deconvolution that removes exactly that inflation. Second, they estimate
#' items by CML for its distribution-free character; here items come from MML,
#' with the empirical latent default preserving the distribution-free spirit -
#' `latent = "normal"` departs from it and is not the faithful choice. Their
#' n.mat of 5000 is now the default here. Whenever n.mat meets or exceeds the
#' table's total number of distinct submatrices, every submatrix is checked
#' EXACTLY ONCE (exhaustive mode - the exact population violation rate, no
#' Monte Carlo noise; typically the case for short tests). n.mat = "all"
#' forces enumeration. Runtime scales linearly in the checks performed;
#' economy configurations remain valid (observed and null share n.mat) but
#' noisier. Coverage is reported in the checks object's means$coverage.
#'
#' @export
cc_bootstrap_null <- function(data, check = "double", n.mat = 5000, B = 100,
                              cutoff = 0.95, alpha = 0.05, ss.lower = 10,
                              null_method = c("conditional_cml", "empirical_mml"),
                              latent = c("empirical", "normal"),
                              person_order = c("complete", "facility", "adjusted"),
                              propagate_item_error = FALSE,
                              mc.cores = 1L, seed = NULL, verbose = TRUE) {
  latent <- match.arg(latent)
  null_method <- match.arg(null_method)
  person_order <- match.arg(person_order)
  if (is.data.frame(data)) data <- as.matrix(data)
  poly <- .is_polytomous(data)
  data <- if (poly || anyNA(data)) .validate_poly(data, allow_na = TRUE)
          else validate_data(data)
  n_obs <- nrow(data); J <- ncol(data)

  if (n_obs < 1000) {
    warning("N = ", n_obs, " is below ~1000; Student & Read (2025) find this ",
            "procedure underpowered at moderate samples. A non-rejection may ",
            "reflect low power rather than interval scalability.")
  }

  # Build the score-by-item count matrices. Polytomous data goes through the
  # adjacent-category recoding of prepare_polytomous(); dichotomous data uses
  # PrepareChecks() directly.
  min_rows <- if (check == "triple") 4 else 3
  prep_fn <- if (poly) {
    function(d) prepare_polytomous(d, ss.lower = ss.lower,
                                   person_order = person_order)
  } else {
    function(d) PrepareChecks(d, ss.lower = ss.lower,
                              person_order = person_order)
  }

  # 1. observed violation rate (mc.cores = 1 and seeded so the random
  #    submatrix sampling is reproducible; mclapply would not control it)
  if (verbose) cat("Computing observed violation rate...\n")
  if (!is.null(seed)) set.seed(seed)
  prep <- prep_fn(data)
  obs_obj <- ConjointChecks(prep$N, prep$n, n.mat = n.mat, check = check,
                            mc.cores = 1L)
  obs <- obs_obj@means$weighted
  # retain the observed violation topography: with observable score bands the
  # cell map reads directly as "cancellation fails for band X on item Y"
  obs_tab <- obs_obj@tab
  obs_counts <- obs_obj@check.counts
  dimnames(obs_tab) <- dimnames(obs_counts) <-
    list(band = rownames(prep$N), item = colnames(prep$N))

  # 2. fit the (partial-credit) Rasch model for a marginal parametric bootstrap:
  #    redraw abilities from N(0, sigma^2) per replicate.
  if (null_method == "conditional_cml") {
    # CML items + patterns conditional on each person's (answered set, total
    # score): no MML fit, no latent estimation, footprints preserved exactly;
    # dichotomous AND polytomous via generalized ESF (conditional_null.R)
    dl_cml <- .cml_fit_general(data)
    sim_fn <- local(function() .conditional_null_general(data, dl_cml))
    fit <- NULL
  } else {
  fit <- suppressWarnings(fit_rm(data, verbose = FALSE))
  sigma <- sqrt(rasch_latent_var(fit))
  rmf <- attr(fit, "rm_fit")
  if (latent == "empirical") {
    # Bock-Aitkin empirical-histogram refit: latent weights re-estimated
    # jointly with the item parameters, so the null reproduces the observed
    # ability distribution (and hence the sum-score group structure)
    rmf <- .rm_empirical_refit(rmf)
  }
  beta <- unlist(rmf$delta_list)
  if (poly) {
    if (propagate_item_error)
      warning("propagate_item_error is not implemented for polytomous data; ignored")
    sim_fn <- function() .simulate_pcm(rmf, n_obs, sigma, latent = latent)
  } else {
    sim_fn <- function() {
      b <- beta
      if (propagate_item_error) {
        # full parametric bootstrap of the item parameters: draw beta* with
        # its sampling variability by refitting a provisional Rasch draw,
        # so the null width reflects item-estimation error (plug-in betas
        # understate it, increasingly with J)
        d0 <- matrix(rbinom(n_obs * J, 1,
              plogis(outer(.rm_draw_theta(rmf, n_obs, sigma, latent), beta, "-"))),
              n_obs, J)
        f0 <- tryCatch(suppressWarnings(fit_rm(d0, verbose = FALSE)),
                       error = function(e) NULL)
        if (!is.null(f0)) b <- unlist(attr(f0, "rm_fit")$delta_list)
      }
      theta <- .rm_draw_theta(rmf, n_obs, sigma, latent)
      matrix(rbinom(n_obs * J, 1, plogis(outer(theta, b, "-"))), n_obs, J)
    }
  }
  }  # end empirical_mml branch

  # 3. simulate B null datasets and compute each violation rate
  if (verbose) cat("Simulating", B, if (poly) "partial-credit" else "Rasch",
                   "null datasets...\n")
  if (!is.null(seed)) set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, B)

  boot_one <- function(b) {
    set.seed(rep_seeds[b])
    d <- (if (null_method == "conditional_cml") sim_fn() else .impose_mask(sim_fn(), data))   # null replicates share the observed missingness
    p <- tryCatch(prep_fn(d), error = function(e) NULL)
    if (is.null(p) || nrow(p$N) < min_rows) {
      return(NA_real_)
    }
    tryCatch(
      ConjointChecks(p$N, p$n, n.mat = n.mat, check = check,
                     mc.cores = 1L)@means$weighted,
      error = function(e) NA_real_)
  }

  raw <- par_lapply(seq_len(B), boot_one, mc.cores)
  null <- vapply(raw, function(z)
    if (is.numeric(z) && length(z) == 1L) z else NA_real_, numeric(1))
  n_failed <- sum(is.na(null))
  null <- null[!is.na(null)]
  if (length(null) == 0L) {
    stop("All ", B, " null simulations failed; cannot calibrate")
  }
  if (n_failed > 0.1 * B) {
    warning(n_failed, " of ", B, " null simulations failed and were dropped")
  }

  percentile <- mean(null < obs)
  # upper-tail p with the (1 + .)/(B + 1) continuity correction (avoids p = 0)
  p_value <- (1 + sum(null >= obs)) / (length(null) + 1)
  structure(list(observed = obs, null = sort(null),
                 percentile = percentile, p_value = p_value,
                 reject = p_value <= alpha, alpha = alpha, cutoff = cutoff,
                 check = check, N = n_obs, J = J, B = length(null),
                 n_failed = n_failed,
                 obs_tab = obs_tab, obs_counts = obs_counts,
                 obs_checks = obs_obj@Checks),
            class = "ccnull")
}

#' Print method for ccnull objects
#'
#' @param x A ccnull object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#' @export
print.ccnull <- function(x, ...) {
  is_kara <- x$check %in% c("omni-KL", "kara-KL")
  if (is_kara) {
    cat("\nOmnibus cancellation-hierarchy test (Karabatsos KL, Rasch null)\n")
    cat("---------------------------------------------\n")
    cat("Statistic        : Karabatsos global KL (banded)\n")
    cat(sprintf("Observed KL      : %.4f\n", x$observed))
  } else {
    cat("\nBootstrapped conjoint-check null (Student & Read, 2025)\n")
    cat("------------------------------------------------------\n")
    cat(sprintf("Axiom            : %s cancellation\n", x$check))
    cat(sprintf("Observed rate    : %.4f\n", x$observed))
  }
  cat(sprintf("N x J            : %d x %d,  %d null datasets\n", x$N, x$J, x$B))
  cat(sprintf("Null (Rasch)     : mean %.4f, 95%%ile %.4f, max %.4f\n",
              mean(x$null), stats::quantile(x$null, 0.95, names = FALSE),
              max(x$null)))
  if (is_kara && !is.null(x$kl_median)) {
    cat(sprintf("Per-cell KL      : median %.3f, Q3 %.3f, 90%% %.3f, max %.3f\n",
                x$kl_median, x$kl_q3, x$kl_p90, x$kl_max))
    if (!is.null(x$ess_min) && is.finite(x$ess_min)) {
      cat(sprintf("Sampler ESS      : min %.0f, median %.0f%s\n", x$ess_min,
                  x$ess_median,
                  if (x$ess_min < 50) "  ** LOW - increase S **" else ""))
    }
  }
  cat(sprintf("Percentile       : %.1f%%  (p = %.3f%s)\n",
              100 * x$percentile, x$p_value,
              if (!is.null(x$p_adjusted)) sprintf(", Holm-adj = %.3f", x$p_adjusted) else ""))
  cat(sprintf("Interval scaling : %s at the %.0f%% cutoff\n",
              if (x$reject) "REJECTED" else "not rejected", 100 * x$cutoff))
  if (x$reject && !is.null(x$obs_tab)) {
    # name the worst cells (min 2 checks so a single draw can't dominate)
    tb <- x$obs_tab; ct <- x$obs_counts
    ok <- !is.na(tb) & !is.na(ct) & ct >= 2
    if (any(ok)) {
      idx <- order(tb[ok], decreasing = TRUE)[seq_len(min(3, sum(ok)))]
      cells <- which(ok, arr.ind = TRUE)[idx, , drop = FALSE]
      lab <- apply(cells, 1L, function(rc) sprintf("band %s x item %s (%.0f%%)",
        if (is.null(rownames(tb))) rc[1] else rownames(tb)[rc[1]],
        if (is.null(colnames(tb))) rc[2] else colnames(tb)[rc[2]],
        100 * tb[rc[1], rc[2]]))
      cat("Worst cells      :", paste(lab, collapse = ", "), "\n")
    }
  }
  invisible(x)
}

#' Plot method for ccnull objects
#'
#' Histogram of the Rasch null distribution with the observed rate marked.
#'
#' @param x A ccnull object.
#' @param ... Passed to [graphics::hist()].
#' @return Invisibly returns x.
#' @export
plot.ccnull <- function(x, ...) {
  is_kara <- x$check %in% c("omni-KL", "kara-KL")
  rng <- range(c(x$null, x$observed))
  graphics::hist(x$null, breaks = "FD", col = "grey85", border = "white",
                 xlim = rng,
                 main = if (is_kara) "Karabatsos KL null vs observed"
                        else "Conjoint-check null vs observed",
                 xlab = if (is_kara) "Global Kullback-Leibler divergence"
                        else "Weighted mean proportion of violations", ...)
  graphics::abline(v = stats::quantile(x$null, x$cutoff, names = FALSE),
                   lty = 2, col = "grey40")
  graphics::abline(v = x$observed, col = "firebrick", lwd = 2)
  graphics::legend("topright", bty = "n",
                   legend = c("observed", sprintf("%.0f%% cutoff", 100 * x$cutoff)),
                   col = c("firebrick", "grey40"), lwd = c(2, 1), lty = c(1, 2))
  invisible(x)
}

#' Sequential hierarchy of calibrated cancellation checks
#'
#' Runs the calibrated cancellation checks in their logical order - single,
#' then double, then triple - stopping at the first level that rejects. The
#' cancellation axioms form a hierarchy: the double-cancellation condition is
#' only a logically distinct requirement once the orderings tested by single
#' cancellation hold, and triple likewise presupposes double. (Each individual
#' [cc_bootstrap_null()] check already imposes the *joint* constraint set up to
#' its level, so a deep check is always well-formed; what it cannot do is say
#' *where* in the hierarchy additivity failed.) Testing sequentially therefore
#' adds attribution: a rejection at the double level, reached only after single
#' passed, is evidence against the genuinely double-cancellation part of
#' additive structure rather than a re-detection of an ordering failure.
#'
#' Each level is calibrated exactly as in [cc_bootstrap_null()]: the observed
#' violation rate is located in a null distribution from data simulated under
#' the Rasch (or, for polytomous data, partial credit) model fitted to `data`,
#' passed through the identical pipeline at the same level.
#'
#' @param data Response matrix (persons x items), dichotomous or polytomous.
#' @param levels Which levels to run, in order (default
#'   `c("single", "double", "triple")`).
#' @param person_order As in [cc_bootstrap_null()].
#' @param n.mat,B,cutoff,ss.lower,latent,mc.cores,seed,verbose As in
#'   [cc_bootstrap_null()]; each level uses `seed`, `seed + 1L`, `seed + 2L`.
#'
#' @return An object of class `cchier`: a list with `levels` (the per-level
#'   `ccnull` objects actually run), `attribution` (`"none"` when every level
#'   passes, otherwise the first level whose HOLM-ADJUSTED p is <= alpha; all
#'   requested levels are always run), `p_raw`/`p_adjusted`, and
#'   `supports_quant` (`TRUE` iff no adjusted p rejects).
#'
#' @seealso [cc_bootstrap_null()], [quant_fit()]
#' @export
cc_bootstrap_hierarchy <- function(data, levels = c("single", "double", "triple"),
                                   n.mat = 5000, B = 100, cutoff = 0.95,
                                   alpha = 0.05, ss.lower = 10,
                                   null_method = c("conditional_cml", "empirical_mml"),
                              latent = c("empirical", "normal"),
                                   person_order = c("complete", "facility", "adjusted"),
                                   mc.cores = 1L, seed = NULL,
                                   verbose = TRUE) {
  person_order <- match.arg(person_order)
  levels <- match.arg(levels, c("single", "double", "triple"),
                      several.ok = TRUE)
  latent <- match.arg(latent)
  null_method <- match.arg(null_method)
  # Run ALL requested levels, then Holm-adjust the family of p-values so the
  # familywise error of "reject if any level rejects" is controlled at alpha.
  # (The former stop-at-first-rejection sequence ran up to three ~alpha tests
  # with no correction - FWER near 1-(1-alpha)^3, visible empirically as
  # ~18% false rejection of true-Rasch data in the TID validation.)
  out <- list()
  for (k in seq_along(levels)) {
    lv <- levels[k]
    if (verbose) cat("Level", k, "of", length(levels), ":", lv,
                     "cancellation...\n")
    out[[lv]] <- cc_bootstrap_null(data, check = lv, n.mat = n.mat, B = B,
                           cutoff = cutoff, alpha = alpha, ss.lower = ss.lower,
                           null_method = null_method,
                           latent = latent, person_order = person_order,
                           mc.cores = mc.cores,
                           seed = if (!is.null(seed)) seed + k - 1L else NULL,
                           verbose = FALSE)
  }
  p_raw <- vapply(out, function(r) r$p_value, numeric(1))
  p_adj <- stats::p.adjust(p_raw, method = "holm")
  for (k in seq_along(out)) {
    out[[k]]$p_adjusted <- p_adj[k]
    out[[k]]$reject <- p_adj[k] <= alpha        # decision on ADJUSTED p
  }
  rejected <- names(out)[p_adj <= alpha]
  attribution <- if (length(rejected)) rejected[1] else "none"
  structure(list(levels = out, attribution = attribution,
                 stopped_at = levels[length(levels)],
                 p_raw = p_raw, p_adjusted = p_adj, alpha = alpha,
                 supports_quant = identical(attribution, "none")),
            class = "cchier")
}

#' Print method for cchier objects
#'
#' @param x A cchier object.
#' @param ... Ignored.
#' @return Invisibly returns x.
#' @export
print.cchier <- function(x, ...) {
  cat("\nSequential calibrated cancellation checks\n")
  cat("------------------------------------------\n")
  for (lv in names(x$levels)) {
    r <- x$levels[[lv]]
    cat(sprintf("  %-7s rate = %.4f  p = %.3f  Holm-adj = %.3f  %s\n",
                lv, r$observed, r$p_value,
                if (!is.null(r$p_adjusted)) r$p_adjusted else NA,
                if (isTRUE(r$reject)) "REJECTED" else "not rejected"))
  }
  if (identical(x$attribution, "none")) {
    cat("\nAll tested levels consistent with additive structure.\n")
  } else {
    cat(sprintf("\nAdditivity fails at the %s-cancellation level%s.\n",
                x$attribution,
                if (x$attribution == "triple") ""
                else " (all levels run; Holm-adjusted decisions)"))
  }
  invisible(x)
}
