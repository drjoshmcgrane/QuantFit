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
#' the more defensible parametric bootstrap of the Rasch model.
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
#' @param ss.lower Minimum sum-score-group size passed to [PrepareChecks()].
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
#' @export
cc_bootstrap_null <- function(data, check = "double", n.mat = 50, B = 100,
                              cutoff = 0.95, ss.lower = 10, mc.cores = 1L,
                              seed = NULL, verbose = TRUE) {
  if (is.data.frame(data)) data <- as.matrix(data)
  poly <- .is_polytomous(data)
  data <- if (poly) .validate_poly(data) else validate_data(data)
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
    function(d) prepare_polytomous(d, ss.lower = ss.lower)
  } else {
    function(d) PrepareChecks(d, ss.lower = ss.lower)
  }

  # 1. observed violation rate (mc.cores = 1 and seeded so the random
  #    submatrix sampling is reproducible; mclapply would not control it)
  if (verbose) cat("Computing observed violation rate...\n")
  if (!is.null(seed)) set.seed(seed)
  obs <- {
    prep <- prep_fn(data)
    ConjointChecks(prep$N, prep$n, n.mat = n.mat, check = check,
                   mc.cores = 1L)@means$weighted
  }

  # 2. fit the (partial-credit) Rasch model for a marginal parametric bootstrap:
  #    redraw abilities from N(0, sigma^2) per replicate.
  fit <- suppressWarnings(fit_rm(data, verbose = FALSE))
  sigma <- sqrt(rasch_latent_var(fit))
  if (poly) {
    mfit <- attr(fit, "mirt_object")
    sim_fn <- function() .simulate_pcm(mfit, n_obs, sigma)
  } else {
    beta <- fit$delta
    sim_fn <- function() {
      theta <- rnorm(n_obs, 0, sigma)
      matrix(rbinom(n_obs * J, 1, plogis(outer(theta, beta, "-"))), n_obs, J)
    }
  }

  # 3. simulate B null datasets and compute each violation rate
  if (verbose) cat("Simulating", B, if (poly) "partial-credit" else "Rasch",
                   "null datasets...\n")
  if (!is.null(seed)) set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, B)

  boot_one <- function(b) {
    set.seed(rep_seeds[b])
    d <- sim_fn()
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
                 reject = percentile >= cutoff, cutoff = cutoff,
                 check = check, N = n_obs, J = J, B = length(null),
                 n_failed = n_failed),
            class = "ccnull")
}

#' Print method for ccnull objects
#'
#' @param x A ccnull object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#' @export
print.ccnull <- function(x, ...) {
  is_kara <- identical(x$check, "kara-KL")
  if (is_kara) {
    cat("\nBootstrapped Karabatsos global-KL null (Rasch)\n")
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
  }
  cat(sprintf("Percentile (p)   : %.1f%%  (p = %.3f)\n",
              100 * x$percentile, x$p_value))
  cat(sprintf("Interval scaling : %s at the %.0f%% cutoff\n",
              if (x$reject) "REJECTED" else "not rejected", 100 * x$cutoff))
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
  is_kara <- identical(x$check, "kara-KL")
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
