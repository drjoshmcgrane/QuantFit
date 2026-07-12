#' Group persons into ability bands for the axiom checks
#'
#' Bins persons into `n_bands` bands by their Rasch ability estimate and
#' aggregates responses into a band-by-item count matrix, with each band's
#' mean ability as the row metric. Used by the Karabatsos route so its
#' additivity test operates on a real ability metric rather than raw
#' (nonlinearly spaced) sum-score groups.
#'
#' @param data Binary response matrix (persons x items).
#' @param n_bands Number of ability bands.
#' @param theta Optional pre-computed person abilities; if `NULL`, a Rasch
#'   model is fitted to obtain them.
#' @return A list with `N`, `n` (band x item counts) and `ability` (band mean
#'   abilities), or `NULL` on failure.
#' @keywords internal
band_by_ability <- function(data, n_bands, theta = NULL) {
  tryCatch({
    if (is.null(theta)) {
      theta <- rm_scores(suppressWarnings(fit_rm(data, verbose = FALSE)))$theta
    }
    br <- stats::quantile(theta, seq(0, 1, length.out = n_bands + 1))
    br[1] <- -Inf; br[length(br)] <- Inf
    br <- unique(br)                       # tied theta -> duplicate breaks
    if (length(br) < 3L) return(NULL)      # need >= 2 usable bands
    grp <- cut(theta, br, labels = FALSE)
    keep <- sort(unique(grp))              # drop any empty band
    J <- ncol(data)
    nb_bands <- length(keep)
    Nb <- nb <- matrix(0, nb_bands, J); ab <- numeric(nb_bands)
    for (k in seq_along(keep)) {
      idx <- grp == keep[k]
      Nb[k, ] <- sum(idx); nb[k, ] <- colSums(data[idx, , drop = FALSE])
      ab[k] <- mean(theta[idx])
    }
    ord <- order(colSums(data))
    list(N = Nb[, ord], n = nb[, ord], ability = ab)
  }, error = function(e) NULL)
}

#' Bootstrapped null distribution for the Karabatsos global KL
#'
#' Calibrates the Karabatsos (2018) global Kullback-Leibler additivity
#' statistic against a null distribution simulated from the Rasch model fitted
#' to the data - the same per-dataset parametric-bootstrap logic that Student
#' & Read (2025) apply to the [ConjointChecks()] violation rate (see
#' [cc_bootstrap_null()]). Karabatsos's fixed KL > 0.01 per-cell criterion is a
#' rule of thumb whose relationship to sampling variability depends on sample
#' size, test length, and item parameters; locating the observed global KL in
#' a Rasch null distribution gives an interpretable, per-dataset percentile
#' p-value instead.
#'
#' @details
#' The observed global KL is computed on an ability-banded matrix (see
#' [band_by_ability()]); `B` datasets are then simulated under the fitted
#' *marginal* Rasch model - abilities redrawn from N(0, sigma^2) per replicate,
#' with the estimated item difficulties - each is re-fitted, re-banded, and
#' passed through [KaraChecks()] at the same `S`, and the observed global KL is
#' located in the resulting null distribution. Observed and null share the
#' identical pipeline and iteration count, so any baseline the pipeline induces
#' on additive data cancels. Interval scaling is rejected when the observed
#' global KL exceeds the `cutoff` percentile of the null.
#'
#' This is computationally heavy - each replicate is a full [KaraChecks()]
#' importance-sampling run - so `B` and `S` default lower than for the
#' cheaper CC bootstrap; increase them for a final analysis.
#'
#' @param data Binary response matrix (persons x items).
#' @param n_bands Number of ability bands (default 6).
#' @param B Number of Rasch-simulated null datasets (default 50).
#' @param cutoff Null percentile above which interval scaling is rejected
#'   (default 0.95).
#' @param S,N_synth Iterations and synthetic datasets for each [KaraChecks()]
#'   run (defaults 10000, 100), used identically for the observed and null.
#' @param mc.cores Cores for the bootstrap (default 1); replicates are seeded
#'   independently so parallel and serial results agree.
#' @param seed Optional integer seed.
#' @param verbose Print progress (default TRUE).
#'
#' @return An object of class `ccnull` (shared with [cc_bootstrap_null()]):
#'   `observed` (global KL), `null`, `percentile`, `p_value`, `reject`, and
#'   settings. Its `check` field is `"kara-KL"`.
#'
#' @references
#' Karabatsos, G. (2018). On Bayesian testing of additive conjoint measurement
#' axioms using synthetic likelihood. \emph{Psychometrika}, 83(2), 321-332.
#'
#' Student, S. R., & Read, W. S. (2025). Applying Bayesian checks of
#' cancellation axioms for interval scaling in limited samples.
#' \emph{Behavior Research Methods}, 57, 305.
#'
#' @seealso [cc_bootstrap_null()], [KaraChecks()], [quant_fit()].
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' theta <- rnorm(1000); beta <- seq(-2, 2, length.out = 20)
#' dat <- matrix(rbinom(1000 * 20, 1, plogis(outer(theta, beta, "-"))), 1000, 20)
#' kara_bootstrap_null(dat, B = 40, mc.cores = 4, seed = 1)
#' }
#' @export
kara_bootstrap_null <- function(data, n_bands = 6L, B = 50, cutoff = 0.95,
                                S = 10000, N_synth = 100, mc.cores = 1L,
                                seed = NULL, verbose = TRUE) {
  data <- validate_data(data)
  n_obs <- nrow(data); J <- ncol(data)

  run_kara <- function(band) {
    if (is.null(band)) return(NULL)
    nr <- nrow(band$N); nc <- ncol(band$N)
    KaraChecks(as.vector(band$N), as.vector(band$n), S = S,
               N_synth = N_synth, mc.cores = 1L, verbose = FALSE,
               testscore = rep(band$ability, nc), item = rep(1:nc, each = nr))
  }
  global_kl <- function(band) {
    kc <- run_kara(band)
    if (is.null(kc)) NA_real_ else kc$global_KL
  }

  # 1. observed global KL + per-cell KL quantiles (banded); seeded so the
  #    KaraChecks sampler is reproducible (it runs at mc.cores = 1 inside)
  if (verbose) cat("Computing observed Karabatsos global KL...\n")
  fit <- suppressWarnings(fit_rm(data, verbose = FALSE))
  theta0 <- rm_scores(fit)$theta
  if (!is.null(seed)) set.seed(seed)
  kc_obs <- run_kara(band_by_ability(data, n_bands, theta = theta0))
  if (is.null(kc_obs)) stop("Observed KaraChecks failed")
  obs <- kc_obs$global_KL
  kl_q <- stats::quantile(as.vector(kc_obs$KL), c(0.5, 0.75, 0.90, 1),
                          names = FALSE)

  # 2. Rasch parameters for a marginal parametric bootstrap (redraw abilities
  #    from N(0, sigma^2) per replicate, not the fixed EAP estimates)
  beta <- fit$delta
  sigma <- sqrt(rasch_latent_var(fit))

  # 3. simulate B Rasch null datasets
  if (verbose) cat("Simulating", B, "Rasch null datasets (KaraChecks each)...\n")
  if (!is.null(seed)) set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, B)
  boot_one <- function(b) {
    set.seed(rep_seeds[b])
    theta <- rnorm(n_obs, 0, sigma)
    d <- matrix(rbinom(n_obs * J, 1, plogis(outer(theta, beta, "-"))), n_obs, J)
    tryCatch(global_kl(band_by_ability(d, n_bands)), error = function(e) NA_real_)
  }
  raw <- par_lapply(seq_len(B), boot_one, mc.cores)
  null <- vapply(raw, function(z)
    if (is.numeric(z) && length(z) == 1L) z else NA_real_, numeric(1))
  n_failed <- sum(is.na(null))
  null <- null[!is.na(null)]
  if (length(null) == 0L) stop("All ", B, " null simulations failed")
  if (n_failed > 0.1 * B) {
    warning(n_failed, " of ", B, " null simulations failed and were dropped")
  }

  percentile <- mean(null < obs)
  p_value <- (1 + sum(null >= obs)) / (length(null) + 1)
  structure(list(observed = obs, null = sort(null),
                 percentile = percentile, p_value = p_value,
                 reject = percentile >= cutoff, cutoff = cutoff,
                 kl_median = kl_q[1], kl_q3 = kl_q[2], kl_p90 = kl_q[3],
                 kl_max = kl_q[4], check = "kara-KL", N = n_obs, J = J,
                 B = length(null), n_failed = n_failed),
            class = "ccnull")
}
