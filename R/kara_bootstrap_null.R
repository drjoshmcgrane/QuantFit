#' Group persons into observable score bands for the axiom checks
#'
#' Bands persons by an observable ordering - the sum score by default - and
#' aggregates responses into a band-by-item count matrix, with each band's
#' mean score as the row metric and items ordered by weighted facility.
#' Band boundaries fall only between distinct score values, so tied persons
#' are never split; when there are more distinct values than `n_bands`,
#' adjacent value-groups are merged toward equal person counts.
#'
#' Conditioning on the observable ordering keeps the additivity test free of
#' any fitted model: under the Rasch null the sum score is sufficient for
#' ability, so nothing is lost when the null is true - and nothing is
#' presupposed when it is false (a fitted-theta banding would condition the
#' axiom test on the very structure under scrutiny).
#'
#' @param data Binary response matrix (persons x items).
#' @param n_bands Maximum number of score bands.
#' @param person_order How to order persons when responses are missing.
#'   `"complete"` (default) uses complete cases only - the assumption-free
#'   frame, since any ordering across different answered item sets imports
#'   extra-ordinal structure. `"facility"` keeps all persons, ordering by
#'   proportion correct among answered items (difficulty-blind).
#'   `"adjusted"` orders by the observed correct count centred and scaled by
#'   the facility-implied expectation over the answered items
#'   (difficulty-aware, but imports a metric commensuration). Complete data
#'   are identical under all three.
#' @return A list with `N`, `n` (band x item observation/correct counts) and
#'   `ability` (band mean scores), or `NULL` on failure.
#' @keywords internal
band_by_score <- function(data, n_bands, person_order = "complete") {
  tryCatch({
    if (person_order == "complete" && anyNA(data)) {
      cc <- stats::complete.cases(data)
      if (sum(cc) < 2L * n_bands) return(NULL)
      data <- data[cc, , drop = FALSE]
    }
    rs <- switch(person_order,
      complete = rowSums(data, na.rm = TRUE),
      facility = rowMeans(data, na.rm = TRUE),
      adjusted = {
        p <- colMeans(data, na.rm = TRUE)          # weighted facility
        obs <- !is.na(data)
        e <- as.vector(obs %*% p)                  # expected correct, answered set
        v <- as.vector(obs %*% (p * (1 - p)))
        (rowSums(data, na.rm = TRUE) - e) / sqrt(pmax(v, 1e-12))
      },
      stop("unknown person_order: ", person_order))
    vals <- sort(unique(rs))
    if (length(vals) < 2L) return(NULL)
    # <= n_bands contiguous value-groups, boundaries only between distinct
    # values (ties never split), greedily balanced person counts
    cnt <- as.integer(table(factor(rs, levels = vals)))
    tgt <- sum(cnt) / min(n_bands, length(vals))
    gid <- integer(length(vals)); g <- 1L; acc <- 0L
    for (k in seq_along(vals)) {
      gid[k] <- g
      acc <- acc + cnt[k]
      if (acc >= tgt && g < n_bands && k < length(vals)) { g <- g + 1L; acc <- 0L }
    }
    grp <- gid[match(rs, vals)]
    keep <- sort(unique(grp))
    J <- ncol(data)
    nb_bands <- length(keep)
    if (nb_bands < 2L) return(NULL)
    Nb <- nb <- matrix(0, nb_bands, J); ab <- numeric(nb_bands)
    for (k in seq_along(keep)) {
      idx <- grp == keep[k]
      sub <- data[idx, , drop = FALSE]
      # per-cell N: respondents in the band who ANSWERED the item (missing
      # responses just lower the cell's observation count)
      Nb[k, ] <- colSums(!is.na(sub))
      nb[k, ] <- colSums(sub, na.rm = TRUE)
      ab[k] <- mean(rs[idx])
    }
    # items ordered by weighted facility (observed-correct / observed-count)
    fac <- colSums(data, na.rm = TRUE) / pmax(colSums(!is.na(data)), 1L)
    ord <- order(fac)
    if (any(Nb == 0)) return(NULL)          # a band with no observations
    list(N = Nb[, ord], n = nb[, ord], ability = ab)
  }, error = function(e) NULL)
}

#' Omnibus cancellation-hierarchy test (Karabatsos KL, bootstrapped null)
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
#' The observed global KL is computed on a sum-score-banded matrix (see
#' [band_by_score()] - conditioning is on the observable score ordering, never
#' on a fitted ability); `B` datasets are then simulated under the fitted
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
#' @param n_bands Maximum number of score bands (default 6).
#' @param B Number of Rasch-simulated null datasets (default 50).
#' @param cutoff Null percentile above which interval scaling is rejected
#'   (default 0.95).
#' @param S,N_synth Iterations and synthetic datasets for each [KaraChecks()]
#'   run (defaults 10000, 100), used identically for the observed and null.
#' @param latent How person abilities are drawn in the null replicates:
#'   `"empirical"` (default) samples from the latent distribution estimated
#'   from the data, `"normal"` draws theta ~ N(0, sigma^2). See
#'   [cc_bootstrap_null()] for the rationale.
#' @param person_order Person ordering used for score banding when responses
#'   are missing; see [band_by_score()]. `"complete"` (default) bands complete
#'   cases only; `"facility"` and `"adjusted"` keep all persons at the cost of
#'   an extra-ordinal commensuration assumption. Complete data are identical
#'   under all three.
#' @param propagate_item_error As in [cc_bootstrap_null()]: bootstrap the item
#'   difficulties per replicate instead of treating plug-ins as exact.
#' @param mc.cores Cores for the bootstrap (default 1); replicates are seeded
#'   independently so parallel and serial results agree.
#' @param seed Optional integer seed.
#' @param verbose Print progress (default TRUE).
#'
#' @return An object of class `ccnull` (shared with [cc_bootstrap_null()]):
#'   `observed` (global KL), `null`, `percentile`, `p_value`, `reject`, and
#'   settings. Its `check` field is `"omni-KL"`.
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
omni_bootstrap_null <- function(data, n_bands = 6L, B = 50, cutoff = 0.95,
                                S = 10000, N_synth = 100,
                                latent = c("empirical", "normal"),
                                person_order = c("complete", "facility", "adjusted"),
                                propagate_item_error = FALSE,
                                mc.cores = 1L, seed = NULL, verbose = TRUE) {
  latent <- match.arg(latent)
  person_order <- match.arg(person_order)
  if (is.data.frame(data)) data <- as.matrix(data)
  poly <- .is_polytomous(data)
  data <- if (poly || anyNA(data)) .validate_poly(data, allow_na = TRUE)
          else validate_data(data)
  n_obs <- nrow(data); J <- ncol(data)

  run_kara <- function(band) {
    if (is.null(band)) return(NULL)
    nr <- nrow(band$N); nc <- ncol(band$N)
    KaraChecks(as.vector(band$N), as.vector(band$n), S = S,
               N_synth = N_synth, mc.cores = 1L, verbose = FALSE,
               testscore = rep(band$ability, nc), item = rep(1:nc, each = nr))
  }
  # KaraChecks on a polytomous dataset via the adjacent-category, total-score
  # conditioned matrix (prepare_polytomous), using score groups as the person
  # factor and sub-items as the item factor.
  run_kara_poly <- function(d) {
    prep <- tryCatch(prepare_polytomous(d, ss.lower = 10, person_order = person_order),
                     error = function(e) NULL)
    if (is.null(prep)) return(NULL)
    nr <- nrow(prep$N); nc <- ncol(prep$N)
    KaraChecks(as.vector(prep$N), as.vector(prep$n), S = S,
               N_synth = N_synth, mc.cores = 1L, verbose = FALSE,
               testscore = rep(as.numeric(rownames(prep$N)), nc),
               item = rep(1:nc, each = nr))
  }
  kara_of <- function(d) {
    if (poly) run_kara_poly(d)
    else run_kara(band_by_score(d, n_bands, person_order = person_order))
  }
  global_kl <- function(d) {
    kc <- kara_of(d)
    if (is.null(kc)) NA_real_ else kc$global_KL
  }

  # 1. observed global KL + per-cell KL quantiles; seeded so the KaraChecks
  #    sampler is reproducible (it runs at mc.cores = 1 inside)
  if (verbose) cat("Computing observed Karabatsos global KL...\n")
  if (!is.null(seed)) set.seed(seed)   # covers the null-generator fit's starts
  fit <- suppressWarnings(fit_rm(data, verbose = FALSE))  # null generator only
  kc_obs <- kara_of(data)
  if (is.null(kc_obs)) stop("Observed KaraChecks failed")
  obs <- kc_obs$global_KL
  kl_q <- stats::quantile(as.vector(kc_obs$KL), c(0.5, 0.75, 0.90, 1),
                          names = FALSE)

  # 2. parameters for a marginal parametric bootstrap (redraw abilities from
  #    N(0, sigma^2) per replicate, not the fixed EAP estimates)
  sigma <- sqrt(rasch_latent_var(fit))
  rmf <- attr(fit, "rm_fit")
  if (latent == "empirical") {
    # Bock-Aitkin empirical-histogram refit (see cc_bootstrap_null())
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
        # parametric bootstrap of item difficulties (see cc_bootstrap_null)
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

  # 3. simulate B null datasets
  if (verbose) cat("Simulating", B, if (poly) "partial-credit" else "Rasch",
                   "null datasets (KaraChecks each)...\n")
  if (!is.null(seed)) set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, B)
  boot_one <- function(b) {
    set.seed(rep_seeds[b])
    tryCatch(global_kl(.impose_mask(sim_fn(), data)), error = function(e) NA_real_)
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
                 kl_max = kl_q[4], check = "omni-KL", N = n_obs, J = J,
                 B = length(null), n_failed = n_failed),
            class = "ccnull")
}

#' @rdname omni_bootstrap_null
#' @export
kara_bootstrap_null <- function(...) omni_bootstrap_null(...)
