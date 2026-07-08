#' Bayesian omnibus test of conjoint measurement axioms (Karabatsos, 2018)
#'
#' Implements the Bayesian omnibus test of the additive conjoint measurement
#' axioms of Karabatsos (2018), based on approximate Bayesian computation
#' with a synthetic likelihood. An importance sampler draws `S` candidate
#' Bernoulli success probabilities per cell from a Beta(`a`, `b`) prior; for
#' each candidate, `N_synth` synthetic data sets are generated and reduced to
#' an axiom-respecting summary statistic (an isotonic regression of the cell
#' proportions on the predictions of a two-stage logistic model, enforcing
#' the ordering implied by the conjoint measurement axioms). A kernel density
#' estimate over the synthetic statistics forms the synthetic likelihood used
#' as importance weights. For each cell, the Kullback-Leibler divergence
#' between the unrestricted posterior mean (from the Beta-Binomial conjugate
#' update) and the axiom-restricted importance-weighted posterior mean is
#' computed; a cell is flagged as violating the axioms when KL > 0.01.
#'
#' @param N Matrix (score groups by items) containing the total number of
#'   responses per cell, or a vector of totals. If a vector is supplied,
#'   `testscore` and `item` must also be supplied.
#' @param n Matrix (score groups by items) containing the number of correct
#'   responses per cell, or a vector of counts.
#' @param S Number of importance sampling iterations (prior draws).
#' @param N_synth Number of synthetic data sets generated per prior draw to
#'   build the synthetic likelihood.
#' @param a,b Shape parameters of the Beta prior on the cell success
#'   probabilities.
#' @param mc.cores The number of cores to parallelize over (forking is used,
#'   so parallelization is unavailable on Windows). When `mc.cores > 1`, set
#'   `RNGkind("L'Ecuyer-CMRG")` (and a seed) before calling for reproducible
#'   results across parallel workers.
#' @param verbose Print progress and summary information.
#' @param testscore Numeric vector giving the score-group (ability) value of
#'   each cell. Required when `N` and `n` are vectors; ignored (and derived
#'   from the matrix layout) when `N` is a matrix.
#' @param item Vector giving the item of each cell. Required when `N` and
#'   `n` are vectors; ignored when `N` is a matrix.
#'
#' @details
#' When `N` is a matrix, its cells are internally vectorized in column-major
#' order, i.e. with `testscore` (score group, the row index) varying fastest
#' within `item` (the column index), exactly as `as.vector(N)` would.
#' Vector inputs are assumed to follow the same ordering when the per-cell
#' results are reshaped into score-group by item matrices: the reshape is
#' only performed when the input was a matrix or when `testscore` and `item`
#' jointly define a complete grid (each score-group/item combination appears
#' exactly once); otherwise the per-cell results are returned as vectors in
#' the input order.
#'
#' @return A list with components:
#' \describe{
#'   \item{`KL`}{Per-cell Kullback-Leibler divergences.}
#'   \item{`global_KL`}{Sum of the per-cell KL divergences (global test
#'     statistic).}
#'   \item{`theta_bar`}{Axiom-restricted (importance-weighted) posterior mean
#'     of each cell probability.}
#'   \item{`theta_0`}{Unrestricted Beta-Binomial posterior mean of each cell
#'     probability.}
#'   \item{`t_obs`}{Observed isotonic summary statistic for each cell.}
#'   \item{`ESS`}{Effective sample size of the importance weights for each
#'     cell.}
#'   \item{`violations`}{Logical matrix flagging cells with KL > 0.01.}
#'   \item{`n_violations`}{Total number of flagged cells.}
#'   \item{`ZMn`}{Importance-weighted posterior mean of the standardized
#'     residuals.}
#'   \item{`parameters`}{The values of `S`, `N_synth`, `a`, and `b` used.}
#' }
#'
#' @references
#' Karabatsos, G. (2018). On Bayesian testing of additive conjoint
#' measurement axioms using synthetic likelihood. \emph{Psychometrika},
#' 83(2), 321-332. \doi{10.1007/s11336-017-9581-x}
#'
#' @seealso [PrepareChecks()], [ConjointChecks()]
#'
#' @examples
#' \dontrun{
#' # simulated Rasch example
#' n.items <- 20
#' n.respondents <- 2000
#' diff <- rnorm(n.items)
#' abil <- rnorm(n.respondents)
#' kern <- outer(abil, diff, "-")
#' pv <- exp(kern)/(1+exp(kern))
#' resp <- ifelse(pv > runif(n.items*n.respondents), 1, 0)
#' tmp <- PrepareChecks(resp)
#' out <- KaraChecks(tmp$N, tmp$n, S = 10000)
#' out$global_KL
#' out$n_violations
#' }
#'
#' @export
KaraChecks <- function(N, n, S = 30000, N_synth = 100, a = 0.5, b = 0.5,
                       mc.cores = parallel::detectCores() - 1, verbose = TRUE,
                       testscore = NULL, item = NULL) {

  if (is.matrix(N)) {
    if (!is.null(testscore) || !is.null(item)) {
      warning("'testscore' and 'item' are ignored when N is a matrix; they are derived from the matrix layout.")
    }
    nr <- nrow(N)
    nc <- ncol(N)

    ability_values <- seq(-(nr - 1) / 2, (nr - 1) / 2, length.out = nr)
    # column-major vectorization: testscore (score group) varies fastest within item
    testscore <- rep(ability_values, nc)
    item <- rep(1:nc, each = nr)

    N_vec <- as.vector(N)
    n_vec <- as.vector(n)
    reshape_output <- TRUE
  } else {
    if (is.null(testscore) || is.null(item)) {
      stop("When N and n are supplied as vectors, both 'testscore' and 'item' must be supplied ",
           "(one value per cell, in the same order as N and n).")
    }
    if (length(testscore) != length(N) || length(item) != length(N)) {
      stop("'testscore' (length ", length(testscore), ") and 'item' (length ", length(item),
           ") must each have the same length as N (length ", length(N), ").")
    }
    if (length(n) != length(N)) {
      stop("'n' (length ", length(n), ") must have the same length as N (length ", length(N), ").")
    }
    N_vec <- N
    n_vec <- n
    nr <- length(unique(testscore))
    nc <- length(unique(item))
    # only reshape per-cell results into nr x nc matrices when testscore/item
    # define a complete grid (each combination appearing exactly once) and the
    # cells are ordered column-major with testscore varying fastest
    complete_grid <- (nr * nc == length(N_vec)) &&
      (anyDuplicated(paste(testscore, item, sep = "\r")) == 0)
    reshape_output <- complete_grid &&
      identical(order(match(item, unique(item)), match(testscore, unique(testscore))),
                seq_along(N_vec))
  }

  IJ <- length(N_vec)

  item_dummies <- model.matrix(~ factor(item) - 1)
  X <- cbind(testscore, item_dummies)

  if (verbose) {
    cat("Running KaraChecks with S =", S, "iterations, N_synth =", N_synth, "\n")
    cat("Using", mc.cores, "cores\n")
  }

  dat <- n_vec / N_vec

  fit1 <- suppressWarnings(glm(cbind(n_vec, N_vec - n_vec) ~ X - 1,
                                family = binomial(link = "logit")))
  xhat <- predict(fit1, type = "response")

  fit2 <- suppressWarnings(glm(cbind(n_vec, N_vec - n_vec) ~ xhat,
                                family = binomial(link = "logit")))
  xhat2 <- predict(fit2, type = "response")

  ty <- lsqisotonic(xhat2, dat, N_vec)

  worker_fn <- function(s_indices) {
    n_batch <- length(s_indices)
    thetas_batch <- matrix(NA, n_batch, IJ)
    ws_batch <- matrix(NA, n_batch, IJ)
    Zs_batch <- matrix(NA, n_batch, IJ)
    TYSTAR <- matrix(NA, N_synth, IJ)

    for (i in seq_along(s_indices)) {
      theta_s <- rbeta(IJ, a, b)
      thetas_batch[i, ] <- theta_s

      for (m in 1:N_synth) {
        rstar <- rbinom(IJ, N_vec, theta_s)

        fit_synth <- tryCatch({
          suppressWarnings(glm(cbind(rstar, N_vec - rstar) ~ X - 1,
                                family = binomial(link = "logit")))
        }, error = function(e) NULL)

        if (!is.null(fit_synth)) {
          xhat_synth <- predict(fit_synth, type = "response")
          TYSTAR[m, ] <- lsqisotonic(xhat_synth, rstar / N_vec, N_vec)
        } else {
          TYSTAR[m, ] <- rstar / N_vec
        }
      }

      for (j in 1:IJ) {
        ws_batch[i, j] <- ksdensity(TYSTAR[, j], ty[j]) + .Machine$double.eps
      }

      Zs_batch[i, ] <- (n_vec - N_vec * theta_s) / sqrt(N_vec * theta_s * (1 - theta_s))
    }

    list(thetas = thetas_batch, ws = ws_batch, Zs = Zs_batch)
  }

  batch_size <- ceiling(S / mc.cores)
  batches <- split(1:S, ceiling(seq_along(1:S) / batch_size))

  if (verbose) cat("Processing", length(batches), "batches...\n")

  if (mc.cores > 1 && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(batches, worker_fn, mc.cores = mc.cores)
  } else {
    results <- lapply(seq_along(batches), function(i) {
      if (verbose) cat("  Batch", i, "of", length(batches), "\n")
      worker_fn(batches[[i]])
    })
  }

  thetas <- do.call(rbind, lapply(results, `[[`, "thetas"))
  ws <- do.call(rbind, lapply(results, `[[`, "ws"))
  Zs <- do.call(rbind, lapply(results, `[[`, "Zs"))

  sumw <- colSums(ws)
  W <- ws / matrix(sumw, S, IJ, byrow = TRUE)

  thetaMn <- colSums(thetas * W)
  ZMn <- colSums(Zs * W)
  ESS <- 1 / colSums(W^2)

  thetau <- (a + n_vec) / (a + b + N_vec)

  thetaMn[thetaMn < 1e-10] <- 1e-10
  thetaMn[thetaMn > 1 - 1e-10] <- 1 - 1e-10
  thetau[thetau < 1e-10] <- 1e-10
  thetau[thetau > 1 - 1e-10] <- 1 - 1e-10

  KL <- thetau * log(thetau / thetaMn) + (1 - thetau) * log((1 - thetau) / (1 - thetaMn))

  violations <- KL > 0.01

  if (reshape_output) {
    KL_mat <- matrix(KL, nr, nc)
    violations_mat <- matrix(violations, nr, nc)
    ESS_mat <- matrix(ESS, nr, nc)
    thetaMn_mat <- matrix(thetaMn, nr, nc)
    thetau_mat <- matrix(thetau, nr, nc)
    ty_mat <- matrix(ty, nr, nc)
    ZMn_mat <- matrix(ZMn, nr, nc)
  } else {
    KL_mat <- KL
    violations_mat <- violations
    ESS_mat <- ESS
    thetaMn_mat <- thetaMn
    thetau_mat <- thetau
    ty_mat <- ty
    ZMn_mat <- ZMn
  }

  global_KL <- sum(KL)

  if (verbose) {
    cat("\nKaraChecks complete.\n")
    cat("Global KL:", round(global_KL, 3), "\n")
    cat("Number of violations (KL > 0.01):", sum(violations), "of", IJ, "cells\n")
    cat("Median ESS:", round(median(ESS), 0), "\n")
  }

  list(
    KL = KL_mat,
    global_KL = global_KL,
    theta_bar = thetaMn_mat,
    theta_0 = thetau_mat,
    t_obs = ty_mat,
    ESS = ESS_mat,
    violations = violations_mat,
    n_violations = sum(violations),
    ZMn = ZMn_mat,
    parameters = list(S = S, N_synth = N_synth, a = a, b = b)
  )
}
