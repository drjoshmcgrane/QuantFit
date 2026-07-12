#' Statistically calibrated log-likelihood model selection
#'
#' @name select_ll
#' @description Parametric-bootstrap likelihood-ratio machinery for comparing
#'   latent structure models whose parameter counts coincide (UN vs MON/IIO/DM)
#'   or that are nested with inequality constraints (LCR in DM), plus an
#'   ordered selection procedure over the six Torres Irribarra & Diakow models.
NULL

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Simulate a dataset from a fitted latent class model
#'
#' Draws class memberships from the fitted class probabilities and then
#' Bernoulli responses from the fitted item probabilities.
#'
#' @param fit A qlfit object with `class_probs` (length C) and `item_probs`
#'   (n_items x C).
#' @param n_obs Number of observations to simulate.
#' @return Binary matrix (n_obs x n_items).
#' @keywords internal
simulate_from_qlfit <- function(fit, n_obs) {
  n_classes <- length(fit$class_probs)
  cls <- sample.int(n_classes, size = n_obs, replace = TRUE,
                    prob = fit$class_probs)
  # item_probs is items x classes; transpose so rows index classes,
  # then expand to an n_obs x n_items probability matrix
  prob <- t(fit$item_probs)[cls, , drop = FALSE]
  matrix(rbinom(length(prob), size = 1L, prob = prob),
         nrow = n_obs, ncol = ncol(prob))
}

#' Apply a function over a list, optionally in parallel
#'
#' Uses [parallel::mclapply()] when `mc.cores > 1` on a non-Windows platform
#' (forking is unavailable on Windows), otherwise falls back to [lapply()].
#' Callers that seed each element independently get identical results in
#' serial and parallel.
#'
#' @param X Vector or list to iterate over.
#' @param FUN Function applied to each element.
#' @param mc.cores Number of cores.
#' @return A list of results.
#' @keywords internal
par_lapply <- function(X, FUN, mc.cores = 1L) {
  if (mc.cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(X, FUN, mc.cores = mc.cores)
  } else {
    lapply(X, FUN)
  }
}

#' Fitted latent variance of a Rasch (RM) qlfit
#'
#' Extracts the estimated latent variance from the underlying mirt object so a
#' marginal parametric bootstrap can redraw abilities from N(0, sigma^2).
#' Falls back to 1 if the variance cannot be recovered.
#'
#' @param fit An RM qlfit (from [fit_rm()]).
#' @return A positive numeric latent variance.
#' @keywords internal
rasch_latent_var <- function(fit) {
  mo <- attr(fit, "mirt_object")
  v <- tryCatch(as.numeric(mirt::coef(mo, simplify = TRUE)$cov),
                error = function(e) NA_real_)
  if (length(v) == 0L || is.na(v[1]) || v[1] <= 0) 1 else v[1]
}

#' Bootstrap test of continuous (RM) vs discrete (LCR) Rasch structure
#'
#' RM and LCR are non-nested and use different estimators (marginal maximum
#' likelihood over a continuous latent density vs the EM algorithm over
#' discrete classes), so a likelihood-ratio bootstrap is not clean. Instead the
#' BIC difference \eqn{BIC_{LCR} - BIC_{RM}} is calibrated against a null
#' simulated from the fitted (continuous) RM: the discrete LCR is preferred only
#' when it fits better than the RM null would produce by chance. This keeps the
#' final step of the hierarchy consistent with the bootstrap tests above.
#'
#' @param data Binary response matrix.
#' @param fit_rm_obj,fit_lcr_obj Fitted RM and LCR models.
#' @param n_classes Number of latent classes for LCR.
#' @param alpha Significance level for choosing LCR over RM (default 0.05).
#' @param B,n_starts,use_cpp,mc.cores,seed As in [select_model_ll()].
#' @return A list with `statistic` (observed BIC difference), `p_value`
#'   (lower-tail, with the `(1 + .)/(B + 1)` correction, `NA` if every
#'   replicate failed), `available` (FALSE when the bootstrap could not be
#'   computed), and `select_lcr` (TRUE when discreteness is significant at
#'   `alpha`). The null is a marginal parametric bootstrap: abilities are
#'   redrawn from N(0, sigma^2) with the RM-estimated latent variance.
#' @keywords internal
rm_vs_lcr_test <- function(data, fit_rm_obj, fit_lcr_obj, n_classes, B = 99,
                           alpha = 0.05, n_starts = 3, use_cpp = TRUE,
                           mc.cores = 1L, seed = NULL) {
  obs <- BIC(fit_lcr_obj) - BIC(fit_rm_obj)     # negative favours LCR
  n <- nrow(data); J <- ncol(data)
  # marginal parametric bootstrap: redraw abilities from N(0, sigma^2)
  beta <- fit_rm_obj$delta
  sigma <- sqrt(rasch_latent_var(fit_rm_obj))
  if (!is.null(seed)) set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, B)
  boot_one <- function(b) {
    set.seed(rep_seeds[b])
    theta <- rnorm(n, 0, sigma)
    d <- matrix(rbinom(n * J, 1, plogis(outer(theta, beta, "-"))), n, J)
    rm_b <- tryCatch(suppressWarnings(fit_rm(d, verbose = FALSE)),
                     error = function(e) NULL)
    lcr_b <- tryCatch(
      suppressWarnings(fit_lcr(d, n_classes, n_starts = n_starts, use_cpp = use_cpp)),
      error = function(e) NULL)
    if (is.null(rm_b) || is.null(lcr_b)) return(NA_real_)
    BIC(lcr_b) - BIC(rm_b)
  }
  raw <- par_lapply(seq_len(B), boot_one, mc.cores)
  null <- vapply(raw, function(z)
    if (is.numeric(z) && length(z) == 1L) z else NA_real_, numeric(1))
  null <- null[!is.na(null)]
  # lower-tail p (LCR favoured) with the (1 + .)/(B + 1) correction
  p_lower <- if (length(null) == 0L) NA_real_
             else (1 + sum(null <= obs)) / (length(null) + 1)
  list(statistic = obs, p_value = p_lower, null = sort(null),
       available = length(null) > 0L, select_lcr = isTRUE(p_lower <= alpha))
}

#' Refit a latent class model of a given type
#'
#' @param model_type One of "UN", "MON", "IIO", "DM", "LCR".
#' @param data Binary response matrix.
#' @param n_classes Number of latent classes.
#' @param n_starts Number of random starts.
#' @param use_cpp Use the compiled EM engine.
#' @return A qlfit object.
#' @keywords internal
refit_model_type <- function(model_type, data, n_classes, n_starts, use_cpp) {
  switch(model_type,
    UN  = fit_un(data, n_classes, n_starts = n_starts, use_cpp = use_cpp),
    MON = fit_mon(data, n_classes, n_starts = n_starts, use_cpp = use_cpp),
    IIO = fit_iio(data, n_classes, n_starts = n_starts, use_cpp = use_cpp),
    DM  = fit_dm(data, n_classes, n_starts = n_starts, use_cpp = use_cpp),
    LCR = fit_lcr(data, n_classes, n_starts = n_starts, use_cpp = use_cpp),
    stop("Bootstrap refitting is not supported for model type '",
         model_type, "'")
  )
}

# ---------------------------------------------------------------------------
# Bootstrap LL-equivalence test
# ---------------------------------------------------------------------------

#' Parametric bootstrap log-likelihood equivalence test
#'
#' Tests whether a constrained latent structure model (e.g. MON, IIO, DM, or
#' LCR) is statistically equivalent to a more general model (typically the
#' unconstrained UN model, or DM when testing LCR) via a parametric bootstrap
#' of the likelihood-ratio statistic.
#'
#' @param data Binary response matrix (persons x items) the models were
#'   fitted to.
#' @param fit_constrained The fitted constrained (null) model: a qlfit object
#'   of type "MON", "IIO", "DM", or "LCR". Bootstrap datasets are simulated
#'   from this fit.
#' @param fit_un The fitted more-general (alternative) model: a qlfit object,
#'   typically of type "UN" (or "DM" when the constrained model is "LCR").
#' @param B Number of bootstrap replicates (default 99).
#' @param n_starts Number of random starts for each bootstrap refit
#'   (default 3; fewer than for the observed-data fits since the bootstrap
#'   truth is close to the initialization heuristics).
#' @param seed Optional integer seed. When supplied the entire procedure,
#'   including all bootstrap refits, is deterministic.
#' @param use_cpp Use the compiled C++ EM engine for the refits
#'   (default TRUE).
#' @param mc.cores Number of cores for the bootstrap refits (default 1).
#'   Values above 1 use [parallel::mclapply()] on non-Windows platforms;
#'   because each replicate is seeded independently, results are identical
#'   to the serial run regardless of `mc.cores`.
#' @param verbose Print progress every 10 replicates (default FALSE; ignored
#'   when `mc.cores > 1`).
#'
#' @return An object of class `qleqtest`: a list with elements
#'   \describe{
#'     \item{statistic}{Observed likelihood-ratio statistic
#'       \eqn{2(\ell_{general} - \ell_{constrained})}, floored at 0.}
#'     \item{p_value}{Bootstrap p-value
#'       \eqn{(1 + \#\{LR^*_b \ge LR_{obs}\}) / (B_{eff} + 1)}.}
#'     \item{B}{Requested number of replicates.}
#'     \item{B_effective}{Replicates that produced a valid statistic.}
#'     \item{n_failed}{Replicates dropped because a refit failed.}
#'     \item{null_distribution}{Sorted vector of bootstrap statistics.}
#'     \item{models}{Character vector `c(constrained_label, general_label)`.}
#'   }
#'
#' @details
#' UN and the ordinally constrained models (MON, IIO, DM) have identical
#' parameter counts: inequality constraints do not reduce the model's
#' dimension, so BIC comparisons reduce to comparing log-likelihoods and the
#' unconstrained model can never lose. The proper test of
#' H0 "the constraints hold" against H1 "unconstrained" uses the
#' likelihood-ratio statistic \eqn{LR = 2(\ell_{UN} - \ell_{constrained})},
#' whose asymptotic null distribution is not a chi-square but a
#' *chi-bar-squared* distribution - a mixture of chi-squares with
#' data-dependent mixing weights that depend on which inequality constraints
#' are active at the truth. Rather than deriving those weights, this function
#' calibrates the test by parametric bootstrap: simulate `B` datasets from
#' the fitted constrained model (the null), refit both models to each,
#' and compare the observed statistic to the resulting null distribution.
#' Each bootstrap statistic is floored at zero, as is the observed statistic
#' (a negative value can only arise from Monte Carlo optimization noise and
#' triggers a warning beyond a small tolerance).
#'
#' The same machinery applies to LCR vs DM, where LCR is nested with
#' *fewer* parameters (equality constraints from the Rasch structure plus
#' the ordering inequalities); simulating from the fitted LCR again gives a
#' correctly calibrated null distribution.
#'
#' Replicates on which either refit fails are dropped; a warning is issued
#' if more than 10 percent are lost. The bootstrap is parallelised over
#' replicates when `mc.cores > 1` (see that argument); each replicate is
#' seeded independently, so the parallel and serial results are identical.
#'
#' @seealso [select_model_ll()] for the full ordered selection procedure,
#'   [select_model_constraint()] for the heuristic diagnostic selector.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 500; J <- 8; C <- 3
#' prob <- plogis(t(apply(matrix(runif(J * C, -4, 4), J, C), 1, sort)))
#' cls <- sample(1:C, n, replace = TRUE)
#' dat <- matrix(rbinom(n * J, 1, t(prob)[cls, ]), n, J)
#' f_un <- fit_un(dat, C, seed = 1)
#' f_mon <- fit_mon(dat, C, seed = 1)
#' ll_equivalence_test(dat, f_mon, f_un, B = 99, seed = 2)
#' }
#' @export
ll_equivalence_test <- function(data, fit_constrained, fit_un,
                                B = 99, n_starts = 3, seed = NULL,
                                use_cpp = TRUE, mc.cores = 1L,
                                verbose = FALSE) {

  if (!inherits(fit_constrained, "qlfit") || !inherits(fit_un, "qlfit")) {
    stop("fit_constrained and fit_un must be qlfit objects")
  }
  supported <- c("UN", "MON", "IIO", "DM", "LCR")
  if (!fit_constrained$model_type %in% supported ||
      !fit_un$model_type %in% supported) {
    stop("ll_equivalence_test supports UN, MON, IIO, DM, and LCR fits ",
         "(RM uses different estimation machinery; compare it by BIC)")
  }

  data <- validate_data(data)
  n_obs <- nrow(data)
  n_classes <- fit_constrained$n_classes
  if (!is.na(fit_un$n_classes) && fit_un$n_classes != n_classes) {
    stop("fit_constrained and fit_un must use the same number of classes")
  }

  type_c <- fit_constrained$model_type
  type_g <- fit_un$model_type

  # Observed statistic, guarded against optimization noise
  lr_obs <- 2 * (fit_un$loglik - fit_constrained$loglik)
  if (lr_obs < -1e-6) {
    warning("Observed LR statistic is negative (",
            formatC(lr_obs, format = "g", digits = 4),
            "): the general model has a lower log-likelihood than the ",
            "constrained model. Treating the statistic as 0; consider ",
            "increasing n_starts for the observed fits.")
  }
  lr_obs <- max(0, lr_obs)

  # Per-replicate seeds so the procedure is deterministic given `seed`
  if (!is.null(seed)) set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, B)

  boot_one <- function(b) {
    set.seed(rep_seeds[b])
    boot_data <- simulate_from_qlfit(fit_constrained, n_obs)

    fit_c_star <- tryCatch(
      suppressWarnings(
        refit_model_type(type_c, boot_data, n_classes, n_starts, use_cpp)),
      error = function(e) NULL)
    fit_g_star <- tryCatch(
      suppressWarnings(
        refit_model_type(type_g, boot_data, n_classes, n_starts, use_cpp)),
      error = function(e) NULL)

    if (verbose && mc.cores == 1L && b %% 10 == 0) {
      cat("  bootstrap replicate", b, "of", B, "\n")
    }

    if (is.null(fit_c_star) || is.null(fit_g_star)) return(NA_real_)
    max(0, 2 * (fit_g_star$loglik - fit_c_star$loglik))
  }

  raw <- par_lapply(seq_len(B), boot_one, mc.cores)
  # coerce any parallel worker failures (try-error, wrong length) to NA
  lr_star <- vapply(raw, function(z) {
    if (is.numeric(z) && length(z) == 1L) z else NA_real_
  }, numeric(1))
  n_failed <- sum(is.na(lr_star))
  lr_star <- lr_star[!is.na(lr_star)]
  b_eff <- length(lr_star)

  if (b_eff == 0L) {
    stop("All ", B, " bootstrap refits failed; cannot compute a p-value")
  }
  if (n_failed > 0.1 * B) {
    warning(n_failed, " of ", B, " bootstrap replicates failed to refit ",
            "and were dropped; the p-value is based on B_effective = ",
            b_eff, " replicates")
  }

  p_value <- (1 + sum(lr_star >= lr_obs)) / (b_eff + 1)

  structure(
    list(
      statistic = lr_obs,
      p_value = p_value,
      B = B,
      B_effective = b_eff,
      n_failed = n_failed,
      null_distribution = sort(lr_star),
      models = c(type_c, type_g)
    ),
    class = "qleqtest"
  )
}

#' Print method for qleqtest objects
#'
#' @param x A qleqtest object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#' @export
print.qleqtest <- function(x, ...) {
  cat("\nParametric bootstrap LL-equivalence test\n")
  cat("-----------------------------------------\n")
  cat("H0: ", x$models[1], " constraints hold   vs   H1: ",
      x$models[2], "\n", sep = "")
  cat(sprintf("LR statistic : %.4f\n", x$statistic))
  cat(sprintf("Bootstrap p  : %.4f\n", x$p_value))
  cat(sprintf("Replicates   : %d effective (of %d requested",
              x$B_effective, x$B))
  if (x$n_failed > 0) cat(sprintf("; %d failed", x$n_failed))
  cat(")\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Ordered selection procedure
# ---------------------------------------------------------------------------

#' Ordered latent structure selection via bootstrap LL-equivalence tests
#'
#' Selects among the six latent structure models of Torres Irribarra &
#' Diakow (UN, MON, IIO, DM, LCR, RM) by walking the hierarchy from the
#' most to the least flexible model, using parametric-bootstrap
#' likelihood-ratio tests where parameter counts coincide or inequality
#' constraints make standard asymptotics invalid, and BIC only where
#' parameter counts genuinely differ.
#'
#' @param data Binary response matrix (persons x items).
#' @param n_classes Number of latent classes for the discrete models. May be
#'   a single integer, or a vector (e.g. `1:6`), in which case the class count
#'   is selected first by BIC with [select_n_classes()] and the structural
#'   comparison is run at the chosen count (the enumeration table is returned
#'   as `n_classes_table`).
#' @param alpha Significance level for the bootstrap tests (default 0.05).
#'   A constrained model is deemed *adequate* when its bootstrap p-value
#'   exceeds `alpha`, i.e. the constraints are rejected when
#'   \eqn{p \le \alpha}. With the p-value estimator
#'   \eqn{(1 + \#\{LR^* \ge LR\})/(B+1)} this is an exact level-\eqn{\alpha}
#'   Monte Carlo test whenever \eqn{\alpha (B + 1)} is an integer
#'   (e.g. `B = 99`, `alpha = 0.05`).
#' @param B Number of bootstrap replicates per test (default 99).
#' @param n_starts Number of random starts for the observed-data fits
#'   (default 5).
#' @param boot_n_starts Number of random starts for each bootstrap refit
#'   (default 3).
#' @param method How the ordinal layer is tested. `"joint"` (default) tests
#'   the doubly-monotone model directly against UN and only falls back to the
#'   single-constraint models on rejection. `"lattice"` instead tests each
#'   constraint edge separately and accepts the deepest model reachable by
#'   non-rejected edges (MON vs UN, IIO vs UN, then the DM increments
#'   DM vs IIO / DM vs MON). The lattice method is more principled - it gives
#'   each constraint family its own targeted test - but in a paired
#'   simulation study (n = 1000, J = 10, C = 3, 30 replicates per generating
#'   model) the two methods were statistically indistinguishable on both
#'   overall recovery and the rate of scale-type errors, while lattice costs
#'   roughly twice the computation (up to five bootstrap tests per data set
#'   rather than one to three). `"joint"` is therefore the default;
#'   `"lattice"` is available for users who prefer the edge-wise procedure.
#' @param seed Optional integer seed; makes the whole procedure
#'   deterministic.
#' @param use_cpp Use the compiled C++ EM engine (default TRUE).
#' @param mc.cores Number of cores for the bootstrap refits within each test
#'   (default 1, passed to [ll_equivalence_test()]). Because replicates are
#'   seeded independently the result is identical to the serial run.
#' @param verbose Print progress messages (default FALSE).
#' @param ... Further arguments passed to the latent class fitting functions
#'   ([fit_un()], [fit_mon()], [fit_iio()], [fit_dm()], [fit_lcr()]).
#'
#' @return An object of class `qlselect_ll`: a list with elements
#'   \describe{
#'     \item{selected}{Label of the selected model.}
#'     \item{interpretation}{Character: "CLASSIFICATORY", "ORDINAL (...)",
#'       or "QUANTITATIVE (...)".}
#'     \item{tests}{Data frame of the bootstrap tests performed
#'       (comparison, statistic, p_value, decision).}
#'     \item{bics}{Named numeric BIC values for LCR and RM (NA when that
#'       step was not reached or a fit failed).}
#'     \item{fits}{Named list of the fitted models (failed fits are NULL).}
#'     \item{alpha, B}{The settings used.}
#'   }
#'
#' @details
#' The procedure walks the paper's hierarchy UN -> MON -> IIO -> DM ->
#' LCR -> RM in three steps:
#'
#' \enumerate{
#'   \item \strong{Ordinal structure.} The way this layer is tested depends
#'     on `method`.
#'
#'     With `method = "joint"` (the default), DM (the most constrained
#'     ordinal model) is tested against UN directly with
#'     [ll_equivalence_test()]. If DM is adequate (p > `alpha`) it becomes
#'     the ordinal candidate and the procedure continues to step 2.
#'     Otherwise IIO and MON are each tested against UN; among those that are
#'     adequate, the one with the \emph{higher log-likelihood} is selected
#'     (IIO and MON are not nested in each other; with equal parameter counts
#'     the higher log-likelihood is the natural tie-break). If neither is
#'     adequate the unconstrained model is selected (CLASSIFICATORY).
#'
#'     With `method = "lattice"`, each edge of the constraint lattice is
#'     tested instead: MON vs UN and IIO vs UN establish whether each single
#'     constraint holds, and DM is reached only when one of them holds
#'     \emph{and} the corresponding increment edge is not rejected - DM vs
#'     IIO tests the added monotonicity given invariant ordering, DM vs MON
#'     tests the added ordering given monotonicity. The deepest model
#'     reachable by a path of non-rejected edges is selected; if both single
#'     constraints hold but the DM increment is rejected from either parent,
#'     the better-fitting single ordinal model is kept. This gives each
#'     constraint family its own targeted test at roughly twice the
#'     computational cost, with recovery performance that matched the joint
#'     method in simulation.
#'   \item \strong{Discrete quantitative structure} (reached only when DM
#'     was adequate). LCR is tested against DM the same way: LCR is nested
#'     in DM with strictly fewer parameters (the Rasch structure imposes
#'     equality constraints on top of the orderings), so
#'     \eqn{LR = 2(\ell_{DM} - \ell_{LCR})} is bootstrapped by simulating
#'     from the fitted LCR and refitting both models. If LCR is adequate it
#'     becomes the quantitative candidate; otherwise DM is selected
#'     (ORDINAL interpretation).
#'   \item \strong{Continuous quantitative structure.} RM (continuous) and
#'     LCR (discrete) are non-nested and estimated with different machinery
#'     (mirt marginal maximum likelihood over a continuous latent density vs
#'     the EM algorithm over discrete classes), so a likelihood-ratio
#'     bootstrap is not clean. Instead the BIC difference
#'     \eqn{BIC_{LCR} - BIC_{RM}} is calibrated against a null simulated from
#'     the fitted RM: the discrete LCR is selected (QUANTITATIVE, discrete)
#'     only when it fits better than the RM null would produce by chance,
#'     otherwise the continuous Rasch model is selected (QUANTITATIVE,
#'     continuous). This keeps the final step consistent with the
#'     bootstrap-calibrated tests above rather than a raw BIC comparison.
#' }
#'
#' \strong{Why bootstrap rather than chi-square?} UN, MON, IIO, and DM all
#' have the same number of free parameters; the ordinal models differ from
#' UN only through inequality constraints. The LR statistic therefore has a
#' chi-bar-squared null distribution (a mixture of chi-squares with
#' data-dependent weights), not a fixed chi-square, and information criteria
#' cannot separate the models at all. The parametric bootstrap calibrates
#' the test without deriving the mixture weights.
#'
#' \strong{Computational cost.} Each bootstrap test refits two models on
#' each of `B` datasets. With the compiled EM engine a single fit takes on
#' the order of tens of milliseconds at moderate n, so the default `B = 99`
#' costs a few seconds per test and the full procedure typically runs 2-3
#' tests. Supplying `seed` makes every fit and every bootstrap draw
#' reproducible.
#'
#' @references
#' Torres Irribarra, D., & Diakow, R. Categorization, Ordering and
#' Quantification: Selecting a Latent Variable Model by Comparing Latent
#' Structures.
#'
#' @seealso [ll_equivalence_test()] for a single comparison,
#'   [select_model_constraint()] for the heuristic diagnostic selector,
#'   [compare_models()] and [successive_comparison()] for pure
#'   information-criterion comparisons.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' theta <- rnorm(500)
#' beta <- seq(-1.5, 1.5, length.out = 8)
#' dat <- matrix(rbinom(500 * 8, 1, plogis(outer(theta, beta, "-"))), 500, 8)
#' sel <- select_model_ll(dat, n_classes = 3, B = 99, seed = 2)
#' print(sel)
#' }
#' @export
select_model_ll <- function(data, n_classes, alpha = 0.05, B = 99,
                            n_starts = 5, boot_n_starts = 3,
                            method = c("joint", "lattice"), seed = NULL,
                            use_cpp = TRUE, mc.cores = 1L, verbose = FALSE,
                            ...) {

  method <- match.arg(method)
  data <- validate_data(data)

  # Optional first stage: if a range of class counts is supplied, select the
  # number of classes by BIC before the structural comparison. The structural
  # models require C >= 2, so the comparison uses the best fittable C >= 2;
  # a globally preferred C = 1 is surfaced as a warning.
  n_classes_table <- NULL
  if (length(n_classes) > 1L) {
    if (verbose) cat("Selecting the number of classes over C =",
                     paste(range(n_classes), collapse = "-"), "...\n")
    nc <- select_n_classes(data, C_range = n_classes, n_starts = n_starts,
                           use_cpp = use_cpp, mc.cores = mc.cores,
                           seed = seed, verbose = verbose, ...)
    n_classes_table <- nc$table
    cand <- nc$table[nc$table$C >= 2L & !is.na(nc$table$BIC), ]
    if (nrow(cand) == 0L) {
      stop("No fittable class count C >= 2 in the requested range")
    }
    n_classes <- cand$C[which.min(cand$BIC)]
    if (!is.na(nc$best_C) && nc$best_C == 1L) {
      warning("A single-class model (C = 1) had the lowest BIC; the latent ",
              "structure comparison is conditional on C = ", n_classes,
              ", but the data show little evidence of multiple classes.")
    }
    if (verbose) cat("Using n_classes =", n_classes,
                     "for the structural comparison\n")
  }

  fit_safe <- function(label, expr) {
    fit <- tryCatch(suppressWarnings(expr), error = function(e) {
      if (verbose) cat("Fit of", label, "failed:", conditionMessage(e), "\n")
      NULL
    })
    fit
  }

  if (verbose) cat("Fitting the six candidate models...\n")
  fits <- list(
    UN  = fit_safe("UN",  fit_un(data, n_classes, n_starts = n_starts,
                                 use_cpp = use_cpp, seed = seed, ...)),
    MON = fit_safe("MON", fit_mon(data, n_classes, n_starts = n_starts,
                                  use_cpp = use_cpp, seed = seed, ...)),
    IIO = fit_safe("IIO", fit_iio(data, n_classes, n_starts = n_starts,
                                  use_cpp = use_cpp, seed = seed, ...)),
    DM  = fit_safe("DM",  fit_dm(data, n_classes, n_starts = n_starts,
                                 use_cpp = use_cpp, seed = seed, ...)),
    LCR = fit_safe("LCR", fit_lcr(data, n_classes, n_starts = n_starts,
                                  use_cpp = use_cpp, seed = seed, ...)),
    RM  = fit_safe("RM",  fit_rm(data, verbose = FALSE))
  )

  if (is.null(fits$UN)) {
    stop("The unconstrained (UN) fit failed; the selection procedure ",
         "requires it as the reference model")
  }

  tests <- data.frame(comparison = character(0), statistic = numeric(0),
                      p_value = numeric(0), decision = character(0),
                      stringsAsFactors = FALSE)
  add_test <- function(tests, comparison, test, adequate) {
    rbind(tests, data.frame(
      comparison = comparison,
      statistic = test$statistic,
      p_value = test$p_value,
      decision = if (adequate) "adequate (p > alpha)" else "rejected (p <= alpha)",
      stringsAsFactors = FALSE))
  }

  run_test <- function(fit_c, fit_g, seed_offset) {
    ll_equivalence_test(
      data, fit_c, fit_g, B = B, n_starts = boot_n_starts,
      seed = if (!is.null(seed)) seed + seed_offset else NULL,
      use_cpp = use_cpp, mc.cores = mc.cores, verbose = verbose)
  }

  bics <- c(LCR = NA_real_, RM = NA_real_)
  selected <- NULL
  interpretation <- NULL

  # Run one edge test, record it, and return whether the constrained model
  # is adequate (p > alpha). NULL fits yield NA (no test performed).
  test_adequate <- function(fit_c, fit_g, label, seed_offset) {
    if (is.null(fit_c) || is.null(fit_g)) return(NA)
    t <- run_test(fit_c, fit_g, seed_offset)
    ok <- t$p_value > alpha
    tests <<- add_test(tests, label, t, ok)
    ok
  }

  # -- Step 1: ordinal structure -------------------------------------------
  # Sets proceed_to_lcr = TRUE when the doubly-monotone model is supported;
  # otherwise it terminates the selection in the ordinal/classificatory layer
  # by setting `selected` and `interpretation`.
  proceed_to_lcr <- FALSE

  if (method == "joint") {
    # Test the most constrained ordinal model (DM) directly against UN;
    # on rejection fall back to the single-constraint models.
    if (verbose) cat("Step 1 (joint): testing DM vs UN...\n")
    dm_adequate <- isTRUE(test_adequate(fits$DM, fits$UN, "DM vs UN", 1000L))

    if (dm_adequate) {
      proceed_to_lcr <- TRUE
    } else {
      candidates <- list()
      if (verbose) cat("Step 1 (joint): testing IIO vs UN and MON vs UN...\n")
      if (isTRUE(test_adequate(fits$IIO, fits$UN, "IIO vs UN", 2000L))) {
        candidates$IIO <- fits$IIO
      }
      if (isTRUE(test_adequate(fits$MON, fits$UN, "MON vs UN", 3000L))) {
        candidates$MON <- fits$MON
      }
      if (length(candidates) == 0L) {
        selected <- "UN"
        interpretation <- "CLASSIFICATORY (no ordinal structure supported)"
      } else {
        lls <- vapply(candidates, function(f) f$loglik, numeric(1))
        selected <- names(candidates)[which.max(lls)]
        interpretation <- if (selected == "IIO") {
          "ORDINAL (invariant item ordering)"
        } else {
          "ORDINAL (class monotonicity)"
        }
      }
    }
  } else {
    # Lattice: test each constraint edge and accept the deepest model
    # reachable by a path of non-rejected edges. The single-constraint
    # edges (MON vs UN, IIO vs UN) are tested first; DM is reached only when
    # one of them holds AND the corresponding increment edge (DM vs IIO for
    # the added monotonicity, or DM vs MON for the added item ordering) is
    # also not rejected. This gives each constraint family its own targeted
    # test instead of diluting the signal across a joint DM-vs-UN test.
    if (verbose) cat("Step 1 (lattice): testing constraint edges...\n")
    mon_ok <- isTRUE(test_adequate(fits$MON, fits$UN, "MON vs UN", 3000L))
    iio_ok <- isTRUE(test_adequate(fits$IIO, fits$UN, "IIO vs UN", 2000L))

    reach_dm <- FALSE
    if (!is.null(fits$DM)) {
      if (iio_ok) {
        # MON increment given IIO
        reach_dm <- isTRUE(test_adequate(fits$DM, fits$IIO, "DM vs IIO", 5000L))
      }
      if (!reach_dm && mon_ok) {
        # IIO increment given MON
        reach_dm <- isTRUE(test_adequate(fits$DM, fits$MON, "DM vs MON", 6000L))
      }
    }

    if (reach_dm) {
      proceed_to_lcr <- TRUE
    } else if (iio_ok && mon_ok) {
      # Both single constraints hold individually but the doubly-monotone
      # increment is rejected from either parent: keep the better-fitting
      # single ordinal model.
      selected <- if (fits$IIO$loglik >= fits$MON$loglik) "IIO" else "MON"
      interpretation <- if (selected == "IIO") {
        "ORDINAL (invariant item ordering)"
      } else {
        "ORDINAL (class monotonicity)"
      }
    } else if (iio_ok) {
      selected <- "IIO"
      interpretation <- "ORDINAL (invariant item ordering)"
    } else if (mon_ok) {
      selected <- "MON"
      interpretation <- "ORDINAL (class monotonicity)"
    } else {
      selected <- "UN"
      interpretation <- "CLASSIFICATORY (no ordinal structure supported)"
    }
  }

  # -- Steps 2-3: quantitative structure (only when DM is supported) --------
  if (proceed_to_lcr) {
    # Step 2: discrete quantitative structure (LCR vs DM)
    if (verbose) cat("Step 2: testing LCR vs DM...\n")
    lcr_adequate <- isTRUE(test_adequate(fits$LCR, fits$DM, "LCR vs DM", 4000L))

    if (!lcr_adequate) {
      selected <- "DM"
      interpretation <- "ORDINAL (double monotonicity)"
    } else {
      # Step 3: continuous (RM) vs discrete (LCR), bootstrap-calibrated on the
      # BIC difference against a null simulated from the fitted RM.
      bics["LCR"] <- BIC(fits$LCR)
      if (!is.null(fits$RM)) bics["RM"] <- BIC(fits$RM)
      if (is.null(fits$RM)) {
        selected <- "LCR"
        interpretation <- "QUANTITATIVE (discrete: latent class Rasch)"
      } else {
        if (verbose) cat("Step 3: testing RM vs LCR (BIC bootstrap)...\n")
        rl <- rm_vs_lcr_test(data, fits$RM, fits$LCR, n_classes, B = B,
                             alpha = alpha, n_starts = boot_n_starts,
                             use_cpp = use_cpp, mc.cores = mc.cores,
                             seed = if (!is.null(seed)) seed + 7000L else NULL)
        tests <- rbind(tests, data.frame(
          comparison = "RM vs LCR", statistic = rl$statistic,
          p_value = rl$p_value,
          decision = if (!isTRUE(rl$available)) "comparison unavailable"
                     else if (rl$select_lcr) "discrete (LCR) supported"
                     else "continuous (RM) preferred",
          stringsAsFactors = FALSE))
        if (!isTRUE(rl$available)) {
          # bootstrap failed; fall back to the raw BIC comparison
          if (bics["RM"] <= bics["LCR"]) {
            selected <- "RM"
            interpretation <- "QUANTITATIVE (continuous: Rasch model)"
          } else {
            selected <- "LCR"
            interpretation <- "QUANTITATIVE (discrete: latent class Rasch)"
          }
        } else if (rl$select_lcr) {
          selected <- "LCR"
          interpretation <- "QUANTITATIVE (discrete: latent class Rasch)"
        } else {
          selected <- "RM"
          interpretation <- "QUANTITATIVE (continuous: Rasch model)"
        }
      }
    }
  }

  structure(
    list(
      selected = selected,
      interpretation = interpretation,
      tests = tests,
      bics = bics,
      fits = fits,
      alpha = alpha,
      B = B,
      method = method,
      n_classes = n_classes,
      n_classes_table = n_classes_table
    ),
    class = "qlselect_ll"
  )
}

#' Print method for qlselect_ll objects
#'
#' @param x A qlselect_ll object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#' @export
print.qlselect_ll <- function(x, ...) {
  cat("\nBootstrap LL-equivalence latent structure selection\n")
  cat("---------------------------------------------------\n")
  cat("Selected model :", x$selected, "\n")
  cat("Interpretation :", x$interpretation, "\n")
  if (!is.null(x$n_classes)) {
    cat("Latent classes :", x$n_classes,
        if (!is.null(x$n_classes_table)) "(selected by BIC)" else "", "\n")
  }
  cat("Method =", if (is.null(x$method)) "joint" else x$method,
      "  Alpha =", x$alpha, "  Bootstrap replicates per test =", x$B, "\n\n")

  if (nrow(x$tests) > 0) {
    cat("Decision path (bootstrap LR tests):\n")
    for (i in seq_len(nrow(x$tests))) {
      cat(sprintf("  %-10s LR = %8.3f   p = %.4f   -> %s\n",
                  x$tests$comparison[i], x$tests$statistic[i],
                  x$tests$p_value[i], x$tests$decision[i]))
    }
  } else {
    cat("No bootstrap tests were performed.\n")
  }

  if (any(!is.na(x$bics))) {
    cat("\nBIC comparison at the quantitative step:\n")
    for (m in names(x$bics)) {
      if (!is.na(x$bics[m])) {
        cat(sprintf("  %-4s BIC = %.2f%s\n", m, x$bics[m],
                    if (x$selected == m) "   (selected)" else ""))
      }
    }
  }
  invisible(x)
}
