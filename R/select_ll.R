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
  if (isTRUE(fit$polytomous)) {
    # item_probs is a list of C x (m_j + 1) category-probability matrices;
    # draw each person's class, then a category per item from that class
    J <- length(fit$item_probs)
    out <- matrix(0L, n_obs, J)
    for (j in seq_len(J)) {
      P <- fit$item_probs[[j]]
      m1 <- ncol(P)
      for (c in seq_len(nrow(P))) {
        who <- which(cls == c)
        if (length(who)) {
          out[who, j] <- sample.int(m1, length(who), replace = TRUE,
                                    prob = P[c, ]) - 1L
        }
      }
    }
    return(out)
  }
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
  rf <- attr(fit, "rm_fit")
  v <- if (is.null(rf)) NA_real_ else rf$sigma^2
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
                           C_range = NULL, boot_reselect = TRUE,
                           mc.cores = 1L, seed = NULL) {
  n <- nrow(data); J <- ncol(data)
  # marginal parametric bootstrap: redraw abilities from N(0, sigma^2).
  # Deliberately NORMAL (not the empirical latent density used by the CC/Kara
  # nulls): this step tests continuous-normal (RM) against discrete (LCR)
  # latent structure, so the latent distribution's shape IS the hypothesis
  # under test, not a nuisance - an empirical density would absorb the very
  # discreteness the test exists to detect.
  poly <- isTRUE(fit_rm_obj$polytomous)
  sigma <- sqrt(rasch_latent_var(fit_rm_obj))
  if (poly) {
    rmf <- attr(fit_rm_obj, "rm_fit")
    sim_fn <- function() .simulate_pcm(rmf, n, sigma)
  } else {
    beta <- fit_rm_obj$delta
    sim_fn <- function() {
      theta <- rnorm(n, 0, sigma)
      matrix(rbinom(n * J, 1, plogis(outer(theta, beta, "-"))), n, J)
    }
  }
  if (!is.null(seed)) set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, B)
  # Post-selection calibration: the OBSERVED statistic uses a class count
  # selected on the observed data (which flatters LCR there); holding that C
  # fixed in the null replicates leaves the null BIC-differences too large and
  # honest RM advantages land in the lower tail (external review; empirically
  # 7/34 true-Rasch datasets mislabelled LCR). With boot_reselect the
  # replicate re-runs the class-count selection over a neighbourhood grid
  # (C_hat +/- 1 within C_range; LCR BIC criterion) so the null carries the
  # selection variability of the two-stage statistic.
  # FULL fixed-range profiling (no data-dependent grid): the statistic is
  # min over C in C_range of BIC(LCR_C) minus BIC(RM) - "best discrete Rasch
  # mixture vs continuous Rasch" - computed identically on observed and null
  # data, so the bootstrap is calibrated without per-replicate re-selection.
  grid <- if (boot_reselect && !is.null(C_range)) {
    sort(unique(pmax(2L, C_range)))
  } else n_classes
  # PROFILED statistic, identical functional on observed and null data:
  # min over the grid of BIC(LCR_C) minus BIC(RM). Profiling only in the
  # null (earlier version) made observed and bootstrap statistics differ -
  # not a calibrated post-selection bootstrap (external review).
  obs_bics <- vapply(grid, function(C) {
    if (C == n_classes) return(BIC(fit_lcr_obj))
    f <- tryCatch(suppressWarnings(fit_lcr(data, C, n_starts = n_starts,
                                           use_cpp = use_cpp)),
                  error = function(e) NULL)
    if (is.null(f)) NA_real_ else BIC(f)
  }, numeric(1))
  obs <- min(obs_bics, na.rm = TRUE) - BIC(fit_rm_obj)  # negative favours LCR
  boot_one <- function(b) {
    set.seed(rep_seeds[b])
    d <- .impose_mask(sim_fn(), data)   # null replicates share the observed missingness
    rm_b <- tryCatch(suppressWarnings(fit_rm(d, verbose = FALSE)),
                     error = function(e) NULL)
    fits <- lapply(grid, function(C) tryCatch(
      suppressWarnings(fit_lcr(d, C, n_starts = n_starts, use_cpp = use_cpp)),
      error = function(e) NULL))
    bics <- vapply(fits, function(f) if (is.null(f)) NA_real_ else BIC(f), numeric(1))
    if (is.null(rm_b) || all(is.na(bics))) return(NA_real_)
    min(bics, na.rm = TRUE) - BIC(rm_b)
  }
  raw <- par_lapply(seq_len(B), boot_one, mc.cores)
  null <- vapply(raw, function(z)
    if (is.numeric(z) && length(z) == 1L) z else NA_real_, numeric(1))
  null <- null[!is.na(null)]
  # lower-tail p (LCR favoured) with the (1 + .)/(B + 1) correction
  p_lower <- if (length(null) == 0L) NA_real_
             else (1 + sum(null <= obs)) / (length(null) + 1)
  list(statistic = obs, p_value = p_lower, null = sort(null),
       grid = grid, obs_bics = obs_bics,
       profiled_C = grid[which.min(obs_bics)],
       profiled_bic = min(obs_bics, na.rm = TRUE),
       B_effective = length(null), B_failed = B - length(null),
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
#' @param reselect_C_range When non-NULL (a class-count range), every null
#'   replicate repeats UN-BIC class selection over this range before both
#'   models are refit (post-selection calibration; used by the quantitative
#'   gate). Default NULL: fixed class count.
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
                                reselect_C_range = NULL,
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

  data <- validate_data_any(data, allow_na = TRUE)
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
    boot_data <- .impose_mask(simulate_from_qlfit(fit_constrained, n_obs), data)

    C_b <- n_classes
    fit_g_star <- NULL
    if (!is.null(reselect_C_range) && length(reselect_C_range) > 1L) {
      # POST-SELECTION CALIBRATION (external review): the observed statistic
      # is two-stage - C selected by UN BIC on the observed data, then the
      # LR computed at that C. Repeat the SAME selection inside every
      # replicate so the null carries the selection variability; holding the
      # observed C fixed made the null too favourable to the general model
      # (the RM->DM leak at the quantitative gate).
      g_fits <- lapply(reselect_C_range, function(C) tryCatch(
        suppressWarnings(refit_model_type(type_g, boot_data, C,
                                          n_starts, use_cpp)),
        error = function(e) NULL))
      g_bics <- vapply(g_fits, function(f)
        if (is.null(f)) NA_real_ else BIC(f), numeric(1))
      usable <- which(is.finite(g_bics) & reselect_C_range >= 2L)
      if (length(usable)) {
        pick <- usable[which.min(g_bics[usable])]
        C_b <- reselect_C_range[pick]
        fit_g_star <- g_fits[[pick]]
      }
    }
    fit_c_star <- tryCatch(
      suppressWarnings(
        refit_model_type(type_c, boot_data, C_b, n_starts, use_cpp)),
      error = function(e) NULL)
    if (is.null(fit_g_star)) fit_g_star <- tryCatch(
      suppressWarnings(
        refit_model_type(type_g, boot_data, C_b, n_starts, use_cpp)),
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

#' Reference distribution of the LR statistic under the GENERAL model
#'
#' Companion to [ll_equivalence_test()]: simulates from the fitted *general*
#' model, refits both models to each replicate, and returns the resulting
#' distribution of the LR statistic. This estimates the LR distribution at the
#' fitted alternative - i.e. the test's power at the estimated departure from
#' the constrained family. Because the models are nested, when the constrained
#' model is true the unrestricted MLE lies essentially inside the constrained
#' set and this distribution coincides with the null (no power, negligible
#' effect size); when the truth is genuinely outside the constrained family it
#' shifts right (noncentral). [select_model_ll()] uses the separation between
#' the two distributions as an effect-size check on rejections.
#'
#' @param data Response matrix.
#' @param fit_constrained,fit_general Fitted qlfit objects (the LR is
#'   \eqn{2(\ell_g - \ell_c)} refitted on data simulated from `fit_general`).
#' @param B,n_starts,seed,use_cpp,mc.cores As in [ll_equivalence_test()].
#' @return Sorted numeric vector of LR draws (length <= B; failures dropped),
#'   or `NULL` if every replicate failed.
#' @keywords internal
ll_general_null <- function(data, fit_constrained, fit_general,
                            B = 99, n_starts = 3, seed = NULL,
                            use_cpp = TRUE, mc.cores = 1L,
                            reselect_C_range = NULL) {
  n_obs <- nrow(data)
  n_classes <- fit_constrained$n_classes
  type_c <- fit_constrained$model_type
  type_g <- fit_general$model_type
  if (!is.null(seed)) set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, B)
  boot_one <- function(b) {
    set.seed(rep_seeds[b])
    boot_data <- .impose_mask(simulate_from_qlfit(fit_general, n_obs), data)
    # the severity comparison must use the SAME statistic as the primary
    # null: when the gate reselects C per replicate, so must this
    if (!is.null(reselect_C_range) && length(reselect_C_range) > 1L) {
      g_fits <- lapply(reselect_C_range, function(C) tryCatch(
        suppressWarnings(refit_model_type(type_g, boot_data, C,
                                          n_starts, use_cpp)),
        error = function(e) NULL))
      g_bics <- vapply(g_fits, function(f)
        if (is.null(f)) NA_real_ else BIC(f), numeric(1))
      usable <- which(is.finite(g_bics) & reselect_C_range >= 2L)
      if (length(usable)) {
        pick <- usable[which.min(g_bics[usable])]
        C_b <- reselect_C_range[pick]
        fit_g_star <- g_fits[[pick]]
        fit_c_star <- tryCatch(suppressWarnings(
          refit_model_type(type_c, boot_data, C_b, n_starts, use_cpp)),
          error = function(e) NULL)
        if (is.null(fit_c_star) || is.null(fit_g_star)) return(NA_real_)
        return(max(0, 2 * (fit_g_star$loglik - fit_c_star$loglik)))
      }
    }
    fit_c_star <- tryCatch(suppressWarnings(
      refit_model_type(type_c, boot_data, n_classes, n_starts, use_cpp)),
      error = function(e) NULL)
    fit_g_star <- tryCatch(suppressWarnings(
      refit_model_type(type_g, boot_data, n_classes, n_starts, use_cpp)),
      error = function(e) NULL)
    if (is.null(fit_c_star) || is.null(fit_g_star)) return(NA_real_)
    max(0, 2 * (fit_g_star$loglik - fit_c_star$loglik))
  }
  raw <- par_lapply(seq_len(B), boot_one, mc.cores)
  lr_star <- vapply(raw, function(z)
    if (is.numeric(z) && length(z) == 1L) z else NA_real_, numeric(1))
  lr_star <- lr_star[!is.na(lr_star)]
  if (length(lr_star) == 0L) return(NULL)
  sort(lr_star)
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
#' @param alpha_quant Significance level for the single quantitative gate
#'   (LCR vs UN), default 0.05 (the conventional level; with `B = 99` this is
#'   an exact level-0.05 Monte Carlo test). This governs the one decision that
#'   demotes a fitted quantitative model to the ordinal layer. True
#'   quantitative models are protected from chance upper-tail rejections by
#'   the estimated-power check (a rejection must also separate from the
#'   general model's reference distribution), so the conventional level does
#'   not carry the compounded false-demotion cost it would in a plain
#'   sequential procedure; lower it (e.g. 0.01) to demote the quantitative
#'   model only on stronger evidence, at the cost of letting more
#'   near-additive ordinal data pass as quantitative.
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
select_model_ll <- function(data, n_classes, alpha = 0.05, alpha_quant = 0.05,
                            B = 99, n_starts = 5, boot_n_starts = 3,
                            method = c("joint", "lattice"), seed = NULL,
                            use_cpp = TRUE, mc.cores = 1L, verbose = FALSE,
                            ...) {

  method <- match.arg(method)
  data <- validate_data_any(data, allow_na = TRUE)

  # Optional first stage: if a range of class counts is supplied, select the
  # number of classes by BIC before the structural comparison. The structural
  # models require C >= 2, so the comparison uses the best fittable C >= 2;
  # a globally preferred C = 1 is surfaced as a warning.
  n_classes_table <- NULL
  orig_C_range <- n_classes                 # full range, for bootstrap re-selection
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

  run_test <- function(fit_c, fit_g, seed_offset, reselect_C_range = NULL) {
    ll_equivalence_test(
      data, fit_c, fit_g, B = B, n_starts = boot_n_starts,
      seed = if (!is.null(seed)) seed + seed_offset else NULL,
      use_cpp = use_cpp, mc.cores = mc.cores,
      reselect_C_range = reselect_C_range, verbose = verbose)
  }

  bics <- c(LCR = NA_real_, RM = NA_real_)
  selected <- NULL
  interpretation <- NULL

  # Run one edge test, record it, and return whether the constrained model
  # is adequate. When the LR test rejects (p_c <= a), the rejection must also
  # pass an estimated-power / effect-size check before the constrained model
  # is demoted: the LR distribution is simulated under the FITTED GENERAL
  # model as well. Because the models are nested, if the constrained model is
  # true the unrestricted MLE lies (essentially) inside the constrained set,
  # so this distribution coincides with the null; if the truth is genuinely
  # outside the constrained family, the fitted general model sits at a
  # detectable distance and the distribution shifts right (noncentral). When
  # the two distributions are NOT separated, the test has no power at the
  # estimated departure - the effect size is negligible and the rejection is
  # an upper-tail draw on a null-sized effect (significance without severity)
  # - so parsimony keeps the constrained model. This eliminates the dominant
  # audit error mode (true constrained models demoted by unlucky tail draws)
  # while leaving genuine violations, where the distributions separate
  # decisively, untouched. Separation criterion: median of the general-model
  # LR distribution exceeds the constrained null's 95th percentile.
  # NULL fits yield NA (no test performed).
  last_test_obj <- NULL
  test_adequate <- function(fit_c, fit_g, label, seed_offset, a = alpha,
                            reselect_C_range = NULL) {
    if (is.null(fit_c) || is.null(fit_g)) return(NA)
    t <- run_test(fit_c, fit_g, seed_offset, reselect_C_range)
    ok <- t$p_value > a
    if (!ok) {
      g_null <- ll_general_null(
        data, fit_c, fit_g, B = B, n_starts = boot_n_starts,
        seed = if (!is.null(seed)) seed + seed_offset + 500L else NULL,
        use_cpp = use_cpp, mc.cores = mc.cores,
        reselect_C_range = reselect_C_range)
      if (!is.null(g_null)) {
        separated <- stats::median(g_null) >
          stats::quantile(t$null_distribution, 0.95, names = FALSE)
        if (!separated) ok <- TRUE   # models indistinguishable -> parsimony
      }
    }
    overridden <- ok && t$p_value <= a
    tests <<- add_test(tests, label, t, ok)
    if (overridden && nrow(tests))
      tests$decision[nrow(tests)] <<- "retained (severity override: no detectable effect)"
    t$label <- label
    t$severity_override <- overridden
    if (exists("g_null", inherits = FALSE)) t$severity_null <- g_null
    last_test_obj <<- t
    ok
  }

  # -- Ordinal / nominal layer ---------------------------------------------
  # Identify the most restrictive ORDINAL model supported (UN, MON, IIO, or DM).
  # This is the fallback used when the data do not support the parametric
  # quantitative model tested at the gate below.
  ordinal_selected <- "UN"
  ordinal_interp <- "CLASSIFICATORY (no ordinal structure supported)"

  if (method == "joint") {
    if (verbose) cat("Ordinal layer (joint): testing DM vs UN...\n")
    if (isTRUE(test_adequate(fits$DM, fits$UN, "DM vs UN", 1000L))) {
      ordinal_selected <- "DM"
      ordinal_interp <- "ORDINAL (double monotonicity)"
    } else {
      candidates <- list()
      if (isTRUE(test_adequate(fits$IIO, fits$UN, "IIO vs UN", 2000L)))
        candidates$IIO <- fits$IIO
      if (isTRUE(test_adequate(fits$MON, fits$UN, "MON vs UN", 3000L)))
        candidates$MON <- fits$MON
      if (length(candidates) > 0L) {
        lls <- vapply(candidates, function(f) f$loglik, numeric(1))
        ordinal_selected <- names(candidates)[which.max(lls)]
        ordinal_interp <- if (ordinal_selected == "IIO")
          "ORDINAL (invariant item ordering)" else "ORDINAL (class monotonicity)"
      }
    }
  } else {
    # Lattice: test each constraint edge and accept the deepest ordinal model
    # reachable by a path of non-rejected edges.
    if (verbose) cat("Ordinal layer (lattice): testing constraint edges...\n")
    mon_ok <- isTRUE(test_adequate(fits$MON, fits$UN, "MON vs UN", 3000L))
    iio_ok <- isTRUE(test_adequate(fits$IIO, fits$UN, "IIO vs UN", 2000L))
    reach_dm <- FALSE
    if (!is.null(fits$DM)) {
      if (iio_ok)
        reach_dm <- isTRUE(test_adequate(fits$DM, fits$IIO, "DM vs IIO", 5000L))
      if (!reach_dm && mon_ok)
        reach_dm <- isTRUE(test_adequate(fits$DM, fits$MON, "DM vs MON", 6000L))
    }
    if (reach_dm) {
      ordinal_selected <- "DM"; ordinal_interp <- "ORDINAL (double monotonicity)"
    } else if (iio_ok && mon_ok) {
      ordinal_selected <- if (fits$IIO$loglik >= fits$MON$loglik) "IIO" else "MON"
      ordinal_interp <- if (ordinal_selected == "IIO")
        "ORDINAL (invariant item ordering)" else "ORDINAL (class monotonicity)"
    } else if (iio_ok) {
      ordinal_selected <- "IIO"; ordinal_interp <- "ORDINAL (invariant item ordering)"
    } else if (mon_ok) {
      ordinal_selected <- "MON"; ordinal_interp <- "ORDINAL (class monotonicity)"
    }
  }

  # -- Quantitative gate ----------------------------------------------------
  # Does the parametric latent-class Rasch model fit as well as the fully
  # unconstrained model? Testing LCR directly against UN - a single gate rather
  # than a sequential DM-then-LCR path - avoids compounding the false-rejection
  # rate, so genuinely quantitative data are not lost to the ordinal layer by
  # chance. A separate alpha_quant governs this one demotion decision
  # (the quantitative model is only overturned on strong evidence).
  proceed_to_lcr <- FALSE
  if (!is.null(fits$LCR)) {
    if (verbose) cat("Quantitative gate: testing LCR vs UN...\n")
    # the gate cannot resolve a level finer than the bootstrap p-value floor
    # 1/(B+1); use the achievable threshold so it can still reject on the
    # strongest evidence (observed LR beyond every null draw) at small B.
    a_q <- alpha_quant
    if (1 / (B + 1) > alpha_quant)
      warning("B = ", B, " cannot resolve alpha_quant = ", alpha_quant,
              " (p-value floor 1/(B+1) = ", signif(1/(B+1), 3),
              "); the gate cannot reject at this B - increase B")
    proceed_to_lcr <- isTRUE(test_adequate(fits$LCR, fits$UN, "LCR vs UN", 4000L,
                                           a = a_q,
                                           reselect_C_range = orig_C_range))
  }

  if (!proceed_to_lcr) {
    selected <- ordinal_selected
    interpretation <- ordinal_interp
  } else {
    {
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
                             C_range = orig_C_range,
                             alpha = alpha, n_starts = boot_n_starts,
                             use_cpp = use_cpp, mc.cores = mc.cores,
                             seed = if (!is.null(seed)) seed + 7000L else NULL)
        # keep the returned object consistent with the profiled decision
        if (isTRUE(rl$available) && !is.na(rl$profiled_C) &&
            rl$profiled_C != n_classes) {
          f2 <- tryCatch(suppressWarnings(fit_lcr(data, rl$profiled_C,
                  n_starts = n_starts, use_cpp = use_cpp)),
                  error = function(e) NULL)
          if (!is.null(f2)) {
            fits$LCR <- f2; bics["LCR"] <- BIC(f2); n_classes <- rl$profiled_C
          }
        }
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
      n_classes_table = n_classes_table,
      rm_vs_lcr = if (exists("rl", inherits = FALSE)) rl else NULL,
      quant_gate = if (exists("last_test_obj", inherits = FALSE))
        last_test_obj else NULL
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
