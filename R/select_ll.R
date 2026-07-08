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
#' @param verbose Print progress every 10 replicates (default FALSE).
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
#' if more than 10 percent are lost. Computation is sequential (`lapply`);
#' users needing parallelism can split `B` across seeds and pool the
#' resulting `null_distribution` vectors.
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
                                use_cpp = TRUE, verbose = FALSE) {

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

    if (verbose && b %% 10 == 0) {
      cat("  bootstrap replicate", b, "of", B, "\n")
    }

    if (is.null(fit_c_star) || is.null(fit_g_star)) return(NA_real_)
    max(0, 2 * (fit_g_star$loglik - fit_c_star$loglik))
  }

  lr_star <- unlist(lapply(seq_len(B), boot_one))
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
#' @param n_classes Number of latent classes for the discrete models.
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
#' @param seed Optional integer seed; makes the whole procedure
#'   deterministic.
#' @param use_cpp Use the compiled C++ EM engine (default TRUE).
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
#'   \item \strong{Ordinal structure.} DM (double monotonicity, the most
#'     constrained ordinal model) is tested against UN with
#'     [ll_equivalence_test()]. If DM is adequate (p > `alpha`), it becomes
#'     the ordinal candidate and the procedure continues to step 2.
#'     Otherwise IIO and MON are each tested against UN; among those that
#'     are adequate, the one with the \emph{higher log-likelihood} is
#'     selected (IIO and MON are not nested in each other, so they cannot
#'     be tested against one another; with equal parameter counts the
#'     higher log-likelihood is the natural tie-break). If neither is
#'     adequate the unconstrained model is selected and the interpretation
#'     is CLASSIFICATORY.
#'   \item \strong{Discrete quantitative structure} (reached only when DM
#'     was adequate). LCR is tested against DM the same way: LCR is nested
#'     in DM with strictly fewer parameters (the Rasch structure imposes
#'     equality constraints on top of the orderings), so
#'     \eqn{LR = 2(\ell_{DM} - \ell_{LCR})} is bootstrapped by simulating
#'     from the fitted LCR and refitting both models. If LCR is adequate it
#'     becomes the quantitative candidate; otherwise DM is selected
#'     (ORDINAL interpretation).
#'   \item \strong{Continuous quantitative structure.} RM is compared with
#'     LCR by BIC. A bootstrap LR across these two is not clean because
#'     they are estimated with different machinery (mirt marginal maximum
#'     likelihood over a continuous latent density vs the EM algorithm over
#'     discrete classes), but BIC is legitimate here because the parameter
#'     counts genuinely differ. If \eqn{BIC_{RM} < BIC_{LCR}} the Rasch
#'     model is selected (QUANTITATIVE, continuous); otherwise LCR
#'     (QUANTITATIVE, discrete).
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
                            n_starts = 5, boot_n_starts = 3, seed = NULL,
                            use_cpp = TRUE, verbose = FALSE, ...) {

  data <- validate_data(data)

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
      use_cpp = use_cpp, verbose = verbose)
  }

  bics <- c(LCR = NA_real_, RM = NA_real_)
  selected <- NULL
  interpretation <- NULL

  # -- Step 1: ordinal structure (DM vs UN, else IIO/MON vs UN) -------------
  dm_adequate <- FALSE
  if (!is.null(fits$DM)) {
    if (verbose) cat("Step 1: testing DM vs UN...\n")
    t_dm <- run_test(fits$DM, fits$UN, 1000L)
    dm_adequate <- t_dm$p_value > alpha
    tests <- add_test(tests, "DM vs UN", t_dm, dm_adequate)
  }

  if (!dm_adequate) {
    candidates <- list()
    if (!is.null(fits$IIO)) {
      if (verbose) cat("Step 1: testing IIO vs UN...\n")
      t_iio <- run_test(fits$IIO, fits$UN, 2000L)
      iio_ok <- t_iio$p_value > alpha
      tests <- add_test(tests, "IIO vs UN", t_iio, iio_ok)
      if (iio_ok) candidates$IIO <- fits$IIO
    }
    if (!is.null(fits$MON)) {
      if (verbose) cat("Step 1: testing MON vs UN...\n")
      t_mon <- run_test(fits$MON, fits$UN, 3000L)
      mon_ok <- t_mon$p_value > alpha
      tests <- add_test(tests, "MON vs UN", t_mon, mon_ok)
      if (mon_ok) candidates$MON <- fits$MON
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
  } else {
    # -- Step 2: discrete quantitative structure (LCR vs DM) ---------------
    lcr_adequate <- FALSE
    if (!is.null(fits$LCR)) {
      if (verbose) cat("Step 2: testing LCR vs DM...\n")
      t_lcr <- run_test(fits$LCR, fits$DM, 4000L)
      lcr_adequate <- t_lcr$p_value > alpha
      tests <- add_test(tests, "LCR vs DM", t_lcr, lcr_adequate)
    }

    if (!lcr_adequate) {
      selected <- "DM"
      interpretation <- "ORDINAL (double monotonicity)"
    } else {
      # -- Step 3: continuous quantitative structure (RM vs LCR by BIC) ----
      bics["LCR"] <- BIC(fits$LCR)
      if (!is.null(fits$RM)) bics["RM"] <- BIC(fits$RM)
      if (!is.na(bics["RM"]) && bics["RM"] < bics["LCR"]) {
        selected <- "RM"
        interpretation <- "QUANTITATIVE (continuous: Rasch model)"
      } else {
        selected <- "LCR"
        interpretation <- "QUANTITATIVE (discrete: latent class Rasch)"
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
      B = B
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
  cat("Alpha =", x$alpha, "  Bootstrap replicates per test =", x$B, "\n\n")

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
