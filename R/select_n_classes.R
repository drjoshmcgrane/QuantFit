#' Select the number of latent classes by information criterion
#'
#' Fits the unconstrained latent class model across a range of class counts
#' and returns an enumeration table together with the class count preferred
#' by the chosen information criterion. This is the standard first stage of a
#' workflow: decide the grain for the ordinal comparison, then decide *what
#' structure* those classes have with [select_model_ll()] or [compare_models()].
#' [select_model_ll()] repeats this UN-BIC enumeration inside its ordinal-edge
#' bootstraps so inference includes the uncertainty introduced by automation.
#' If [select_model_ll()] reaches DM, its DM-to-LCR bridge uses the fixed
#' support-point count required by the LCR/Rasch equivalence result; the range
#' here is then used to profile LCR for the final grain comparison.
#'
#' @param data Binary response matrix (persons x items).
#' @param C_range Integer vector of class counts to consider (default
#'   `1:6`). `C = 1` is the single-class independence baseline (a homogeneous
#'   population with no latent heterogeneity); `C >= 2` are fitted with
#'   [fit_un()].
#' @param n_starts Number of random starts for each unconstrained fit
#'   (default 5).
#' @param criterion Information criterion used to pick the class count,
#'   `"BIC"` (default) or `"AIC"`.
#' @param use_cpp Use the compiled C++ EM engine (default TRUE).
#' @param mc.cores Number of cores. The class counts are fitted concurrently
#'   with [parallel::mclapply()] when `mc.cores > 1` on a non-Windows
#'   platform; each fit is seeded independently so the result is identical to
#'   the serial run.
#' @param seed Optional integer seed; fit for class count `C` uses
#'   `seed + C`, making the whole enumeration reproducible.
#' @param verbose Print progress (default FALSE).
#' @param ... Further arguments passed to [fit_un()].
#'
#' @return An object of class `qlnclasses`: a list with
#'   \describe{
#'     \item{table}{Data frame with one row per class count and columns
#'       `C`, `loglik`, `n_par`, `AIC`, `BIC`, `entropy` (normalised
#'       classification entropy, `NA` for `C = 1`), and `converged`.}
#'     \item{best_C}{The class count minimising `criterion`.}
#'     \item{criterion}{The criterion used.}
#'     \item{fits}{Named list of the fitted `qlfit` objects (`NULL` for
#'       `C = 1` and for any fit that failed), keyed by class count.}
#'   }
#'
#' @details
#' The unconstrained model is used for enumeration because it is the least
#' restrictive of the latent class models and the natural reference for "how
#' many groups are there". Its parameter count is `(C - 1) + C * J`, so BIC
#' genuinely rewards parsimony across class counts (unlike the *structural*
#' comparison at a fixed `C`, where UN, MON, IIO and DM share a parameter
#' count). The single-class model (`C = 1`) has each item at its sample
#' proportion and `J` parameters; comparing it against `C >= 2` tests whether
#' there is any latent heterogeneity at all.
#'
#' A low `best_C` (especially `C = 1`) is itself informative: it indicates the
#' data support little or no class structure, and any subsequent structural
#' comparison should be read in that light.
#'
#' @seealso [select_model_ll()] (which accepts a range for its `n_classes`
#'   argument and calls this function internally), [fit_un()],
#'   [entropy_r2()].
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' theta <- rnorm(600)
#' beta <- seq(-1.5, 1.5, length.out = 8)
#' dat <- matrix(rbinom(600 * 8, 1, plogis(outer(theta, beta, "-"))), 600, 8)
#' nc <- select_n_classes(dat, C_range = 1:6, mc.cores = 2, seed = 1)
#' nc
#' }
#' @export
select_n_classes <- function(data, C_range = 1:6, n_starts = 5,
                             criterion = c("BIC", "AIC"), use_cpp = TRUE,
                             mc.cores = 1L, seed = NULL, verbose = FALSE, ...) {
  criterion <- match.arg(criterion)
  data <- validate_data(data)
  n_obs <- nrow(data)
  C_range <- sort(unique(as.integer(C_range)))
  if (length(C_range) == 0L || any(C_range < 1L)) {
    stop("C_range must be a vector of integers >= 1")
  }

  fit_one_C <- function(C) {
    if (verbose && mc.cores == 1L) cat("Fitting C =", C, "...\n")
    if (C == 1L) {
      # single-class independence baseline: each item at its sample proportion
      p <- pmin(pmax(colMeans(data), 1e-10), 1 - 1e-10)
      s1 <- colSums(data)
      ll <- sum(s1 * log(p) + (n_obs - s1) * log(1 - p))
      npar <- ncol(data)
      return(list(C = 1L, loglik = ll, n_par = as.numeric(npar),
                  AIC = -2 * ll + 2 * npar,
                  BIC = -2 * ll + npar * log(n_obs),
                  entropy = NA_real_, converged = TRUE, fit = NULL))
    }
    f <- tryCatch(
      suppressWarnings(fit_un(data, C, n_starts = n_starts, use_cpp = use_cpp,
                              seed = if (!is.null(seed)) seed + C else NULL,
                              ...)),
      error = function(e) NULL)
    if (is.null(f)) {
      return(list(C = C, loglik = NA_real_, n_par = NA_real_, AIC = NA_real_,
                  BIC = NA_real_, entropy = NA_real_, converged = FALSE,
                  fit = NULL))
    }
    list(C = C, loglik = f$loglik, n_par = as.numeric(f$n_par),
         AIC = AIC(f), BIC = BIC(f),
         entropy = tryCatch(entropy_r2(f), error = function(e) NA_real_),
         converged = isTRUE(f$convergence), fit = f)
  }

  res <- par_lapply(C_range, fit_one_C, mc.cores)

  tab <- data.frame(
    C = vapply(res, function(x) as.integer(x$C), integer(1)),
    loglik = vapply(res, function(x) x$loglik, numeric(1)),
    n_par = vapply(res, function(x) x$n_par, numeric(1)),
    AIC = vapply(res, function(x) x$AIC, numeric(1)),
    BIC = vapply(res, function(x) x$BIC, numeric(1)),
    entropy = vapply(res, function(x) x$entropy, numeric(1)),
    converged = vapply(res, function(x) x$converged, logical(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )

  crit_vals <- tab[[criterion]]
  if (all(is.na(crit_vals))) {
    stop("All class-enumeration fits failed; cannot select a class count")
  }
  best_C <- tab$C[which.min(crit_vals)]

  fits <- stats::setNames(lapply(res, function(x) x$fit), tab$C)

  structure(list(table = tab, best_C = best_C, criterion = criterion,
                 fits = fits),
            class = "qlnclasses")
}

#' Print method for qlnclasses objects
#'
#' @param x A qlnclasses object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#' @export
print.qlnclasses <- function(x, ...) {
  cat("\nLatent class enumeration (unconstrained model)\n")
  cat("----------------------------------------------\n")
  tab <- x$table
  disp <- data.frame(
    C = tab$C,
    logLik = round(tab$loglik, 1),
    n_par = tab$n_par,
    AIC = round(tab$AIC, 1),
    BIC = round(tab$BIC, 1),
    entropy = ifelse(is.na(tab$entropy), "-", formatC(tab$entropy, digits = 3, format = "f")),
    conv = ifelse(tab$converged, "yes", "NO"),
    stringsAsFactors = FALSE
  )
  disp$best <- ifelse(tab$C == x$best_C, paste0("  <- min ", x$criterion), "")
  print(disp, row.names = FALSE)
  cat("\nSelected C =", x$best_C, "by", x$criterion, "\n")
  invisible(x)
}
