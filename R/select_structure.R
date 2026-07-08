#' Constraint-based latent structure selection
#'
#' Selects among the six latent structure models by checking which ordinal
#' constraints the *unconstrained* latent class estimates satisfy, then using
#' BIC only where parameter counts genuinely differ (UN/MON/IIO/DM vs LCR/RM).
#'
#' Because UN, MON, IIO, and DM have identical parameter counts, information
#' criteria reduce to comparing log-likelihoods, and the unconstrained model
#' can never lose such a comparison. This function instead follows the logic
#' of Torres Irribarra & Diakow: fit the unconstrained model, order its
#' classes by mean response probability, and diagnose which constraints
#' (class monotonicity, invariant item ordering, double monotonicity, Rasch
#' additivity) the estimates satisfy within `tolerance`. Quantitative models
#' (LCR/RM) are then adopted only when the Rasch-structure diagnostic holds
#' and BIC supports them.
#'
#' @param data Binary response matrix (persons x items).
#' @param n_classes Number of latent classes for the discrete models.
#' @param n_starts Number of random starts for the latent class fits.
#' @param tolerance Maximum tolerated proportion of constraint violations
#'   for the MON/IIO diagnostics (default 0.15).
#' @param rasch_tolerance Maximum tolerated coefficient of variation of
#'   inter-class logit differences for the Rasch-additivity diagnostic
#'   (default 0.25).
#' @param seed Optional seed passed to the model fits for reproducibility.
#' @param ... Further arguments passed to [fit_un()], [fit_lcr()].
#'
#' @return An object of class `qlselect`: a list with elements
#'   `selected` (model label), `interpretation` (character),
#'   `constraints` (logical diagnostics for mon/iio/dm/rasch),
#'   `bics` (named numeric for UN/LCR/RM), and `fits` (the fitted UN, LCR,
#'   and RM models).
#'
#' @references
#' Torres Irribarra, D., & Diakow, R. Categorization, Ordering and
#' Quantification: Selecting a Latent Variable Model by Comparing Latent
#' Structures.
#'
#' @seealso [compare_models()], [successive_comparison()] for the
#'   information-criterion-based procedures.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' theta <- rnorm(500)
#' beta <- seq(-1.5, 1.5, length.out = 8)
#' dat <- matrix(rbinom(500 * 8, 1, plogis(outer(theta, beta, "-"))), 500, 8)
#' sel <- select_model_constraint(dat, n_classes = 3)
#' print(sel)
#' }
#' @export
select_model_constraint <- function(data, n_classes, n_starts = 5,
                                    tolerance = 0.15, rasch_tolerance = 0.25,
                                    seed = NULL, ...) {

  un_fit <- tryCatch(
    suppressWarnings(fit_un(data, n_classes, n_starts = n_starts,
                            seed = seed, ...)),
    error = function(e) NULL)

  if (is.null(un_fit)) {
    return(structure(
      list(selected = NA_character_, interpretation = "FIT_FAILED",
           constraints = list(mon = NA, iio = NA, dm = NA, rasch = NA),
           bics = c(UN = NA_real_, LCR = NA_real_, RM = NA_real_),
           fits = list()),
      class = "qlselect"))
  }

  # classes x items, as the diagnostics expect
  prob <- t(un_fit$item_probs)

  mon_ok <- check_mon_constraint(prob, tolerance)
  iio_ok <- check_iio_constraint(prob, tolerance)
  dm_ok <- mon_ok && iio_ok
  rasch_ok <- check_rasch_constraint(prob, rasch_tolerance)

  rm_fit <- tryCatch(suppressWarnings(fit_rm(data, verbose = FALSE)),
                     error = function(e) NULL)
  lcr_fit <- tryCatch(
    suppressWarnings(fit_lcr(data, n_classes, n_starts = n_starts,
                             seed = seed, ...)),
    error = function(e) NULL)

  un_bic <- BIC(un_fit)
  rm_bic <- if (!is.null(rm_fit)) BIC(rm_fit) else Inf
  lcr_bic <- if (!is.null(lcr_fit)) BIC(lcr_fit) else Inf
  best_quant_bic <- min(rm_bic, lcr_bic)

  if (dm_ok && rasch_ok) {
    if (best_quant_bic < un_bic) {
      selected <- if (rm_bic < lcr_bic) "RM" else "LCR"
      interpretation <- "QUANTITATIVE (DM + Rasch structure)"
    } else {
      selected <- "DM"
      interpretation <- "ORDINAL (double monotonicity with Rasch-like structure)"
    }
  } else if (dm_ok) {
    # DM without Rasch additivity: adopt a quantitative model only on a
    # strong BIC margin
    if (best_quant_bic < un_bic - 20) {
      selected <- if (rm_bic < lcr_bic) "RM" else "LCR"
      interpretation <- "QUANTITATIVE (DM + BIC preference)"
    } else {
      selected <- "DM"
      interpretation <- "ORDINAL (double monotonicity)"
    }
  } else if (rasch_ok && best_quant_bic < un_bic) {
    selected <- if (rm_bic < lcr_bic) "RM" else "LCR"
    interpretation <- "QUANTITATIVE (Rasch structure detected)"
  } else if (iio_ok) {
    selected <- "IIO"
    interpretation <- "ORDINAL (invariant item ordering)"
  } else if (mon_ok) {
    selected <- "MON"
    interpretation <- "ORDINAL (class monotonicity)"
  } else {
    selected <- "UN"
    interpretation <- "CLASSIFICATORY (no ordinal structure)"
  }

  structure(
    list(selected = selected,
         interpretation = interpretation,
         constraints = list(mon = mon_ok, iio = iio_ok, dm = dm_ok,
                            rasch = rasch_ok),
         bics = c(UN = un_bic, LCR = lcr_bic, RM = rm_bic),
         fits = list(UN = un_fit, LCR = lcr_fit, RM = rm_fit)),
    class = "qlselect")
}

#' @export
print.qlselect <- function(x, ...) {
  cat("Constraint-based latent structure selection\n")
  cat("-------------------------------------------\n")
  if (is.na(x$selected)) {
    cat("Model fitting failed; no selection made.\n")
    return(invisible(x))
  }
  cat("Selected model :", x$selected, "\n")
  cat("Interpretation :", x$interpretation, "\n\n")
  cat("Constraint diagnostics on UN estimates:\n")
  cat(sprintf("  MON   (class monotonicity)      : %s\n", x$constraints$mon))
  cat(sprintf("  IIO   (invariant item ordering) : %s\n", x$constraints$iio))
  cat(sprintf("  DM    (double monotonicity)     : %s\n", x$constraints$dm))
  cat(sprintf("  Rasch (additive logit structure): %s\n", x$constraints$rasch))
  cat("\nBIC: ", paste(names(x$bics),
                       formatC(x$bics, format = "f", digits = 1),
                       sep = " = ", collapse = ", "), "\n")
  invisible(x)
}

#' Check class monotonicity of a probability matrix
#'
#' Orders classes by mean response probability and reports whether the
#' proportion of adjacent-class monotonicity violations is within tolerance.
#'
#' @param prob Probability matrix (classes x items).
#' @param tolerance Maximum tolerated violation proportion.
#' @return Logical.
#' @keywords internal
check_mon_constraint <- function(prob, tolerance = 0.15) {
  class_order <- order(rowMeans(prob))
  prob_ordered <- prob[class_order, , drop = FALSE]

  n_classes <- nrow(prob_ordered)
  n_items <- ncol(prob_ordered)
  if (n_classes < 2) return(TRUE)

  violations <- 0L
  comparisons <- 0L
  for (j in seq_len(n_items)) {
    for (c in seq_len(n_classes - 1)) {
      comparisons <- comparisons + 1L
      if (prob_ordered[c + 1, j] < prob_ordered[c, j]) {
        violations <- violations + 1L
      }
    }
  }
  violations / comparisons <= tolerance
}

#' Check invariant item ordering of a probability matrix
#'
#' Reports whether the proportion of item-pair order reversals across class
#' pairs (ignoring differences smaller than 0.05) is within tolerance.
#'
#' @inheritParams check_mon_constraint
#' @return Logical.
#' @keywords internal
check_iio_constraint <- function(prob, tolerance = 0.15) {
  n_classes <- nrow(prob)
  n_items <- ncol(prob)
  if (n_classes < 2 || n_items < 2) return(TRUE)

  violations <- 0L
  comparisons <- 0L
  for (c1 in seq_len(n_classes - 1)) {
    for (c2 in (c1 + 1):n_classes) {
      for (i in seq_len(n_items - 1)) {
        for (j in (i + 1):n_items) {
          comparisons <- comparisons + 1L
          diff_c1 <- prob[c1, i] - prob[c1, j]
          diff_c2 <- prob[c2, i] - prob[c2, j]
          if (sign(diff_c1) != sign(diff_c2) &&
              abs(diff_c1) > 0.05 && abs(diff_c2) > 0.05) {
            violations <- violations + 1L
          }
        }
      }
    }
  }
  if (comparisons == 0L) return(TRUE)
  violations / comparisons <= tolerance
}

#' Check Rasch (additive logit) structure of a probability matrix
#'
#' Under a Rasch parameterization, logit(P) is additive: inter-class logit
#' differences are constant across items. This reports whether the mean
#' coefficient of variation of those differences is below tolerance.
#'
#' @inheritParams check_mon_constraint
#' @param tolerance Maximum tolerated mean coefficient of variation.
#' @return Logical.
#' @keywords internal
check_rasch_constraint <- function(prob, tolerance = 0.25) {
  prob_clipped <- pmax(pmin(prob, 0.999), 0.001)
  logit <- log(prob_clipped / (1 - prob_clipped))

  class_order <- order(rowMeans(logit))
  logit_ordered <- logit[class_order, , drop = FALSE]
  n_classes <- nrow(logit_ordered)
  if (n_classes < 2) return(FALSE)

  cv_values <- numeric(0)
  for (c in seq_len(n_classes - 1)) {
    row_diffs <- logit_ordered[c + 1, ] - logit_ordered[c, ]
    if (sd(row_diffs) > 0 && abs(mean(row_diffs)) > 0.1) {
      cv_values <- c(cv_values, sd(row_diffs) / abs(mean(row_diffs)))
    }
  }
  if (length(cv_values) == 0) return(FALSE)
  mean(cv_values) < tolerance
}
