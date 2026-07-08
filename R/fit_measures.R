#' Fit Measures for QuantFit Models
#'
#' @name fit_measures
#' @description Functions for computing model fit statistics including
#'   log-likelihood, AIC, BIC, SABIC, and related measures
NULL

#' Extract fit measures from a qlfit object
#'
#' @param object A qlfit object
#' @param ... Additional arguments (ignored)
#'
#' @return Data frame with fit measures
#'
#' @examples
#' \dontrun{
#' fit <- fit_un(data, n_classes = 3)
#' fit_measures(fit)
#' }
#'
#' @export
fit_measures <- function(object, ...) {
  UseMethod("fit_measures")
}

#' @export
fit_measures.qlfit <- function(object, ...) {
  data.frame(
    model = object$model_type,
    loglik = object$loglik,
    n_par = object$n_par,
    n_obs = object$n_obs,
    AIC = AIC.qlfit(object),
    BIC = BIC.qlfit(object),
    SABIC = sabic(object),
    converged = object$convergence,
    stringsAsFactors = FALSE
  )
}

#' Compute G-squared (likelihood ratio chi-square)
#'
#' G² = 2 * sum(O * log(O/E)) where O is observed and E is expected frequency
#' (E = n * P(pattern) under the fitted model; patterns with O = 0 contribute
#' zero to the sum).
#'
#' @param object A qlfit object
#' @param data Original data matrix (needed for computation)
#'
#' @return G-squared value with degrees of freedom
#'
#' @details
#' Degrees of freedom are computed against the full multinomial over all
#' \eqn{2^I} response patterns: \eqn{df = 2^I - 1 - n_{par}}. Note that when
#' many cells are empty (sparse tables, common for larger I), the chi-square
#' reference distribution is questionable and the p-value should be treated
#' with caution.
#'
#' @export
g_squared <- function(object, data) {
  UseMethod("g_squared")
}

#' @export
g_squared.qlfit <- function(object, data) {
  if (object$model_type == "RM") {
    warning("G-squared not implemented for RM model")
    return(NA)
  }

  # Get observed pattern frequencies
  patterns <- pattern_frequencies(data)

  # Compute expected frequencies under the model: E = n * P(pattern)
  n_obs <- nrow(data)
  n_items <- ncol(data)
  expected_probs <- compute_expected_frequencies(
    object$item_probs,
    object$class_probs,
    patterns$pattern
  )
  expected <- expected_probs * n_obs

  # Compute G-squared (patterns with O = 0 contribute 0)
  observed <- patterns$frequency
  valid <- observed > 0 & expected > 0
  g2 <- 2 * sum(observed[valid] * log(observed[valid] / expected[valid]))

  # Degrees of freedom against the full multinomial over 2^I patterns.
  # Caveat: with many empty cells the chi-square approximation is poor.
  df <- 2^n_items - 1 - object$n_par

  result <- list(
    statistic = g2,
    df = df,
    p_value = if (df > 0) pchisq(g2, df, lower.tail = FALSE) else NA_real_
  )
  class(result) <- "g_squared"

  result
}

#' Compute expected pattern frequencies
#'
#' @param item_probs Item probability matrix (I x C)
#' @param class_probs Class probability vector
#' @param patterns Character vector of response patterns
#'
#' @return Vector of expected frequencies
#' @keywords internal
compute_expected_frequencies <- function(item_probs, class_probs, patterns) {
  n_patterns <- length(patterns)
  n_classes <- length(class_probs)
  n_obs <- 1  # Will be scaled later

  expected <- numeric(n_patterns)

  for (p in 1:n_patterns) {
    pattern_vec <- pattern_to_vector(patterns[p])

    # P(pattern) = sum_c P(c) * prod_i P(x_i | c)
    prob <- 0
    for (c in 1:n_classes) {
      class_prob <- class_probs[c]
      item_prob <- prod(
        item_probs[, c]^pattern_vec * (1 - item_probs[, c])^(1 - pattern_vec)
      )
      prob <- prob + class_prob * item_prob
    }

    expected[p] <- prob
  }

  # Returns pattern probabilities; the caller scales by sample size
  expected
}

#' Compute Pearson chi-square
#'
#' X² = sum((O - E)² / E)
#'
#' @param object A qlfit object
#' @param data Original data matrix
#'
#' @return Chi-square value with degrees of freedom
#' @export
pearson_chisq <- function(object, data) {
  UseMethod("pearson_chisq")
}

#' @export
pearson_chisq.qlfit <- function(object, data) {
  if (object$model_type == "RM") {
    warning("Pearson chi-square not implemented for RM model")
    return(NA)
  }

  # Get observed pattern frequencies
  patterns <- pattern_frequencies(data)

  # Compute expected frequencies under the model
  n_obs <- nrow(data)
  expected_probs <- compute_expected_frequencies(
    object$item_probs,
    object$class_probs,
    patterns$pattern
  )
  expected <- expected_probs * n_obs

  # Compute chi-square over observed patterns, plus the contribution of the
  # unobserved patterns: for O = 0 the term (O - E)^2 / E equals E, and the
  # expected frequencies over all 2^I patterns sum to n, so the unobserved
  # contribution is n - sum(E over observed patterns).
  observed <- patterns$frequency
  valid <- expected > 0
  x2 <- sum((observed[valid] - expected[valid])^2 / expected[valid]) +
    max(0, n_obs - sum(expected[valid]))

  # Degrees of freedom against the full multinomial over 2^I patterns
  # (same sparseness caveat as g_squared)
  n_items <- ncol(data)
  df <- 2^n_items - 1 - object$n_par

  result <- list(
    statistic = x2,
    df = df,
    p_value = if (df > 0) pchisq(x2, df, lower.tail = FALSE) else NA_real_
  )
  class(result) <- "pearson_chisq"

  result
}

#' Likelihood ratio test between two models
#'
#' @param model0 Restricted model (qlfit object)
#' @param model1 Full model (qlfit object)
#'
#' @return List with test statistic, df, and p-value
#'
#' @details
#' This test is appropriate when model0 is nested within model1.
#' For QuantFit models, this applies when comparing RM to LCR (with many classes)
#' but NOT when comparing UN to constrained models (MON, IIO, DM): those share
#' the same parameter count and differ only by inequality constraints, so the
#' LR statistic follows a chi-bar-squared (mixture of chi-squares)
#' distribution, not a chi-square with 0 df. In that case the p-value is
#' returned as NA with an explanatory message.
#'
#' When comparing LCR or RM to other models, note that parameters may lie on
#' the boundary of the parameter space, in which case the chi-square
#' approximation is conservative.
#'
#' @export
lr_test <- function(model0, model1) {
  if (!inherits(model0, "qlfit") || !inherits(model1, "qlfit")) {
    stop("Both models must be qlfit objects")
  }

  # Check that model1 has more (or equal) parameters
  if (model1$n_par < model0$n_par) {
    # Swap
    temp <- model0
    model0 <- model1
    model1 <- temp
    message("Models swapped: model0 is now the restricted model")
  }

  # Compute LR statistic
  lr_stat <- -2 * (model0$loglik - model1$loglik)
  df <- model1$n_par - model0$n_par

  note <- NULL
  if (df == 0) {
    # Same parameter count: models differ by inequality constraints only
    # (e.g. UN vs MON/IIO/DM). The LR statistic follows a chi-bar-squared
    # mixture, not chi-square with 0 df, so no standard p-value exists.
    note <- paste(
      "The models have the same number of parameters and differ only by",
      "inequality constraints. The LR statistic follows a chi-bar-squared",
      "(mixture of chi-squares) distribution; the standard chi-square test",
      "does not apply. Returning NA p-value. Consider information criteria",
      "or a parametric bootstrap instead."
    )
    warning(note, call. = FALSE)
    p_value <- NA_real_
  } else {
    # Boundary-condition caveat for LCR/RM comparisons
    if (any(c(model0$model_type, model1$model_type) %in% c("LCR", "RM"))) {
      warning("Comparison involves LCR/RM: parameters may lie on the ",
              "boundary of the parameter space, so the chi-square p-value ",
              "may be conservative.", call. = FALSE)
    }
    p_value <- pchisq(lr_stat, df, lower.tail = FALSE)
  }

  result <- list(
    statistic = lr_stat,
    df = df,
    p_value = p_value,
    model0 = model0$model_type,
    model1 = model1$model_type,
    note = note
  )
  class(result) <- "lr_test"

  result
}

#' Print method for lr_test
#'
#' @param x An lr_test object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.lr_test <- function(x, ...) {
  cat("\nLikelihood Ratio Test\n")
  cat("---------------------\n")
  cat("Restricted model:", x$model0, "\n")
  cat("Full model:", x$model1, "\n")
  cat("\nLR statistic:", round(x$statistic, 4), "\n")
  cat("Degrees of freedom:", x$df, "\n")
  if (is.na(x$p_value)) {
    cat("P-value: NA\n")
    if (!is.null(x$note)) cat("Note:", x$note, "\n")
  } else {
    cat("P-value:", format.pval(x$p_value), "\n")
  }

  invisible(x)
}

#' Compute entropy-based classification quality
#'
#' Entropy R² measures how well the model classifies observations
#'
#' @param object A qlfit object
#'
#' @return Entropy R² value (0 to 1, higher is better)
#' @export
entropy_r2 <- function(object) {
  UseMethod("entropy_r2")
}

#' @export
entropy_r2.qlfit <- function(object) {
  if (is.null(object$posteriors)) {
    warning("Posteriors not available")
    return(NA)
  }

  n_obs <- nrow(object$posteriors)
  n_classes <- ncol(object$posteriors)

  # Entropy of posterior probabilities
  # E = -sum(p * log(p)) / (N * log(K))
  posteriors <- bound_probs(object$posteriors)  # Avoid log(0)
  entropy <- -sum(posteriors * log(posteriors))
  max_entropy <- n_obs * log(n_classes)

  # Entropy R²
  1 - entropy / max_entropy
}

#' Compute classification accuracy
#'
#' Modal assignment accuracy based on posterior probabilities
#'
#' @param object A qlfit object
#'
#' @return Mean posterior probability for modal class assignments
#' @export
classification_accuracy <- function(object) {
  UseMethod("classification_accuracy")
}

#' @export
classification_accuracy.qlfit <- function(object) {
  if (is.null(object$posteriors)) {
    warning("Posteriors not available")
    return(NA)
  }

  # For each observation, get the maximum posterior probability
  max_posteriors <- apply(object$posteriors, 1, max)

  mean(max_posteriors)
}

#' Get modal class assignments
#'
#' @param object A qlfit object
#'
#' @return Vector of class assignments (1 to n_classes)
#' @export
class_assignments <- function(object) {
  UseMethod("class_assignments")
}

#' @export
class_assignments.qlfit <- function(object) {
  if (is.null(object$posteriors)) {
    stop("Posteriors not available")
  }

  apply(object$posteriors, 1, which.max)
}

#' Compute relative fit indices
#'
#' Computes CFI, TLI, and RMSEA (approximations for LCA)
#'
#' @param object A qlfit object
#' @param data Original data matrix
#'
#' @return List with fit indices
#' @export
relative_fit <- function(object, data) {
  UseMethod("relative_fit")
}

#' @export
relative_fit.qlfit <- function(object, data) {
  # Get G-squared for the model
  g2_model <- g_squared(object, data)

  if (length(g2_model) == 1 && is.na(g2_model)) {
    return(list(CFI = NA, TLI = NA, RMSEA = NA))
  }

  n_items <- ncol(data)
  n_obs <- nrow(data)

  # Baseline: a genuine 1-class (independence) model. Its MLE is the vector
  # of item means; compute its G-squared directly against the observed
  # pattern frequencies.
  item_means <- bound_probs(colMeans(data))
  patterns <- pattern_frequencies(data)
  null_probs <- compute_expected_frequencies(
    matrix(item_means, ncol = 1),  # I x 1 item probability matrix
    1,                              # single class with probability 1
    patterns$pattern
  )
  expected_null <- null_probs * n_obs
  observed <- patterns$frequency
  valid <- observed > 0 & expected_null > 0
  g2_null <- 2 * sum(observed[valid] * log(observed[valid] / expected_null[valid]))
  df_null <- 2^n_items - 1 - n_items

  df_m <- g2_model$df
  num <- max(g2_model$statistic - df_m, 0)
  den <- max(g2_model$statistic - df_m, g2_null - df_null, 0)

  # CFI
  cfi <- if (den > 0) 1 - num / den else 1

  # TLI
  tli <- if (df_null > 0 && df_m > 0 && (g2_null / df_null) != 1) {
    ((g2_null / df_null) - (g2_model$statistic / df_m)) /
      ((g2_null / df_null) - 1)
  } else {
    NA_real_
  }

  # RMSEA
  rmsea <- if (df_m > 0) {
    sqrt(max(0, (g2_model$statistic - df_m) / (df_m * (n_obs - 1))))
  } else {
    NA_real_
  }

  list(
    CFI = cfi,
    TLI = tli,
    RMSEA = rmsea,
    G2 = g2_model$statistic,
    df = g2_model$df,
    p_value = g2_model$p_value,
    G2_null = g2_null,
    df_null = df_null
  )
}

#' Compare fit of multiple models
#'
#' @param ... qlfit objects to compare
#' @param measures Which measures to include: "all", "ic" (information criteria only),
#'   or a character vector of specific measures
#'
#' @return Data frame with fit measures for all models
#' @export
compare_fit <- function(..., measures = "ic") {
  models <- list(...)

  if (length(models) == 1 && is.list(models[[1]])) {
    models <- models[[1]]
  }

  # Extract fit measures from each model
  fit_list <- lapply(models, function(m) {
    if (!inherits(m, "qlfit")) {
      stop("All inputs must be qlfit objects")
    }
    fit_measures(m)
  })

  # Combine into data frame
  result <- do.call(rbind, fit_list)

  # Add model names if available
  if (!is.null(names(models))) {
    result$model_name <- names(models)
  }

  # Filter measures if requested
  if (measures == "ic") {
    result <- result[, c("model", "loglik", "n_par", "AIC", "BIC", "SABIC", "converged")]
  }

  # Sort by BIC
  result <- result[order(result$BIC), ]
  rownames(result) <- NULL

  result
}
