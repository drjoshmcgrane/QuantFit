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
#'
#' @param object A qlfit object
#' @param data Original data matrix (needed for computation)
#'
#' @return G-squared value with degrees of freedom
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

  # Compute expected frequencies under the model
  n_obs <- nrow(data)
  expected <- compute_expected_frequencies(
    object$item_probs,
    object$class_probs,
    patterns$pattern
  )

  # Compute G-squared
  observed <- patterns$frequency
  # Avoid log(0) issues
  valid <- observed > 0 & expected > 0
  g2 <- 2 * sum(observed[valid] * log(observed[valid] / expected[valid]))

  # Degrees of freedom: number of patterns - 1 - number of parameters
  n_patterns <- nrow(patterns)
  df <- n_patterns - 1 - object$n_par

  result <- list(
    statistic = g2,
    df = df,
    p_value = pchisq(g2, df, lower.tail = FALSE)
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

  # Scale by sample size (extract from patterns via sum)
  total_observed <- sum(sapply(patterns, function(p) {
    # This needs the frequency info which we don't have here
    # So the caller should handle scaling
    1
  }))

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

  # Compute chi-square
  observed <- patterns$frequency
  valid <- expected > 0
  x2 <- sum((observed[valid] - expected[valid])^2 / expected[valid])

  # Degrees of freedom
  n_patterns <- nrow(patterns)
  df <- n_patterns - 1 - object$n_par

  result <- list(
    statistic = x2,
    df = df,
    p_value = pchisq(x2, df, lower.tail = FALSE)
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
#' but NOT when comparing UN to constrained models (same parameters, different constraints).
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

  if (df <= 0) {
    warning("Degrees of freedom <= 0. Models may not be nested.")
    df <- abs(df)
  }

  p_value <- pchisq(lr_stat, df, lower.tail = FALSE)

  result <- list(
    statistic = lr_stat,
    df = df,
    p_value = p_value,
    model0 = model0$model_type,
    model1 = model1$model_type
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
  cat("P-value:", format.pval(x$p_value), "\n")

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

  if (is.na(g2_model$statistic)) {
    return(list(CFI = NA, TLI = NA, RMSEA = NA))
  }

  # Baseline model: independence (each item has single probability)
  n_items <- ncol(data)
  n_obs <- nrow(data)

  # Fit independence model (1-class LCA)
  # Log-likelihood under independence
  item_means <- colMeans(data)
  ll_indep <- sum(apply(data, 1, function(x) {
    sum(x * log(bound_probs(item_means)) +
        (1 - x) * log(bound_probs(1 - item_means)))
  }))

  # G-squared for independence model
  g2_null <- -2 * (ll_indep - object$loglik) +
             2 * (object$n_par - n_items)  # Approximate

  # CFI
  cfi <- max(0, 1 - (g2_model$statistic / g2_model$df) /
               (g2_null / max(1, g2_model$df + object$n_par - n_items)))

  # TLI
  tli <- max(0, ((g2_null / (g2_model$df + object$n_par - n_items)) -
                 (g2_model$statistic / g2_model$df)) /
               ((g2_null / (g2_model$df + object$n_par - n_items)) - 1))

  # RMSEA
  rmsea <- sqrt(max(0, (g2_model$statistic - g2_model$df) /
                       (g2_model$df * (n_obs - 1))))

  list(
    CFI = cfi,
    TLI = tli,
    RMSEA = rmsea,
    G2 = g2_model$statistic,
    df = g2_model$df,
    p_value = g2_model$p_value
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
