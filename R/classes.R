#' S3 Class Definitions and Methods for QuantFit
#'
#' @name qlfit-class
#' @description S3 class for fitted QuantFit models
NULL

#' Create a qlfit object
#'
#' @param model_type Character string indicating model type
#'   ("UN", "MON", "IIO", "DM", "LCR", "RM")
#' @param item_probs Matrix of item probabilities (items x classes)
#' @param class_probs Vector of class probabilities
#' @param posteriors Matrix of posterior class memberships (n x classes)
#' @param loglik Log-likelihood value
#' @param n_par Number of estimated parameters
#' @param n_obs Number of observations
#' @param n_items Number of items
#' @param n_classes Number of latent classes
#' @param convergence Logical indicating convergence
#' @param iterations Number of EM iterations
#' @param call Original function call
#' @param theta Class locations (for LCR/RM models)
#' @param delta Item difficulties (for LCR/RM models)
#' @param item_order Item ordering (for IIO/DM models)
#' @param constraints List of active constraints
#' @param se Standard errors (if computed)
#'
#' @return Object of class "qlfit"
#' @keywords internal
new_qlfit <- function(model_type,
                      item_probs,
                      class_probs,
                      posteriors,
                      loglik,
                      n_par,
                      n_obs,
                      n_items,
                      n_classes,
                      convergence,
                      iterations,
                      call = NULL,
                      theta = NULL,
                      delta = NULL,
                      item_order = NULL,
                      constraints = NULL,
                      se = NULL) {

  structure(
    list(
      model_type = model_type,
      item_probs = item_probs,
      class_probs = class_probs,
      posteriors = posteriors,
      loglik = loglik,
      n_par = n_par,
      n_obs = n_obs,
      n_items = n_items,
      n_classes = n_classes,
      convergence = convergence,
      iterations = iterations,
      call = call,
      theta = theta,
      delta = delta,
      item_order = item_order,
      constraints = constraints,
      se = se
    ),
    class = "qlfit"
  )
}

#' Print method for qlfit objects
#'
#' @param x A qlfit object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns x
#' @export
print.qlfit <- function(x, ...) {
  model_names <- c(
    UN = "Unconstrained Latent Class",
    MON = "Class Monotonicity",
    IIO = "Invariant Item Ordering",
    DM = "Double Monotonicity",
    LCR = "Latent Class Rasch",
    RM = "Rasch Model"
  )

  cat("\n")
  cat("QuantFit Model:", model_names[x$model_type], "(", x$model_type, ")\n")
  cat(rep("-", 50), "\n", sep = "")
  cat("Number of observations:", x$n_obs, "\n")
  cat("Number of items:", x$n_items, "\n")

  if (x$model_type != "RM") {
    cat("Number of classes:", x$n_classes, "\n")
  }

  cat("\n")
  cat("Log-likelihood:", format(x$loglik, nsmall = 2), "\n")
  cat("Number of parameters:", x$n_par, "\n")
  cat("AIC:", format(AIC.qlfit(x), nsmall = 2), "\n")
  cat("BIC:", format(BIC.qlfit(x), nsmall = 2), "\n")

  cat("\n")
  if (x$convergence) {
    cat("Converged in", x$iterations, "iterations\n")
  } else {
    cat("WARNING: Did not converge (", x$iterations, " iterations)\n")
  }

  invisible(x)
}

#' Summary method for qlfit objects
#'
#' @param object A qlfit object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns a summary list
#' @export
summary.qlfit <- function(object, ...) {
  model_names <- c(
    UN = "Unconstrained Latent Class",
    MON = "Class Monotonicity",
    IIO = "Invariant Item Ordering",
    DM = "Double Monotonicity",
    LCR = "Latent Class Rasch",
    RM = "Rasch Model"
  )

  cat("\n")
  cat("QuantFit Model Summary\n")
  cat("=========================\n\n")

  cat("Model:", model_names[object$model_type], "(", object$model_type, ")\n\n")

  # Fit statistics
  cat("Fit Statistics:\n")
  cat("  Log-likelihood:", format(object$loglik, nsmall = 2), "\n")
  cat("  Parameters:", object$n_par, "\n")
  cat("  AIC:", format(AIC.qlfit(object), nsmall = 2), "\n")
  cat("  BIC:", format(BIC.qlfit(object), nsmall = 2), "\n")
  cat("  SABIC:", format(sabic(object), nsmall = 2), "\n\n")

  # Class proportions (for LCA models)
  if (object$model_type != "RM") {
    cat("Class Proportions:\n")
    class_props <- setNames(object$class_probs,
                            paste0("Class ", seq_along(object$class_probs)))
    print(round(class_props, 4))
    cat("\n")

    # Item probabilities
    cat("Item Response Probabilities:\n")
    item_probs_df <- as.data.frame(object$item_probs)
    colnames(item_probs_df) <- paste0("Class ", 1:object$n_classes)
    rownames(item_probs_df) <- paste0("Item ", 1:object$n_items)
    print(round(item_probs_df, 4))
  }

  # Rasch parameters (for LCR/RM)
  if (object$model_type %in% c("LCR", "RM")) {
    cat("\nItem Difficulties (delta):\n")
    delta_vec <- setNames(object$delta, paste0("Item ", seq_along(object$delta)))
    print(round(delta_vec, 4))

    if (object$model_type == "LCR") {
      cat("\nClass Locations (theta):\n")
      theta_vec <- setNames(object$theta,
                            paste0("Class ", seq_along(object$theta)))
      print(round(theta_vec, 4))
    }
  }

  # Item ordering (for IIO/DM)
  if (!is.null(object$item_order)) {
    cat("\nEstimated Item Order (easiest to hardest):\n")
    cat(object$item_order, "\n")
  }

  # Convergence info
  cat("\n")
  if (object$convergence) {
    cat("Converged in", object$iterations, "iterations\n")
  } else {
    cat("WARNING: Did not converge after", object$iterations, "iterations\n")
  }

  invisible(list(
    model_type = object$model_type,
    loglik = object$loglik,
    n_par = object$n_par,
    aic = AIC.qlfit(object),
    bic = BIC.qlfit(object),
    class_probs = object$class_probs,
    item_probs = object$item_probs
  ))
}

#' Extract coefficients from qlfit object
#'
#' @param object A qlfit object
#' @param type Type of coefficients: "probs" (default), "logit", or "rasch"
#' @param ... Additional arguments (ignored)
#'
#' @return Matrix or list of coefficients
#' @export
coef.qlfit <- function(object, type = c("probs", "logit", "rasch"), ...) {
  type <- match.arg(type)

  if (type == "probs") {
    result <- list(
      class_probs = object$class_probs,
      item_probs = object$item_probs
    )
  } else if (type == "logit") {
    result <- list(
      class_probs = object$class_probs,
      item_logits = qlogis(object$item_probs)
    )
  } else if (type == "rasch") {
    if (is.null(object$theta) || is.null(object$delta)) {
      stop("Rasch parameters only available for LCR and RM models")
    }
    result <- list(
      theta = object$theta,
      delta = object$delta,
      class_probs = object$class_probs
    )
  }

  result
}

#' Log-likelihood for qlfit object
#'
#' @param object A qlfit object
#' @param ... Additional arguments (ignored)
#'
#' @return Log-likelihood value with attributes
#' @export
logLik.qlfit <- function(object, ...) {
  ll <- object$loglik
  attr(ll, "df") <- object$n_par
  attr(ll, "nobs") <- object$n_obs
  class(ll) <- "logLik"
  ll
}

#' AIC for qlfit object
#'
#' @param object A qlfit object
#' @param ... Additional arguments (ignored)
#' @param k Penalty parameter (default 2)
#'
#' @return AIC value
#' @export
AIC.qlfit <- function(object, ..., k = 2) {
  -2 * object$loglik + k * object$n_par
}

#' BIC for qlfit object
#'
#' @param object A qlfit object
#' @param ... Additional arguments (ignored)
#'
#' @return BIC value
#' @export
BIC.qlfit <- function(object, ...) {
  -2 * object$loglik + object$n_par * log(object$n_obs)
}

#' Sample-adjusted BIC (SABIC) for qlfit object
#'
#' @param object A qlfit object
#'
#' @return SABIC value
#' @export
sabic <- function(object) {
  UseMethod("sabic")
}

#' @export
sabic.qlfit <- function(object) {
  n_star <- (object$n_obs + 2) / 24
  -2 * object$loglik + object$n_par * log(n_star)
}

#' Plot method for qlfit objects
#'
#' @param x A qlfit object
#' @param type Type of plot: "irf" (default) or "class"
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Invisibly returns x
#' @export
plot.qlfit <- function(x, type = c("irf", "class"), ...) {
  type <- match.arg(type)

  if (type == "irf") {
    plot_irfs(x, ...)
  } else if (type == "class") {
    plot_class_profiles(x, ...)
  }

  invisible(x)
}


# ============================================================================
# qlcompare class for model comparison results
# ============================================================================

#' Create a qlcompare object
#'
#' @param fits List of qlfit objects
#' @param comparison_table Data frame with comparison statistics
#' @param best_model Name of best model by BIC
#' @param call Original function call
#'
#' @return Object of class "qlcompare"
#' @keywords internal
new_qlcompare <- function(fits, comparison_table, best_model, call = NULL) {
  structure(
    list(
      fits = fits,
      comparison_table = comparison_table,
      best_model = best_model,
      call = call
    ),
    class = "qlcompare"
  )
}

#' Print method for qlcompare objects
#'
#' @param x A qlcompare object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns x
#' @export
print.qlcompare <- function(x, ...) {
  cat("\nQuantFit Model Comparison\n")
  cat("============================\n\n")

  # Print comparison table
  print(x$comparison_table, row.names = FALSE)

  cat("\n")
  cat("Best model by BIC:", x$best_model, "\n")

  invisible(x)
}

#' Summary method for qlcompare objects
#'
#' @param object A qlcompare object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the comparison table
#' @export
summary.qlcompare <- function(object, ...) {
  cat("\nQuantFit Model Comparison Summary\n")
  cat("====================================\n\n")

  # Enhanced comparison table with delta values
  tbl <- object$comparison_table
  tbl$delta_AIC <- tbl$AIC - min(tbl$AIC)
  tbl$delta_BIC <- tbl$BIC - min(tbl$BIC)

  print(tbl, row.names = FALSE)

  cat("\n")
  cat("Interpretation guide:\n")
  cat("  - Lower AIC/BIC indicates better fit\n")
  cat("  - delta_BIC < 2: Models essentially equivalent\n")
  cat("  - delta_BIC 2-6: Positive evidence for better model\n")
  cat("  - delta_BIC 6-10: Strong evidence\n")
  cat("  - delta_BIC > 10: Very strong evidence\n")
  cat("\n")
  cat("Best model by BIC:", object$best_model, "\n")

  invisible(tbl)
}

#' Plot method for qlcompare objects
#'
#' @param x A qlcompare object
#' @param ... Additional arguments passed to plot_comparison
#'
#' @return Invisibly returns x
#' @export
plot.qlcompare <- function(x, ...) {
  plot_comparison(x, ...)
  invisible(x)
}

#' Extract a fitted model from comparison
#'
#' @param object A qlcompare object
#' @param model Model name to extract
#'
#' @return A qlfit object
#' @export
get_model <- function(object, model) {
  UseMethod("get_model")
}

#' @export
get_model.qlcompare <- function(object, model) {
  if (!model %in% names(object$fits)) {
    stop("Model '", model, "' not found. Available models: ",
         paste(names(object$fits), collapse = ", "))
  }
  object$fits[[model]]
}
