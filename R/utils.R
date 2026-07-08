#' Utility Functions for QuantFit
#'
#' @name utils
#' @description Helper functions for data validation, initialization, and common operations
NULL

#' Validate input data for latent class models
#'
#' @param data Matrix or data frame of binary responses (0/1)
#' @param allow_na Logical, whether to allow missing values
#'
#' @return Validated data matrix
#' @keywords internal
validate_data <- function(data, allow_na = FALSE) {
  # Convert to matrix
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  if (!is.matrix(data)) {
    stop("Data must be a matrix or data frame")
  }

  # Check for binary values
  unique_vals <- unique(as.vector(data[!is.na(data)]))
  if (!all(unique_vals %in% c(0, 1))) {
    stop("Data must contain only binary values (0 and 1)")
  }

  # Check for missing values
  if (!allow_na && any(is.na(data))) {
    stop("Data contains missing values. Set allow_na = TRUE to allow")
  }

  # Ensure numeric
  storage.mode(data) <- "double"

  data
}

#' Initialize class probabilities
#'
#' @param n_classes Number of latent classes
#' @param method Initialization method: "uniform", "random", or "dirichlet"
#' @param seed Random seed (optional)
#'
#' @return Vector of class probabilities summing to 1
#' @keywords internal
init_class_probs <- function(n_classes, method = c("uniform", "random", "dirichlet"),
                             seed = NULL) {
  method <- match.arg(method)

  if (!is.null(seed)) set.seed(seed)

  if (method == "uniform") {
    probs <- rep(1 / n_classes, n_classes)
  } else if (method == "random") {
    probs <- runif(n_classes)
    probs <- probs / sum(probs)
  } else if (method == "dirichlet") {
    # Simple Dirichlet draw with alpha = 1 (uniform)
    probs <- rgamma(n_classes, shape = 1, rate = 1)
    probs <- probs / sum(probs)
  }

  probs
}

#' Initialize item probabilities for LCA
#'
#' @param data Data matrix
#' @param n_classes Number of latent classes
#' @param method Initialization method: "random", "kmeans", or "quantiles"
#' @param seed Random seed (optional)
#'
#' @return Matrix of item probabilities (n_items x n_classes)
#' @keywords internal
init_item_probs <- function(data, n_classes,
                            method = c("random", "kmeans", "quantiles"),
                            seed = NULL) {
  method <- match.arg(method)

  if (!is.null(seed)) set.seed(seed)

  n_items <- ncol(data)
  n_obs <- nrow(data)

  if (method == "random") {
    # Random probabilities between 0.1 and 0.9
    probs <- matrix(runif(n_items * n_classes, 0.1, 0.9),
                    nrow = n_items, ncol = n_classes)
  } else if (method == "kmeans") {
    # Use k-means clustering for initialization
    km <- tryCatch({
      kmeans(data, centers = n_classes, nstart = 5)
    }, error = function(e) NULL)

    if (is.null(km)) {
      # Fall back to random if kmeans fails
      return(init_item_probs(data, n_classes, "random", seed))
    }

    # Compute class-specific item means
    probs <- matrix(0, nrow = n_items, ncol = n_classes)
    for (c in 1:n_classes) {
      class_data <- data[km$cluster == c, , drop = FALSE]
      if (nrow(class_data) > 0) {
        probs[, c] <- colMeans(class_data, na.rm = TRUE)
      } else {
        probs[, c] <- runif(n_items, 0.3, 0.7)
      }
    }

    # Ensure probabilities are bounded away from 0 and 1
    probs <- pmax(pmin(probs, 0.99), 0.01)

  } else if (method == "quantiles") {
    # Initialize based on total score quantiles
    total_scores <- rowSums(data, na.rm = TRUE)
    quantile_cuts <- quantile(total_scores, probs = seq(0, 1, length.out = n_classes + 1))

    probs <- matrix(0, nrow = n_items, ncol = n_classes)
    for (c in 1:n_classes) {
      if (c == 1) {
        in_class <- total_scores <= quantile_cuts[c + 1]
      } else if (c == n_classes) {
        in_class <- total_scores > quantile_cuts[c]
      } else {
        in_class <- total_scores > quantile_cuts[c] & total_scores <= quantile_cuts[c + 1]
      }

      if (sum(in_class) > 0) {
        probs[, c] <- colMeans(data[in_class, , drop = FALSE], na.rm = TRUE)
      } else {
        probs[, c] <- (c - 0.5) / n_classes
      }
    }

    # Add jitter so multiple random starts actually differ (the quantile
    # split itself is deterministic given the data)
    probs <- probs + matrix(runif(n_items * n_classes, -0.05, 0.05),
                            nrow = n_items, ncol = n_classes)

    # Ensure probabilities are bounded
    probs <- pmax(pmin(probs, 0.99), 0.01)
  }

  probs
}

#' Initialize item probabilities with monotonicity constraints
#'
#' @param data Data matrix
#' @param n_classes Number of latent classes
#' @param seed Random seed (optional)
#'
#' @return Matrix of item probabilities satisfying class monotonicity
#' @keywords internal
init_item_probs_monotonic <- function(data, n_classes, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_items <- ncol(data)

  # Start with quantile-based initialization
  probs <- init_item_probs(data, n_classes, "quantiles", seed)

  # Ensure monotonicity: sort each row
  for (i in 1:n_items) {
    probs[i, ] <- sort(probs[i, ])
  }

  probs
}

#' Compute log-sum-exp in a numerically stable way
#'
#' @param x Numeric vector
#'
#' @return log(sum(exp(x)))
#' @keywords internal
log_sum_exp <- function(x) {
  max_x <- max(x)
  if (is.infinite(max_x)) return(max_x)
  max_x + log(sum(exp(x - max_x)))
}

#' Compute log-sum-exp for matrix rows
#'
#' @param x Numeric matrix
#'
#' @return Vector of log(sum(exp(x))) for each row
#' @keywords internal
row_log_sum_exp <- function(x) {
  apply(x, 1, log_sum_exp)
}

#' Softmax function
#'
#' @param x Numeric vector
#'
#' @return Vector of probabilities summing to 1
#' @keywords internal
softmax <- function(x) {
  exp_x <- exp(x - max(x))
  exp_x / sum(exp_x)
}

#' Inverse logit (logistic) function
#'
#' @param x Numeric vector
#'
#' @return Vector of probabilities
#' @keywords internal
inv_logit <- function(x) {

  1 / (1 + exp(-x))
}

#' Logit function
#'
#' @param p Probability vector (values in (0,1))
#'
#' @return Vector of log-odds
#' @keywords internal
logit <- function(p) {
  log(p / (1 - p))
}

#' Bound probabilities away from 0 and 1
#'
#' @param p Probability vector or matrix
#' @param eps Minimum distance from 0 and 1 (default 1e-10)
#'
#' @return Bounded probabilities
#' @keywords internal
bound_probs <- function(p, eps = 1e-10) {
  pmax(pmin(p, 1 - eps), eps)
}

#' Compute response pattern frequencies
#'
#' @param data Binary data matrix
#'
#' @return Data frame with patterns and frequencies
#' @keywords internal
pattern_frequencies <- function(data) {
  # Convert each row to string pattern
  patterns <- apply(data, 1, paste, collapse = "")

  # Count frequencies
  freq_table <- table(patterns)

  data.frame(
    pattern = names(freq_table),
    frequency = as.integer(freq_table),
    stringsAsFactors = FALSE
  )
}

#' Convert pattern string back to numeric vector
#'
#' @param pattern Character string of 0s and 1s
#'
#' @return Numeric vector
#' @keywords internal
pattern_to_vector <- function(pattern) {
  as.numeric(strsplit(pattern, "")[[1]])
}

#' Check if model converged
#'
#' @param ll_history Vector of log-likelihood values
#' @param tol Convergence tolerance
#' @param min_iter Minimum iterations before checking
#'
#' @return Logical indicating convergence
#' @keywords internal
check_convergence <- function(ll_history, tol = 1e-6, min_iter = 3) {
  n <- length(ll_history)
  if (n < min_iter) return(FALSE)

  # Relative change in log-likelihood
  change <- ll_history[n] - ll_history[n - 1]
  rel_change <- change / abs(ll_history[n - 1])

  # Monotonicity guard: a meaningful decrease violates the (generalized)
  # EM ascent property and must not be declared convergence
  if (rel_change < -tol) {
    warning("Log-likelihood decreased by ", format(-change),
            " between EM iterations (generalized EM ascent property ",
            "violated); not declaring convergence.", call. = FALSE)
    return(FALSE)
  }

  abs(rel_change) < tol
}

#' Estimate item ordering from data
#'
#' @param data Binary data matrix
#'
#' @return Vector of item indices ordered from easiest to hardest
#' @keywords internal
estimate_item_order <- function(data) {
  # Order items by proportion correct (easiest to hardest)
  item_means <- colMeans(data, na.rm = TRUE)
  order(item_means, decreasing = TRUE)
}

#' Count number of free parameters
#'
#' @param model_type Model type code
#' @param n_items Number of items
#' @param n_classes Number of classes
#'
#' @return Number of free parameters
#' @keywords internal
count_parameters <- function(model_type, n_items, n_classes) {
  switch(model_type,
    UN = {
      # Class probs: (C-1) + Item probs: I * C
      (n_classes - 1) + n_items * n_classes
    },
    MON = {
      # Same as UN (constraints don't reduce parameters, just constrain them)
      (n_classes - 1) + n_items * n_classes
    },
    IIO = {
      # Same as UN
      (n_classes - 1) + n_items * n_classes
    },
    DM = {
      # Same as UN
      (n_classes - 1) + n_items * n_classes
    },
    LCR = {
      # The implementation leaves all C thetas free and uses mean(delta) = 0
      # as the single identification constraint:
      # (C-1) class probs + C theta + (I-1) delta = 2C + I - 2
      (n_classes - 1) + n_classes + (n_items - 1)
    },
    RM = {
      # mirt Rasch parameterization: I item intercepts + latent variance
      # (latent mean fixed at 0) = I + 1 estimated parameters
      n_items + 1
    },
    stop("Unknown model type: ", model_type)
  )
}

#' Multiple random starts wrapper
#'
#' @param fit_fn Fitting function to call
#' @param n_starts Number of random starts
#' @param seed Base random seed
#' @param verbose Print progress
#' @param ... Arguments passed to fit_fn
#'
#' @return Best fit across all starts
#' @keywords internal
multiple_starts <- function(fit_fn, n_starts = 10, seed = NULL, verbose = FALSE, ...) {
  if (!is.null(seed)) set.seed(seed)

  best_fit <- NULL
  best_ll <- -Inf

  for (i in 1:n_starts) {
    if (verbose) cat("Start", i, "of", n_starts, "... ")

    fit <- tryCatch({
      fit_fn(seed = seed + i, ...)
    }, error = function(e) {
      if (verbose) cat("failed:", e$message, "\n")
      NULL
    })

    if (!is.null(fit) && fit$loglik > best_ll) {
      best_fit <- fit
      best_ll <- fit$loglik
      if (verbose) cat("LL =", round(best_ll, 2), "(new best)\n")
    } else if (verbose && !is.null(fit)) {
      cat("LL =", round(fit$loglik, 2), "\n")
    }
  }

  if (is.null(best_fit)) {
    stop("All random starts failed")
  }

  best_fit
}

#' Print progress bar
#'
#' @param iteration Current iteration
#' @param total Total iterations
#' @param width Width of progress bar
#'
#' @keywords internal
progress_bar <- function(iteration, total, width = 50) {
  pct <- iteration / total
  filled <- round(pct * width)
  bar <- paste0(
    "\r[",
    paste(rep("=", filled), collapse = ""),
    paste(rep(" ", width - filled), collapse = ""),
    "] ",
    sprintf("%3.0f%%", pct * 100)
  )
  cat(bar)
  if (iteration == total) cat("\n")
}
