#' Fit Invariant Item Ordering Model (IIO)
#'
#' @name fit_iio
#' @description Fits an ordered latent class model with invariant item ordering constraints.
#'   This model constrains item probabilities to maintain the same ordering across all
#'   latent classes.
NULL

#' Fit Invariant Item Ordering Model
#'
#' @param data A matrix or data frame of binary responses (0/1) with subjects in rows
#'   and items in columns.
#' @param n_classes Number of latent classes (integer >= 2)
#' @param item_order Optional vector of item indices specifying the order from easiest
#'   to hardest. If NULL (default), the order is estimated from the data using
#'   overall item proportions.
#' @param n_starts Number of random starting values to try (default 10)
#' @param max_iter Maximum number of EM iterations per start (default 500)
#' @param tol Convergence tolerance for log-likelihood (default 1e-6)
#' @param method M-step method: "pava" (default) uses the exact constrained
#'   M-step via isotonic regression per class; "optimizer" uses a general
#'   constrained optimizer
#' @param optimizer Optimizer for the constrained M-step when
#'   method = "optimizer": "alabama" (default) or "nloptr"
#' @param seed Random seed for reproducibility (optional)
#' @param verbose Print progress messages (default FALSE)
#'
#' @return A qlfit object containing:
#' \describe{
#'   \item{model_type}{"IIO" for invariant item ordering}
#'   \item{item_probs}{Matrix of item probabilities (n_items x n_classes) satisfying
#'     item ordering constraints}
#'   \item{class_probs}{Vector of class probabilities}
#'   \item{posteriors}{Matrix of posterior class memberships (n_obs x n_classes)}
#'   \item{loglik}{Maximized log-likelihood}
#'   \item{n_par}{Number of estimated parameters}
#'   \item{item_order}{Vector of item indices from easiest to hardest}
#'   \item{convergence}{Logical indicating successful convergence}
#'   \item{constraints}{Active constraint specification}
#' }
#'
#' @details
#' The invariant item ordering model adds the constraint that for each class c
#' and for the specified item ordering i < i':
#' \deqn{P(X_i = 1 | c) \geq P(X_{i'} = 1 | c)}
#'
#' This means easier items (higher probability of success) maintain their relative
#' ordering across all latent classes. This is related to the concept of
#' "invariant item ordering" from Mokken scale analysis.
#'
#' If \code{item_order} is not specified, it is estimated from the data by ordering
#' items from highest to lowest overall proportion correct.
#'
#' @examples
#' \dontrun{
#' # Generate data with item ordering
#' set.seed(123)
#' n <- 500
#' data <- matrix(rbinom(n * 10, 1, 0.5), nrow = n)
#'
#' # Fit with estimated item order
#' fit <- fit_iio(data, n_classes = 3)
#' print(fit)
#' cat("Estimated item order:", fit$item_order, "\n")
#'
#' # Fit with specified item order
#' fit2 <- fit_iio(data, n_classes = 3, item_order = 1:10)
#' }
#'
#' @seealso \code{\link{fit_un}} for unconstrained model,
#'   \code{\link{fit_mon}} for class monotonicity,
#'   \code{\link{fit_dm}} for double monotonicity
#'
#' @export
fit_iio <- function(data, n_classes,
                    item_order = NULL,
                    n_starts = 10,
                    max_iter = 500,
                    tol = 1e-6,
                    method = c("pava", "optimizer"),
                    optimizer = c("alabama", "nloptr"),
                    seed = NULL,
                    verbose = FALSE) {

  # Capture call
  call <- match.call()

  # Validate inputs
  data <- validate_data(data)
  method <- match.arg(method)
  optimizer <- match.arg(optimizer)

  n_obs <- nrow(data)
  n_items <- ncol(data)

  if (!is.numeric(n_classes) || n_classes < 2 || n_classes != round(n_classes)) {
    stop("n_classes must be an integer >= 2")
  }

  # Estimate or validate item order
  if (is.null(item_order)) {
    item_order <- estimate_item_order(data)
    if (verbose) {
      cat("Estimated item order (easiest to hardest):", item_order, "\n")
    }
  } else {
    # Validate provided order
    if (length(item_order) != n_items) {
      stop("item_order must have length equal to number of items (", n_items, ")")
    }
    if (!all(sort(item_order) == 1:n_items)) {
      stop("item_order must contain each integer from 1 to ", n_items, " exactly once")
    }
  }

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Build constraint function
  constraint_fn <- build_item_ordering_fn(n_items, n_classes, item_order)
  constraints_spec <- specify_constraints(item_ordering = TRUE, item_order = item_order)

  # Run multiple starts
  best_fit <- NULL
  best_ll <- -Inf

  for (start in 1:n_starts) {
    if (verbose) {
      cat("Start", start, "of", n_starts, "... ")
    }

    # Initialize with item ordering satisfied
    start_seed <- if (!is.null(seed)) seed + start else NULL

    init_probs <- init_item_probs_iio(data, n_classes, item_order, start_seed)
    init_class_probs <- init_class_probs(n_classes, "random", start_seed)

    # Run constrained EM
    fit <- tryCatch({
      em_constrained(
        data = data,
        n_classes = n_classes,
        constraints_spec = constraints_spec,
        item_order = item_order,
        constraint_fn = constraint_fn,
        init_probs = init_probs,
        init_class_probs = init_class_probs,
        max_iter = max_iter,
        tol = tol,
        method = method,
        optimizer = optimizer,
        verbose = FALSE
      )
    }, error = function(e) {
      if (verbose) cat("failed: ", e$message, "\n")
      NULL
    })

    # Project onto the constraint space only if the solution is infeasible
    # (within numerical tolerance), using the weighted projection with the
    # final E-step weights; then recompute the log-likelihood
    if (!is.null(fit)) {
      feas <- check_constraints(fit$item_probs, constraints_spec, item_order)
      if (!feas$satisfied) {
        class_counts <- colSums(fit$posteriors)
        fit$item_probs <- project_constraints_weighted(
          fit$item_probs, constraints_spec, item_order,
          class_weights = class_counts
        )
        e_final <- e_step(data, fit$item_probs, fit$class_probs)
        fit$posteriors <- e_final$posteriors
        fit$loglik <- e_final$loglik
      }

      if (verbose) {
        cat("LL =", round(fit$loglik, 2))
        if (fit$loglik > best_ll) cat(" (new best)")
        cat("\n")
      }

      if (fit$loglik > best_ll) {
        best_fit <- fit
        best_ll <- fit$loglik
      }
    }
  }

  if (is.null(best_fit)) {
    stop("All random starts failed. Try different initialization or check data.")
  }

  # Verify constraints are satisfied
  constraint_check <- check_constraints(best_fit$item_probs, constraints_spec, item_order)
  if (!constraint_check$satisfied) {
    warning("Final solution has ", constraint_check$n_violations,
            " constraint violations. Consider re-fitting.")
  }

  # Count parameters
  n_par <- count_parameters("IIO", n_items, n_classes)

  # Create qlfit object
  result <- new_qlfit(
    model_type = "IIO",
    item_probs = best_fit$item_probs,
    class_probs = best_fit$class_probs,
    posteriors = best_fit$posteriors,
    loglik = best_fit$loglik,
    n_par = n_par,
    n_obs = n_obs,
    n_items = n_items,
    n_classes = n_classes,
    convergence = best_fit$converged,
    iterations = best_fit$iterations,
    call = call,
    theta = NULL,
    delta = NULL,
    item_order = item_order,
    constraints = constraints_spec,
    se = NULL
  )

  # Flag collapsed classes so users know the effective number of classes
  result$degenerate <- isTRUE(best_fit$degenerate)

  result
}

#' Initialize item probabilities satisfying IIO constraints
#'
#' @param data Data matrix
#' @param n_classes Number of classes
#' @param item_order Item order (easiest to hardest)
#' @param seed Random seed
#'
#' @return Item probability matrix satisfying IIO constraints
#' @keywords internal
init_item_probs_iio <- function(data, n_classes, item_order, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_items <- ncol(data)

  # Start with quantile-based initialization
  probs <- init_item_probs(data, n_classes, "quantiles", seed)

  # Ensure item ordering is satisfied within each class
  for (c in 1:n_classes) {
    # Get probabilities in item order
    ordered_probs <- probs[item_order, c]

    # Apply PAVA to make them decreasing (easier items have higher prob)
    ordered_probs <- pava_decreasing(ordered_probs)

    # Put back in original item positions
    probs[item_order, c] <- ordered_probs
  }

  # Ensure probabilities are bounded
  probs <- bound_probs(probs)

  probs
}

#' Fit IIO model using projection method
#'
#' @param data Data matrix
#' @param n_classes Number of classes
#' @param item_order Item order
#' @param n_starts Number of random starts
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param seed Random seed
#' @param verbose Print progress
#'
#' @return qlfit object
#' @keywords internal
fit_iio_projection <- function(data, n_classes,
                               item_order = NULL,
                               n_starts = 10,
                               max_iter = 500,
                               tol = 1e-6,
                               seed = NULL,
                               verbose = FALSE) {

  call <- match.call()
  data <- validate_data(data)

  n_obs <- nrow(data)
  n_items <- ncol(data)

  if (is.null(item_order)) {
    item_order <- estimate_item_order(data)
  }

  if (!is.null(seed)) set.seed(seed)

  constraints_spec <- specify_constraints(item_ordering = TRUE, item_order = item_order)

  best_fit <- NULL
  best_ll <- -Inf

  for (start in 1:n_starts) {
    if (verbose) cat("Start", start, "of", n_starts, "... ")

    start_seed <- if (!is.null(seed)) seed + start else NULL
    init_probs <- init_item_probs_iio(data, n_classes, item_order, start_seed)
    init_class_probs <- init_class_probs(n_classes, "random", start_seed)

    # Modified EM with projection
    item_probs <- init_probs
    class_probs <- init_class_probs

    ll_history <- numeric(max_iter)
    converged <- FALSE

    for (iter in 1:max_iter) {
      # E-step
      e_result <- e_step(data, item_probs, class_probs)
      posteriors <- e_result$posteriors
      ll_history[iter] <- e_result$loglik

      # Check convergence
      if (iter > 1 && check_convergence(ll_history[1:iter], tol)) {
        converged <- TRUE
        break
      }

      # Exact constrained M-step (isotonic regression per class)
      m_result <- m_step_exact(data, posteriors, constraints_spec, item_order)
      item_probs <- m_result$item_probs
      class_probs <- m_result$class_probs
    }

    ll_history <- ll_history[1:iter]

    # Consistency on max_iter exit: recompute E-step for final parameters
    if (!converged) {
      e_result <- e_step(data, item_probs, class_probs)
      posteriors <- e_result$posteriors
      ll_history <- c(ll_history, e_result$loglik)
    }

    final_ll <- ll_history[length(ll_history)]

    if (verbose) {
      cat("LL =", round(final_ll, 2))
      if (final_ll > best_ll) cat(" (new best)")
      cat("\n")
    }

    if (final_ll > best_ll) {
      best_fit <- list(
        item_probs = item_probs,
        class_probs = class_probs,
        posteriors = posteriors,
        loglik = final_ll,
        converged = converged,
        iterations = iter
      )
      best_ll <- final_ll
    }
  }

  if (is.null(best_fit)) {
    stop("All random starts failed.")
  }

  n_par <- count_parameters("IIO", n_items, n_classes)

  new_qlfit(
    model_type = "IIO",
    item_probs = best_fit$item_probs,
    class_probs = best_fit$class_probs,
    posteriors = best_fit$posteriors,
    loglik = best_fit$loglik,
    n_par = n_par,
    n_obs = n_obs,
    n_items = n_items,
    n_classes = n_classes,
    convergence = best_fit$converged,
    iterations = best_fit$iterations,
    call = call,
    theta = NULL,
    delta = NULL,
    item_order = item_order,
    constraints = constraints_spec,
    se = NULL
  )
}
