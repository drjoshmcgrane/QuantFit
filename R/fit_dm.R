#' Fit Double Monotonicity Model (DM)
#'
#' @name fit_dm
#' @description Fits an ordered latent class model with both class monotonicity
#'   and invariant item ordering constraints. This model represents the strongest
#'   ordinal interpretation, related to the Mokken scale model.
NULL

#' Fit Double Monotonicity Model
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
#' @param use_cpp Use the compiled C++ EM engine (default TRUE). Set to FALSE
#'   to run the pure-R reference implementation; both paths produce
#'   numerically equivalent results.
#' @param seed Random seed for reproducibility (optional)
#' @param verbose Print progress messages (default FALSE)
#'
#' @return A qlfit object containing:
#' \describe{
#'   \item{model_type}{"DM" for double monotonicity}
#'   \item{item_probs}{Matrix of item probabilities (n_items x n_classes) satisfying
#'     both monotonicity constraints}
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
#' The double monotonicity model combines both constraints:
#'
#' \strong{Class Monotonicity (MON)}: For each item i:
#' \deqn{P(X_i = 1 | c) \leq P(X_i = 1 | c')  \text{ for } c < c'}
#'
#' \strong{Invariant Item Ordering (IIO)}: For each class c and items i < i':
#' \deqn{P(X_i = 1 | c) \geq P(X_{i'} = 1 | c)}
#'
#' This combination implies that the item response functions do not cross,
#' which is a key property of Mokken scales and suggests unidimensionality.
#'
#' The double monotonicity model is more restrictive than either MON or IIO alone,
#' and when it fits well, provides stronger evidence for an ordinal interpretation
#' of the latent variable.
#'
#' Estimation uses a projected EM algorithm whose M-step is exact: the
#' constrained maximizer of the expected complete-data log-likelihood equals
#' the 2D weighted isotonic regression of the unconstrained M-step means
#' (with expected class counts as weights), computed via Dykstra's
#' alternating-projections algorithm.
#'
#' @examples
#' \dontrun{
#' # Generate data with double monotonicity structure
#' set.seed(123)
#' n <- 500
#' data <- matrix(rbinom(n * 10, 1, 0.5), nrow = n)
#'
#' # Fit double monotonicity model
#' fit <- fit_dm(data, n_classes = 3)
#' print(fit)
#'
#' # Check constraints
#' # Class monotonicity: each row is non-decreasing
#' all(apply(fit$item_probs, 1, function(x) all(diff(x) >= 0)))
#'
#' # Item ordering: each column maintains item order
#' item_order <- fit$item_order
#' all(sapply(1:3, function(c) {
#'   all(diff(fit$item_probs[item_order, c]) <= 0)
#' }))
#' }
#'
#' @seealso \code{\link{fit_mon}} for class monotonicity only,
#'   \code{\link{fit_iio}} for item ordering only,
#'   \code{\link{fit_lcr}} for Latent Class Rasch model
#'
#' @export
fit_dm <- function(data, n_classes,
                   item_order = NULL,
                   n_starts = 10,
                   max_iter = 500,
                   tol = 1e-6,
                   use_cpp = TRUE,
                   seed = NULL,
                   verbose = FALSE) {

  # Capture call
  call <- match.call()

  # Dispatch polytomous data to the multinomial EM engine
  if (is.data.frame(data)) data <- as.matrix(data)
  if (.is_polytomous(data)) {
    return(fit_lca_poly(data, n_classes, "DM", n_starts = n_starts,
                        max_iter = max_iter, tol = tol, seed = seed,
                        item_order = item_order, use_cpp = use_cpp, call = call))
  }

  # Validate inputs
  data <- validate_data(data)

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

  constraints_spec <- specify_constraints(
    class_monotonicity = TRUE,
    item_ordering = TRUE,
    item_order = item_order
  )

  # Projected EM with the exact Dykstra M-step
  best_fit <- NULL
  best_ll <- -Inf

  for (start in 1:n_starts) {
    if (verbose) {
      cat("Start", start, "of", n_starts, "... ")
    }

    start_seed <- if (!is.null(seed)) seed + start else NULL
    init_probs <- init_item_probs_dm(data, n_classes, item_order, start_seed)
    init_class_probs <- init_class_probs(n_classes, "random", start_seed)

    fit <- tryCatch({
      em_constrained(
        data = data,
        n_classes = n_classes,
        constraints_spec = constraints_spec,
        item_order = item_order,
        init_probs = init_probs,
        init_class_probs = init_class_probs,
        max_iter = max_iter,
        tol = tol,
        method = "pava",
        use_cpp = use_cpp,
        verbose = FALSE
      )
    }, error = function(e) {
      if (verbose) cat("failed: ", e$message, "\n")
      NULL
    })

    if (!is.null(fit)) {
      # Project only if infeasible (numerical tolerance), using the weighted
      # projection with the final E-step weights; recompute LL afterwards
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
  n_par <- count_parameters("DM", n_items, n_classes)

  # Create qlfit object
  result <- new_qlfit(
    model_type = "DM",
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

#' Initialize item probabilities satisfying double monotonicity
#'
#' @param data Data matrix
#' @param n_classes Number of classes
#' @param item_order Item order (easiest to hardest)
#' @param seed Random seed
#'
#' @return Item probability matrix satisfying both DM constraints
#' @keywords internal
init_item_probs_dm <- function(data, n_classes, item_order, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_items <- ncol(data)

  # Create probabilities that satisfy both constraints by construction
  # Use item means to set relative item difficulties
  item_means <- colMeans(data, na.rm = TRUE)

  # Create grid of probabilities
  # For class c (1 to n_classes), use base probability that increases
  # For item i, shift based on item mean

  probs <- matrix(0, nrow = n_items, ncol = n_classes)

  for (c in 1:n_classes) {
    # Base probability for this class
    base_prob <- (c - 0.5) / n_classes

    for (i in 1:n_items) {
      # Item-specific adjustment
      item_adj <- (item_means[i] - 0.5) * 0.3

      # Final probability
      probs[i, c] <- base_prob + item_adj

      # Add small random noise
      probs[i, c] <- probs[i, c] + runif(1, -0.05, 0.05)
    }
  }

  # Ensure bounded
  probs <- pmax(pmin(probs, 0.95), 0.05)

  # Project onto double monotonicity constraints
  probs <- project_dm_constraints(probs, item_order)

  probs
}

#' Project onto double monotonicity constraint space
#'
#' Exact (weighted) L2 projection onto the intersection of the class
#' monotonicity and item ordering constraint sets, computed via Dykstra's
#' alternating-projections algorithm (see
#' \code{\link{dykstra_dm_projection}}).
#'
#' @param item_probs Item probability matrix (I x C)
#' @param item_order Item order vector
#' @param class_weights Optional class weights (e.g. expected class counts);
#'   NULL means unit weights
#' @param max_iter Maximum projection cycles
#' @param tol Convergence tolerance
#'
#' @return Projected item probability matrix
#' @keywords internal
project_dm_constraints <- function(item_probs, item_order,
                                   class_weights = NULL,
                                   max_iter = 500, tol = 1e-10) {
  probs <- dykstra_dm_projection(item_probs, item_order,
                                 class_weights = class_weights,
                                 tol = tol, max_cycles = max_iter)
  bound_probs(probs)
}

#' Fit DM model using projection method
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
fit_dm_projection <- function(data, n_classes,
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

  constraints_spec <- specify_constraints(
    class_monotonicity = TRUE,
    item_ordering = TRUE,
    item_order = item_order
  )

  best_fit <- NULL
  best_ll <- -Inf

  for (start in 1:n_starts) {
    if (verbose) cat("Start", start, "of", n_starts, "... ")

    start_seed <- if (!is.null(seed)) seed + start else NULL
    init_probs <- init_item_probs_dm(data, n_classes, item_order, start_seed)
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

      # Exact constrained M-step (Dykstra with expected class counts)
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

  n_par <- count_parameters("DM", n_items, n_classes)

  new_qlfit(
    model_type = "DM",
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
