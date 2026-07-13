#' Fit Class Monotonicity Model (MON)
#'
#' @name fit_mon
#' @description Fits an ordered latent class model with class monotonicity constraints.
#'   This model constrains item probabilities to be non-decreasing across classes,
#'   implying an ordinal interpretation of the latent variable.
NULL

#' Fit Class Monotonicity Model
#'
#' @param data A matrix or data frame of binary responses (0/1) with subjects in rows
#'   and items in columns.
#' @param n_classes Number of latent classes (integer >= 2)
#' @param n_starts Number of random starting values to try (default 10)
#' @param max_iter Maximum number of EM iterations per start (default 500)
#' @param tol Convergence tolerance for log-likelihood (default 1e-6)
#' @param method M-step method: "pava" (default) uses the exact constrained
#'   M-step via weighted isotonic regression (PAVA with expected class counts
#'   as weights); "optimizer" uses a general constrained optimizer
#' @param optimizer Optimizer for the constrained M-step when
#'   method = "optimizer": "alabama" (default) or "nloptr"
#' @param use_cpp Use the compiled C++ EM engine (default TRUE); only applies
#'   to method = "pava". Set to FALSE to run the pure-R reference
#'   implementation; both paths produce numerically equivalent results.
#' @param seed Random seed for reproducibility (optional)
#' @param verbose Print progress messages (default FALSE)
#'
#' @return A qlfit object containing:
#' \describe{
#'   \item{model_type}{"MON" for class monotonicity}
#'   \item{item_probs}{Matrix of item probabilities (n_items x n_classes) satisfying
#'     monotonicity constraints}
#'   \item{class_probs}{Vector of class probabilities}
#'   \item{posteriors}{Matrix of posterior class memberships (n_obs x n_classes)}
#'   \item{loglik}{Maximized log-likelihood}
#'   \item{n_par}{Number of estimated parameters}
#'   \item{convergence}{Logical indicating successful convergence}
#'   \item{constraints}{Active constraint specification}
#' }
#'
#' @details
#' The class monotonicity model adds the constraint that for each item i:
#' \deqn{P(X_i = 1 | c) \leq P(X_i = 1 | c')  \text{ for } c < c'}
#'
#' This means item success probabilities are non-decreasing across classes,
#' establishing an ordinal ordering of the latent classes.
#'
#' The model is estimated via a constrained EM algorithm. By default the
#' M-step is solved exactly: the constrained maximizer of the expected
#' complete-data log-likelihood equals the weighted isotonic regression
#' (PAVA) of the unconstrained M-step means with the expected class counts
#' as weights (Robertson, Wright & Dykstra). Alternatively,
#' method = "optimizer" uses the \code{alabama} or \code{nloptr} package
#' for general constrained optimization.
#'
#' @examples
#' \dontrun{
#' # Generate ordinal data
#' set.seed(123)
#' n <- 500
#' # Data with natural class ordering
#' data <- matrix(rbinom(n * 10, 1, 0.5), nrow = n)
#'
#' # Fit class monotonicity model
#' fit <- fit_mon(data, n_classes = 3)
#' print(fit)
#'
#' # Check that monotonicity holds
#' all(apply(fit$item_probs, 1, function(x) all(diff(x) >= 0)))
#' }
#'
#' @seealso \code{\link{fit_un}} for unconstrained model,
#'   \code{\link{fit_iio}} for item ordering constraints,
#'   \code{\link{fit_dm}} for double monotonicity
#'
#' @export
fit_mon <- function(data, n_classes,
                    n_starts = 10,
                    max_iter = 500,
                    tol = 1e-6,
                    method = c("pava", "optimizer"),
                    optimizer = c("alabama", "nloptr"),
                    use_cpp = TRUE,
                    seed = NULL,
                    verbose = FALSE) {

  # Capture call
  call <- match.call()

  # Dispatch polytomous data to the multinomial EM engine
  if (is.data.frame(data)) data <- as.matrix(data)
  if (.needs_poly_engine(data)) {
    return(fit_lca_poly(data, n_classes, "MON", n_starts = n_starts,
                        max_iter = max_iter, tol = tol, seed = seed, use_cpp = use_cpp,
                        call = call))
  }

  # Validate inputs
  data <- validate_data(data)
  method <- match.arg(method)
  optimizer <- match.arg(optimizer)

  n_obs <- nrow(data)
  n_items <- ncol(data)

  if (!is.numeric(n_classes) || n_classes < 2 || n_classes != round(n_classes)) {
    stop("n_classes must be an integer >= 2")
  }

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Build constraint function
  constraint_fn <- build_class_monotonicity_fn(n_items, n_classes)
  constraints_spec <- specify_constraints(class_monotonicity = TRUE)

  # Run multiple starts
  best_fit <- NULL
  best_ll <- -Inf

  for (start in 1:n_starts) {
    if (verbose) {
      cat("Start", start, "of", n_starts, "... ")
    }

    # Initialize with monotonic starting values
    start_seed <- if (!is.null(seed)) seed + start else NULL

    init_probs <- init_item_probs_monotonic(data, n_classes, start_seed)
    init_class_probs <- init_class_probs(n_classes, "random", start_seed)

    # Run constrained EM
    fit <- tryCatch({
      em_constrained(
        data = data,
        n_classes = n_classes,
        constraints_spec = constraints_spec,
        constraint_fn = constraint_fn,
        init_probs = init_probs,
        init_class_probs = init_class_probs,
        max_iter = max_iter,
        tol = tol,
        method = method,
        optimizer = optimizer,
        use_cpp = use_cpp,
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
      feas <- check_constraints(fit$item_probs, constraints_spec)
      if (!feas$satisfied) {
        class_counts <- colSums(fit$posteriors)
        fit$item_probs <- project_constraints_weighted(
          fit$item_probs, constraints_spec, class_weights = class_counts
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
  constraint_check <- check_constraints(best_fit$item_probs, constraints_spec)
  if (!constraint_check$satisfied) {
    warning("Final solution has ", constraint_check$n_violations,
            " constraint violations. Consider re-fitting.")
  }

  # Count parameters (same as UN - constraints don't reduce parameters)
  n_par <- count_parameters("MON", n_items, n_classes)

  # Create qlfit object
  result <- new_qlfit(
    model_type = "MON",
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
    item_order = NULL,
    constraints = constraints_spec,
    se = NULL
  )

  # Flag collapsed classes so users know the effective number of classes
  result$degenerate <- isTRUE(best_fit$degenerate)

  result
}

#' Fit MON model using projection method
#'
#' Alternative fitting approach that uses the exact weighted-PAVA M-step
#' (projection of the unconstrained M-step means in the weighted L2 metric)
#' after each E-step. Equivalent to the default path of \code{fit_mon}.
#'
#' @param data Data matrix
#' @param n_classes Number of classes
#' @param n_starts Number of random starts
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param seed Random seed
#' @param verbose Print progress
#'
#' @return qlfit object
#' @keywords internal
fit_mon_projection <- function(data, n_classes,
                               n_starts = 10,
                               max_iter = 500,
                               tol = 1e-6,
                               seed = NULL,
                               verbose = FALSE) {

  call <- match.call()
  data <- validate_data(data)

  n_obs <- nrow(data)
  n_items <- ncol(data)

  if (!is.null(seed)) set.seed(seed)

  constraints_spec <- specify_constraints(class_monotonicity = TRUE)

  best_fit <- NULL
  best_ll <- -Inf

  for (start in 1:n_starts) {
    if (verbose) cat("Start", start, "of", n_starts, "... ")

    start_seed <- if (!is.null(seed)) seed + start else NULL
    init_probs <- init_item_probs_monotonic(data, n_classes, start_seed)
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

      # Exact constrained M-step (weighted PAVA)
      m_result <- m_step_exact(data, posteriors, constraints_spec)
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

  n_par <- count_parameters("MON", n_items, n_classes)

  new_qlfit(
    model_type = "MON",
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
    item_order = NULL,
    constraints = constraints_spec,
    se = NULL
  )
}
