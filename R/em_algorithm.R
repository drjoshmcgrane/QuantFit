#' Core EM Algorithm for Latent Class Analysis
#'
#' @name em_algorithm
#' @description Core EM algorithm implementation with support for constrained optimization
NULL

#' Standard EM algorithm for Latent Class Analysis
#'
#' @param data Binary data matrix (n x I)
#' @param n_classes Number of latent classes
#' @param init_probs Initial item probabilities (I x C matrix)
#' @param init_class_probs Initial class probabilities (vector of length C)
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance
#' @param verbose Print progress
#'
#' @return List with estimated parameters and diagnostics
#' @keywords internal
em_lca <- function(data, n_classes,
                   init_probs = NULL,
                   init_class_probs = NULL,
                   max_iter = 1000,
                   tol = 1e-6,
                   verbose = FALSE) {

  n_obs <- nrow(data)
  n_items <- ncol(data)

  # Initialize parameters
  if (is.null(init_probs)) {
    item_probs <- init_item_probs(data, n_classes, "quantiles")
  } else {
    item_probs <- init_probs
  }

  if (is.null(init_class_probs)) {
    class_probs <- init_class_probs(n_classes, "uniform")
  } else {
    class_probs <- init_class_probs
  }

  # EM iterations
  ll_history <- numeric(max_iter)
  converged <- FALSE

  for (iter in 1:max_iter) {
    # E-step: Compute posteriors
    e_result <- e_step(data, item_probs, class_probs)
    posteriors <- e_result$posteriors
    ll_history[iter] <- e_result$loglik

    if (verbose && iter %% 10 == 0) {
      cat("Iteration", iter, "- Log-likelihood:", round(ll_history[iter], 4), "\n")
    }

    # Check convergence
    if (iter > 1 && check_convergence(ll_history[1:iter], tol)) {
      converged <- TRUE
      break
    }

    # M-step: Update parameters
    m_result <- m_step(data, posteriors)
    item_probs <- m_result$item_probs
    class_probs <- m_result$class_probs
  }

  # Trim history
  ll_history <- ll_history[1:iter]

  # If we exited on max_iter, the parameters were updated after the last
  # E-step; run a final E-step so loglik/posteriors match returned parameters
  if (!converged) {
    e_result <- e_step(data, item_probs, class_probs)
    posteriors <- e_result$posteriors
    ll_history <- c(ll_history, e_result$loglik)
  }

  degenerate <- any(colSums(posteriors) < 1)
  if (degenerate) {
    warning("One or more latent classes collapsed (expected class count < 1). ",
            "The effective number of classes is smaller than requested.",
            call. = FALSE)
  }

  list(
    item_probs = item_probs,
    class_probs = class_probs,
    posteriors = posteriors,
    loglik = ll_history[length(ll_history)],
    ll_history = ll_history,
    converged = converged,
    iterations = iter,
    degenerate = degenerate
  )
}

#' E-step: Compute posterior class memberships
#'
#' @param data Binary data matrix (n x I)
#' @param item_probs Item probability matrix (I x C)
#' @param class_probs Class probability vector (length C)
#'
#' @return List with posteriors and log-likelihood
#' @keywords internal
e_step <- function(data, item_probs, class_probs) {
  n_obs <- nrow(data)
  n_classes <- length(class_probs)

  # Compute log-likelihood contribution for each observation and class
  log_lik_mat <- matrix(0, nrow = n_obs, ncol = n_classes)

  for (c in 1:n_classes) {
    # Log probability of data given class c
    # P(x_i | c) = prod_j [ p_jc^x_ij * (1-p_jc)^(1-x_ij) ]
    log_p <- log(bound_probs(item_probs[, c]))
    log_1mp <- log(bound_probs(1 - item_probs[, c]))

    # For each observation: sum of log probabilities
    log_lik_mat[, c] <- data %*% log_p + (1 - data) %*% log_1mp + log(class_probs[c])
  }

  # Compute posteriors using log-sum-exp for numerical stability
  log_row_sums <- row_log_sum_exp(log_lik_mat)
  log_posteriors <- log_lik_mat - log_row_sums
  posteriors <- exp(log_posteriors)

  # Compute total log-likelihood
  loglik <- sum(log_row_sums)

  list(
    posteriors = posteriors,
    loglik = loglik
  )
}

#' M-step: Update parameters (unconstrained)
#'
#' @param data Binary data matrix (n x I)
#' @param posteriors Posterior class membership matrix (n x C)
#'
#' @return List with updated item_probs and class_probs
#' @keywords internal
m_step <- function(data, posteriors) {
  n_obs <- nrow(data)
  n_items <- ncol(data)
  n_classes <- ncol(posteriors)

  # Update class probabilities
  class_probs <- colSums(posteriors) / n_obs

  # Update item probabilities
  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)
  # Flag classes whose expected count has collapsed below ~1 observation:
  # their item probabilities are essentially unidentified and the effective
  # number of classes is smaller than requested.
  degenerate <- any(colSums(posteriors) < 1)
  for (c in 1:n_classes) {
    # Weighted mean of data
    weights <- posteriors[, c]
    if (sum(weights) > 0) {
      item_probs[, c] <- colSums(data * weights) / sum(weights)
    } else {
      item_probs[, c] <- 0.5
    }
  }

  # Bound probabilities
  item_probs <- bound_probs(item_probs)

  list(
    item_probs = item_probs,
    class_probs = class_probs,
    degenerate = degenerate
  )
}

#' Exact constrained M-step via weighted isotonic regression
#'
#' The constrained M-step maximizes
#' \eqn{Q(p) = \sum_c N_c [\bar{p}_{ic} \log p_{ic} +
#' (1-\bar{p}_{ic}) \log(1-p_{ic})]} where \eqn{N_c} are the expected class
#' counts and \eqn{\bar{p}_{ic}} the unconstrained weighted means. By the
#' exponential-family isotonic MLE theorem (Robertson, Wright & Dykstra),
#' the constrained maximizer equals the weighted L2 isotonic regression of
#' \eqn{\bar{p}} with weights \eqn{N_c}, computed by
#' \code{\link{project_constraints_weighted}}. This M-step is exact, so the
#' EM retains its monotone ascent property.
#'
#' @param data Binary data matrix
#' @param posteriors Posterior class memberships (n x C)
#' @param constraints_spec A ql_constraints object
#' @param item_order Item order vector (needed for item ordering constraints)
#'
#' @return List with updated item_probs, class_probs, degenerate flag
#' @keywords internal
m_step_exact <- function(data, posteriors, constraints_spec, item_order = NULL) {
  m_result <- m_step(data, posteriors)
  class_counts <- colSums(posteriors)

  item_probs <- project_constraints_weighted(
    m_result$item_probs, constraints_spec, item_order,
    class_weights = class_counts
  )

  list(
    item_probs = item_probs,
    class_probs = m_result$class_probs,
    degenerate = m_result$degenerate
  )
}

#' EM algorithm with inequality constraints
#'
#' @param data Binary data matrix (n x I)
#' @param n_classes Number of latent classes
#' @param constraints_spec A ql_constraints object describing the constraints
#' @param item_order Item order vector (needed for item ordering constraints)
#' @param constraint_fn Function that returns inequality constraints
#'   (should be >= 0); only needed for method = "optimizer"
#' @param constraint_grad Gradient of constraint function (optional)
#' @param init_probs Initial item probabilities (I x C matrix)
#' @param init_class_probs Initial class probabilities (vector of length C)
#' @param max_iter Maximum number of EM iterations
#' @param tol Convergence tolerance
#' @param method M-step method: "pava" (default; exact weighted isotonic
#'   regression M-step) or "optimizer" (general constrained optimizer)
#' @param optimizer Which optimizer to use when method = "optimizer":
#'   "alabama" or "nloptr"
#' @param verbose Print progress
#'
#' @return List with estimated parameters and diagnostics
#' @keywords internal
em_constrained <- function(data, n_classes,
                           constraints_spec,
                           item_order = NULL,
                           constraint_fn = NULL,
                           constraint_grad = NULL,
                           init_probs = NULL,
                           init_class_probs = NULL,
                           max_iter = 500,
                           tol = 1e-6,
                           method = c("pava", "optimizer"),
                           optimizer = c("alabama", "nloptr"),
                           verbose = FALSE) {

  method <- match.arg(method)
  optimizer <- match.arg(optimizer)

  n_obs <- nrow(data)
  n_items <- ncol(data)

  if (method == "optimizer" && is.null(constraint_fn)) {
    constraint_fn <- build_constraint_fn(constraints_spec, n_items, n_classes,
                                         item_order)
  }

  # Initialize parameters
  if (is.null(init_probs)) {
    item_probs <- init_item_probs(data, n_classes, "quantiles")
  } else {
    item_probs <- init_probs
  }

  if (is.null(init_class_probs)) {
    class_probs <- init_class_probs(n_classes, "uniform")
  } else {
    class_probs <- init_class_probs
  }

  # EM iterations
  ll_history <- numeric(max_iter)
  converged <- FALSE
  optimizer_warned <- FALSE

  for (iter in 1:max_iter) {
    # E-step: Compute posteriors (same as unconstrained)
    e_result <- e_step(data, item_probs, class_probs)
    posteriors <- e_result$posteriors
    ll_history[iter] <- e_result$loglik

    if (verbose && iter %% 10 == 0) {
      cat("Iteration", iter, "- Log-likelihood:", round(ll_history[iter], 4), "\n")
    }

    # Check convergence
    if (iter > 1 && check_convergence(ll_history[1:iter], tol)) {
      converged <- TRUE
      break
    }

    # M-step
    if (method == "pava") {
      # Exact constrained M-step via weighted isotonic regression
      m_result <- m_step_exact(data, posteriors, constraints_spec, item_order)
    } else {
      m_result <- m_step_constrained(
        data, posteriors, item_probs, class_probs,
        constraint_fn, constraint_grad, optimizer,
        constraints_spec = constraints_spec, item_order = item_order
      )
      if (isTRUE(m_result$optimizer_failed) && !optimizer_warned) {
        warning("Constrained optimizer failed in the M-step; ",
                "falling back to the exact weighted-PAVA M-step.",
                call. = FALSE)
        optimizer_warned <- TRUE
      }
    }
    item_probs <- m_result$item_probs
    class_probs <- m_result$class_probs
  }

  # Trim history
  ll_history <- ll_history[1:iter]

  # If we exited on max_iter, the parameters were updated after the last
  # E-step; run a final E-step so loglik/posteriors match returned parameters
  if (!converged) {
    e_result <- e_step(data, item_probs, class_probs)
    posteriors <- e_result$posteriors
    ll_history <- c(ll_history, e_result$loglik)
  }

  degenerate <- any(colSums(posteriors) < 1)
  if (degenerate) {
    warning("One or more latent classes collapsed (expected class count < 1). ",
            "The effective number of classes is smaller than requested.",
            call. = FALSE)
  }

  list(
    item_probs = item_probs,
    class_probs = class_probs,
    posteriors = posteriors,
    loglik = ll_history[length(ll_history)],
    ll_history = ll_history,
    converged = converged,
    iterations = iter,
    degenerate = degenerate
  )
}

#' Constrained M-step using optimization
#'
#' @param data Binary data matrix
#' @param posteriors Posterior class memberships
#' @param current_probs Current item probabilities (starting values)
#' @param current_class_probs Current class probabilities
#' @param constraint_fn Constraint function
#' @param constraint_grad Gradient of constraints
#' @param optimizer Optimizer choice
#' @param constraints_spec ql_constraints object used for the exact
#'   weighted-PAVA fallback if the optimizer fails
#' @param item_order Item order vector for the fallback (if needed)
#'
#' @return List with updated parameters and an optimizer_failed flag
#' @keywords internal
m_step_constrained <- function(data, posteriors, current_probs, current_class_probs,
                               constraint_fn, constraint_grad, optimizer,
                               constraints_spec = NULL, item_order = NULL) {

  n_items <- ncol(data)
  n_classes <- ncol(posteriors)
  n_obs <- nrow(data)

  # Update class probabilities (no constraints on these typically)
  class_probs <- colSums(posteriors) / n_obs

  # Convert item_probs to parameter vector (logit scale for unconstrained optimization)
  # Work on probability scale directly but use bounds
  par_init <- as.vector(current_probs)

  # Define objective: negative expected complete data log-likelihood for item probs
  objective_fn <- function(par) {
    item_probs <- matrix(par, nrow = n_items, ncol = n_classes)
    item_probs <- bound_probs(item_probs)

    # Expected complete data log-likelihood for item parameters
    ll <- 0
    for (c in 1:n_classes) {
      weights <- posteriors[, c]
      log_p <- log(item_probs[, c])
      log_1mp <- log(1 - item_probs[, c])

      # Contribution from this class
      ll <- ll + sum(weights * (data %*% log_p + (1 - data) %*% log_1mp))
    }

    -ll  # Return negative for minimization
  }

  # Gradient of objective (analytical for speed)
  objective_grad <- function(par) {
    item_probs <- matrix(par, nrow = n_items, ncol = n_classes)
    item_probs <- bound_probs(item_probs)

    grad <- matrix(0, nrow = n_items, ncol = n_classes)
    for (c in 1:n_classes) {
      weights <- posteriors[, c]
      # d/dp [-sum(w * (x*log(p) + (1-x)*log(1-p)))]
      # = -sum(w * (x/p - (1-x)/(1-p)))
      for (i in 1:n_items) {
        grad[i, c] <- -sum(weights * (data[, i] / item_probs[i, c] -
                                      (1 - data[, i]) / (1 - item_probs[i, c])))
      }
    }

    as.vector(grad)
  }

  # Wrapper for constraint function expecting vector input
  hin <- function(par) {
    item_probs <- matrix(par, nrow = n_items, ncol = n_classes)
    constraint_fn(item_probs)
  }

  # Bounds: probabilities between eps and 1-eps
  eps <- 1e-6
  lower <- rep(eps, length(par_init))
  upper <- rep(1 - eps, length(par_init))

  # Run constrained optimization
  optimizer_failed <- FALSE

  if (optimizer == "alabama") {
    result <- tryCatch({
      alabama::constrOptim.nl(
        par = par_init,
        fn = objective_fn,
        gr = objective_grad,
        hin = hin,
        control.outer = list(trace = FALSE, eps = 1e-6),
        control.optim = list(maxit = 100)
      )
    }, error = function(e) NULL)

    if (is.null(result)) {
      optimizer_failed <- TRUE
    } else {
      new_par <- result$par
    }

  } else if (optimizer == "nloptr") {
    result <- tryCatch({
      nloptr::nloptr(
        x0 = par_init,
        eval_f = objective_fn,
        eval_grad_f = objective_grad,
        eval_g_ineq = function(x) -hin(x),  # nloptr wants g(x) <= 0
        lb = lower,
        ub = upper,
        opts = list(
          algorithm = "NLOPT_LD_SLSQP",
          xtol_rel = 1e-6,
          maxeval = 100,
          print_level = 0
        )
      )
    }, error = function(e) NULL)

    if (is.null(result) || result$status < 0) {
      optimizer_failed <- TRUE
    } else {
      new_par <- result$solution
    }
  }

  if (optimizer_failed) {
    # Fall back to the exact weighted-PAVA/Dykstra M-step rather than
    # silently returning the unchanged parameters (which stalled the EM)
    if (!is.null(constraints_spec)) {
      m_exact <- m_step_exact(data, posteriors, constraints_spec, item_order)
      return(list(
        item_probs = m_exact$item_probs,
        class_probs = class_probs,
        optimizer_failed = TRUE
      ))
    }
    # No constraint specification available: keep current values but signal
    new_par <- par_init
  }

  # Convert back to matrix
  item_probs <- matrix(new_par, nrow = n_items, ncol = n_classes)
  item_probs <- bound_probs(item_probs)

  list(
    item_probs = item_probs,
    class_probs = class_probs,
    optimizer_failed = optimizer_failed
  )
}

#' EM algorithm for Latent Class Rasch model
#'
#' @param data Binary data matrix (n x I)
#' @param n_classes Number of latent classes
#' @param init_theta Initial class locations (vector of length C)
#' @param init_delta Initial item difficulties (vector of length I)
#' @param init_class_probs Initial class probabilities (vector of length C)
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance
#' @param verbose Print progress
#'
#' @return List with estimated parameters and diagnostics
#' @keywords internal
em_lcr <- function(data, n_classes,
                   init_theta = NULL,
                   init_delta = NULL,
                   init_class_probs = NULL,
                   max_iter = 500,
                   tol = 1e-6,
                   verbose = FALSE) {

  n_obs <- nrow(data)
  n_items <- ncol(data)

  # Initialize parameters
  if (is.null(init_theta)) {
    # Spread theta values across reasonable range
    theta <- seq(-2, 2, length.out = n_classes)
  } else {
    theta <- init_theta
  }

  if (is.null(init_delta)) {
    # Initialize based on item proportions
    item_means <- colMeans(data, na.rm = TRUE)
    delta <- -qlogis(bound_probs(item_means))
  } else {
    delta <- init_delta
  }

  if (is.null(init_class_probs)) {
    class_probs <- init_class_probs(n_classes, "uniform")
  } else {
    class_probs <- init_class_probs
  }

  # Identification constraint: mean(theta) = 0 or theta[1] = 0
  # We'll use mean(delta) = 0
  delta <- delta - mean(delta)

  # EM iterations
  ll_history <- numeric(max_iter)
  converged <- FALSE

  for (iter in 1:max_iter) {
    # Compute item probabilities from Rasch parameters
    item_probs <- compute_rasch_probs(theta, delta)

    # E-step
    e_result <- e_step(data, item_probs, class_probs)
    posteriors <- e_result$posteriors
    ll_history[iter] <- e_result$loglik

    if (verbose && iter %% 10 == 0) {
      cat("Iteration", iter, "- Log-likelihood:", round(ll_history[iter], 4), "\n")
    }

    # Check convergence
    if (iter > 1 && check_convergence(ll_history[1:iter], tol)) {
      converged <- TRUE
      break
    }

    # M-step: Update theta, delta, and class_probs
    m_result <- m_step_rasch(data, posteriors, theta, delta, class_probs)
    theta <- m_result$theta
    delta <- m_result$delta
    class_probs <- m_result$class_probs
  }

  # Final item probabilities
  item_probs <- compute_rasch_probs(theta, delta)

  # Trim history
  ll_history <- ll_history[1:iter]

  # If we exited on max_iter, the parameters were updated after the last
  # E-step; run a final E-step so loglik/posteriors match returned parameters
  if (!converged) {
    e_result <- e_step(data, item_probs, class_probs)
    posteriors <- e_result$posteriors
    ll_history <- c(ll_history, e_result$loglik)
  }

  degenerate <- any(colSums(posteriors) < 1)
  if (degenerate) {
    warning("One or more latent classes collapsed (expected class count < 1). ",
            "The effective number of classes is smaller than requested.",
            call. = FALSE)
  }

  list(
    theta = theta,
    delta = delta,
    item_probs = item_probs,
    class_probs = class_probs,
    posteriors = posteriors,
    loglik = ll_history[length(ll_history)],
    ll_history = ll_history,
    converged = converged,
    iterations = iter,
    degenerate = degenerate
  )
}

#' Compute item probabilities from Rasch parameters
#'
#' @param theta Class locations (vector of length C)
#' @param delta Item difficulties (vector of length I)
#'
#' @return Item probability matrix (I x C)
#' @keywords internal
compute_rasch_probs <- function(theta, delta) {
  n_items <- length(delta)
  n_classes <- length(theta)

  probs <- matrix(0, nrow = n_items, ncol = n_classes)
  for (c in 1:n_classes) {
    probs[, c] <- inv_logit(theta[c] - delta)
  }

  probs
}

#' M-step for Latent Class Rasch model
#'
#' @param data Binary data matrix
#' @param posteriors Posterior class memberships
#' @param theta Current theta values
#' @param delta Current delta values
#' @param class_probs Current class probabilities
#'
#' @return List with updated parameters
#' @keywords internal
m_step_rasch <- function(data, posteriors, theta, delta, class_probs) {
  n_obs <- nrow(data)
  n_items <- ncol(data)
  n_classes <- ncol(posteriors)

  # Update class probabilities
  class_probs <- colSums(posteriors) / n_obs

  # Free parameters: theta (all C) and delta_2..delta_J, with
  # delta_1 = -sum(delta_2..delta_J) implied by the mean(delta) = 0
  # identification constraint. Because the current delta already satisfies
  # mean(delta) = 0, the objective evaluates exactly the current Q at
  # par_init, preserving the generalized-EM ascent property.
  par_init <- c(theta, delta[-1])

  # Objective: negative expected complete data log-likelihood
  objective <- function(par) {
    theta_new <- par[1:n_classes]
    delta_free <- par[(n_classes + 1):length(par)]
    delta_new <- c(-sum(delta_free), delta_free)  # mean(delta) = 0 implied

    item_probs <- compute_rasch_probs(theta_new, delta_new)
    item_probs <- bound_probs(item_probs)

    ll <- 0
    for (c in 1:n_classes) {
      weights <- posteriors[, c]
      log_p <- log(item_probs[, c])
      log_1mp <- log(1 - item_probs[, c])
      ll <- ll + sum(weights * (data %*% log_p + (1 - data) %*% log_1mp))
    }

    -ll
  }

  # Optimize
  result <- optim(
    par = par_init,
    fn = objective,
    method = "BFGS",
    control = list(maxit = 50)
  )

  # Guard the generalized-EM ascent property: only accept the update if it
  # does not worsen the expected complete-data log-likelihood
  if (result$value <= objective(par_init) + 1e-10) {
    theta <- result$par[1:n_classes]
    delta_free <- result$par[(n_classes + 1):length(result$par)]
    delta <- c(-sum(delta_free), delta_free)
  }

  list(
    theta = theta,
    delta = delta,
    class_probs = class_probs
  )
}

#' Compute observed data log-likelihood
#'
#' @param data Binary data matrix
#' @param item_probs Item probability matrix (I x C)
#' @param class_probs Class probability vector
#'
#' @return Log-likelihood value
#' @keywords internal
compute_loglik <- function(data, item_probs, class_probs) {
  e_step(data, item_probs, class_probs)$loglik
}
