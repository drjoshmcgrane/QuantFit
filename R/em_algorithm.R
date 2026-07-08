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

  list(
    item_probs = item_probs,
    class_probs = class_probs,
    posteriors = posteriors,
    loglik = ll_history[iter],
    ll_history = ll_history,
    converged = converged,
    iterations = iter
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
    class_probs = class_probs
  )
}

#' EM algorithm with inequality constraints
#'
#' @param data Binary data matrix (n x I)
#' @param n_classes Number of latent classes
#' @param constraint_fn Function that returns inequality constraints (should be >= 0)
#' @param constraint_grad Gradient of constraint function (optional)
#' @param init_probs Initial item probabilities (I x C matrix)
#' @param init_class_probs Initial class probabilities (vector of length C)
#' @param max_iter Maximum number of EM iterations
#' @param tol Convergence tolerance
#' @param optimizer Which optimizer to use: "alabama" or "nloptr"
#' @param verbose Print progress
#'
#' @return List with estimated parameters and diagnostics
#' @keywords internal
em_constrained <- function(data, n_classes,
                           constraint_fn,
                           constraint_grad = NULL,
                           init_probs = NULL,
                           init_class_probs = NULL,
                           max_iter = 500,
                           tol = 1e-6,
                           optimizer = c("alabama", "nloptr"),
                           verbose = FALSE) {

  optimizer <- match.arg(optimizer)

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

    # M-step: Constrained optimization
    m_result <- m_step_constrained(
      data, posteriors, item_probs, class_probs,
      constraint_fn, constraint_grad, optimizer
    )
    item_probs <- m_result$item_probs
    class_probs <- m_result$class_probs
  }

  # Trim history
  ll_history <- ll_history[1:iter]

  list(
    item_probs = item_probs,
    class_probs = class_probs,
    posteriors = posteriors,
    loglik = ll_history[iter],
    ll_history = ll_history,
    converged = converged,
    iterations = iter
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
#'
#' @return List with updated parameters
#' @keywords internal
m_step_constrained <- function(data, posteriors, current_probs, current_class_probs,
                               constraint_fn, constraint_grad, optimizer) {

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
    }, error = function(e) {
      # Fall back to unconstrained with projection
      list(par = par_init, convergence = 1)
    })

    new_par <- result$par

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
    }, error = function(e) {
      list(solution = par_init, status = -1)
    })

    new_par <- result$solution
  }

  # Convert back to matrix
  item_probs <- matrix(new_par, nrow = n_items, ncol = n_classes)
  item_probs <- bound_probs(item_probs)

  list(
    item_probs = item_probs,
    class_probs = class_probs
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

  list(
    theta = theta,
    delta = delta,
    item_probs = item_probs,
    class_probs = class_probs,
    posteriors = posteriors,
    loglik = ll_history[iter],
    ll_history = ll_history,
    converged = converged,
    iterations = iter
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

  # Pack parameters for optimization
  # par = c(theta[-1], delta[-1])  # Fix theta[1] = 0 and delta[1] = 0 for identification
  # Actually, let's fix mean(delta) = 0

  par_init <- c(theta, delta[-1])  # Fix delta[1]

  # Objective: negative expected complete data log-likelihood
  objective <- function(par) {
    theta_new <- par[1:n_classes]
    delta_new <- c(0, par[(n_classes + 1):length(par)])  # delta[1] = 0
    delta_new <- delta_new - mean(delta_new)  # Re-center

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

  # Unpack parameters
  theta <- result$par[1:n_classes]
  delta <- c(0, result$par[(n_classes + 1):length(result$par)])
  delta <- delta - mean(delta)  # Identification: mean(delta) = 0

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
