#' Fit Latent Class Rasch Model (LCR)
#'
#' @name fit_lcr
#' @description Fits a Latent Class Rasch model where item response probabilities
#'   follow the Rasch parameterization with discrete latent classes. This model
#'   represents a quantitative interpretation with discrete ability levels.
NULL

#' Fit Latent Class Rasch Model
#'
#' @param data A matrix or data frame of binary responses (0/1) with subjects in rows
#'   and items in columns.
#' @param n_classes Number of latent classes (integer >= 2)
#' @param n_starts Number of random starting values to try (default 10)
#' @param max_iter Maximum number of EM iterations per start (default 500)
#' @param tol Convergence tolerance for log-likelihood (default 1e-6)
#' @param use_cpp Use the compiled C++ E-step and M-step objective
#'   (default TRUE); the BFGS optimization itself always runs via
#'   \code{stats::optim} in R so both paths follow the identical optimizer
#'   trajectory. Set to FALSE for the pure-R reference implementation.
#' @param seed Random seed for reproducibility (optional)
#' @param verbose Print progress messages (default FALSE)
#'
#' @return A qlfit object containing:
#' \describe{
#'   \item{model_type}{"LCR" for Latent Class Rasch}
#'   \item{item_probs}{Matrix of item probabilities (n_items x n_classes) derived from
#'     Rasch parameters}
#'   \item{class_probs}{Vector of class probabilities}
#'   \item{theta}{Vector of class ability locations}
#'   \item{delta}{Vector of item difficulty parameters}
#'   \item{posteriors}{Matrix of posterior class memberships (n_obs x n_classes)}
#'   \item{loglik}{Maximized log-likelihood}
#'   \item{n_par}{Number of estimated parameters}
#'   \item{convergence}{Logical indicating successful convergence}
#' }
#'
#' @details
#' The Latent Class Rasch model parameterizes item response probabilities as:
#' \deqn{P(X_{ij} = 1 | c) = \frac{1}{1 + \exp(-(\theta_c - \delta_j))}}
#'
#' where:
#' \itemize{
#'   \item \eqn{\theta_c} is the ability level for class c
#'   \item \eqn{\delta_j} is the difficulty of item j
#' }
#'
#' This parameterization is more restrictive than the general LCA models (UN, MON, IIO, DM)
#' because item probabilities are determined by the difference between class ability
#' and item difficulty. The model has fewer free parameters
#' (\eqn{2C + I - 2} in total):
#' \itemize{
#'   \item \eqn{C-1} class proportion parameters
#'   \item \eqn{C} class location parameters (all \eqn{\theta_c} free)
#'   \item \eqn{I-1} item difficulty parameters (the single identification
#'     constraint is \eqn{\sum_j \delta_j = 0})
#' }
#'
#' The LCR model implies interval-level measurement when it fits well, as the
#' distance between class locations \eqn{\theta_c} can be interpreted on the same
#' scale as item difficulties.
#'
#' @examples
#' \dontrun{
#' # Generate Rasch-structured data
#' set.seed(123)
#' n <- 500
#' n_items <- 10
#' theta_true <- rnorm(n)
#' delta_true <- seq(-1, 1, length.out = n_items)
#'
#' data <- matrix(0, n, n_items)
#' for (i in 1:n) {
#'   for (j in 1:n_items) {
#'     prob <- 1 / (1 + exp(-(theta_true[i] - delta_true[j])))
#'     data[i, j] <- rbinom(1, 1, prob)
#'   }
#' }
#'
#' # Fit LCR model
#' fit <- fit_lcr(data, n_classes = 5)
#' print(fit)
#' summary(fit)
#'
#' # Compare with Double Monotonicity
#' dm_fit <- fit_dm(data, n_classes = 5)
#' compare_fit(fit, dm_fit)
#' }
#'
#' @seealso \code{\link{fit_rm}} for continuous Rasch model,
#'   \code{\link{fit_dm}} for double monotonicity model
#'
#' @export
fit_lcr <- function(data, n_classes,
                    n_starts = 10,
                    max_iter = 500,
                    tol = 1e-6,
                    use_cpp = TRUE,
                    seed = NULL,
                    verbose = FALSE) {

  # Capture call
  call <- match.call()

  # Dispatch polytomous data to the partial-credit latent-class engine
  if (is.data.frame(data)) data <- as.matrix(data)
  if (.is_polytomous(data)) {
    return(fit_lcr_poly(data, n_classes, n_starts = n_starts,
                        max_iter = max_iter, tol = tol, seed = seed,
                        call = call))
  }

  # Validate inputs
  data <- validate_data(data)

  n_obs <- nrow(data)
  n_items <- ncol(data)

  if (!is.numeric(n_classes) || n_classes < 2 || n_classes != round(n_classes)) {
    stop("n_classes must be an integer >= 2")
  }

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Run multiple starts
  best_fit <- NULL
  best_ll <- -Inf

  for (start in 1:n_starts) {
    if (verbose) {
      cat("Start", start, "of", n_starts, "... ")
    }

    # Initialize parameters
    start_seed <- if (!is.null(seed)) seed + start else NULL

    init_result <- init_lcr_params(data, n_classes, start_seed)

    # Run LCR-specific EM
    fit <- tryCatch({
      em_lcr(
        data = data,
        n_classes = n_classes,
        init_theta = init_result$theta,
        init_delta = init_result$delta,
        init_class_probs = init_result$class_probs,
        max_iter = max_iter,
        tol = tol,
        use_cpp = use_cpp,
        verbose = FALSE
      )
    }, error = function(e) {
      if (verbose) cat("failed: ", e$message, "\n")
      NULL
    })

    if (!is.null(fit)) {
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

  # Order classes by theta (lowest to highest ability)
  class_order <- order(best_fit$theta)
  best_fit$theta <- best_fit$theta[class_order]
  best_fit$class_probs <- best_fit$class_probs[class_order]
  best_fit$item_probs <- best_fit$item_probs[, class_order, drop = FALSE]
  best_fit$posteriors <- best_fit$posteriors[, class_order, drop = FALSE]

  # Count parameters:
  # (C-1) class probs + C theta (all free) + (I-1) delta (mean(delta) = 0
  # is the single identification constraint) = 2C + I - 2
  n_par <- count_parameters("LCR", n_items, n_classes)

  # Create qlfit object
  result <- new_qlfit(
    model_type = "LCR",
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
    theta = best_fit$theta,
    delta = best_fit$delta,
    item_order = NULL,
    constraints = NULL,
    se = NULL
  )

  # Flag collapsed classes so users know the effective number of classes
  result$degenerate <- isTRUE(best_fit$degenerate)

  result
}

#' Initialize LCR parameters
#'
#' @param data Data matrix
#' @param n_classes Number of classes
#' @param seed Random seed
#'
#' @return List with initial theta, delta, and class_probs
#' @keywords internal
init_lcr_params <- function(data, n_classes, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_items <- ncol(data)
  n_obs <- nrow(data)

  # Initialize delta from item proportions
  item_means <- colMeans(data, na.rm = TRUE)
  item_means <- pmax(pmin(item_means, 0.99), 0.01)  # Bound
  delta <- -qlogis(item_means)
  delta <- delta - mean(delta)  # Center for identification

  # Initialize theta based on total scores
  total_scores <- rowSums(data, na.rm = TRUE)

  # Use quantiles to define class boundaries
  quantile_probs <- seq(0, 1, length.out = n_classes + 1)
  score_cuts <- quantile(total_scores, probs = quantile_probs)

  theta <- numeric(n_classes)
  for (c in 1:n_classes) {
    if (c == 1) {
      in_class <- total_scores <= score_cuts[2]
    } else if (c == n_classes) {
      in_class <- total_scores > score_cuts[c]
    } else {
      in_class <- total_scores > score_cuts[c] & total_scores <= score_cuts[c + 1]
    }

    if (sum(in_class) > 0) {
      # Use mean score to initialize theta
      mean_score <- mean(total_scores[in_class])
      # Convert to logit scale (roughly)
      theta[c] <- qlogis(bound_probs(mean_score / n_items))
    } else {
      theta[c] <- qlogis((c - 0.5) / n_classes)
    }
  }

  # Add some random variation
  theta <- theta + runif(n_classes, -0.2, 0.2)

  # Initialize class probabilities
  class_probs <- rep(1 / n_classes, n_classes)

  list(
    theta = theta,
    delta = delta,
    class_probs = class_probs
  )
}

#' Fit LCR model using joint MLE
#'
#' Alternative fitting approach using direct likelihood optimization
#' rather than EM.
#'
#' @param data Data matrix
#' @param n_classes Number of classes
#' @param n_starts Number of random starts
#' @param seed Random seed
#' @param verbose Print progress
#'
#' @return qlfit object
#' @keywords internal
fit_lcr_mle <- function(data, n_classes,
                        n_starts = 10,
                        seed = NULL,
                        verbose = FALSE) {

  call <- match.call()
  data <- validate_data(data)

  n_obs <- nrow(data)
  n_items <- ncol(data)

  if (!is.null(seed)) set.seed(seed)

  best_fit <- NULL
  best_ll <- -Inf

  for (start in 1:n_starts) {
    if (verbose) cat("Start", start, "of", n_starts, "... ")

    start_seed <- if (!is.null(seed)) seed + start else NULL
    init_result <- init_lcr_params(data, n_classes, start_seed)

    # Pack parameters for optimization
    # theta: n_classes values (but we'll fix one for identification)
    # delta: n_items values (but mean = 0)
    # class_probs: use softmax parameterization

    par_init <- c(
      init_result$theta,
      init_result$delta[-1],  # Fix delta[1]
      log(init_result$class_probs[-n_classes] / init_result$class_probs[n_classes])  # Softmax
    )

    # Negative log-likelihood function
    neg_loglik <- function(par) {
      theta <- par[1:n_classes]
      delta <- c(0, par[(n_classes + 1):(n_classes + n_items - 1)])
      delta <- delta - mean(delta)  # Re-center

      # Class probs via softmax
      log_ratios <- par[(n_classes + n_items):length(par)]
      class_probs <- c(exp(log_ratios), 1)
      class_probs <- class_probs / sum(class_probs)

      # Compute item probs
      item_probs <- compute_rasch_probs(theta, delta)
      item_probs <- bound_probs(item_probs)

      # Log-likelihood
      ll <- 0
      for (i in 1:n_obs) {
        lik_i <- 0
        for (c in 1:n_classes) {
          lik_ic <- class_probs[c] *
            prod(item_probs[, c]^data[i, ] * (1 - item_probs[, c])^(1 - data[i, ]))
          lik_i <- lik_i + lik_ic
        }
        ll <- ll + log(max(lik_i, 1e-300))
      }

      -ll
    }

    # Optimize
    result <- tryCatch({
      optim(
        par = par_init,
        fn = neg_loglik,
        method = "BFGS",
        control = list(maxit = 500)
      )
    }, error = function(e) {
      if (verbose) cat("failed: ", e$message, "\n")
      NULL
    })

    if (!is.null(result) && result$convergence == 0) {
      if (verbose) {
        cat("LL =", round(-result$value, 2))
        if (-result$value > best_ll) cat(" (new best)")
        cat("\n")
      }

      if (-result$value > best_ll) {
        # Unpack parameters
        theta <- result$par[1:n_classes]
        delta <- c(0, result$par[(n_classes + 1):(n_classes + n_items - 1)])
        delta <- delta - mean(delta)

        log_ratios <- result$par[(n_classes + n_items):length(result$par)]
        class_probs <- c(exp(log_ratios), 1)
        class_probs <- class_probs / sum(class_probs)

        item_probs <- compute_rasch_probs(theta, delta)

        # Compute posteriors
        e_result <- e_step(data, item_probs, class_probs)

        best_fit <- list(
          theta = theta,
          delta = delta,
          item_probs = item_probs,
          class_probs = class_probs,
          posteriors = e_result$posteriors,
          loglik = -result$value,
          converged = TRUE,
          iterations = NA
        )
        best_ll <- -result$value
      }
    }
  }

  if (is.null(best_fit)) {
    stop("All random starts failed.")
  }

  n_par <- count_parameters("LCR", n_items, n_classes)

  new_qlfit(
    model_type = "LCR",
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
    theta = best_fit$theta,
    delta = best_fit$delta,
    item_order = NULL,
    constraints = NULL,
    se = NULL
  )
}

#' Extract Rasch-scale scores from LCR model
#'
#' @param object A qlfit object from fit_lcr
#' @param type Type of score: "eap" (expected a posteriori) or "modal"
#'
#' @return Vector of estimated person scores on the theta scale
#' @export
lcr_scores <- function(object, type = c("eap", "modal")) {
  if (object$model_type != "LCR") {
    stop("lcr_scores only applies to LCR models")
  }

  type <- match.arg(type)

  if (type == "eap") {
    # Expected a posteriori: weighted sum of theta values
    scores <- object$posteriors %*% object$theta
    as.vector(scores)
  } else if (type == "modal") {
    # Modal: theta for the most likely class
    modal_class <- apply(object$posteriors, 1, which.max)
    object$theta[modal_class]
  }
}
