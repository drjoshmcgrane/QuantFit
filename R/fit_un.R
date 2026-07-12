#' Fit Unconstrained Latent Class Model (UN)
#'
#' @name fit_un
#' @description Fits an unconstrained latent class analysis model using the EM algorithm.
#'   This is the baseline model in the Torres Irribarra & Diakow framework with no
#'   ordering constraints on class or item probabilities.
NULL

#' Fit Unconstrained Latent Class Model
#'
#' @param data A matrix or data frame of binary responses (0/1) with subjects in rows
#'   and items in columns.
#' @param n_classes Number of latent classes (integer >= 2)
#' @param n_starts Number of random starting values to try (default 10)
#' @param max_iter Maximum number of EM iterations per start (default 1000)
#' @param tol Convergence tolerance for log-likelihood (default 1e-6)
#' @param init_method Initialization method: "quantiles" (default), "kmeans", or "random"
#' @param use_cpp Use the compiled C++ EM engine (default TRUE). Set to FALSE
#'   to run the pure-R reference implementation; both paths produce
#'   numerically equivalent results.
#' @param seed Random seed for reproducibility (optional)
#' @param verbose Print progress messages (default FALSE)
#'
#' @return A qlfit object containing:
#' \describe{
#'   \item{model_type}{"UN" for unconstrained latent class}
#'   \item{item_probs}{Matrix of item probabilities (n_items x n_classes)}
#'   \item{class_probs}{Vector of class probabilities}
#'   \item{posteriors}{Matrix of posterior class memberships (n_obs x n_classes)}
#'   \item{loglik}{Maximized log-likelihood}
#'   \item{n_par}{Number of estimated parameters}
#'   \item{convergence}{Logical indicating successful convergence}
#'   \item{iterations}{Number of EM iterations}
#' }
#'
#' @details
#' The unconstrained latent class model assumes that:
#' \itemize{
#'   \item Subjects belong to one of K latent classes
#'   \item Items are conditionally independent given class membership
#'   \item No ordering constraints on item probabilities across classes
#' }
#'
#' The model is estimated via the EM algorithm with multiple random starts
#' to avoid local maxima. The solution with the highest log-likelihood is returned.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 500
#' data <- matrix(rbinom(n * 10, 1, 0.5), nrow = n)
#'
#' # Fit 3-class model
#' fit <- fit_un(data, n_classes = 3)
#' print(fit)
#' summary(fit)
#'
#' # Get class assignments
#' assignments <- class_assignments(fit)
#' }
#'
#' @export
fit_un <- function(data, n_classes,
                   n_starts = 10,
                   max_iter = 1000,
                   tol = 1e-6,
                   init_method = c("quantiles", "kmeans", "random"),
                   use_cpp = TRUE,
                   seed = NULL,
                   verbose = FALSE) {

  # Capture call
  call <- match.call()

  # Dispatch polytomous data to the multinomial EM engine
  if (is.data.frame(data)) data <- as.matrix(data)
  if (.is_polytomous(data)) {
    return(fit_lca_poly(data, n_classes, "UN", n_starts = n_starts,
                        max_iter = max_iter, tol = tol, seed = seed, use_cpp = use_cpp,
                        call = call))
  }

  # Validate inputs
  data <- validate_data(data)
  init_method <- match.arg(init_method)

  n_obs <- nrow(data)
  n_items <- ncol(data)

  if (!is.numeric(n_classes) || n_classes < 2 || n_classes != round(n_classes)) {
    stop("n_classes must be an integer >= 2")
  }

  if (n_classes > n_obs / 2) {
    warning("Number of classes is large relative to sample size")
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

    init_probs <- tryCatch({
      init_item_probs(data, n_classes, init_method, start_seed)
    }, error = function(e) {
      init_item_probs(data, n_classes, "random", start_seed)
    })

    init_class_probs <- init_class_probs(n_classes, "random", start_seed)

    # Run EM
    fit <- tryCatch({
      em_lca(
        data = data,
        n_classes = n_classes,
        init_probs = init_probs,
        init_class_probs = init_class_probs,
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

  # Order classes by prevalence (largest first)
  class_order <- order(best_fit$class_probs, decreasing = TRUE)
  best_fit$class_probs <- best_fit$class_probs[class_order]
  best_fit$item_probs <- best_fit$item_probs[, class_order, drop = FALSE]
  best_fit$posteriors <- best_fit$posteriors[, class_order, drop = FALSE]

  # Count parameters
  n_par <- count_parameters("UN", n_items, n_classes)

  # Create qlfit object
  result <- new_qlfit(
    model_type = "UN",
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
    constraints = NULL,
    se = NULL
  )

  result
}

#' Fit UN model with specific starting values
#'
#' @param data Data matrix
#' @param n_classes Number of classes
#' @param init_probs Starting item probabilities
#' @param init_class_probs Starting class probabilities
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param use_cpp Use the compiled C++ EM engine (default TRUE)
#' @param verbose Print progress
#'
#' @return qlfit object
#' @keywords internal
fit_un_single <- function(data, n_classes,
                          init_probs,
                          init_class_probs,
                          max_iter = 1000,
                          tol = 1e-6,
                          use_cpp = TRUE,
                          verbose = FALSE) {

  data <- validate_data(data)
  n_obs <- nrow(data)
  n_items <- ncol(data)

  # Run EM
  fit <- em_lca(
    data = data,
    n_classes = n_classes,
    init_probs = init_probs,
    init_class_probs = init_class_probs,
    max_iter = max_iter,
    tol = tol,
    use_cpp = use_cpp,
    verbose = verbose
  )

  # Count parameters
  n_par <- count_parameters("UN", n_items, n_classes)

  # Create qlfit object
  new_qlfit(
    model_type = "UN",
    item_probs = fit$item_probs,
    class_probs = fit$class_probs,
    posteriors = fit$posteriors,
    loglik = fit$loglik,
    n_par = n_par,
    n_obs = n_obs,
    n_items = n_items,
    n_classes = n_classes,
    convergence = fit$converged,
    iterations = fit$iterations,
    call = NULL,
    theta = NULL,
    delta = NULL,
    item_order = NULL,
    constraints = NULL,
    se = NULL
  )
}
