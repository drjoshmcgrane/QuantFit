#' Fit Rasch Model (RM)
#'
#' @name fit_rm
#' @description Fits a standard Rasch model with continuous latent trait using
#'   the mirt package. This represents the full quantitative interpretation
#'   of the latent variable.
NULL

#' Fit Rasch Model
#'
#' @param data A matrix or data frame of binary responses (0/1) with subjects in rows
#'   and items in columns.
#' @param method Estimation method: "EM" (default), "MHRM", or "QMCEM"
#' @param quadpts Number of quadrature points for numerical integration (default 61)
#' @param verbose Print progress messages (default FALSE)
#' @param ... Additional arguments passed to \code{\link[mirt]{mirt}}
#'
#' @return A qlfit object containing:
#' \describe{
#'   \item{model_type}{"RM" for Rasch Model}
#'   \item{item_probs}{Matrix of item probabilities at quadrature points (for compatibility)}
#'   \item{class_probs}{Not used for RM (set to NULL)}
#'   \item{delta}{Vector of item difficulty parameters}
#'   \item{posteriors}{Not directly available (set to NULL)}
#'   \item{loglik}{Maximized log-likelihood}
#'   \item{n_par}{Number of estimated parameters (I item intercepts plus the
#'     latent variance, as counted by \code{mirt::extract.mirt(fit, "nest")})}
#'   \item{convergence}{Logical indicating successful convergence}
#'   \item{mirt_object}{The underlying mirt object for additional methods}
#' }
#'
#' @details
#' The Rasch model assumes:
#' \deqn{P(X_{ij} = 1 | \theta_i) = \frac{1}{1 + \exp(-(\theta_i - \delta_j))}}
#'
#' where:
#' \itemize{
#'   \item \eqn{\theta_i} is the continuous latent trait for person i
#'   \item \eqn{\delta_j} is the difficulty of item j
#' }
#'
#' Unlike the Latent Class Rasch model (LCR), the continuous Rasch model assumes
#' the latent trait is continuously distributed (typically normal). This
#' represents the strongest quantitative interpretation.
#'
#' The model is estimated using the \code{mirt} package with itemtype = "Rasch".
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
#' # Fit Rasch model
#' fit <- fit_rm(data)
#' print(fit)
#' summary(fit)
#'
#' # Get person ability estimates
#' theta_est <- rm_scores(fit)
#' plot(theta_true, theta_est)
#' }
#'
#' @seealso \code{\link{fit_lcr}} for Latent Class Rasch model
#'
#' @export
fit_rm <- function(data, method = c("EM", "MHRM", "QMCEM"),
                   quadpts = 61, verbose = FALSE, ...) {

  # Check for mirt
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' is required for fit_rm(). Please install it.")
  }

  # Capture call
  call <- match.call()

  # Validate inputs
  data <- validate_data(data)
  method <- match.arg(method)

  n_obs <- nrow(data)
  n_items <- ncol(data)

  # Fit Rasch model using mirt
  mirt_fit <- tryCatch({
    mirt::mirt(
      data = as.data.frame(data),
      model = 1,
      itemtype = "Rasch",
      method = method,
      quadpts = quadpts,
      verbose = verbose,
      ...
    )
  }, error = function(e) {
    stop("mirt estimation failed: ", e$message)
  })

  # Extract parameters
  coefs <- mirt::coef(mirt_fit, simplify = TRUE)
  delta <- -coefs$items[, "d"]  # mirt uses d = -delta convention
  names(delta) <- rownames(coefs$items)

  # Note: delta is NOT re-centered. mirt's scale (latent mean fixed at 0)
  # is already identified, and re-centering would make delta inconsistent
  # with rm_scores() thetas and item_probs.

  # Extract log-likelihood
  loglik <- mirt::extract.mirt(mirt_fit, "logLik")

  # Check convergence
  converged <- mirt::extract.mirt(mirt_fit, "converged")
  iterations <- mirt::extract.mirt(mirt_fit, "iterations")

  # Compute item probabilities at quadrature points (for compatibility with other models)
  # This gives a matrix similar to LCA item_probs but at continuous theta points
  theta_grid <- seq(-4, 4, length.out = 21)
  item_probs <- matrix(0, nrow = n_items, ncol = length(theta_grid))
  for (t in seq_along(theta_grid)) {
    item_probs[, t] <- inv_logit(theta_grid[t] - delta)
  }

  # Number of parameters: mirt's Rasch estimates I item intercepts plus the
  # latent variance (latent mean fixed at 0) = I + 1. Use mirt's own count.
  n_par <- tryCatch(
    as.integer(mirt::extract.mirt(mirt_fit, "nest")),
    error = function(e) n_items + 1L
  )

  # Create qlfit object
  result <- new_qlfit(
    model_type = "RM",
    item_probs = item_probs,
    class_probs = NULL,
    posteriors = NULL,
    loglik = as.numeric(loglik),
    n_par = n_par,
    n_obs = n_obs,
    n_items = n_items,
    n_classes = NA,
    convergence = converged,
    iterations = iterations,
    call = call,
    theta = theta_grid,  # Grid points, not class locations
    delta = delta,
    item_order = NULL,
    constraints = NULL,
    se = NULL
  )

  # Store mirt object for additional functionality
  attr(result, "mirt_object") <- mirt_fit

  result
}

#' Extract Rasch scores from RM model
#'
#' @param object A qlfit object from fit_rm
#' @param type Type of score: "EAP" (expected a posteriori), "MAP" (maximum a posteriori),
#'   "ML" (maximum likelihood), or "WLE" (weighted likelihood)
#'
#' @return Data frame with person scores and standard errors
#' @export
rm_scores <- function(object, type = c("EAP", "MAP", "ML", "WLE")) {
  if (object$model_type != "RM") {
    stop("rm_scores only applies to RM models")
  }

  type <- match.arg(type)

  # Get mirt object
  mirt_fit <- attr(object, "mirt_object")
  if (is.null(mirt_fit)) {
    stop("mirt object not found. Re-fit the model using fit_rm().")
  }

  # Extract scores using mirt
  scores <- mirt::fscores(mirt_fit, method = type, full.scores = TRUE)

  data.frame(
    theta = scores[, 1],
    se = if (ncol(scores) > 1) scores[, 2] else NA
  )
}

#' Item fit statistics for Rasch model
#'
#' @param object A qlfit object from fit_rm
#'
#' @return Data frame with item fit statistics
#' @export
rm_itemfit <- function(object) {
  if (object$model_type != "RM") {
    stop("rm_itemfit only applies to RM models")
  }

  mirt_fit <- attr(object, "mirt_object")
  if (is.null(mirt_fit)) {
    stop("mirt object not found. Re-fit the model using fit_rm().")
  }

  mirt::itemfit(mirt_fit)
}

#' Person fit statistics for Rasch model
#'
#' @param object A qlfit object from fit_rm
#'
#' @return Data frame with person fit statistics
#' @export
rm_personfit <- function(object) {
  if (object$model_type != "RM") {
    stop("rm_personfit only applies to RM models")
  }

  mirt_fit <- attr(object, "mirt_object")
  if (is.null(mirt_fit)) {
    stop("mirt object not found. Re-fit the model using fit_rm().")
  }

  mirt::personfit(mirt_fit)
}

#' Item information function for Rasch model
#'
#' @param object A qlfit object from fit_rm
#' @param theta Vector of theta values at which to compute information
#'
#' @return Matrix of item information values (items x theta)
#' @export
rm_item_info <- function(object, theta = seq(-4, 4, by = 0.1)) {
  if (object$model_type != "RM") {
    stop("rm_item_info only applies to RM models")
  }

  # Item information for Rasch model: I(theta) = P(theta) * Q(theta)
  delta <- object$delta
  n_items <- length(delta)

  info <- matrix(0, nrow = n_items, ncol = length(theta))
  for (i in 1:n_items) {
    p <- inv_logit(theta - delta[i])
    info[i, ] <- p * (1 - p)
  }

  rownames(info) <- names(delta)
  colnames(info) <- round(theta, 2)

  info
}

#' Test information function for Rasch model
#'
#' @param object A qlfit object from fit_rm
#' @param theta Vector of theta values
#'
#' @return Vector of test information values
#' @export
rm_test_info <- function(object, theta = seq(-4, 4, by = 0.1)) {
  item_info <- rm_item_info(object, theta)
  colSums(item_info)
}

#' Compare Rasch model with LCR approximation
#'
#' Fits both RM and LCR models and compares how well LCR approximates
#' the continuous Rasch model.
#'
#' @param data Data matrix
#' @param n_classes_lcr Number of classes for LCR model (default 21)
#' @param verbose Print progress
#'
#' @return List with both fits and comparison statistics
#' @export
compare_rm_lcr <- function(data, n_classes_lcr = 21, verbose = FALSE) {
  # Fit both models
  rm_fit <- fit_rm(data, verbose = verbose)
  lcr_fit <- fit_lcr(data, n_classes = n_classes_lcr, verbose = verbose)

  # Compare item parameters
  delta_rm <- rm_fit$delta
  delta_lcr <- lcr_fit$delta

  # Scale LCR to match RM (they should be on similar scale)
  delta_diff <- delta_rm - delta_lcr

  # Compare fit statistics
  comparison <- data.frame(
    model = c("RM", "LCR"),
    loglik = c(rm_fit$loglik, lcr_fit$loglik),
    n_par = c(rm_fit$n_par, lcr_fit$n_par),
    AIC = c(AIC(rm_fit), AIC(lcr_fit)),
    BIC = c(BIC(rm_fit), BIC(lcr_fit))
  )

  list(
    rm = rm_fit,
    lcr = lcr_fit,
    comparison = comparison,
    delta_comparison = data.frame(
      item = 1:length(delta_rm),
      delta_rm = delta_rm,
      delta_lcr = delta_lcr,
      difference = delta_diff
    )
  )
}
