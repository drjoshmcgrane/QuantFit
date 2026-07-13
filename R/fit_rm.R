#' Fit Rasch Model (RM)
#'
#' @name fit_rm
#' @description Fits a Rasch / partial-credit model with a continuous latent
#'   trait by marginal maximum likelihood over Gauss-Hermite quadrature, using
#'   the package's own EM engine (no external IRT package). This represents the
#'   full quantitative interpretation of the latent variable. Dichotomous data
#'   give the Rasch model; polytomous data give the partial credit model.
NULL

#' Fit Rasch Model
#'
#' @param data A matrix or data frame of item responses with subjects in rows
#'   and items in columns. Dichotomous (0/1) or polytomous (0..m) scoring.
#' @param quadpts Number of Gauss-Hermite quadrature points (default 61)
#' @param use_cpp Use the compiled C++ engine for the E-step (default TRUE)
#' @param seed Optional random seed for the starting values
#' @param verbose Print progress messages (default FALSE)
#' @param ... Ignored (accepted for backward compatibility)
#'
#' @return A qlfit object containing:
#' \describe{
#'   \item{model_type}{"RM" for Rasch Model}
#'   \item{item_probs}{Expected item-score curves on a theta grid (for
#'     dichotomous items, the usual response probabilities)}
#'   \item{class_probs}{Not used for RM (set to NULL)}
#'   \item{delta}{Vector of item step (difficulty) parameters}
#'   \item{loglik}{Maximized marginal log-likelihood}
#'   \item{n_par}{Number of estimated parameters (total category steps plus the
#'     latent variance)}
#'   \item{convergence}{Logical indicating successful convergence}
#' }
#' The fitted quadrature model is stored in \code{attr(fit, "rm_fit")}.
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
#' the latent trait is continuously distributed (normal). This represents the
#' strongest quantitative interpretation.
#'
#' The model is estimated by marginal maximum likelihood: person locations
#' \eqn{\theta \sim N(0, \sigma^2)} are integrated out by Gauss-Hermite
#' quadrature and the quadrature nodes play the role of latent classes in the
#' same EM engine used for the latent-class models, so no external IRT package
#' is required. Validated against \code{mirt} to a log-likelihood difference
#' below 0.02.
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
fit_rm <- function(data, quadpts = 61L, use_cpp = TRUE, seed = NULL,
                   verbose = FALSE, ...) {

  # Capture call
  call <- match.call()

  # Estimate the Rasch / partial-credit model with the package's own marginal
  # maximum-likelihood engine (Gauss-Hermite quadrature EM), for both
  # dichotomous and polytomous data. No external IRT package is required.
  if (is.data.frame(data)) data <- as.matrix(data)
  poly <- .is_polytomous(data)
  # Missing responses are allowed: the quadrature E-step uses the masked
  # likelihood (each person contributes only their observed cells), valid
  # under MAR.
  data <- if (anyNA(data)) .validate_poly(data, allow_na = TRUE)
          else validate_data_any(data)
  n_obs <- nrow(data); n_items <- ncol(data)

  mml <- em_rasch_mml(data, n_quad = quadpts, use_cpp = use_cpp, seed = seed)
  cat_counts <- mml$cat_counts
  item_names <- colnames(data)
  if (is.null(item_names)) item_names <- paste0("Item", seq_len(n_items))
  idx <- split(seq_along(mml$delta), rep(seq_len(n_items), cat_counts))
  delta_list <- lapply(idx, function(ii) mml$delta[ii])

  # expected item-score curves E(X_j | theta) on a grid (item_probs slot); for
  # dichotomous items this is the usual P(X = 1 | theta)
  theta_grid <- seq(-4, 4, length.out = 21L)
  item_probs <- t(vapply(seq_len(n_items), function(j) {
    P <- cpp_pcm_probs(theta_grid, delta_list[[j]])
    as.numeric(P %*% (0:cat_counts[j]))
  }, numeric(length(theta_grid))))
  rownames(item_probs) <- item_names

  result <- new_qlfit(
    model_type = "RM", item_probs = item_probs, class_probs = NULL,
    posteriors = NULL, loglik = mml$loglik, n_par = mml$n_par, n_obs = n_obs,
    n_items = n_items, n_classes = NA, convergence = mml$converged,
    iterations = mml$iterations, call = call, theta = theta_grid,
    delta = stats::setNames(mml$delta, NULL), item_order = NULL,
    constraints = NULL, se = NULL)
  if (poly) result$polytomous <- TRUE
  result$cat_counts <- cat_counts
  attr(result, "rm_fit") <- mml
  return(result)
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
  rf <- attr(object, "rm_fit")
  if (is.null(rf)) stop("RM fit not found. Re-fit the model using fit_rm().")

  nodes <- rf$nodes; post <- rf$posteriors      # posteriors: n x Q over nodes
  if (type == "EAP") {
    eap <- as.numeric(post %*% nodes)
    v <- as.numeric(post %*% (nodes^2)) - eap^2
    return(data.frame(theta = eap, se = sqrt(pmax(v, 0))))
  }
  if (type == "MAP") {
    return(data.frame(theta = nodes[max.col(post, ties.method = "first")],
                      se = NA_real_))
  }

  # ML / WLE via the sufficient total score (valid for the Rasch/PCM family):
  # solve the test-score equation for each distinct observed total score.
  # With missing responses the total is not comparable across persons (it
  # depends on which items were answered), so the score-to-theta mapping is
  # not defined; use EAP/MAP, which come from the masked likelihood.
  if (anyNA(rf$data)) {
    stop("ML/WLE scores are score-sufficiency based and require complete ",
         "data; use type = \"EAP\" (or \"MAP\") with missing responses.")
  }
  dl <- rf$delta_list; sumM <- sum(rf$cat_counts)
  TS <- function(th) sum(vapply(dl, function(dj) .pcm_moments(th, dj)$E, numeric(1)))
  Inf_fn <- function(th) sum(vapply(dl, function(dj) .pcm_moments(th, dj)$V, numeric(1)))
  warm <- function(th) sum(vapply(dl, function(dj) .pcm_moments(th, dj)$M3, numeric(1)))
  solve_theta <- function(target) {
    if (target <= 0) return(if (type == "ML") -Inf else NA_real_)
    if (target >= sumM) return(if (type == "ML") Inf else NA_real_)
    f <- if (type == "ML") function(th) TS(th) - target
         else function(th) TS(th) + warm(th) / (2 * Inf_fn(th)) - target
    tryCatch(stats::uniroot(f, c(-10, 10))$root, error = function(e) NA_real_)
  }
  us <- sort(unique(rf$scores))
  tmap <- vapply(us, solve_theta, numeric(1)); names(tmap) <- as.character(us)
  th <- tmap[as.character(rf$scores)]
  se <- vapply(th, function(x) if (is.finite(x)) 1 / sqrt(Inf_fn(x)) else NA_real_,
               numeric(1))
  data.frame(theta = as.numeric(th), se = as.numeric(se))
}

#' Item fit statistics for Rasch / partial-credit model
#'
#' @param object A qlfit object from fit_rm
#'
#' @return Data frame with outfit and infit mean-square statistics per item.
#' @export
rm_itemfit <- function(object) {
  if (object$model_type != "RM") stop("rm_itemfit only applies to RM models")
  rf <- attr(object, "rm_fit")
  if (is.null(rf)) stop("RM fit not found. Re-fit the model using fit_rm().")
  data <- rf$data; dl <- rf$delta_list
  th <- rm_scores(object, "EAP")$theta
  J <- ncol(data)
  outfit <- infit <- numeric(J)
  for (j in seq_len(J)) {
    m <- .pcm_moments(th, dl[[j]])
    resid2 <- (data[, j] - m$E)^2
    outfit[j] <- mean(resid2 / m$V)              # unweighted mean-square
    infit[j] <- sum(resid2) / sum(m$V)           # information-weighted
  }
  data.frame(item = rownames(object$item_probs),
             outfit_MSQ = outfit, infit_MSQ = infit,
             stringsAsFactors = FALSE)
}

#' Person fit statistics for Rasch / partial-credit model
#'
#' @param object A qlfit object from fit_rm
#'
#' @return Data frame with outfit and infit mean-square statistics per person.
#' @export
rm_personfit <- function(object) {
  if (object$model_type != "RM") stop("rm_personfit only applies to RM models")
  rf <- attr(object, "rm_fit")
  if (is.null(rf)) stop("RM fit not found. Re-fit the model using fit_rm().")
  data <- rf$data; dl <- rf$delta_list
  th <- rm_scores(object, "EAP")$theta
  n <- nrow(data); J <- ncol(data)
  E <- V <- matrix(0, n, J)
  for (j in seq_len(J)) { m <- .pcm_moments(th, dl[[j]]); E[, j] <- m$E; V[, j] <- m$V }
  resid2 <- (data - E)^2
  data.frame(outfit_MSQ = rowMeans(resid2 / V),
             infit_MSQ = rowSums(resid2) / rowSums(V))
}

#' Item information function for Rasch / partial-credit model
#'
#' @param object A qlfit object from fit_rm
#' @param theta Vector of theta values at which to compute information
#'
#' @return Matrix of item information values (items x theta)
#' @export
rm_item_info <- function(object, theta = seq(-4, 4, by = 0.1)) {
  if (object$model_type != "RM") stop("rm_item_info only applies to RM models")
  rf <- attr(object, "rm_fit")
  if (is.null(rf)) stop("RM fit not found. Re-fit the model using fit_rm().")
  # partial-credit item information is the response variance at theta
  info <- t(vapply(rf$delta_list, function(dj) .pcm_moments(theta, dj)$V,
                   numeric(length(theta))))
  rownames(info) <- rownames(object$item_probs)
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
