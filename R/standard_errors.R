#' Standard Errors for QuantFit Models
#'
#' @name standard_errors
#' @description Observed-information (numerical Hessian) and nonparametric
#'   bootstrap standard errors for fitted QuantFit models.
NULL

# ============================================================================
# Internal helpers
# ============================================================================

#' Softmax with the last element as reference category
#'
#' @param alpha Vector of C-1 multinomial-logit parameters
#' @return Probability vector of length C (softmax of c(alpha, 0))
#' @keywords internal
#' @noRd
softmax_ref <- function(alpha) {
  z <- c(alpha, 0)
  z <- z - max(z)
  ez <- exp(z)
  ez / sum(ez)
}

#' Enumerate all permutations of 1:k
#'
#' @param k Integer (k! grows fast; used only for k <= 6)
#' @return Matrix with k! rows, each row a permutation of 1:k
#' @keywords internal
#' @noRd
all_permutations <- function(k) {
  k <- as.integer(k)
  if (k == 1L) return(matrix(1L, 1L, 1L))
  sub <- all_permutations(k - 1L)
  out <- vector("list", k)
  for (i in seq_len(k)) {
    others <- setdiff(seq_len(k), i)
    out[[i]] <- cbind(i, matrix(others[sub], nrow = nrow(sub), ncol = k - 1L))
  }
  res <- do.call(rbind, out)
  dimnames(res) <- NULL
  res
}

#' Moore-Penrose pseudo-inverse via SVD
#'
#' @param m Square matrix
#' @param tol Relative tolerance for singular values
#' @return Pseudo-inverse of m
#' @keywords internal
#' @noRd
pseudo_inverse <- function(m, tol = 1e-10) {
  s <- svd(m)
  pos <- s$d > tol * max(s$d)
  if (!any(pos)) return(matrix(NA_real_, nrow(m), ncol(m)))
  s$v[, pos, drop = FALSE] %*%
    diag(1 / s$d[pos], sum(pos)) %*%
    t(s$u[, pos, drop = FALSE])
}

#' Invert an observed information matrix, with pseudo-inverse fallback
#'
#' @param info Observed information matrix (-Hessian of the log-likelihood)
#' @return Variance-covariance matrix
#' @keywords internal
#' @noRd
invert_information <- function(info) {
  V <- tryCatch(solve(info), error = function(e) NULL)
  if (is.null(V) || anyNA(V)) {
    warning("Observed information matrix is singular or ill-conditioned; ",
            "using an SVD-based pseudo-inverse. Standard errors may be ",
            "unreliable.", call. = FALSE)
    V <- pseudo_inverse(info)
  }
  V
}

#' Convert variances to standard errors (NA for non-positive variances)
#'
#' @param v Vector of variances
#' @return Vector of standard errors
#' @keywords internal
#' @noRd
se_from_var <- function(v) {
  ifelse(is.finite(v) & v > 0, sqrt(pmax(v, 0)), NA_real_)
}

#' Variance-covariance of class probabilities via the softmax Jacobian
#'
#' @param V_alpha Vcov of the C-1 multinomial-logit parameters
#' @param cp Fitted class probabilities (length C)
#' @return C x C variance-covariance matrix of the class probabilities
#' @keywords internal
#' @noRd
class_prob_vcov <- function(V_alpha, cp) {
  C <- length(cp)
  Jm <- matrix(0, C, C - 1L)
  for (cc in seq_len(C)) {
    for (k in seq_len(C - 1L)) {
      Jm[cc, k] <- cp[cc] * ((cc == k) - cp[k])
    }
  }
  Jm %*% V_alpha %*% t(Jm)
}

#' Observed-data log-likelihood of an LCA model as a function of eta
#'
#' eta = (multinomial logits of class probs (C-1 values), logits of item
#' probs (J*C values, column-major)).
#'
#' @keywords internal
#' @noRd
lca_loglik_eta <- function(eta, data, n_items, n_classes, use_cpp) {
  n_alpha <- n_classes - 1L
  cp <- softmax_ref(eta[seq_len(n_alpha)])
  ip <- matrix(stats::plogis(eta[-seq_len(n_alpha)]), n_items, n_classes)
  if (use_cpp) {
    cpp_e_step(data, ip, cp)$loglik
  } else {
    e_step(data, ip, cp)$loglik
  }
}

#' Observed-data log-likelihood of the LCR model as a function of eta
#'
#' eta = (theta (C values), delta_2..delta_J (J-1 values, with
#' delta_1 = -sum implied), multinomial logits of class probs (C-1 values)).
#'
#' @keywords internal
#' @noRd
lcr_loglik_eta <- function(eta, data, n_items, n_classes, use_cpp) {
  theta <- eta[seq_len(n_classes)]
  dfree <- eta[n_classes + seq_len(n_items - 1L)]
  delta <- c(-sum(dfree), dfree)
  cp <- softmax_ref(eta[n_classes + n_items - 1L + seq_len(n_classes - 1L)])
  ip <- bound_probs(compute_rasch_probs(theta, delta))
  if (use_cpp) {
    cpp_e_step(data, ip, cp)$loglik
  } else {
    e_step(data, ip, cp)$loglik
  }
}

#' Detect active inequality constraints in a constrained LCA fit
#'
#' A constraint is active when two adjacent fitted probabilities (adjacent
#' classes within an item for MON, adjacent items along the item order within
#' a class for IIO, both for DM) are equal within tol.
#'
#' @param object A qlfit object
#' @param tol Equality tolerance (default 1e-6)
#' @return List with n_active (count) and cells (J x C logical matrix flagging
#'   item-probability cells involved in an active constraint)
#' @keywords internal
#' @noRd
detect_active_constraints <- function(object, tol = 1e-6) {
  J <- object$n_items
  C <- object$n_classes
  ip <- object$item_probs
  cells <- matrix(FALSE, J, C)
  n_active <- 0L

  if (object$model_type %in% c("MON", "DM")) {
    for (j in seq_len(J)) {
      act <- which(abs(diff(ip[j, ])) < tol)
      if (length(act)) {
        cells[j, act] <- TRUE
        cells[j, act + 1L] <- TRUE
        n_active <- n_active + length(act)
      }
    }
  }

  if (object$model_type %in% c("IIO", "DM")) {
    ord <- object$item_order
    if (is.null(ord)) ord <- seq_len(J)
    for (cc in seq_len(C)) {
      act <- which(abs(diff(ip[ord, cc])) < tol)
      if (length(act)) {
        cells[ord[act], cc] <- TRUE
        cells[ord[act + 1L], cc] <- TRUE
        n_active <- n_active + length(act)
      }
    }
  }

  list(n_active = n_active, cells = cells)
}

#' Recover the data used to produce a qlfit object
#'
#' @param object A qlfit object
#' @param data Data supplied by the user (may be NULL)
#' @param caller_env Environment in which to evaluate the stored call
#' @return Validated binary data matrix
#' @keywords internal
#' @noRd
resolve_se_data <- function(object, data, caller_env) {
  if (!is.null(data)) return(validate_data(data))

  d <- attr(object, "data")
  if (!is.null(d)) return(validate_data(d))

  if (object$model_type == "RM") {
    rf <- attr(object, "rm_fit")
    if (!is.null(rf) && !is.null(rf$data)) return(validate_data_any(rf$data))
  }

  if (!is.null(object$call) && !is.null(object$call$data)) {
    d <- tryCatch(eval(object$call$data, envir = caller_env),
                  error = function(e) NULL)
    if (!is.null(d)) {
      d <- tryCatch(validate_data(d), error = function(e) NULL)
      if (!is.null(d)) return(d)
    }
  }

  stop("Could not recover the original data from the fit. ",
       "Pass it explicitly via compute_se(object, data = ...).",
       call. = FALSE)
}

# ============================================================================
# Hessian-based standard errors
# ============================================================================

#' Hessian SEs for LCA-type models (UN, MON, IIO, DM)
#'
#' @keywords internal
#' @noRd
se_hessian_lca <- function(object, data, use_cpp, verbose) {
  J <- object$n_items
  C <- object$n_classes
  ip <- object$item_probs
  cp <- pmax(object$class_probs, 1e-12)

  eps_boundary <- 1e-6
  boundary <- ip < eps_boundary | ip > 1 - eps_boundary
  active <- detect_active_constraints(object)

  eta <- c(log(cp[-C] / cp[C]), stats::qlogis(bound_probs(ip)))
  if (verbose) {
    message("Computing numerical Hessian of the observed-data log-likelihood (",
            length(eta), " parameters)...")
  }

  H <- numDeriv::hessian(function(e) lca_loglik_eta(e, data, J, C, use_cpp),
                         eta)
  V <- invert_information(-H)

  n_alpha <- C - 1L
  idx_alpha <- seq_len(n_alpha)
  V_alpha <- V[idx_alpha, idx_alpha, drop = FALSE]
  se_class <- se_from_var(diag(class_prob_vcov(V_alpha, object$class_probs)))

  # Delta method on the logit scale: SE(p) = SE(eta) * p * (1 - p)
  se_eta_items <- se_from_var(diag(V)[-idx_alpha])
  se_items <- matrix(se_eta_items, J, C) * ip * (1 - ip)

  if (any(boundary)) {
    se_items[boundary] <- NA_real_
    warning(sum(boundary), " item probabilit",
            if (sum(boundary) == 1L) "y is" else "ies are",
            " at the 0/1 boundary; the logit transform is degenerate there ",
            "and the corresponding standard errors are set to NA.",
            call. = FALSE)
  }

  if (active$n_active > 0L) {
    se_items[active$cells] <- NA_real_
    warning(active$n_active, " constraint",
            if (active$n_active == 1L) " is" else "s are",
            " active at the solution; hessian SEs are unreliable for the ",
            "affected parameters (set to NA), consider method = 'bootstrap'.",
            call. = FALSE)
  }

  list(
    se = list(item_probs = se_items, class_probs = se_class),
    vcov_eta = V,
    eta_index = list(
      alpha = idx_alpha,
      item_probs = matrix(n_alpha + seq_len(J * C), J, C)
    ),
    active = active,
    boundary = boundary
  )
}

#' Hessian SEs for the Latent Class Rasch model
#'
#' @keywords internal
#' @noRd
se_hessian_lcr <- function(object, data, use_cpp, verbose) {
  J <- object$n_items
  C <- object$n_classes
  cp <- pmax(object$class_probs, 1e-12)

  eta <- c(object$theta, object$delta[-1], log(cp[-C] / cp[C]))
  if (verbose) {
    message("Computing numerical Hessian of the observed-data log-likelihood (",
            length(eta), " parameters)...")
  }

  H <- numDeriv::hessian(function(e) lcr_loglik_eta(e, data, J, C, use_cpp),
                         eta)
  V <- invert_information(-H)

  idx_theta <- seq_len(C)
  idx_dfree <- C + seq_len(J - 1L)
  idx_alpha <- C + J - 1L + seq_len(C - 1L)

  se_theta <- se_from_var(diag(V)[idx_theta])

  # delta_1 = -sum(delta_2..delta_J), so Var(delta_1) = 1' V_dfree 1
  V_dfree <- V[idx_dfree, idx_dfree, drop = FALSE]
  se_delta <- c(se_from_var(sum(V_dfree)), se_from_var(diag(V_dfree)))
  names(se_delta) <- names(object$delta)

  V_alpha <- V[idx_alpha, idx_alpha, drop = FALSE]
  se_class <- se_from_var(diag(class_prob_vcov(V_alpha, object$class_probs)))

  list(
    se = list(theta = se_theta, delta = se_delta, class_probs = se_class),
    vcov_eta = V,
    eta_index = list(theta = idx_theta, delta_free = idx_dfree,
                     alpha = idx_alpha)
  )
}

#' Hessian SEs for the Rasch / partial-credit model
#'
#' Observed-information standard errors from the numeric Hessian of the marginal
#' (Gauss-Hermite) log-likelihood evaluated at the fitted step parameters and
#' log latent SD. No external IRT package is used.
#'
#' @keywords internal
#' @noRd
se_hessian_rm <- function(object, verbose) {
  rf <- attr(object, "rm_fit")
  if (is.null(rf)) {
    stop("RM fit not found on the object. Re-fit the model using fit_rm().",
         call. = FALSE)
  }
  data <- rf$data; cat_counts <- rf$cat_counts; sumM <- sum(cat_counts)
  w <- rf$weights
  z0 <- rf$nodes / rf$sigma                       # base N(0,1) quadrature nodes
  idx <- split(seq_len(sumM), rep(seq_len(ncol(data)), cat_counts))

  negll <- function(par) {
    delta <- par[seq_len(sumM)]; sigma <- exp(par[sumM + 1L])
    dl <- lapply(idx, function(ii) delta[ii])
    ip <- compute_pcm_probs(sigma * z0, dl, use_cpp = TRUE)
    -poly_estep(data, ip, w, use_cpp = TRUE)$loglik
  }
  par_hat <- c(rf$delta, log(rf$sigma))
  H <- numDeriv::hessian(negll, par_hat)
  V <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(V)) {
    warning("Rasch information matrix is singular; standard errors are NA.",
            call. = FALSE)
    se_all <- rep(NA_real_, length(par_hat))
  } else {
    se_all <- sqrt(pmax(diag(V), 0))
  }
  se_delta <- se_all[seq_len(sumM)]
  step <- unlist(lapply(seq_along(cat_counts), function(j) seq_len(cat_counts[j])))
  item <- rep(seq_along(cat_counts), cat_counts)
  names(se_delta) <- if (sumM == length(cat_counts)) rownames(object$item_probs)
                     else sprintf("%s.s%d", rownames(object$item_probs)[item], step)
  list(se = list(delta = se_delta, log_sigma = se_all[sumM + 1L]))
}

# ============================================================================
# Bootstrap standard errors
# ============================================================================

#' Refit the same model on (resampled) data
#'
#' @keywords internal
#' @noRd
refit_model <- function(object, bd, n_starts, seed_b, use_cpp) {
  switch(object$model_type,
    UN = fit_un(bd, n_classes = object$n_classes, n_starts = n_starts,
                use_cpp = use_cpp, seed = seed_b),
    MON = fit_mon(bd, n_classes = object$n_classes, n_starts = n_starts,
                  use_cpp = use_cpp, seed = seed_b),
    IIO = fit_iio(bd, n_classes = object$n_classes,
                  item_order = object$item_order, n_starts = n_starts,
                  use_cpp = use_cpp, seed = seed_b),
    DM = fit_dm(bd, n_classes = object$n_classes,
                item_order = object$item_order, n_starts = n_starts,
                use_cpp = use_cpp, seed = seed_b),
    LCR = fit_lcr(bd, n_classes = object$n_classes, n_starts = n_starts,
                  use_cpp = use_cpp, seed = seed_b),
    RM = fit_rm(bd, verbose = FALSE),
    stop("Unknown model type: ", object$model_type)
  )
}

#' Find the class permutation minimizing total item-probability distance
#'
#' @keywords internal
#' @noRd
best_permutation <- function(ip, ref, perms) {
  best <- 1L
  best_val <- Inf
  for (r in seq_len(nrow(perms))) {
    v <- sum(abs(ip[, perms[r, ], drop = FALSE] - ref))
    if (v < best_val) {
      best_val <- v
      best <- r
    }
  }
  perms[best, ]
}

#' Nonparametric bootstrap standard errors
#'
#' @keywords internal
#' @noRd
se_bootstrap <- function(object, data, B, n_starts, seed, use_cpp, verbose) {
  model <- object$model_type
  n <- nrow(data)
  is_lca <- model %in% c("UN", "MON", "IIO", "DM")

  if (is_lca && object$n_classes > 6L) {
    stop("Bootstrap label alignment enumerates class permutations and ",
         "currently supports at most 6 classes (got ", object$n_classes,
         "). Sorry - please reduce the number of classes or use ",
         "method = 'hessian'.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)
  # Draw all resampling indices up front so results are deterministic given
  # the seed regardless of how many refits fail
  idx_mat <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n, ncol = B)

  perms <- if (is_lca) all_permutations(object$n_classes) else NULL

  ip_list <- list()
  cp_list <- list()
  th_list <- list()
  de_list <- list()
  n_fail <- 0L

  for (b in seq_len(B)) {
    bd <- data[idx_mat[, b], , drop = FALSE]
    seed_b <- if (!is.null(seed)) seed + b else NULL

    fit_b <- tryCatch(
      suppressWarnings(suppressMessages(
        refit_model(object, bd, n_starts, seed_b, use_cpp)
      )),
      error = function(e) NULL
    )

    if (is.null(fit_b)) {
      n_fail <- n_fail + 1L
      next
    }

    if (is_lca) {
      # Resolve label switching: align classes to the original fit
      pm <- best_permutation(fit_b$item_probs, object$item_probs, perms)
      ip_list[[length(ip_list) + 1L]] <- fit_b$item_probs[, pm, drop = FALSE]
      cp_list[[length(cp_list) + 1L]] <- fit_b$class_probs[pm]
    } else if (model == "LCR") {
      # fit_lcr already orders classes by theta; re-order defensively
      ord <- order(fit_b$theta)
      th_list[[length(th_list) + 1L]] <- fit_b$theta[ord]
      cp_list[[length(cp_list) + 1L]] <- fit_b$class_probs[ord]
      de_list[[length(de_list) + 1L]] <- as.numeric(fit_b$delta)
    } else {  # RM
      de_list[[length(de_list) + 1L]] <- as.numeric(fit_b$delta)
    }

    if (verbose && b %% 10L == 0L) {
      message("Bootstrap replicate ", b, " of ", B)
    }
  }

  B_eff <- B - n_fail
  if (B_eff < 2L) {
    stop("Fewer than 2 bootstrap refits succeeded (", B_eff, " of ", B,
         "); cannot compute bootstrap standard errors.", call. = FALSE)
  }
  if (n_fail > 0.1 * B) {
    warning(n_fail, " of ", B, " bootstrap refits failed (more than 10%); ",
            "bootstrap standard errors may be unreliable.", call. = FALSE)
  }

  se <- list()
  if (length(ip_list)) {
    arr <- array(unlist(ip_list),
                 dim = c(object$n_items, object$n_classes, length(ip_list)))
    se$item_probs <- apply(arr, c(1, 2), stats::sd)
  }
  if (length(cp_list)) {
    se$class_probs <- apply(do.call(rbind, cp_list), 2, stats::sd)
  }
  if (length(th_list)) {
    se$theta <- apply(do.call(rbind, th_list), 2, stats::sd)
  }
  if (length(de_list)) {
    se$delta <- apply(do.call(rbind, de_list), 2, stats::sd)
    if (!is.null(names(object$delta))) names(se$delta) <- names(object$delta)
  }

  list(se = se, B_effective = B_eff, n_failed = n_fail)
}

# ============================================================================
# Exported interface
# ============================================================================

#' Compute Standard Errors for a Fitted QuantFit Model
#'
#' Computes standard errors for the parameters of a fitted QuantFit model,
#' either from the observed information matrix (numerical Hessian of the
#' observed-data log-likelihood) or by a nonparametric bootstrap.
#'
#' @param object A qlfit object from \code{\link{fit_un}},
#'   \code{\link{fit_mon}}, \code{\link{fit_iio}}, \code{\link{fit_dm}},
#'   \code{\link{fit_lcr}}, or \code{\link{fit_rm}}.
#' @param method Either \code{"hessian"} (default; observed-information SEs)
#'   or \code{"bootstrap"} (nonparametric bootstrap over rows of the data).
#' @param B Number of bootstrap replicates (bootstrap only, default 200).
#' @param n_starts Number of random starts for each bootstrap refit
#'   (bootstrap only, default 3).
#' @param seed Random seed; results are deterministic given the seed.
#' @param use_cpp Use the compiled C++ likelihood/EM machinery (default TRUE).
#' @param verbose Print progress messages (default FALSE).
#' @param data The original data matrix used to obtain \code{object}. Fitted
#'   objects do not store the data, so it is recovered from (in order): this
#'   argument, a \code{"data"} attribute of the fit, the stored mirt object
#'   (RM only), or by re-evaluating the data argument of the stored call. If
#'   none of these succeed, an error asks for the data explicitly.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"qlse"}: a list with elements
#' \describe{
#'   \item{model}{Model type ("UN", "MON", "IIO", "DM", "LCR", "RM")}
#'   \item{method}{"hessian" or "bootstrap"}
#'   \item{se}{List of standard errors matching the parameter structure:
#'     \code{item_probs} (items x classes matrix) and \code{class_probs} for
#'     the LCA models; \code{theta}, \code{delta}, and \code{class_probs} for
#'     LCR; \code{delta} for RM}
#'   \item{estimates}{The corresponding point estimates from the fit}
#'   \item{vcov_eta, eta_index}{(hessian only) variance-covariance matrix of
#'     the unconstrained working parameters and the index map into it}
#'   \item{active_constraints}{(MON/IIO/DM, hessian only) count and location
#'     of active inequality constraints}
#'   \item{boundary_cells}{(LCA models, hessian only) logical matrix flagging
#'     item probabilities on the 0/1 boundary}
#'   \item{B, B_effective, n_failed}{(bootstrap only) requested replicates,
#'     successful replicates, and failures}
#' }
#'
#' @details
#' \strong{Hessian method.} The observed-data log-likelihood is
#' reparameterized in terms of unconstrained working parameters
#' (multinomial-logit class probabilities and logit item probabilities for
#' the LCA models; \code{theta}, free \code{delta}, and multinomial-logit
#' class probabilities for LCR), differentiated numerically with
#' \pkg{numDeriv}, and inverted to obtain the variance-covariance matrix.
#' Standard errors are mapped back to the probability scale by the delta
#' method (full softmax Jacobian for the class probabilities). For RM the
#' SEs are taken from \pkg{mirt} (re-fitting with \code{SE = TRUE} when the
#' stored fit does not carry them); \code{delta = -d}, so the intercept SEs
#' apply unchanged.
#'
#' For MON/IIO/DM the inequality constraints do not change the interior
#' information, but when constraints are \emph{active} (adjacent fitted
#' probabilities equal within 1e-6) the MLE lies on the boundary of the
#' parameter space and Hessian-based SEs are invalid there: affected cells
#' are set to NA with a warning. Item probabilities on the 0/1 boundary also
#' get NA SEs (the logit transform is degenerate there).
#'
#' \strong{Bootstrap method.} Rows of the data are resampled with
#' replacement B times and the same model is refit on each replicate. Label
#' switching is resolved by aligning each replicate's classes to the
#' original fit via the permutation minimizing the total absolute difference
#' in item probabilities (LCR classes are aligned by their theta order). The
#' SE is the standard deviation of each parameter across successful
#' replicates.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' data <- matrix(rbinom(500 * 6, 1, 0.5), nrow = 500)
#' fit <- fit_un(data, n_classes = 2, seed = 1)
#'
#' # Observed-information standard errors
#' se_h <- compute_se(fit, data = data)
#' print(se_h)
#'
#' # Bootstrap standard errors
#' se_b <- compute_se(fit, method = "bootstrap", B = 200, seed = 2,
#'                    data = data)
#' print(se_b)
#' }
#'
#' @export
compute_se <- function(object,
                       method = c("hessian", "bootstrap"),
                       B = 200,
                       n_starts = 3,
                       seed = NULL,
                       use_cpp = TRUE,
                       verbose = FALSE,
                       data = NULL,
                       ...) {

  if (!inherits(object, "qlfit")) {
    stop("object must be a qlfit object (from fit_un, fit_mon, fit_iio, ",
         "fit_dm, fit_lcr, or fit_rm)")
  }
  method <- match.arg(method)
  model <- object$model_type
  caller_env <- parent.frame()

  # The RM hessian path gets everything from the stored mirt object;
  # every other path needs the original data
  if (!(method == "hessian" && model == "RM")) {
    data <- resolve_se_data(object, data, caller_env)
  }

  res <- if (method == "hessian") {
    switch(model,
      UN = ,
      MON = ,
      IIO = ,
      DM = se_hessian_lca(object, data, use_cpp, verbose),
      LCR = se_hessian_lcr(object, data, use_cpp, verbose),
      RM = se_hessian_rm(object, verbose),
      stop("Unknown model type: ", model)
    )
  } else {
    se_bootstrap(object, data, B, n_starts, seed, use_cpp, verbose)
  }

  estimates <- switch(model,
    UN = ,
    MON = ,
    IIO = ,
    DM = list(item_probs = object$item_probs,
              class_probs = object$class_probs),
    LCR = list(theta = object$theta,
               delta = object$delta,
               class_probs = object$class_probs),
    RM = list(delta = object$delta)
  )

  structure(
    list(
      model = model,
      method = method,
      se = res$se,
      estimates = estimates,
      n_obs = object$n_obs,
      n_items = object$n_items,
      n_classes = object$n_classes,
      vcov_eta = res$vcov_eta,
      eta_index = res$eta_index,
      active_constraints = res$active,
      boundary_cells = res$boundary,
      B = if (method == "bootstrap") B else NULL,
      B_effective = res$B_effective,
      n_failed = res$n_failed
    ),
    class = "qlse"
  )
}

#' Print method for qlse objects
#'
#' @param x A qlse object from \code{\link{compute_se}}
#' @param digits Number of digits to display (default 4)
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns x
#' @export
print.qlse <- function(x, digits = 4, ...) {
  fmt <- function(est, se) {
    paste0(formatC(est, digits = digits, format = "f"),
           " (", ifelse(is.na(se), "NA",
                        formatC(se, digits = digits, format = "f")), ")")
  }

  cat("\nStandard errors for", x$model, "model\n")
  if (x$method == "hessian") {
    cat("Method: observed information (numerical Hessian)\n")
  } else {
    cat("Method: nonparametric bootstrap (", x$B_effective, " of ", x$B,
        " replicates successful)\n", sep = "")
  }
  cat(rep("-", 55), "\n", sep = "")

  if (!is.null(x$se$class_probs)) {
    cat("\nClass probabilities [estimate (SE)]:\n")
    v <- fmt(x$estimates$class_probs, x$se$class_probs)
    names(v) <- paste0("Class ", seq_along(v))
    print(v, quote = FALSE)
  }

  if (!is.null(x$se$theta)) {
    cat("\nClass locations theta [estimate (SE)]:\n")
    v <- fmt(x$estimates$theta, x$se$theta)
    names(v) <- paste0("Class ", seq_along(v))
    print(v, quote = FALSE)
  }

  if (!is.null(x$se$delta)) {
    cat("\nItem difficulties delta [estimate (SE)]:\n")
    v <- fmt(x$estimates$delta, x$se$delta)
    nm <- names(x$estimates$delta)
    names(v) <- if (!is.null(nm)) nm else paste0("Item ", seq_along(v))
    print(v, quote = FALSE)
  }

  if (!is.null(x$se$item_probs)) {
    cat("\nItem response probabilities [estimate (SE)]:\n")
    m <- matrix(fmt(x$estimates$item_probs, x$se$item_probs),
                nrow = nrow(x$se$item_probs))
    dimnames(m) <- list(paste0("Item ", seq_len(nrow(m))),
                        paste0("Class ", seq_len(ncol(m))))
    print(m, quote = FALSE)
  }

  if (!is.null(x$active_constraints) && x$active_constraints$n_active > 0L) {
    cat("\nNote: ", x$active_constraints$n_active,
        " active inequality constraint(s); SEs for the affected cells are ",
        "NA (consider method = 'bootstrap').\n", sep = "")
  }
  if (!is.null(x$boundary_cells) && any(x$boundary_cells)) {
    cat("Note: ", sum(x$boundary_cells),
        " item probabilit(y/ies) on the 0/1 boundary; SEs set to NA.\n",
        sep = "")
  }
  if (!is.null(x$n_failed) && x$n_failed > 0L) {
    cat("\nNote: ", x$n_failed, " bootstrap refit(s) failed.\n", sep = "")
  }

  invisible(x)
}
