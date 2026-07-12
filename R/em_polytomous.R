#' Polytomous latent-class EM engine
#'
#' Internal machinery for fitting the latent-structure models to *polytomous*
#' item response data, where each item is scored with consecutive integer
#' categories \eqn{0, 1, \ldots, m_j}. Each item-by-class combination is
#' described by a full category-probability vector on the simplex, generalising
#' the single success probability of the dichotomous engine. The unconstrained,
#' class-monotone, invariant-item-ordering, and double-monotone M-steps are all
#' handled here; the dichotomous case (every \eqn{m_j = 1}) reduces to the
#' Bernoulli engine in \code{em_algorithm.R}.
#'
#' @name em_polytomous
#' @keywords internal
NULL

# Is this response matrix polytomous (any score above 1)?
.is_polytomous <- function(data) {
  vals <- unique(as.vector(data[!is.na(data)]))
  any(vals > 1)
}

# Validate either dichotomous or polytomous data, dispatching by content.
validate_data_any <- function(data) {
  if (is.data.frame(data)) data <- as.matrix(data)
  if (.is_polytomous(data)) .validate_poly(data) else validate_data(data)
}

# Number of category steps m_j per item (categories run 0..m_j).
.item_cat_counts <- function(data) apply(data, 2L, max)

# Validate polytomous responses: complete, non-negative integers, and each item
# uses consecutive categories starting at 0.
.validate_poly <- function(data) {
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data)) stop("Data must be a matrix or data frame")
  if (any(is.na(data))) stop("Polytomous fitting requires complete data (no NA).")
  if (any(data != round(data)) || any(data < 0)) {
    stop("Polytomous responses must be non-negative integer category scores ",
         "starting at 0.")
  }
  storage.mode(data) <- "integer"
  for (j in seq_len(ncol(data))) {
    m <- max(data[, j])
    if (m < 1L) stop("Item ", j, " uses a single category; every item must ",
                     "use at least two categories.")
    used <- sort(unique(data[, j]))
    if (!identical(used, 0:m)) {
      stop("Item ", j, " must use consecutive integer categories 0..", m,
           "; collapse or recode empty categories first.")
    }
  }
  data
}

# Polytomous parameter counts (reduce to the dichotomous counts when sumM = J).
count_parameters_poly <- function(model_type, cat_counts, n_classes) {
  sumM <- sum(cat_counts)
  switch(model_type,
    UN = ,
    MON = ,
    IIO = ,
    DM = (n_classes - 1) + n_classes * sumM,
    LCR = (n_classes - 1) + n_classes + (sumM - 1),
    RM  = sumM + 1,
    stop("Unknown model type: ", model_type)
  )
}

.row_lse <- function(M) {
  mx <- apply(M, 1L, max)
  mx + log(rowSums(exp(M - mx)))
}
.bound <- function(p, eps = 1e-12) pmin(pmax(p, eps), 1)

#' Polytomous E-step
#' @param data integer matrix (n x J), entries 0..m_j
#' @param item_probs list length J; `item_probs[[j]]` is C x (m_j + 1), rows sum 1
#' @param class_probs length C
#' @keywords internal
poly_estep <- function(data, item_probs, class_probs) {
  n <- nrow(data); J <- ncol(data); C <- length(class_probs)
  ll <- matrix(log(class_probs), n, C, byrow = TRUE)
  for (j in seq_len(J)) {
    lp <- log(.bound(item_probs[[j]]))          # C x (m+1)
    idx <- data[, j] + 1L
    ll <- ll + t(lp[, idx, drop = FALSE])        # n x C
  }
  lrs <- .row_lse(ll)
  list(posteriors = exp(ll - lrs), loglik = sum(lrs))
}

# Expected category counts: list length J of C x (m_j + 1) matrices.
poly_expected_counts <- function(data, posteriors, cat_counts) {
  J <- ncol(data)
  lapply(seq_len(J), function(j) {
    m <- cat_counts[j]
    ec <- matrix(0, ncol(posteriors), m + 1L)
    for (x in 0:m) ec[, x + 1L] <- colSums(posteriors * (data[, j] == x))
    ec
  })
}

# Unconstrained M-step: weighted category proportions per item and class.
poly_mstep_un <- function(ec) lapply(ec, function(e) e / rowSums(e))

# Return the M-step function (expected counts, warm start) -> item_probs list.
# UN is closed form; MON/IIO/DM solve a linearly-constrained weighted
# multinomial MLE (see .poly_solve_constrained).
poly_mstep_for <- function(model_type, cat_counts, item_order = NULL) {
  switch(model_type,
    UN  = function(ec, warm) poly_mstep_un(ec),
    MON = function(ec, warm)
            .poly_solve_constrained(ec, cat_counts, mon = TRUE, warm = warm),
    IIO = function(ec, warm)
            .poly_solve_constrained(ec, cat_counts, iio = TRUE,
                                    item_order = item_order, warm = warm),
    DM  = function(ec, warm)
            .poly_solve_constrained(ec, cat_counts, mon = TRUE, iio = TRUE,
                                    item_order = item_order, warm = warm),
    stop("Unknown polytomous model type: ", model_type)
  )
}

# ---- Linearly-constrained weighted multinomial M-step -------------------
# Maximise sum_{j,c,x} ec_j[c,x] log P_j[c,x] subject to: each (item,class) row
# on the simplex; MON = cumulative P(X>=x|c) non-decreasing across classes (per
# item, per threshold); IIO = expected item scores E(X_j|c) ordered across items
# (per class) in item_order. All constraints are linear in P, so SLSQP with the
# analytic gradient (-ec/P) is fast and reliable.

.poly_pack   <- function(ip) unlist(lapply(ip, function(P) as.vector(t(P))))
.poly_unpack <- function(p, cc, C) {
  out <- vector("list", length(cc)); pos <- 0L
  for (j in seq_along(cc)) {
    K <- cc[j] + 1L
    out[[j]] <- matrix(p[(pos + 1L):(pos + C * K)], C, K, byrow = TRUE)
    pos <- pos + C * K
  }
  out
}
.poly_offset <- function(cc, C) c(0L, cumsum((cc + 1L) * C))

.poly_Aeq <- function(cc, C) {
  nc <- sum((cc + 1L) * C); off <- .poly_offset(cc, C); rows <- list(); r <- 0L
  for (j in seq_along(cc)) {
    K <- cc[j] + 1L
    for (cl in seq_len(C)) {
      v <- numeric(nc); v[off[j] + (cl - 1L) * K + seq_len(K)] <- 1
      rows[[r <- r + 1L]] <- v
    }
  }
  do.call(rbind, rows)
}
.poly_Bmon <- function(cc, C) {
  if (C < 2) return(NULL)
  nc <- sum((cc + 1L) * C); off <- .poly_offset(cc, C); rows <- list(); r <- 0L
  for (j in seq_along(cc)) {
    K <- cc[j] + 1L
    for (x in seq_len(K - 1L)) for (cl in seq_len(C - 1L)) {
      v <- numeric(nc); cats <- (x + 1L):K
      v[off[j] + (cl - 1L) * K + cats] <-  1   # G[c,x]
      v[off[j] + cl * K + cats]        <- -1   # -G[c+1,x]  (<=0 => increasing)
      rows[[r <- r + 1L]] <- v
    }
  }
  if (r == 0) NULL else do.call(rbind, rows)
}
.poly_Biio <- function(cc, C, order) {
  nc <- sum((cc + 1L) * C); off <- .poly_offset(cc, C); rows <- list(); r <- 0L
  for (cl in seq_len(C)) for (k in seq_len(length(order) - 1L)) {
    a <- order[k]; b <- order[k + 1L]; v <- numeric(nc)
    v[off[a] + (cl - 1L) * (cc[a] + 1L) + seq_len(cc[a] + 1L)] <-  (0:cc[a])
    v[off[b] + (cl - 1L) * (cc[b] + 1L) + seq_len(cc[b] + 1L)] <- -(0:cc[b])
    rows[[r <- r + 1L]] <- v
  }
  if (r == 0) NULL else do.call(rbind, rows)
}

.poly_solve_constrained <- function(ec, cc, mon = FALSE, iio = FALSE,
                                    item_order = NULL, warm = NULL) {
  C <- nrow(ec[[1]])
  w <- .poly_pack(ec)
  p0 <- if (is.null(warm)) .poly_pack(poly_mstep_un(ec)) else .poly_pack(warm)
  p0 <- pmin(pmax(p0, 1e-6), 1)
  Aeq <- .poly_Aeq(cc, C)
  B <- rbind(if (mon) .poly_Bmon(cc, C),
             if (iio) .poly_Biio(cc, C, item_order))
  opts <- list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8, maxeval = 500L)
  res <- nloptr::nloptr(
    x0 = p0,
    eval_f = function(p) -sum(w * log(pmax(p, 1e-12))),
    eval_grad_f = function(p) -w / pmax(p, 1e-12),
    eval_g_eq = function(p) as.numeric(Aeq %*% p - 1),
    eval_jac_g_eq = function(p) Aeq,
    eval_g_ineq = if (is.null(B)) NULL else function(p) as.numeric(B %*% p),
    eval_jac_g_ineq = if (is.null(B)) NULL else function(p) B,
    lb = rep(0, length(p0)), ub = rep(1, length(p0)), opts = opts)
  P <- .poly_unpack(res$solution, cc, C)
  lapply(P, function(m) { m <- pmax(m, 0); m / rowSums(m) })
}

# Estimate a fixed item order (easiest -> hardest) from marginal expected
# scores; used to orient the IIO / DM constraints.
.poly_item_order <- function(data, cat_counts) {
  escore <- vapply(seq_len(ncol(data)), function(j) mean(data[, j]), numeric(1))
  order(escore)
}

# ---- Latent-class partial credit model (LCR, polytomous) ----------------
# Category probabilities under the PCM: for class location theta_c and item
# step parameters delta_j (length m_j),
#   P(X_j = x | c) prop exp( sum_{k<=x} (theta_c - delta_jk) ),  x = 0..m_j.
compute_pcm_probs <- function(theta, delta_list) {
  lapply(delta_list, function(dj) {
    eta <- outer(theta, dj, "-")               # C x m
    num <- cbind(0, t(apply(eta, 1L, cumsum)))  # C x (m+1)
    num <- num - apply(num, 1L, max)
    e <- exp(num); e / rowSums(e)
  })
}

# EM for the polytomous LCR. Identification: all C class locations free, the
# pooled item steps sum to zero (2C + sumM - 2 parameters with class probs).
em_lcr_poly <- function(data, n_classes, max_iter = 500L, tol = 1e-6,
                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(data); J <- ncol(data); C <- n_classes
  cat_counts <- .item_cat_counts(data)
  sumM <- sum(cat_counts)
  idx <- split(seq_len(sumM), rep(seq_len(J), cat_counts))  # delta positions/item
  unpack_delta <- function(dfree) {
    d <- c(dfree, -sum(dfree))                  # sum-zero identification
    lapply(idx, function(ii) d[ii])
  }
  # starts
  theta <- sort(stats::rnorm(C, 0, 1))
  dfree <- stats::rnorm(sumM - 1L, 0, 0.5)
  class_probs <- rep(1 / C, C)

  ll_old <- -Inf; conv <- FALSE; it <- 0L; post <- NULL
  for (it in seq_len(max_iter)) {
    delta_list <- unpack_delta(dfree)
    item_probs <- compute_pcm_probs(theta, delta_list)
    es <- poly_estep(data, item_probs, class_probs)
    post <- es$posteriors
    if (is.finite(ll_old) && abs(es$loglik - ll_old) < tol * (abs(ll_old) + tol)) {
      conv <- TRUE; break
    }
    ll_old <- es$loglik
    class_probs <- colSums(post) / n
    ec <- poly_expected_counts(data, post, cat_counts)
    # M-step: maximise expected complete-data PCM log-likelihood over (theta, dfree)
    negQ <- function(par) {
      th <- par[seq_len(C)]; df <- par[(C + 1L):(C + sumM - 1L)]
      ip <- compute_pcm_probs(th, unpack_delta(df))
      -sum(vapply(seq_len(J), function(j) sum(ec[[j]] * log(.bound(ip[[j]]))), numeric(1)))
    }
    opt <- stats::optim(c(theta, dfree), negQ, method = "BFGS",
                        control = list(maxit = 100L))
    theta <- opt$par[seq_len(C)]; dfree <- opt$par[(C + 1L):(C + sumM - 1L)]
  }
  delta_list <- unpack_delta(dfree)
  item_probs <- compute_pcm_probs(theta, delta_list)
  es <- poly_estep(data, item_probs, class_probs)
  # order classes by location
  ord <- order(theta)
  list(item_probs = lapply(item_probs, function(P) P[ord, , drop = FALSE]),
       class_probs = class_probs[ord], posteriors = es$posteriors[, ord, drop = FALSE],
       theta = theta[ord], delta = unlist(delta_list), loglik = es$loglik,
       iterations = it, converged = conv, cat_counts = cat_counts)
}

# High-level polytomous LCR fit with multiple starts.
fit_lcr_poly <- function(data, n_classes, n_starts = 10L, max_iter = 500L,
                         tol = 1e-6, seed = NULL, call = NULL) {
  data <- .validate_poly(data)
  best <- NULL; best_ll <- -Inf
  for (s in seq_len(n_starts)) {
    ss <- if (!is.null(seed)) seed + s else NULL
    fit <- tryCatch(em_lcr_poly(data, n_classes, max_iter, tol, ss),
                    error = function(e) NULL)
    if (!is.null(fit) && fit$loglik > best_ll) { best <- fit; best_ll <- fit$loglik }
  }
  if (is.null(best)) stop("All polytomous LCR starts failed.")
  res <- .build_poly_qlfit(best, "LCR", data, call = call)
  res$theta <- best$theta; res$delta <- best$delta
  res
}

# ---- Rasch (partial credit) model via mirt, polytomous ------------------
.fit_rm_poly <- function(data, method = "EM", quadpts = 61, verbose = FALSE,
                         call = NULL, ...) {
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' is required for fit_rm(). Please install it.")
  }
  data <- .validate_poly(data)
  n_obs <- nrow(data); n_items <- ncol(data)
  cat_counts <- .item_cat_counts(data)
  mfit <- tryCatch(
    mirt::mirt(as.data.frame(data), model = 1, itemtype = "Rasch",
               method = method, quadpts = quadpts, verbose = verbose, ...),
    error = function(e) stop("mirt PCM estimation failed: ", e$message))
  loglik <- as.numeric(mirt::extract.mirt(mfit, "logLik"))
  n_par <- tryCatch(as.integer(mirt::extract.mirt(mfit, "nest")),
                    error = function(e) sum(cat_counts) + 1L)
  # expected item-score curves at a theta grid (n_items x grid)
  theta_grid <- seq(-4, 4, length.out = 21)
  escore <- t(vapply(seq_len(n_items), function(j)
    as.numeric(mirt::expected.item(mirt::extract.item(mfit, j), theta_grid)),
    numeric(length(theta_grid))))
  res <- new_qlfit(
    model_type = "RM", item_probs = escore, class_probs = NULL,
    posteriors = NULL, loglik = loglik, n_par = n_par, n_obs = n_obs,
    n_items = n_items, n_classes = NA,
    convergence = mirt::extract.mirt(mfit, "converged"),
    iterations = mirt::extract.mirt(mfit, "iterations"),
    call = call, theta = theta_grid, delta = NULL, item_order = NULL,
    constraints = NULL, se = NULL)
  res$polytomous <- TRUE
  res$cat_counts <- cat_counts
  attr(res, "mirt_object") <- mfit
  res
}

# Simulate polytomous responses from a fitted mirt partial-credit model, drawing
# person locations theta ~ N(0, sigma^2) (marginal parametric bootstrap). Used
# by the polytomous CC / Kara bootstrap nulls.
.simulate_pcm <- function(mfit, n_obs, sigma, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  co <- mirt::coef(mfit, simplify = TRUE)$items
  dcols <- grep("^d[0-9]+$", colnames(co), value = TRUE)
  D <- co[, dcols, drop = FALSE]                 # items x (maxcat+1), NA-padded
  J <- nrow(D)
  theta <- stats::rnorm(n_obs, 0, sigma)
  out <- matrix(0L, n_obs, J)
  for (j in seq_len(J)) {
    dj <- D[j, ]; dj <- dj[!is.na(dj)]            # d0..d_mj (mirt: d0 = 0)
    m <- length(dj) - 1L
    num <- outer(theta, 0:m) + matrix(dj, n_obs, m + 1L, byrow = TRUE)
    num <- num - apply(num, 1L, max)
    P <- exp(num); P <- P / rowSums(P)
    cdf <- t(apply(P, 1L, cumsum))
    out[, j] <- rowSums(stats::runif(n_obs) > cdf)
  }
  out
}

# One EM run from given starting values, for a supplied M-step function that
# maps expected counts -> item_probs list.
.poly_em_run <- function(data, item_probs, class_probs, cat_counts, mstep_fn,
                         max_iter = 1000L, tol = 1e-6) {
  ll_old <- -Inf; conv <- FALSE; it <- 0L
  post <- NULL
  for (it in seq_len(max_iter)) {
    es <- poly_estep(data, item_probs, class_probs)
    post <- es$posteriors
    if (is.finite(ll_old) &&
        abs(es$loglik - ll_old) < tol * (abs(ll_old) + tol)) {
      conv <- TRUE; break
    }
    ll_old <- es$loglik
    class_probs <- colSums(post) / nrow(data)
    ec <- poly_expected_counts(data, post, cat_counts)
    item_probs <- mstep_fn(ec, item_probs)   # warm-start from current probs
  }
  es <- poly_estep(data, item_probs, class_probs)
  list(item_probs = item_probs, class_probs = class_probs,
       posteriors = es$posteriors, loglik = es$loglik,
       iterations = it, converged = conv)
}

# Random simplex starting values per item and class.
.poly_init <- function(cat_counts, n_classes, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  lapply(cat_counts, function(m) {
    P <- matrix(stats::runif(n_classes * (m + 1L), 0.5, 1.5), n_classes, m + 1L)
    P / rowSums(P)
  })
}

#' Fit a polytomous latent-structure model by EM with multiple starts
#'
#' @param data integer response matrix (validated)
#' @param n_classes number of latent classes
#' @param model_type one of "UN","MON","IIO","DM"
#' @param n_starts random starts
#' @param max_iter maximum EM iterations per start
#' @param tol relative log-likelihood convergence tolerance
#' @param seed optional base random seed (each start uses seed + start index)
#' @param item_order optional fixed item order for IIO/DM
#' @keywords internal
poly_lca_fit <- function(data, n_classes, model_type = "UN",
                         n_starts = 10L, max_iter = 1000L, tol = 1e-6,
                         seed = NULL, item_order = NULL) {
  cat_counts <- .item_cat_counts(data)
  if (model_type %in% c("IIO", "DM") && is.null(item_order)) {
    item_order <- .poly_item_order(data, cat_counts)
  }
  mstep_fn <- poly_mstep_for(model_type, cat_counts, item_order)

  best <- NULL; best_ll <- -Inf
  for (s in seq_len(n_starts)) {
    ss <- if (!is.null(seed)) seed + s else NULL
    ip <- .poly_init(cat_counts, n_classes, ss)
    cp <- init_class_probs(n_classes, "random", ss)
    fit <- tryCatch(
      .poly_em_run(data, ip, cp, cat_counts, mstep_fn, max_iter, tol),
      error = function(e) NULL)
    if (!is.null(fit) && fit$loglik > best_ll) { best <- fit; best_ll <- fit$loglik }
  }
  if (is.null(best)) stop("All polytomous EM starts failed.")

  best$cat_counts <- cat_counts
  best$item_order <- item_order
  .poly_order_classes(best, model_type)
}

# Order classes so the object is interpretable. UN orders by prevalence
# (largest first, matching the dichotomous engine); the ordered models order by
# mean category level (low -> high) so class index tracks the latent order.
.poly_order_classes <- function(fit, model_type) {
  if (model_type == "UN") {
    ord <- order(fit$class_probs, decreasing = TRUE)
  } else {
    lvl <- .poly_class_levels(fit$item_probs, fit$cat_counts)
    ord <- order(lvl)
  }
  fit$class_probs <- fit$class_probs[ord]
  fit$posteriors  <- fit$posteriors[, ord, drop = FALSE]
  fit$item_probs  <- lapply(fit$item_probs, function(P) P[ord, , drop = FALSE])
  fit
}

# Mean overall category level per class (summed expected item scores).
.poly_class_levels <- function(item_probs, cat_counts) {
  C <- nrow(item_probs[[1]])
  lvl <- numeric(C)
  for (j in seq_along(item_probs)) {
    lvl <- lvl + as.numeric(item_probs[[j]] %*% (0:cat_counts[j]))
  }
  lvl
}

# Assemble a qlfit from a polytomous EM fit. item_probs is a list of
# C x (m_j + 1) category-probability matrices; the polytomous flag and per-item
# category counts are attached so downstream methods can dispatch.
.build_poly_qlfit <- function(fit, model_type, data, call = NULL,
                              item_order = NULL) {
  n_par <- count_parameters_poly(model_type, fit$cat_counts,
                                  length(fit$class_probs))
  res <- new_qlfit(
    model_type = model_type,
    item_probs = fit$item_probs,
    class_probs = fit$class_probs,
    posteriors = fit$posteriors,
    loglik = fit$loglik,
    n_par = n_par,
    n_obs = nrow(data),
    n_items = ncol(data),
    n_classes = length(fit$class_probs),
    convergence = fit$converged,
    iterations = fit$iterations,
    call = call,
    theta = NULL, delta = NULL,
    item_order = item_order, constraints = NULL, se = NULL
  )
  res$polytomous <- TRUE
  res$cat_counts <- fit$cat_counts
  res
}

# High-level entry: fit any of the LCA-family models, polytomous data.
fit_lca_poly <- function(data, n_classes, model_type, n_starts = 10L,
                         max_iter = 1000L, tol = 1e-6, seed = NULL,
                         item_order = NULL, call = NULL) {
  data <- .validate_poly(data)
  if (!is.numeric(n_classes) || n_classes < 2 || n_classes != round(n_classes)) {
    stop("n_classes must be an integer >= 2")
  }
  fit <- poly_lca_fit(data, n_classes, model_type, n_starts = n_starts,
                      max_iter = max_iter, tol = tol, seed = seed,
                      item_order = item_order)
  .build_poly_qlfit(fit, model_type, data, call = call,
                    item_order = fit$item_order)
}
