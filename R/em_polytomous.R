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
poly_estep <- function(data, item_probs, class_probs, use_cpp = TRUE) {
  if (use_cpp) return(cpp_poly_estep(data, item_probs, class_probs))
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
poly_expected_counts <- function(data, posteriors, cat_counts, use_cpp = TRUE) {
  if (use_cpp) {
    return(cpp_poly_expected_counts(data, posteriors, as.integer(cat_counts)))
  }
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
  if (model_type == "UN") return(function(ec, warm) poly_mstep_un(ec))
  cache <- new.env(parent = emptyenv())              # constant constraint matrices

  # MON alone is separable across items (each item's stochastic-ordering
  # constraint touches only that item), so it is solved as J small independent
  # problems - far cheaper than one joint problem. IIO and DM couple items
  # (via the cross-item expected-score ordering) and are solved jointly.
  if (model_type == "MON") {
    return(function(ec, warm) {
      if (is.null(cache$built)) {
        C <- nrow(ec[[1]])
        cache$mats <- lapply(cat_counts, function(m) {
          list(Aeq = .poly_Aeq(m, C), B = .poly_Bmon(m, C))
        })
        cache$built <- TRUE
      }
      .poly_solve_mon_separable(ec, warm, cache$mats)
    })
  }

  mon <- model_type == "DM"
  function(ec, warm) {
    if (is.null(cache$built)) {
      C <- nrow(ec[[1]])
      cache$Aeq <- .poly_Aeq(cat_counts, C)
      cache$B <- rbind(if (mon) .poly_Bmon(cat_counts, C),
                       .poly_Biio(cat_counts, C, item_order))
      # per-item MON matrices for the DM shortcut
      if (mon) cache$mon_mats <- lapply(cat_counts, function(m)
        list(Aeq = .poly_Aeq(m, C), B = .poly_Bmon(m, C)))
      cache$built <- TRUE
    }
    # Exact shortcut: DM = MON n IIO, IIO = UN n IIO. Solve the cheaper relaxation
    # (per-item MON for DM, closed-form UN for IIO); if it already satisfies the
    # cross-item ordering it is feasible for - and maximises over a superset of -
    # the full region, so it IS the constrained maximiser. Only fall back to the
    # coupled joint solve when the ordering is violated.
    relax <- if (mon) .poly_solve_mon_separable(ec, warm, cache$mon_mats)
             else poly_mstep_un(ec)
    if (.poly_iio_ok(relax, cat_counts, item_order)) return(relax)
    .poly_solve_constrained(ec, cat_counts, mon = mon, iio = TRUE,
                            item_order = item_order, warm = warm,
                            Aeq = cache$Aeq, B = cache$B)
  }
}

# Does an item_probs list satisfy invariant item ordering: within every class,
# expected item scores are non-decreasing across items in item_order?
.poly_iio_ok <- function(item_probs, cat_counts, item_order, tol = 1e-8) {
  C <- nrow(item_probs[[1]])
  escore <- vapply(seq_along(item_probs), function(j)
    as.numeric(item_probs[[j]] %*% (0:cat_counts[j])), numeric(C))  # C x J
  ordered <- escore[, item_order, drop = FALSE]
  all(apply(ordered, 1L, function(r) all(diff(r) >= -tol)))
}

# Separable per-item MON solve: each item is an independent small SLSQP problem.
.poly_solve_mon_separable <- function(ec, warm, mats) {
  lapply(seq_along(ec), function(j) {
    e <- ec[[j]]; C <- nrow(e); K <- ncol(e)
    w <- as.vector(t(e))
    p0 <- if (is.null(warm)) as.vector(t(e / rowSums(e))) else as.vector(t(warm[[j]]))
    p0 <- pmin(pmax(p0, 1e-6), 1)
    sol <- cpp_poly_mstep_solve(w, p0, mats[[j]]$Aeq, mats[[j]]$B,
                                xtol_rel = 1e-8, maxeval = 500L)
    P <- matrix(sol, C, K, byrow = TRUE); P <- pmax(P, 0); P / rowSums(P)
  })
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
                                    item_order = NULL, warm = NULL,
                                    Aeq = NULL, B = NULL, use_cpp = TRUE) {
  C <- nrow(ec[[1]])
  w <- .poly_pack(ec)
  p0 <- if (is.null(warm)) .poly_pack(poly_mstep_un(ec)) else .poly_pack(warm)
  p0 <- pmin(pmax(p0, 1e-6), 1)
  if (is.null(Aeq)) Aeq <- .poly_Aeq(cc, C)
  if (is.null(B)) {
    B <- rbind(if (mon) .poly_Bmon(cc, C),
               if (iio) .poly_Biio(cc, C, item_order))
  }
  if (use_cpp) {
    # same SLSQP algorithm as below, run entirely in C++ (no R callbacks)
    Bm <- if (is.null(B)) matrix(0, 0, length(p0)) else B
    sol <- cpp_poly_mstep_solve(w, p0, Aeq, Bm, xtol_rel = 1e-8, maxeval = 500L)
  } else {
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
    sol <- res$solution
  }
  P <- .poly_unpack(sol, cc, C)
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
compute_pcm_probs <- function(theta, delta_list, use_cpp = TRUE) {
  if (use_cpp) return(lapply(delta_list, function(dj) cpp_pcm_probs(theta, dj)))
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
  Rmats <- lapply(cat_counts, function(m) {
    K <- m + 1L; R <- matrix(0, K, K); R[lower.tri(R, diag = TRUE)] <- 1; R
  })
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
    # analytic gradient over (theta_1..C, dfree_1..sumM-1); dfree maps to
    # sum-zero delta, so d/d(dfree_i) = gd_i - gd_last
    grQ <- function(par) {
      th <- par[seq_len(C)]; df <- par[(C + 1L):(C + sumM - 1L)]
      ip <- compute_pcm_probs(th, unpack_delta(df))
      g_theta <- numeric(C); gd <- numeric(sumM)
      for (j in seq_len(J)) {
        P <- ip[[j]]; k <- 0:(ncol(P) - 1L)
        ecj <- ec[[j]]; n_c <- rowSums(ecj)
        g_theta <- g_theta - (as.numeric(ecj %*% k) - n_c * as.numeric(P %*% k))
        resid <- (ecj %*% Rmats[[j]]) - (P %*% Rmats[[j]]) * n_c
        gd[idx[[j]]] <- colSums(resid)[-1L]
      }
      c(g_theta, gd[-sumM] - gd[sumM])
    }
    opt <- stats::optim(c(theta, dfree), negQ, gr = grQ, method = "BFGS",
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

# Expected score, variance (= item information) and third central moment of a
# partial-credit item at locations theta (each a vector over theta).
.pcm_moments <- function(theta, delta_j) {
  P <- cpp_pcm_probs(theta, delta_j)             # length(theta) x (m+1)
  k <- 0:(ncol(P) - 1L)
  e <- as.numeric(P %*% k)
  v <- as.numeric(P %*% (k^2)) - e^2
  m3 <- as.numeric(P %*% (k^3)) - 3 * e * as.numeric(P %*% (k^2)) + 2 * e^3
  list(E = e, V = v, M3 = m3)
}

# ---- Rasch / partial-credit model by marginal ML (our own engine) -------
# Gauss-Hermite nodes and weights for a standard normal N(0,1) (weights sum
# to 1), via the Golub-Welsch eigendecomposition of the Jacobi matrix.
.gauss_hermite <- function(n) {
  i <- seq_len(n - 1L); b <- sqrt(i / 2)
  Jm <- matrix(0, n, n)
  Jm[cbind(i, i + 1L)] <- b; Jm[cbind(i + 1L, i)] <- b
  eg <- eigen(Jm, symmetric = TRUE)
  ord <- order(eg$values)
  x <- eg$values[ord]; v1 <- eg$vectors[1, ord]
  list(nodes = sqrt(2) * x, weights = v1^2)   # transformed to N(0,1); sum(w)=1
}

# Marginal-ML partial-credit / Rasch model. Person locations theta ~ N(0,
# sigma^2) are integrated out by Gauss-Hermite quadrature; the quadrature nodes
# play the role of latent classes in the polytomous E-step, so the same engine
# fits the latent-class models and this continuous model. Estimates the item
# step parameters delta (pooled sum-zero) and the latent SD sigma. Dichotomous
# data (every m_j = 1) gives the ordinary Rasch model.
em_rasch_mml <- function(data, n_quad = 61L, max_iter = 500L, tol = 1e-7,
                         use_cpp = TRUE, seed = NULL) {
  data <- as.matrix(data); storage.mode(data) <- "integer"
  n <- nrow(data); J <- ncol(data)
  cat_counts <- .item_cat_counts(data); sumM <- sum(cat_counts)
  gh <- .gauss_hermite(n_quad); z <- gh$nodes; w <- gh$weights
  idx <- split(seq_len(sumM), rep(seq_len(J), cat_counts))
  # The N(0, sigma^2) nodes fix the latent mean at 0, so all sumM step
  # parameters are free (location is already identified); estimate delta + sigma.
  unpack <- function(dfree) lapply(idx, function(ii) dfree[ii])
  if (!is.null(seed)) set.seed(seed)
  # right-cumulative matrices per item: (P %*% Rmat)[, b] = P(X >= b - 1)
  Rmats <- lapply(cat_counts, function(m) {
    K <- m + 1L; R <- matrix(0, K, K); R[lower.tri(R, diag = TRUE)] <- 1; R
  })
  sigma <- 1; dfree <- stats::rnorm(sumM, 0, 0.3)
  ll_old <- -Inf; conv <- FALSE; it <- 0L
  for (it in seq_len(max_iter)) {
    ip <- compute_pcm_probs(sigma * z, unpack(dfree), use_cpp)
    es <- poly_estep(data, ip, w, use_cpp)
    if (is.finite(ll_old) && abs(es$loglik - ll_old) < tol * (abs(ll_old) + tol)) {
      conv <- TRUE; break
    }
    ll_old <- es$loglik
    ec <- poly_expected_counts(data, es$posteriors, cat_counts, use_cpp)
    negQ <- function(par) {
      s <- exp(par[1L]); ip2 <- compute_pcm_probs(s * z, unpack(par[-1L]), use_cpp)
      -sum(vapply(seq_len(J), function(j)
        sum(ec[[j]] * log(.bound(ip2[[j]]))), numeric(1)))
    }
    # analytic gradient: expected-count residuals of the partial-credit model
    grQ <- function(par) {
      s <- exp(par[1L]); ip2 <- compute_pcm_probs(s * z, unpack(par[-1L]), use_cpp)
      g_delta <- numeric(sumM); g_theta <- numeric(length(z))
      for (j in seq_len(J)) {
        P <- ip2[[j]]; K <- ncol(P); k <- 0:(K - 1L)
        ecj <- ec[[j]]; n_c <- rowSums(ecj)
        g_theta <- g_theta + (as.numeric(ecj %*% k) - n_c * as.numeric(P %*% k))
        # d(negQ)/d(delta_jl) = sum_c [ N(X>=l) - n_c P(X>=l) ], l = 1..K-1
        resid <- (ecj %*% Rmats[[j]]) - (P %*% Rmats[[j]]) * n_c
        g_delta[idx[[j]]] <- colSums(resid)[-1L]
      }
      c(-s * sum(z * g_theta), g_delta)
    }
    opt <- stats::optim(c(log(sigma), dfree), negQ, gr = grQ, method = "BFGS",
                        control = list(maxit = 50L))
    sigma <- exp(opt$par[1L]); dfree <- opt$par[-1L]
  }
  ip <- compute_pcm_probs(sigma * z, unpack(dfree), use_cpp)
  es <- poly_estep(data, ip, w, use_cpp)
  list(loglik = es$loglik, sigma = sigma, delta = unlist(unpack(dfree)),
       delta_list = unpack(dfree), nodes = sigma * z, weights = w,
       posteriors = es$posteriors, cat_counts = cat_counts, n_par = sumM + 1L,
       data = data, scores = rowSums(data), iterations = it, converged = conv)
}

# Re-estimate a fitted Rasch / partial-credit model with a FREE latent
# distribution: the Bock-Aitkin empirical-histogram (semi-nonparametric MML)
# refit. Nodes are fixed on an equally-spaced grid; the node weights are
# re-estimated as posterior masses jointly with the item step parameters by
# EM. Unlike a single posterior-mass pass at the normal-prior solution (which
# is smeared by measurement error toward the prior), the joint iteration
# deconvolves, so the recovered latent distribution reproduces bimodal,
# skewed, or censored ability distributions. Returns a modified rm_fit whose
# nodes / weights / delta / posteriors describe the empirical-latent solution.
.rm_empirical_refit <- function(rm_fit, n_nodes = 41L, max_iter = 200L,
                                tol = 1e-7, use_cpp = TRUE) {
  data <- rm_fit$data
  cat_counts <- rm_fit$cat_counts
  J <- ncol(data); sumM <- sum(cat_counts)
  idx <- split(seq_len(sumM), rep(seq_len(J), cat_counts))
  unpack <- function(d) lapply(idx, function(ii) d[ii])
  Rmats <- lapply(cat_counts, function(m) {
    K <- m + 1L; R <- matrix(0, K, K); R[lower.tri(R, diag = TRUE)] <- 1; R
  })
  z <- seq(-4, 4, length.out = n_nodes) * rm_fit$sigma
  w <- rep(1 / n_nodes, n_nodes)
  dfree <- rm_fit$delta
  ll_old <- -Inf
  for (it in seq_len(max_iter)) {
    ip <- compute_pcm_probs(z, unpack(dfree), use_cpp)
    es <- poly_estep(data, ip, w, use_cpp)
    if (is.finite(ll_old) &&
        abs(es$loglik - ll_old) < tol * (abs(ll_old) + tol)) break
    ll_old <- es$loglik
    w <- pmax(colMeans(es$posteriors), 1e-10); w <- w / sum(w)
    ec <- poly_expected_counts(data, es$posteriors, cat_counts, use_cpp)
    negQ <- function(par) {
      ip2 <- compute_pcm_probs(z, unpack(par), use_cpp)
      -sum(vapply(seq_len(J), function(j)
        sum(ec[[j]] * log(.bound(ip2[[j]]))), numeric(1)))
    }
    grQ <- function(par) {
      ip2 <- compute_pcm_probs(z, unpack(par), use_cpp)
      gd <- numeric(sumM)
      for (j in seq_len(J)) {
        P <- ip2[[j]]; ecj <- ec[[j]]; n_c <- rowSums(ecj)
        resid <- (ecj %*% Rmats[[j]]) - (P %*% Rmats[[j]]) * n_c
        gd[idx[[j]]] <- colSums(resid)[-1L]
      }
      gd
    }
    opt <- stats::optim(dfree, negQ, gr = grQ, method = "BFGS",
                        control = list(maxit = 50L))
    dfree <- opt$par
  }
  es <- poly_estep(data, compute_pcm_probs(z, unpack(dfree), use_cpp), w,
                   use_cpp)
  out <- rm_fit
  out$nodes <- z; out$weights <- w; out$posteriors <- es$posteriors
  out$delta <- dfree; out$delta_list <- unpack(dfree)
  out$empirical_latent <- TRUE
  out
}

# Draw person locations for a marginal parametric bootstrap from a fitted
# Rasch / partial-credit model.
#   latent = "normal":    theta ~ N(0, sigma^2), the parametric assumption.
#   latent = "empirical": theta from the empirical-histogram latent
#     distribution (see .rm_empirical_refit - the rm_fit passed in should be
#     the refit object so the weights are deconvolved), sampled
#     histogram-style with uniform jitter between node midpoints. This
#     reproduces the observed ability distribution (bimodal, skewed, censored,
#     ...) so the null matches the data's sum-score group structure, and any
#     observed-vs-null difference is attributable to non-additivity rather
#     than to population shape. Additive conjoint structure itself is
#     distribution-free, so the latent shape is a nuisance here.
.rm_draw_theta <- function(rm_fit, n, sigma, latent = "empirical") {
  if (identical(latent, "normal") || is.null(rm_fit) ||
      is.null(rm_fit$weights)) {
    return(stats::rnorm(n, 0, sigma))
  }
  z <- rm_fit$nodes
  w <- if (isTRUE(rm_fit$empirical_latent)) rm_fit$weights
       else colMeans(rm_fit$posteriors)
  w <- pmax(w, 0); w <- w / sum(w)
  Q <- length(z)
  mid <- c(z[1] - (z[2] - z[1]) / 2,
           (z[-1] + z[-Q]) / 2,
           z[Q] + (z[Q] - z[Q - 1]) / 2)          # histogram bin edges
  idx <- sample.int(Q, n, replace = TRUE, prob = w)
  stats::runif(n, mid[idx], mid[idx + 1L])
}

# Simulate polytomous responses from a fitted partial-credit model (our own
# em_rasch_mml fit) for a marginal parametric bootstrap. Person locations are
# supplied via `theta` or drawn according to `latent` (see .rm_draw_theta).
# Used by the polytomous CC / Kara bootstrap nulls.
.simulate_pcm <- function(rm_fit, n_obs, sigma, seed = NULL, theta = NULL,
                          latent = "normal") {
  if (!is.null(seed)) set.seed(seed)
  dl <- rm_fit$delta_list                        # per-item step parameters
  J <- length(dl)
  if (is.null(theta)) theta <- .rm_draw_theta(rm_fit, n_obs, sigma, latent)
  out <- matrix(0L, n_obs, J)
  for (j in seq_len(J)) {
    P <- cpp_pcm_probs(theta, dl[[j]])           # n_obs x (m+1)
    cdf <- t(apply(P, 1L, cumsum))
    out[, j] <- rowSums(stats::runif(n_obs) > cdf)
  }
  out
}

# One EM run from given starting values, for a supplied M-step function that
# maps expected counts -> item_probs list.
.poly_em_run <- function(data, item_probs, class_probs, cat_counts, mstep_fn,
                         max_iter = 1000L, tol = 1e-6, use_cpp = TRUE) {
  ll_old <- -Inf; conv <- FALSE; it <- 0L
  post <- NULL
  for (it in seq_len(max_iter)) {
    es <- poly_estep(data, item_probs, class_probs, use_cpp = use_cpp)
    post <- es$posteriors
    if (is.finite(ll_old) &&
        abs(es$loglik - ll_old) < tol * (abs(ll_old) + tol)) {
      conv <- TRUE; break
    }
    ll_old <- es$loglik
    class_probs <- colSums(post) / nrow(data)
    ec <- poly_expected_counts(data, post, cat_counts, use_cpp = use_cpp)
    item_probs <- mstep_fn(ec, item_probs)   # warm-start from current probs
  }
  es <- poly_estep(data, item_probs, class_probs, use_cpp = use_cpp)
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
                         seed = NULL, item_order = NULL, use_cpp = TRUE) {
  storage.mode(data) <- "integer"          # C++ E-step expects an integer matrix
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
      .poly_em_run(data, ip, cp, cat_counts, mstep_fn, max_iter, tol, use_cpp),
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
                         item_order = NULL, call = NULL, use_cpp = TRUE) {
  data <- .validate_poly(data)
  if (!is.numeric(n_classes) || n_classes < 2 || n_classes != round(n_classes)) {
    stop("n_classes must be an integer >= 2")
  }
  fit <- poly_lca_fit(data, n_classes, model_type, n_starts = n_starts,
                      max_iter = max_iter, tol = tol, seed = seed,
                      item_order = item_order, use_cpp = use_cpp)
  .build_poly_qlfit(fit, model_type, data, call = call,
                    item_order = fit$item_order)
}
