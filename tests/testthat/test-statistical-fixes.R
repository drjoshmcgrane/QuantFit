# Tests for the statistical fixes: weighted PAVA, exact constrained M-steps,
# Dykstra 2D projection, parameter counts, G-squared, and LR test behavior

# ---------------------------------------------------------------------------
# Weighted PAVA
# ---------------------------------------------------------------------------

test_that("pava_increasing matches stats::isoreg for unit weights", {
  set.seed(42)
  for (rep in 1:5) {
    x <- rnorm(25)
    expect_equal(pava_increasing(x), isoreg(x)$yf, tolerance = 1e-12)
  }

  # Ties and already-monotone inputs
  x <- c(1, 1, 1, 2, 3)
  expect_equal(pava_increasing(x), x)
  x <- c(3, 2, 1)
  expect_equal(pava_increasing(x), rep(2, 3))
})

test_that("weighted pava_increasing solves the weighted isotonic QP", {
  skip_if_not_installed("nloptr")

  set.seed(123)
  for (rep in 1:3) {
    n <- 15
    x <- rnorm(n)
    w <- runif(n, 0.1, 10)

    y_pava <- pava_increasing(x, w)

    # Monotone
    expect_true(all(diff(y_pava) >= -1e-12))

    # Brute-force QP reference via SLSQP
    obj <- function(y) list(objective = sum(w * (y - x)^2),
                            gradient = 2 * w * (y - x))
    jac <- matrix(0, n - 1, n)
    for (i in 1:(n - 1)) { jac[i, i] <- 1; jac[i, i + 1] <- -1 }
    con <- function(y) list(constraints = -diff(y), jacobian = jac)

    res <- nloptr::nloptr(
      x0 = sort(x), eval_f = obj, eval_g_ineq = con,
      opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-12,
                  maxeval = 5000)
    )

    expect_equal(sum(w * (y_pava - x)^2), res$objective, tolerance = 1e-8)
    # 1e-3 not 1e-6: SLSQP's solution position converges more loosely than
    # its objective (flat directions), and testthat's tolerance is relative,
    # so elements near zero magnify solver noise (seen at ~6e-4 on Linux CI).
    # PAVA is the exact analytic answer; the QP is only a cross-check.
    expect_equal(y_pava, res$solution, tolerance = 1e-3)
  }
})

test_that("weighted PAVA does not double-count weights (regression test)", {
  # With the old bug (both pooled elements kept combined weight), pooling
  # x = c(2, 1, 0) with unit weights gave a wrong second merge.
  x <- c(2, 1, 0)
  expect_equal(pava_increasing(x), rep(1, 3))  # mean of all three

  # Weighted case with a known analytic answer:
  # x = (4, 0), w = (1, 3) -> pooled value = (4*1 + 0*3)/4 = 1
  expect_equal(pava_increasing(c(4, 0), c(1, 3)), c(1, 1))

  # Three elements where the merge cascades:
  # x = (3, 2, -1), w = (1, 1, 2): pool (2,-1) -> 0 with weight 3;
  # then 3 > 0 -> pool -> (3*1 + 0*3)/4 = 0.75
  expect_equal(pava_increasing(c(3, 2, -1), c(1, 1, 2)), rep(0.75, 3))
})

test_that("pava_decreasing is the reverse-apply-reverse of pava_increasing", {
  set.seed(7)
  x <- rnorm(12)
  w <- runif(12, 0.5, 3)
  expect_equal(pava_decreasing(x, w), rev(pava_increasing(rev(x), rev(w))))
  expect_true(all(diff(pava_decreasing(x, w)) <= 1e-12))
})

# ---------------------------------------------------------------------------
# Exact MON M-step
# ---------------------------------------------------------------------------

test_that("MON M-step via weighted PAVA is the exact Q maximizer", {
  skip_if_not_installed("nloptr")

  set.seed(11)
  n <- 60
  n_items <- 4
  n_classes <- 3
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  # Random valid posteriors
  posteriors <- matrix(rgamma(n * n_classes, 1), n, n_classes)
  posteriors <- posteriors / rowSums(posteriors)

  constraints_spec <- specify_constraints(class_monotonicity = TRUE)
  m_exact <- m_step_exact(data, posteriors, constraints_spec)

  # Expected complete-data log-likelihood for item parameters
  Q <- function(item_probs) {
    ll <- 0
    for (c in 1:n_classes) {
      w <- posteriors[, c]
      p <- bound_probs(item_probs[, c])
      ll <- ll + sum(w * (data %*% log(p) + (1 - data) %*% log(1 - p)))
    }
    ll
  }

  # Reference: maximize Q directly under monotonicity with SLSQP
  n_par <- n_items * n_classes
  x0 <- as.vector(t(apply(matrix(0.5, n_items, n_classes), 1, function(r)
    seq(0.3, 0.7, length.out = n_classes))))
  obj <- function(par) -Q(matrix(par, n_items, n_classes))
  con <- function(par) {
    m <- matrix(par, n_items, n_classes)
    -as.vector(t(m[, -1, drop = FALSE] - m[, -n_classes, drop = FALSE]))
  }
  res <- nloptr::nloptr(
    x0 = x0, eval_f = obj, eval_g_ineq = con,
    lb = rep(1e-6, n_par), ub = rep(1 - 1e-6, n_par),
    opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-10,
                maxeval = 50000)
  )

  q_pava <- Q(m_exact$item_probs)
  q_opt <- -res$objective

  # Solution satisfies constraints
  for (i in 1:n_items) {
    expect_true(all(diff(m_exact$item_probs[i, ]) >= -1e-8))
  }

  # PAVA M-step should be at least as good as the numerical optimizer
  expect_gte(q_pava, q_opt - 1e-4)
})

# ---------------------------------------------------------------------------
# Dykstra 2D projection for DM
# ---------------------------------------------------------------------------

test_that("Dykstra projection satisfies both constraints and beats naive
           alternating projection in Q", {
  set.seed(21)
  n_items <- 5
  n_classes <- 3
  item_order <- c(2, 4, 1, 5, 3)
  class_counts <- c(120, 15, 260)  # unequal weights matter here

  p_bar <- matrix(runif(n_items * n_classes, 0.05, 0.95), n_items, n_classes)

  p_dyk <- dykstra_dm_projection(p_bar, item_order, class_counts)

  # Both constraint sets satisfied
  for (i in 1:n_items) expect_true(all(diff(p_dyk[i, ]) >= -1e-8))
  for (c in 1:n_classes) expect_true(all(diff(p_dyk[item_order, c]) <= 1e-8))

  # Old naive alternating projection (no correction vectors, unweighted)
  naive <- p_bar
  for (iter in 1:100) {
    old <- naive
    for (i in 1:n_items) naive[i, ] <- pava_increasing(naive[i, ])
    for (c in 1:n_classes) {
      naive[item_order, c] <- pava_decreasing(naive[item_order, c])
    }
    if (max(abs(naive - old)) < 1e-10) break
  }

  # Q(p) = sum_c N_c sum_i [pbar log p + (1 - pbar) log(1 - p)]
  Qfun <- function(p) {
    p <- bound_probs(p)
    sum(sweep(p_bar * log(p) + (1 - p_bar) * log(1 - p), 2, class_counts, "*"))
  }

  expect_gte(Qfun(p_dyk), Qfun(naive) - 1e-8)

  # Weighted L2 distance: Dykstra should also be closer to p_bar
  d2 <- function(p) sum(sweep((p - p_bar)^2, 2, class_counts, "*"))
  expect_lte(d2(p_dyk), d2(naive) + 1e-8)
})

test_that("Dykstra with unit weights matches direct QP projection", {
  skip_if_not_installed("nloptr")

  set.seed(31)
  n_items <- 3
  n_classes <- 3
  item_order <- 1:3
  p_bar <- matrix(runif(9, 0.1, 0.9), n_items, n_classes)

  p_dyk <- dykstra_dm_projection(p_bar, item_order)

  # Reference QP: min ||p - p_bar||^2 s.t. rows increasing, cols (in order)
  # decreasing
  obj <- function(par) sum((par - as.vector(p_bar))^2)
  con <- function(par) {
    m <- matrix(par, n_items, n_classes)
    c(
      -as.vector(m[, -1] - m[, -n_classes]),          # rows increasing
      as.vector(m[item_order, , drop = FALSE][-1, ] -
                m[item_order, , drop = FALSE][-n_items, ])  # cols decreasing
    )
  }
  res <- nloptr::nloptr(
    x0 = rep(0.5, 9), eval_f = obj, eval_g_ineq = con,
    lb = rep(0, 9), ub = rep(1, 9),
    opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-12,
                maxeval = 100000)
  )

  expect_equal(sum((as.vector(p_dyk) - as.vector(p_bar))^2), res$objective,
               tolerance = 1e-5)
})

# ---------------------------------------------------------------------------
# Parameter counts
# ---------------------------------------------------------------------------

test_that("RM n_par equals mirt's estimated parameter count (nest)", {
  skip_if_not_installed("mirt")

  set.seed(41)
  n <- 300
  n_items <- 6
  theta <- rnorm(n)
  delta <- seq(-1, 1, length.out = n_items)
  data <- matrix(0, n, n_items)
  for (j in 1:n_items) data[, j] <- rbinom(n, 1, plogis(theta - delta[j]))

  fit <- fit_rm(data, verbose = FALSE)
  mirt_fit <- attr(fit, "mirt_object")

  expect_equal(fit$n_par, as.integer(mirt::extract.mirt(mirt_fit, "nest")))
  expect_equal(fit$n_par, n_items + 1)

  # delta must NOT be re-centered: it should match mirt's -d exactly
  mirt_d <- -mirt::coef(mirt_fit, simplify = TRUE)$items[, "d"]
  expect_equal(unname(fit$delta), unname(mirt_d), tolerance = 1e-10)
})

test_that("LCR n_par equals 2C + J - 2", {
  set.seed(51)
  n <- 150
  n_items <- 5
  n_classes <- 3
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit <- fit_lcr(data, n_classes = n_classes, n_starts = 2, seed = 1)
  expect_equal(fit$n_par, 2 * n_classes + n_items - 2)
  expect_equal(count_parameters("LCR", n_items, n_classes),
               2 * n_classes + n_items - 2)
})

# ---------------------------------------------------------------------------
# Nested-model log-likelihood ordering
# ---------------------------------------------------------------------------

test_that("log-likelihoods respect nesting: UN >= MON, IIO >= DM", {
  set.seed(61)
  # Data with genuine double-monotone structure (Rasch-like)
  n <- 400
  n_items <- 6
  n_classes <- 3
  theta_c <- c(-1.5, 0, 1.5)
  delta_j <- seq(-1.2, 1.2, length.out = n_items)
  cls <- sample(1:n_classes, n, replace = TRUE)
  data <- matrix(0, n, n_items)
  for (j in 1:n_items) data[, j] <- rbinom(n, 1, plogis(theta_c[cls] - delta_j[j]))

  item_order <- estimate_item_order(data)

  fit_un_r <- fit_un(data, n_classes, n_starts = 5, seed = 100)
  fit_mon_r <- fit_mon(data, n_classes, n_starts = 5, seed = 100)
  fit_iio_r <- fit_iio(data, n_classes, item_order = item_order,
                       n_starts = 5, seed = 100)
  fit_dm_r <- fit_dm(data, n_classes, item_order = item_order,
                     n_starts = 5, seed = 100)

  tol <- 1e-2
  expect_gte(fit_un_r$loglik, fit_mon_r$loglik - tol)
  expect_gte(fit_un_r$loglik, fit_iio_r$loglik - tol)
  expect_gte(fit_mon_r$loglik, fit_dm_r$loglik - tol)
  expect_gte(fit_iio_r$loglik, fit_dm_r$loglik - tol)

  # On well-fitting monotone data the MON fit should be close to UN
  expect_lte(fit_un_r$loglik - fit_mon_r$loglik, 5)
})

# ---------------------------------------------------------------------------
# G-squared on a fully enumerated example
# ---------------------------------------------------------------------------

test_that("g_squared matches hand computation on a tiny example", {
  # 2 items, 8 observations: patterns 00 x2, 01 x1, 10 x1, 11 x4
  data <- rbind(
    matrix(rep(c(0, 0), 2), ncol = 2, byrow = TRUE),
    c(0, 1),
    c(1, 0),
    matrix(rep(c(1, 1), 4), ncol = 2, byrow = TRUE)
  )

  # Single-class model with p = (0.5, 0.5): every pattern has prob 0.25,
  # expected frequency E = 8 * 0.25 = 2
  object <- structure(
    list(
      model_type = "UN",
      item_probs = matrix(c(0.5, 0.5), nrow = 2, ncol = 1),
      class_probs = 1,
      n_par = 2
    ),
    class = "qlfit"
  )

  res <- g_squared(object, data)

  # G2 = 2 * [2 log(2/2) + 1 log(1/2) + 1 log(1/2) + 4 log(4/2)] = 4 log 2
  expect_equal(res$statistic, 4 * log(2), tolerance = 1e-10)

  # df = 2^2 - 1 - 2 = 1
  expect_equal(res$df, 1)
  expect_equal(res$p_value, pchisq(4 * log(2), 1, lower.tail = FALSE),
               tolerance = 1e-10)
})

test_that("pearson_chisq accounts for unobserved patterns and full df", {
  # Same setup; X^2 = sum (O-E)^2/E over all 4 cells (all observed here)
  data <- rbind(
    matrix(rep(c(0, 0), 2), ncol = 2, byrow = TRUE),
    c(0, 1),
    c(1, 0),
    matrix(rep(c(1, 1), 4), ncol = 2, byrow = TRUE)
  )
  object <- structure(
    list(
      model_type = "UN",
      item_probs = matrix(c(0.5, 0.5), nrow = 2, ncol = 1),
      class_probs = 1,
      n_par = 2
    ),
    class = "qlfit"
  )

  res <- pearson_chisq(object, data)
  # E = 2 in each of 4 cells: (0 + 1 + 1 + 4)/2 = 3
  expect_equal(res$statistic, 3, tolerance = 1e-10)
  expect_equal(res$df, 1)
})

# ---------------------------------------------------------------------------
# LR test with df = 0
# ---------------------------------------------------------------------------

test_that("lr_test returns NA p-value with message when df = 0", {
  make_fit <- function(type, ll, k) {
    structure(list(model_type = type, loglik = ll, n_par = k, n_obs = 100),
              class = "qlfit")
  }

  m_un <- make_fit("UN", -500, 17)
  m_mon <- make_fit("MON", -502, 17)

  expect_warning(res <- lr_test(m_mon, m_un), "chi-bar-squared")
  expect_true(is.na(res$p_value))
  expect_equal(res$df, 0)
  expect_true(!is.null(res$note))

  # df > 0 still gives a standard p-value
  m_lcr <- make_fit("LCR", -510, 9)
  suppressWarnings(res2 <- lr_test(m_lcr, m_un))
  expect_false(is.na(res2$p_value))
  expect_equal(res2$df, 8)
})

# ---------------------------------------------------------------------------
# Convergence guard
# ---------------------------------------------------------------------------

test_that("check_convergence warns on a meaningful log-likelihood decrease", {
  ll <- c(-1000, -900, -950)  # big decrease at the end
  expect_warning(res <- check_convergence(ll, tol = 1e-6),
                 "decreased")
  expect_false(res)

  # Small increase within tolerance -> converged, no warning
  ll2 <- c(-1000, -900, -900 + 1e-9)
  expect_silent(res2 <- check_convergence(ll2, tol = 1e-6))
  expect_true(res2)
})

# ---------------------------------------------------------------------------
# Quantile initialization jitter
# ---------------------------------------------------------------------------

test_that("quantile initialization differs across seeds (jitter)", {
  set.seed(71)
  data <- matrix(rbinom(100 * 5, 1, 0.5), nrow = 100)

  p1 <- init_item_probs(data, 3, "quantiles", seed = 1)
  p2 <- init_item_probs(data, 3, "quantiles", seed = 2)
  expect_gt(max(abs(p1 - p2)), 0)
})
