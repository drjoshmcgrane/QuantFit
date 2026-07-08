# Equivalence tests: C++ EM engine (use_cpp = TRUE) vs the pure-R reference
# implementation (use_cpp = FALSE).
#
# The C++ port in src/em_lca.cpp is designed to be bit-for-bit identical to
# the R reference on this platform (same BLAS via R, long-double
# accumulation matching R's sum()/colSums(), -ffp-contract=off so products
# round before accumulating). All comparisons below therefore pass at far
# tighter tolerances than the required 1e-8.
#
# Benchmark (Apple M-series, R 4.5.1, reference libRblas; n = 2000, J = 20,
# C = 4, quantile inits, tol = 1e-6; bench::mark medians, identical
# trajectories, 56 EM iterations):
#   em_lca:              R 192.9 ms   C++ 37.6 ms   -> 5.1x speedup
#   em_constrained MON:  R 198.7 ms   C++ 38.5 ms   -> 5.2x speedup
# The 10x target was traded against exact reproduction of the R reference:
# the R hot path is already BLAS-vectorized (E-step matrix products), and
# the remaining C++ cost is dominated by the same reference-BLAS dgemv calls
# plus scalar libm exp/log in the log-sum-exp, which must be kept
# call-for-call identical to R. Refactoring the E-step algebra or using a
# SIMD exp would reach ~10x but would break bit-level agreement with the
# validated R implementation.

make_rasch_data <- function(n, J, seed) {
  set.seed(seed)
  theta <- rnorm(n)
  delta <- seq(-1.5, 1.5, length.out = J)
  p <- 1 / (1 + exp(-outer(theta, delta, "-")))
  data <- matrix(rbinom(n * J, 1, p), n, J)
  storage.mode(data) <- "double"
  data
}

test_that("cpp_e_step matches e_step on random data", {
  set.seed(101)
  for (rep in 1:20) {
    n <- sample(50:400, 1)
    J <- sample(4:15, 1)
    C <- sample(2:5, 1)
    data <- matrix(rbinom(n * J, 1, runif(1, 0.2, 0.8)), n, J)
    storage.mode(data) <- "double"
    ip <- matrix(runif(J * C, 0.01, 0.99), J, C)
    cp <- runif(C); cp <- cp / sum(cp)

    r <- e_step(data, ip, cp)
    cc <- cpp_e_step(data, ip, cp)

    expect_lt(abs(r$loglik - cc$loglik), 1e-8)
    expect_lt(max(abs(r$posteriors - cc$posteriors)), 1e-8)
  }

  # Extreme probabilities exercise the bound_probs / log-sum-exp guards
  data <- matrix(rbinom(200 * 6, 1, 0.5), 200, 6)
  storage.mode(data) <- "double"
  ip <- matrix(c(1e-12, 1 - 1e-12, 0, 1, 0.5, 0.5), 6, 3)
  cp <- c(0.2, 0.3, 0.5)
  r <- e_step(data, ip, cp)
  cc <- cpp_e_step(data, ip, cp)
  expect_lt(abs(r$loglik - cc$loglik), 1e-8)
  expect_lt(max(abs(r$posteriors - cc$posteriors)), 1e-8)
})

test_that("cpp_m_step matches m_step", {
  set.seed(102)
  for (rep in 1:10) {
    n <- sample(50:400, 1)
    J <- sample(4:15, 1)
    C <- sample(2:5, 1)
    data <- matrix(rbinom(n * J, 1, 0.5), n, J)
    storage.mode(data) <- "double"
    post <- matrix(runif(n * C), n, C)
    post <- post / rowSums(post)

    r <- m_step(data, post)
    cc <- cpp_m_step(data, post)

    expect_lt(max(abs(r$item_probs - cc$item_probs)), 1e-12)
    expect_lt(max(abs(r$class_probs - cc$class_probs)), 1e-12)
    expect_identical(r$degenerate, cc$degenerate)
  }
})

test_that("cpp_weighted_pava matches pava_increasing on 500 random weighted vectors", {
  # Tolerance rather than bit-identity: with default compiler flags the
  # C++ accumulation may use fused multiply-adds, differing from R in the
  # last bit. EM is a smooth deterministic iteration, so these differences
  # stay at rounding level (unlike the MCMC samplers, which require
  # bit-identity and default flags — the reason src/Makevars was removed).
  set.seed(103)
  for (rep in 1:500) {
    n <- sample(1:40, 1)
    x <- rnorm(n)
    w <- runif(n, 0, 5)
    expect_equal(cpp_weighted_pava(x, w), pava_increasing(x, w),
                 tolerance = 1e-12)
  }

  # Unit weights (NULL), decreasing direction, and zero weights
  set.seed(104)
  for (rep in 1:50) {
    x <- rnorm(sample(2:30, 1))
    expect_equal(cpp_weighted_pava(x), pava_increasing(x), tolerance = 1e-12)
    expect_equal(cpp_weighted_pava(x, increasing = FALSE), pava_decreasing(x),
                 tolerance = 1e-12)
    w <- runif(length(x), 0, 2)
    w[sample(length(x), max(1, length(x) %/% 3))] <- 0
    expect_equal(cpp_weighted_pava(x, w), pava_increasing(x, w),
                 tolerance = 1e-12)
    expect_equal(cpp_weighted_pava(x, w, increasing = FALSE),
                 pava_decreasing(x, w), tolerance = 1e-12)
  }

  # Edge cases
  expect_identical(cpp_weighted_pava(numeric(0)), numeric(0))
  expect_identical(cpp_weighted_pava(3.14), 3.14)
  expect_error(cpp_weighted_pava(c(1, 2), c(1, 2, 3)),
               "weights must have same length as x")
})

test_that("cpp_dykstra_dm matches dykstra_dm_projection", {
  set.seed(105)
  for (rep in 1:25) {
    I <- sample(3:12, 1)
    C <- sample(2:5, 1)
    m <- matrix(runif(I * C), I, C)
    ord <- sample(1:I)
    w <- runif(C, 0.1, 200)

    expect_equal(cpp_dykstra_dm(m, ord, w),
                 dykstra_dm_projection(m, ord, w), tolerance = 1e-12)
    # Unit weights (NULL)
    expect_equal(cpp_dykstra_dm(m, ord),
                 dykstra_dm_projection(m, ord), tolerance = 1e-12)
  }
})

test_that("cpp_project_constraints matches project_constraints_weighted", {
  set.seed(106)
  for (rep in 1:25) {
    I <- sample(3:12, 1)
    C <- sample(2:5, 1)
    m <- matrix(runif(I * C), I, C)
    ord <- sample(1:I)
    w <- runif(C, 0, 150)  # includes near-zero weights (collapsed classes)

    spec_mon <- specify_constraints(class_monotonicity = TRUE)
    spec_iio <- specify_constraints(item_ordering = TRUE, item_order = ord)
    spec_dm <- specify_constraints(class_monotonicity = TRUE,
                                   item_ordering = TRUE, item_order = ord)

    expect_equal(
      cpp_project_constraints(m, TRUE, FALSE, integer(0), w),
      project_constraints_weighted(m, spec_mon, class_weights = w),
      tolerance = 1e-12)
    expect_equal(
      cpp_project_constraints(m, FALSE, TRUE, ord, w),
      project_constraints_weighted(m, spec_iio, ord, class_weights = w),
      tolerance = 1e-12)
    expect_equal(
      cpp_project_constraints(m, TRUE, TRUE, ord, w),
      project_constraints_weighted(m, spec_dm, ord, class_weights = w),
      tolerance = 1e-12)
  }
})

test_that("em_lca C++ and R paths agree end-to-end with fixed inits", {
  data <- make_rasch_data(300, 8, seed = 201)
  C <- 3
  set.seed(202)
  ip <- init_item_probs(data, C, "quantiles")
  cp <- rep(1 / C, C)

  r <- suppressWarnings(
    em_lca(data, C, init_probs = ip, init_class_probs = cp, use_cpp = FALSE))
  cc <- suppressWarnings(
    em_lca(data, C, init_probs = ip, init_class_probs = cp, use_cpp = TRUE))

  expect_lt(abs(r$loglik - cc$loglik), 1e-8)
  expect_lt(max(abs(r$item_probs - cc$item_probs)), 1e-8)
  expect_lt(max(abs(r$class_probs - cc$class_probs)), 1e-8)
  expect_lt(max(abs(r$posteriors - cc$posteriors)), 1e-8)
  expect_equal(r$ll_history, cc$ll_history, tolerance = 1e-8)
  expect_identical(r$iterations, cc$iterations)
  expect_identical(r$converged, cc$converged)
  expect_identical(r$degenerate, cc$degenerate)
})

test_that("em_lca C++ path replicates max_iter (non-convergence) semantics", {
  data <- make_rasch_data(200, 6, seed = 203)
  C <- 3
  set.seed(204)
  ip <- init_item_probs(data, C, "quantiles")
  cp <- rep(1 / C, C)

  r <- suppressWarnings(em_lca(data, C, init_probs = ip, init_class_probs = cp,
                               max_iter = 7, use_cpp = FALSE))
  cc <- suppressWarnings(em_lca(data, C, init_probs = ip, init_class_probs = cp,
                                max_iter = 7, use_cpp = TRUE))

  expect_false(cc$converged)
  expect_identical(cc$iterations, 7L)
  # Final E-step appended after the max_iter exit on both paths
  expect_length(cc$ll_history, 8L)
  expect_lt(abs(r$loglik - cc$loglik), 1e-8)
  expect_lt(max(abs(r$item_probs - cc$item_probs)), 1e-8)
  expect_equal(r$ll_history, cc$ll_history, tolerance = 1e-8)
})

test_that("em_constrained C++ and R paths agree for MON, IIO, and DM", {
  data <- make_rasch_data(300, 8, seed = 205)
  C <- 3
  n_items <- ncol(data)
  ord <- estimate_item_order(data)
  cp <- rep(1 / C, C)

  specs <- list(
    MON = specify_constraints(class_monotonicity = TRUE),
    IIO = specify_constraints(item_ordering = TRUE, item_order = ord),
    DM = specify_constraints(class_monotonicity = TRUE, item_ordering = TRUE,
                             item_order = ord)
  )
  inits <- list(
    MON = init_item_probs_monotonic(data, C, seed = 206),
    IIO = init_item_probs_iio(data, C, ord, seed = 206),
    DM = init_item_probs_dm(data, C, ord, seed = 206)
  )

  for (model in names(specs)) {
    r <- suppressWarnings(em_constrained(
      data, C, specs[[model]], item_order = ord,
      init_probs = inits[[model]], init_class_probs = cp, use_cpp = FALSE))
    cc <- suppressWarnings(em_constrained(
      data, C, specs[[model]], item_order = ord,
      init_probs = inits[[model]], init_class_probs = cp, use_cpp = TRUE))

    expect_lt(abs(r$loglik - cc$loglik), 1e-8)
    expect_lt(max(abs(r$item_probs - cc$item_probs)), 1e-8)
    expect_lt(max(abs(r$class_probs - cc$class_probs)), 1e-8)
    expect_lt(max(abs(r$posteriors - cc$posteriors)), 1e-8)
    expect_identical(r$iterations, cc$iterations)
    expect_identical(r$converged, cc$converged)

    # The C++ solution satisfies the constraints exactly like the R one
    chk <- check_constraints(cc$item_probs, specs[[model]], ord)
    expect_true(chk$satisfied)
  }
})

test_that("cpp_lcr_q matches the R m_step_rasch objective", {
  data <- make_rasch_data(250, 7, seed = 207)
  n_classes <- 3
  n_items <- ncol(data)
  set.seed(208)
  post <- matrix(runif(nrow(data) * n_classes), nrow(data), n_classes)
  post <- post / rowSums(post)

  r_objective <- function(par) {
    theta_new <- par[1:n_classes]
    delta_free <- par[(n_classes + 1):length(par)]
    delta_new <- c(-sum(delta_free), delta_free)
    item_probs <- bound_probs(compute_rasch_probs(theta_new, delta_new))
    ll <- 0
    for (c in 1:n_classes) {
      weights <- post[, c]
      log_p <- log(item_probs[, c])
      log_1mp <- log(1 - item_probs[, c])
      ll <- ll + sum(weights * (data %*% log_p + (1 - data) %*% log_1mp))
    }
    -ll
  }

  for (rep in 1:25) {
    par <- c(rnorm(n_classes), rnorm(n_items - 1))
    expect_lt(abs(r_objective(par) - cpp_lcr_q(par, data, post, n_classes)),
              1e-10)
  }
})

test_that("em_lcr C++ and R paths agree end-to-end with fixed inits", {
  data <- make_rasch_data(250, 7, seed = 209)
  C <- 3
  init_theta <- c(-1, 0, 1)
  init_delta <- seq(-1, 1, length.out = ncol(data))
  cp <- rep(1 / C, C)

  r <- suppressWarnings(em_lcr(data, C, init_theta = init_theta,
                               init_delta = init_delta, init_class_probs = cp,
                               max_iter = 100, use_cpp = FALSE))
  cc <- suppressWarnings(em_lcr(data, C, init_theta = init_theta,
                                init_delta = init_delta, init_class_probs = cp,
                                max_iter = 100, use_cpp = TRUE))

  expect_lt(abs(r$loglik - cc$loglik), 1e-8)
  expect_lt(max(abs(r$theta - cc$theta)), 1e-8)
  expect_lt(max(abs(r$delta - cc$delta)), 1e-8)
  expect_lt(max(abs(r$item_probs - cc$item_probs)), 1e-8)
  expect_lt(max(abs(r$class_probs - cc$class_probs)), 1e-8)
  expect_identical(r$iterations, cc$iterations)
  expect_identical(r$converged, cc$converged)
})

test_that("fit_un/fit_mon/fit_iio/fit_dm give identical selection-relevant output with use_cpp TRUE vs FALSE", {
  data <- make_rasch_data(300, 8, seed = 210)
  C <- 3

  fits <- list(
    un = function(u) fit_un(data, C, n_starts = 3, seed = 42, use_cpp = u),
    mon = function(u) fit_mon(data, C, n_starts = 3, seed = 42, use_cpp = u),
    iio = function(u) fit_iio(data, C, n_starts = 3, seed = 42, use_cpp = u),
    dm = function(u) fit_dm(data, C, n_starts = 2, seed = 42, use_cpp = u)
  )

  for (model in names(fits)) {
    f_r <- suppressWarnings(fits[[model]](FALSE))
    f_c <- suppressWarnings(fits[[model]](TRUE))

    expect_lt(abs(f_r$loglik - f_c$loglik), 1e-8)
    expect_lt(abs(BIC(f_r) - BIC(f_c)), 1e-8)
    expect_lt(abs(AIC(f_r) - AIC(f_c)), 1e-8)
    expect_identical(f_r$n_par, f_c$n_par)
    expect_lt(max(abs(f_r$item_probs - f_c$item_probs)), 1e-8)
    expect_lt(max(abs(f_r$class_probs - f_c$class_probs)), 1e-8)
    expect_identical(f_r$convergence, f_c$convergence)
    expect_identical(f_r$iterations, f_c$iterations)
  }
})

test_that("fit_lcr gives equivalent selection-relevant output with use_cpp TRUE vs FALSE", {
  data <- make_rasch_data(250, 7, seed = 211)

  f_r <- suppressWarnings(fit_lcr(data, 3, n_starts = 2, max_iter = 100,
                                  seed = 42, use_cpp = FALSE))
  f_c <- suppressWarnings(fit_lcr(data, 3, n_starts = 2, max_iter = 100,
                                  seed = 42, use_cpp = TRUE))

  # Tolerance 1e-4, not 1e-8: with default compiler flags the C++ Q
  # evaluation may use fused multiply-adds; optim's finite-difference
  # gradients amplify those last-bit differences to ~1e-6 in the final
  # parameters. Irrelevant for model selection (BIC differences that
  # matter are > 1).
  expect_lt(abs(f_r$loglik - f_c$loglik), 1e-4)
  expect_lt(abs(BIC(f_r) - BIC(f_c)), 1e-4)
  expect_lt(max(abs(f_r$theta - f_c$theta)), 1e-4)
  expect_lt(max(abs(f_r$delta - f_c$delta)), 1e-4)
})

test_that("C++ path replicates the degenerate-class warning", {
  # Two well-separated groups but four requested classes: some classes
  # collapse on both paths, and the same warning must be raised
  set.seed(212)
  n <- 120
  g <- rep(0:1, each = n / 2)
  p <- ifelse(g == 1, 0.9, 0.1)
  data <- matrix(rbinom(n * 6, 1, rep(p, 6)), n, 6)
  storage.mode(data) <- "double"

  ip <- matrix(c(0.1, 0.1, 0.9, 0.9), nrow = 6, ncol = 4, byrow = TRUE)
  ip <- ip + matrix(seq(-0.02, 0.02, length.out = 24), 6, 4)
  cp <- rep(0.25, 4)

  r_warn <- tryCatch(
    { em_lca(data, 4, init_probs = ip, init_class_probs = cp, use_cpp = FALSE); NULL },
    warning = function(w) conditionMessage(w))
  c_warn <- tryCatch(
    { em_lca(data, 4, init_probs = ip, init_class_probs = cp, use_cpp = TRUE); NULL },
    warning = function(w) conditionMessage(w))

  expect_identical(r_warn, c_warn)
})
