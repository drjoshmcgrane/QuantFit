# Tests for Latent Class Rasch Model (LCR)

test_that("fit_lcr returns valid qlfit object", {
  set.seed(123)
  n <- 200
  n_items <- 6
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit <- fit_lcr(data, n_classes = 3, n_starts = 3, verbose = FALSE)

  expect_s3_class(fit, "qlfit")
  expect_equal(fit$model_type, "LCR")
  expect_equal(fit$n_obs, n)
  expect_equal(fit$n_items, n_items)
  expect_equal(fit$n_classes, 3)
})

test_that("fit_lcr returns theta and delta parameters", {
  set.seed(456)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  fit <- fit_lcr(data, n_classes = 4, n_starts = 3)

  expect_true(!is.null(fit$theta))
  expect_true(!is.null(fit$delta))
  expect_equal(length(fit$theta), 4)
  expect_equal(length(fit$delta), 5)

  # Delta should be centered
  expect_true(abs(mean(fit$delta)) < 0.01)
})

test_that("fit_lcr item_probs follow Rasch parameterization", {
  set.seed(789)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  fit <- fit_lcr(data, n_classes = 3, n_starts = 3)

  # Verify item_probs match Rasch formula
  expected_probs <- compute_rasch_probs(fit$theta, fit$delta)
  expect_equal(fit$item_probs, expected_probs, tolerance = 1e-6)
})

test_that("fit_lcr has correct number of parameters", {
  set.seed(101)
  n <- 200
  n_items <- 6
  n_classes <- 4
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit <- fit_lcr(data, n_classes = n_classes, n_starts = 3)

  # All C thetas are free; mean(delta) = 0 is the single identification
  # constraint: n_par = (C-1) + C + (I-1) = 2C + I - 2 = 8 + 6 - 2 = 12
  expected_npar <- 2 * n_classes + n_items - 2
  expect_equal(fit$n_par, expected_npar)
})

test_that("fit_lcr theta values are ordered", {
  set.seed(202)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  fit <- fit_lcr(data, n_classes = 4, n_starts = 5)

  # Classes should be ordered by theta
  expect_true(all(diff(fit$theta) >= 0))
})

test_that("lcr_scores returns valid scores", {
  set.seed(303)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  fit <- fit_lcr(data, n_classes = 4, n_starts = 3)

  # EAP scores
  scores_eap <- lcr_scores(fit, type = "eap")
  expect_equal(length(scores_eap), n)
  expect_true(all(scores_eap >= min(fit$theta) & scores_eap <= max(fit$theta)))

  # Modal scores
  scores_modal <- lcr_scores(fit, type = "modal")
  expect_equal(length(scores_modal), n)
  expect_true(all(scores_modal %in% fit$theta))
})

test_that("fit_lcr recovers Rasch-structured data", {
  set.seed(404)

  # Generate true Rasch data
  n <- 500
  n_items <- 8

  true_theta <- rnorm(n)
  true_delta <- seq(-1.5, 1.5, length.out = n_items)

  data <- matrix(0, n, n_items)
  for (i in 1:n) {
    for (j in 1:n_items) {
      prob <- 1 / (1 + exp(-(true_theta[i] - true_delta[j])))
      data[i, j] <- rbinom(1, 1, prob)
    }
  }

  fit <- fit_lcr(data, n_classes = 7, n_starts = 10, seed = 505)

  # Delta should correlate well with true_delta
  cor_delta <- cor(fit$delta, true_delta - mean(true_delta))
  expect_true(abs(cor_delta) > 0.9)
})

test_that("compute_rasch_probs works correctly", {
  theta <- c(-1, 0, 1)
  delta <- c(-0.5, 0, 0.5)

  probs <- compute_rasch_probs(theta, delta)

  expect_equal(dim(probs), c(3, 3))

  # Check specific values
  # P(correct | theta=0, delta=0) should be 0.5
  expect_equal(probs[2, 2], 0.5)

  # When theta > delta, prob should be > 0.5
  expect_true(probs[1, 3] > 0.5)  # theta=1, delta=-0.5
  expect_true(probs[3, 1] < 0.5)  # theta=-1, delta=0.5
})

test_that("fit_lcr has fewer parameters than UN with same n_classes", {
  set.seed(606)
  n <- 200
  n_items <- 8
  n_classes <- 4
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit_un_result <- fit_un(data, n_classes = n_classes, n_starts = 3)
  fit_lcr_result <- fit_lcr(data, n_classes = n_classes, n_starts = 3)

  # LCR should have fewer parameters
  # UN: (C-1) + I*C = 3 + 32 = 35
  # LCR: 2C + I - 2 = 8 + 8 - 2 = 14
  expect_true(fit_lcr_result$n_par < fit_un_result$n_par)
})

test_that("em_lcr convergence improves likelihood", {
  set.seed(707)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  init <- init_lcr_params(data, 3)

  result <- em_lcr(data, 3, init$theta, init$delta, init$class_probs,
                   max_iter = 100, verbose = FALSE)

  # Log-likelihood should be non-decreasing
  ll_diff <- diff(result$ll_history)
  expect_true(all(ll_diff >= -1e-6))
})

test_that("init_lcr_params produces valid initialization", {
  set.seed(808)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  init <- init_lcr_params(data, n_classes = 4)

  expect_equal(length(init$theta), 4)
  expect_equal(length(init$delta), 5)
  expect_equal(length(init$class_probs), 4)
  expect_equal(sum(init$class_probs), 1)

  # Delta should be centered
  expect_true(abs(mean(init$delta)) < 0.01)
})
