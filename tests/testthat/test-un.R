# Tests for Unconstrained Latent Class Model (UN)

test_that("fit_un returns valid qlfit object", {
  set.seed(123)
  # Generate simple LCA data
  n <- 200
  n_items <- 6
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit <- fit_un(data, n_classes = 2, n_starts = 3, verbose = FALSE)

  expect_s3_class(fit, "qlfit")
  expect_equal(fit$model_type, "UN")
  expect_equal(fit$n_obs, n)
  expect_equal(fit$n_items, n_items)
  expect_equal(fit$n_classes, 2)
})

test_that("fit_un item probabilities are bounded", {
  set.seed(456)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  fit <- fit_un(data, n_classes = 3, n_starts = 3)

  expect_true(all(fit$item_probs > 0))
  expect_true(all(fit$item_probs < 1))
})

test_that("fit_un class probabilities sum to 1", {
  set.seed(789)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  fit <- fit_un(data, n_classes = 3, n_starts = 3)

  expect_equal(sum(fit$class_probs), 1, tolerance = 1e-6)
})

test_that("fit_un posteriors are valid", {
  set.seed(101)
  n <- 100
  data <- matrix(rbinom(n * 4, 1, 0.5), nrow = n)

  fit <- fit_un(data, n_classes = 2, n_starts = 2)

  # Check dimensions

  expect_equal(dim(fit$posteriors), c(n, 2))

  # Check rows sum to 1
  row_sums <- rowSums(fit$posteriors)
  expect_true(all(abs(row_sums - 1) < 1e-6))

  # Check all values in [0, 1]
  expect_true(all(fit$posteriors >= 0 & fit$posteriors <= 1))
})

test_that("fit_un recovers known parameters", {
  set.seed(202)

  # Generate data from known 2-class model
  n <- 500
  n_items <- 5

  true_class_probs <- c(0.4, 0.6)
  true_item_probs <- matrix(c(
    0.2, 0.8,  # Item 1
    0.3, 0.7,  # Item 2
    0.4, 0.6,  # Item 3
    0.25, 0.75,  # Item 4
    0.35, 0.65   # Item 5
  ), nrow = 5, byrow = TRUE)

  # Generate data
  class_assignments <- sample(1:2, n, replace = TRUE, prob = true_class_probs)
  data <- matrix(0, n, n_items)
  for (i in 1:n) {
    c <- class_assignments[i]
    for (j in 1:n_items) {
      data[i, j] <- rbinom(1, 1, true_item_probs[j, c])
    }
  }

  fit <- fit_un(data, n_classes = 2, n_starts = 10, seed = 303)

  # Check that recovery is reasonable (allowing for label switching)
  # Either fit matches true or is reversed
  if (fit$class_probs[1] < fit$class_probs[2]) {
    # Classes might be in same order as true
    expect_true(
      cor(as.vector(fit$item_probs), as.vector(true_item_probs)) > 0.8 ||
      cor(as.vector(fit$item_probs[, 2:1]), as.vector(true_item_probs)) > 0.8
    )
  }
})

test_that("fit_un log-likelihood increases or stays same across iterations", {
  set.seed(404)
  n <- 100
  data <- matrix(rbinom(n * 4, 1, 0.5), nrow = n)

  # Get single run to examine history
  init_probs <- init_item_probs(data, 2, "random")
  init_class <- init_class_probs(2, "uniform")

  result <- em_lca(data, 2, init_probs, init_class, max_iter = 50)

  # Log-likelihood should be non-decreasing
  ll_diff <- diff(result$ll_history)
  expect_true(all(ll_diff >= -1e-6))  # Allow tiny numerical errors
})

test_that("fit_un S3 methods work", {
  set.seed(505)
  n <- 100
  data <- matrix(rbinom(n * 4, 1, 0.5), nrow = n)

  fit <- fit_un(data, n_classes = 2, n_starts = 2)

  # print method
  expect_output(print(fit), "Unconstrained")

  # coef method
  coefs <- coef(fit)
  expect_true("class_probs" %in% names(coefs))
  expect_true("item_probs" %in% names(coefs))

  # logLik method
  ll <- logLik(fit)
  expect_true(is.finite(ll))
  expect_equal(attr(ll, "df"), fit$n_par)

  # AIC/BIC methods
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
  expect_true(BIC(fit) > AIC(fit))  # BIC penalizes more for n > e^2
})

test_that("fit_un handles edge cases", {
  set.seed(606)

  # Minimum viable data
  data <- matrix(c(0, 0, 1, 1, 0, 1, 1, 0), nrow = 4)
  expect_error(fit_un(data, n_classes = 5))  # Too many classes

  # All same responses
  data_same <- matrix(1, nrow = 20, ncol = 3)
  # This should still run, though results may be degenerate
  fit_same <- fit_un(data_same, n_classes = 2, n_starts = 2)
  expect_s3_class(fit_same, "qlfit")
})
