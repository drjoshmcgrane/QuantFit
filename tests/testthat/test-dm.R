# Tests for Double Monotonicity Model (DM)

test_that("fit_dm returns valid qlfit object", {
  set.seed(123)
  n <- 200
  n_items <- 5
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit <- fit_dm(data, n_classes = 3, n_starts = 3, verbose = FALSE)

  expect_s3_class(fit, "qlfit")
  expect_equal(fit$model_type, "DM")
  expect_equal(fit$n_obs, n)
  expect_equal(fit$n_items, n_items)
  expect_equal(fit$n_classes, 3)
})

test_that("fit_dm satisfies both monotonicity constraints", {
  set.seed(456)
  n <- 300
  n_items <- 6
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit <- fit_dm(data, n_classes = 4, n_starts = 5)

  # Check class monotonicity (each row non-decreasing)
  for (i in 1:n_items) {
    diffs <- diff(fit$item_probs[i, ])
    expect_true(all(diffs >= -1e-6),
                info = paste("Item", i, "violates class monotonicity"))
  }

  # Check item ordering (within each class, ordered items decreasing)
  item_order <- fit$item_order
  for (c in 1:4) {
    ordered_probs <- fit$item_probs[item_order, c]
    diffs <- diff(ordered_probs)
    expect_true(all(diffs <= 1e-6),
                info = paste("Class", c, "violates item ordering"))
  }
})

test_that("fit_dm returns item_order", {
  set.seed(789)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  fit <- fit_dm(data, n_classes = 3, n_starts = 3)

  expect_true(!is.null(fit$item_order))
  expect_equal(length(fit$item_order), 5)
  expect_equal(sort(fit$item_order), 1:5)
})

test_that("fit_dm has both constraints in specification", {
  set.seed(101)
  n <- 100
  data <- matrix(rbinom(n * 4, 1, 0.5), nrow = n)

  fit <- fit_dm(data, n_classes = 2, n_starts = 2)

  expect_true(!is.null(fit$constraints))
  expect_s3_class(fit$constraints, "ql_constraints")
  expect_true(fit$constraints$class_monotonicity)
  expect_true(fit$constraints$item_ordering)
})

test_that("fit_dm projection method works", {
  set.seed(202)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  fit <- fit_dm_projection(data, n_classes = 3, n_starts = 3)

  expect_s3_class(fit, "qlfit")
  expect_equal(fit$model_type, "DM")

  # Verify constraints
  for (i in 1:5) {
    expect_true(all(diff(fit$item_probs[i, ]) >= -1e-6))
  }
})

test_that("project_dm_constraints works correctly", {
  set.seed(303)

  # Create random probabilities that violate constraints
  probs <- matrix(runif(12), nrow = 4, ncol = 3)
  item_order <- c(1, 3, 2, 4)

  projected <- project_dm_constraints(probs, item_order)

  # Check class monotonicity
  for (i in 1:4) {
    expect_true(all(diff(projected[i, ]) >= -1e-6))
  }

  # Check item ordering
  for (c in 1:3) {
    ordered_probs <- projected[item_order, c]
    expect_true(all(diff(ordered_probs) <= 1e-6))
  }

  # Check bounded
  expect_true(all(projected > 0 & projected < 1))
})

test_that("build_double_monotonicity_fn combines constraints", {
  item_order <- c(1, 2, 3)
  constraint_fn <- build_double_monotonicity_fn(3, 3, item_order)

  # Test with doubly monotonic probabilities
  good_probs <- matrix(c(
    0.7, 0.8, 0.9,  # Item 1 (easiest)
    0.4, 0.5, 0.6,  # Item 2
    0.1, 0.2, 0.3   # Item 3 (hardest)
  ), nrow = 3, byrow = TRUE)

  constraints <- constraint_fn(good_probs)
  expect_true(all(constraints >= -1e-6))

  # Number of constraints: 3*(3-1) for class mon + 3*(3-1) for item ord = 12
  expect_equal(length(constraints), 12)
})

test_that("fit_dm log-likelihood is not better than fit_mon or fit_iio", {
  set.seed(404)
  n <- 200
  n_items <- 5
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  item_order <- estimate_item_order(data)

  fit_dm_result <- fit_dm(data, n_classes = 3, item_order = item_order, n_starts = 5)
  fit_mon_result <- fit_mon(data, n_classes = 3, n_starts = 5)
  fit_iio_result <- fit_iio(data, n_classes = 3, item_order = item_order, n_starts = 5)

  # DM is more constrained, so LL should be <= both MON and IIO
  expect_true(fit_dm_result$loglik <= fit_mon_result$loglik + 1e-2)
  expect_true(fit_dm_result$loglik <= fit_iio_result$loglik + 1e-2)
})

test_that("fit_dm converges on doubly monotonic data", {
  set.seed(505)

  # Generate truly doubly monotonic data
  n <- 300
  n_items <- 5
  n_classes <- 3

  true_class_probs <- c(0.3, 0.4, 0.3)

  # Item order: 1 is easiest, 5 is hardest
  # Class order: 1 is lowest, 3 is highest
  true_item_probs <- matrix(c(
    # Class 1, Class 2, Class 3
    0.5, 0.7, 0.9,   # Item 1 (easiest)
    0.4, 0.6, 0.8,   # Item 2
    0.3, 0.5, 0.7,   # Item 3
    0.2, 0.4, 0.6,   # Item 4
    0.1, 0.3, 0.5    # Item 5 (hardest)
  ), nrow = 5, byrow = TRUE)

  # Generate data
  class_assignments <- sample(1:n_classes, n, replace = TRUE, prob = true_class_probs)
  data <- matrix(0, n, n_items)
  for (i in 1:n) {
    c <- class_assignments[i]
    for (j in 1:n_items) {
      data[i, j] <- rbinom(1, 1, true_item_probs[j, c])
    }
  }

  fit <- fit_dm(data, n_classes = 3, item_order = 1:5, n_starts = 10, seed = 606)

  expect_true(fit$convergence)
  expect_true(fit$loglik > -1000)
})

test_that("init_item_probs_dm produces valid doubly monotonic initialization", {
  set.seed(707)
  n <- 200
  n_items <- 5
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  item_order <- 1:5
  init_probs <- init_item_probs_dm(data, n_classes = 3, item_order)

  expect_equal(dim(init_probs), c(5, 3))
  expect_true(all(init_probs > 0 & init_probs < 1))

  # Check class monotonicity
  for (i in 1:5) {
    expect_true(all(diff(init_probs[i, ]) >= -1e-6))
  }

  # Check item ordering
  for (c in 1:3) {
    ordered_probs <- init_probs[item_order, c]
    expect_true(all(diff(ordered_probs) <= 1e-6))
  }
})
