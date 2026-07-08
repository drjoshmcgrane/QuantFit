# Tests for Class Monotonicity Model (MON)

test_that("fit_mon returns valid qlfit object", {
  set.seed(123)
  n <- 200
  n_items <- 5
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit <- fit_mon(data, n_classes = 3, n_starts = 3, verbose = FALSE)

  expect_s3_class(fit, "qlfit")
  expect_equal(fit$model_type, "MON")
  expect_equal(fit$n_obs, n)
  expect_equal(fit$n_items, n_items)
  expect_equal(fit$n_classes, 3)
})

test_that("fit_mon satisfies monotonicity constraints", {
  set.seed(456)
  n <- 300
  n_items <- 6
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit <- fit_mon(data, n_classes = 4, n_starts = 5)

  # Check that each row (item) is non-decreasing across classes
  for (i in 1:n_items) {
    diffs <- diff(fit$item_probs[i, ])
    expect_true(all(diffs >= -1e-6),
                info = paste("Item", i, "violates monotonicity"))
  }
})

test_that("fit_mon projection method works", {
  set.seed(789)
  n <- 200
  n_items <- 5
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit <- fit_mon_projection(data, n_classes = 3, n_starts = 3)

  expect_s3_class(fit, "qlfit")
  expect_equal(fit$model_type, "MON")

  # Check monotonicity
  for (i in 1:n_items) {
    expect_true(all(diff(fit$item_probs[i, ]) >= -1e-6))
  }
})

test_that("fit_mon has constraints attribute", {
  set.seed(101)
  n <- 100
  data <- matrix(rbinom(n * 4, 1, 0.5), nrow = n)

  fit <- fit_mon(data, n_classes = 2, n_starts = 2)

  expect_true(!is.null(fit$constraints))
  expect_s3_class(fit$constraints, "ql_constraints")
  expect_true(fit$constraints$class_monotonicity)
  expect_false(fit$constraints$item_ordering)
})

test_that("fit_mon log-likelihood is not better than fit_un", {
  set.seed(202)
  n <- 200
  n_items <- 5
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit_unconstrained <- fit_un(data, n_classes = 3, n_starts = 5, seed = 1)
  fit_monotonic <- fit_mon(data, n_classes = 3, n_starts = 5, seed = 1)

  # MON should have same or lower log-likelihood (more constrained)
  expect_true(fit_monotonic$loglik <= fit_unconstrained$loglik + 1e-3)
})

test_that("fit_mon converges on monotonic data", {
  set.seed(303)

  # Generate data with true monotonic structure
  n <- 300
  n_items <- 5
  n_classes <- 3

  # True class probs
  true_class_probs <- c(0.3, 0.4, 0.3)

  # True item probs (monotonically increasing across classes)
  true_item_probs <- matrix(c(
    0.1, 0.4, 0.7,
    0.2, 0.5, 0.8,
    0.15, 0.45, 0.75,
    0.25, 0.55, 0.85,
    0.3, 0.6, 0.9
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

  fit <- fit_mon(data, n_classes = 3, n_starts = 10, seed = 404)

  expect_true(fit$convergence)

  # Should have reasonable log-likelihood
  expect_true(fit$loglik > -1000)
})

test_that("build_class_monotonicity_fn works correctly", {
  constraint_fn <- build_class_monotonicity_fn(3, 4)  # 3 items, 4 classes

  # Test with monotonic probabilities (should all be >= 0)
  monotonic_probs <- matrix(c(
    0.1, 0.3, 0.5, 0.7,
    0.2, 0.4, 0.6, 0.8,
    0.15, 0.35, 0.55, 0.75
  ), nrow = 3, byrow = TRUE)

  constraints <- constraint_fn(monotonic_probs)
  expect_true(all(constraints >= 0))

  # Test with non-monotonic probabilities (should have some < 0)
  non_monotonic <- matrix(c(
    0.5, 0.3, 0.4, 0.7,  # Violation at position 2
    0.2, 0.4, 0.6, 0.8,
    0.15, 0.35, 0.55, 0.75
  ), nrow = 3, byrow = TRUE)

  constraints_bad <- constraint_fn(non_monotonic)
  expect_true(any(constraints_bad < 0))
})

test_that("check_constraints detects violations", {
  constraints <- specify_constraints(class_monotonicity = TRUE)

  # Monotonic case
  good_probs <- matrix(c(0.2, 0.5, 0.8, 0.3, 0.6, 0.9), nrow = 2, byrow = TRUE)
  result_good <- check_constraints(good_probs, constraints)
  expect_true(result_good$satisfied)
  expect_equal(result_good$n_violations, 0)

  # Non-monotonic case
  bad_probs <- matrix(c(0.8, 0.5, 0.2, 0.3, 0.6, 0.9), nrow = 2, byrow = TRUE)
  result_bad <- check_constraints(bad_probs, constraints)
  expect_false(result_bad$satisfied)
  expect_true(result_bad$n_violations > 0)
})
