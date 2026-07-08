# Tests for Invariant Item Ordering Model (IIO)

test_that("fit_iio returns valid qlfit object", {
  set.seed(123)
  n <- 200
  n_items <- 5
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit <- fit_iio(data, n_classes = 3, n_starts = 3, verbose = FALSE)

  expect_s3_class(fit, "qlfit")
  expect_equal(fit$model_type, "IIO")
  expect_equal(fit$n_obs, n)
  expect_equal(fit$n_items, n_items)
  expect_equal(fit$n_classes, 3)
})

test_that("fit_iio returns item_order", {
  set.seed(456)
  n <- 200
  data <- matrix(rbinom(n * 6, 1, 0.5), nrow = n)

  fit <- fit_iio(data, n_classes = 3, n_starts = 3)

  expect_true(!is.null(fit$item_order))
  expect_equal(length(fit$item_order), 6)
  expect_equal(sort(fit$item_order), 1:6)
})

test_that("fit_iio respects provided item_order", {
  set.seed(789)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  custom_order <- c(5, 3, 1, 2, 4)
  fit <- fit_iio(data, n_classes = 3, item_order = custom_order, n_starts = 3)

  expect_equal(fit$item_order, custom_order)
})

test_that("fit_iio satisfies item ordering constraints", {
  set.seed(101)
  n <- 300
  n_items <- 6
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  fit <- fit_iio(data, n_classes = 4, n_starts = 5)

  # Check that within each class, items follow the ordering
  # (easier items have higher probs)
  item_order <- fit$item_order

  for (c in 1:4) {
    ordered_probs <- fit$item_probs[item_order, c]
    diffs <- diff(ordered_probs)
    expect_true(all(diffs <= 1e-6),
                info = paste("Class", c, "violates item ordering"))
  }
})

test_that("fit_iio projection method works", {
  set.seed(202)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  fit <- fit_iio_projection(data, n_classes = 3, n_starts = 3)

  expect_s3_class(fit, "qlfit")
  expect_equal(fit$model_type, "IIO")
})

test_that("fit_iio has constraints attribute", {
  set.seed(303)
  n <- 100
  data <- matrix(rbinom(n * 4, 1, 0.5), nrow = n)

  fit <- fit_iio(data, n_classes = 2, n_starts = 2)

  expect_true(!is.null(fit$constraints))
  expect_s3_class(fit$constraints, "ql_constraints")
  expect_false(fit$constraints$class_monotonicity)
  expect_true(fit$constraints$item_ordering)
})

test_that("build_item_ordering_fn works correctly", {
  item_order <- c(1, 3, 2)  # Item 1 easiest, then 3, then 2
  constraint_fn <- build_item_ordering_fn(3, 2, item_order)

  # Test with proper ordering (should all be >= 0)
  good_probs <- matrix(c(
    0.8, 0.9,  # Item 1 (easiest)
    0.4, 0.5,  # Item 2 (hardest)
    0.6, 0.7   # Item 3 (middle)
  ), nrow = 3, byrow = TRUE)

  constraints <- constraint_fn(good_probs)
  expect_true(all(constraints >= 0))

  # Test with violated ordering (should have some < 0)
  bad_probs <- matrix(c(
    0.3, 0.4,  # Item 1 supposed to be easiest but low
    0.4, 0.5,  # Item 2
    0.6, 0.7   # Item 3
  ), nrow = 3, byrow = TRUE)

  constraints_bad <- constraint_fn(bad_probs)
  expect_true(any(constraints_bad < 0))
})

test_that("fit_iio validates item_order input", {
  set.seed(404)
  n <- 100
  data <- matrix(rbinom(n * 4, 1, 0.5), nrow = n)

  # Wrong length
  expect_error(fit_iio(data, n_classes = 2, item_order = c(1, 2, 3)))

  # Invalid values
  expect_error(fit_iio(data, n_classes = 2, item_order = c(1, 2, 3, 5)))
})

test_that("init_item_probs_iio produces valid initialization", {
  set.seed(505)
  n <- 200
  n_items <- 5
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  item_order <- estimate_item_order(data)
  init_probs <- init_item_probs_iio(data, n_classes = 3, item_order)

  expect_equal(dim(init_probs), c(5, 3))
  expect_true(all(init_probs > 0 & init_probs < 1))

  # Check item ordering within each class
  for (c in 1:3) {
    ordered_probs <- init_probs[item_order, c]
    expect_true(all(diff(ordered_probs) <= 0))
  }
})
