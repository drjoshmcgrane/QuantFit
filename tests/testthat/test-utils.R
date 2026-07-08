# Tests for utility functions

test_that("validate_data works correctly", {
  # Valid data
  data <- matrix(c(0, 1, 1, 0, 1, 0), nrow = 2)
  expect_silent(validate_data(data))

  # Data frame input
  df <- data.frame(x1 = c(0, 1), x2 = c(1, 0))
  result <- validate_data(df)
  expect_true(is.matrix(result))

  # Invalid values
  bad_data <- matrix(c(0, 1, 2), nrow = 1)
  expect_error(validate_data(bad_data), "binary")

  # NA handling
  na_data <- matrix(c(0, 1, NA, 1), nrow = 2)
  expect_error(validate_data(na_data), "missing")
  expect_silent(validate_data(na_data, allow_na = TRUE))
})

test_that("init_class_probs returns valid probabilities", {
  probs <- init_class_probs(4, "uniform")
  expect_equal(length(probs), 4)
  expect_equal(sum(probs), 1)
  expect_true(all(probs > 0))

  probs_random <- init_class_probs(3, "random", seed = 123)
  expect_equal(sum(probs_random), 1)
})

test_that("init_item_probs returns valid matrix", {
  set.seed(42)
  data <- matrix(rbinom(100, 1, 0.5), nrow = 20, ncol = 5)

  probs <- init_item_probs(data, n_classes = 3, method = "random")
  expect_equal(dim(probs), c(5, 3))
  expect_true(all(probs > 0 & probs < 1))

  probs_q <- init_item_probs(data, n_classes = 3, method = "quantiles")
  expect_equal(dim(probs_q), c(5, 3))
})

test_that("log_sum_exp is numerically stable", {
  # Normal case
  x <- c(1, 2, 3)
  expect_equal(log_sum_exp(x), log(sum(exp(x))))

  # Large values (would overflow without stability)
  x_large <- c(1000, 1001, 1002)
  result <- log_sum_exp(x_large)
  expect_false(is.infinite(result))
  expect_true(result > 1000)

  # Small values
  x_small <- c(-1000, -1001, -1002)
  result_small <- log_sum_exp(x_small)
  expect_false(is.infinite(result_small))
})

test_that("softmax returns valid probabilities", {
  x <- c(1, 2, 3)
  probs <- softmax(x)
  expect_equal(sum(probs), 1)
  expect_true(all(probs > 0))
  expect_true(probs[3] > probs[2] && probs[2] > probs[1])
})

test_that("inv_logit and logit are inverses", {
  p <- c(0.1, 0.5, 0.9)
  expect_equal(inv_logit(logit(p)), p, tolerance = 1e-10)

  x <- c(-2, 0, 2)
  expect_equal(logit(inv_logit(x)), x, tolerance = 1e-10)
})

test_that("bound_probs keeps values in (0, 1)", {
  p <- c(-0.1, 0, 0.5, 1, 1.1)
  bounded <- bound_probs(p)
  expect_true(all(bounded > 0 & bounded < 1))
})

test_that("pattern_frequencies computes correctly", {
  data <- matrix(c(0, 0, 1, 1, 0, 1), nrow = 3, byrow = TRUE)
  freq <- pattern_frequencies(data)

  expect_true("pattern" %in% names(freq))
  expect_true("frequency" %in% names(freq))
  expect_equal(sum(freq$frequency), 3)
})

test_that("estimate_item_order works correctly", {
  # Create data with clear ordering (item 1 easiest, item 3 hardest)
  data <- matrix(c(
    1, 1, 0,
    1, 0, 0,
    1, 1, 1,
    1, 1, 0,
    1, 0, 0
  ), nrow = 5, byrow = TRUE)

  order <- estimate_item_order(data)
  expect_equal(order[1], 1)  # Easiest (highest mean)
  expect_equal(order[3], 3)  # Hardest (lowest mean)
})

test_that("count_parameters is correct for each model", {
  n_items <- 5
  n_classes <- 3

  expect_equal(count_parameters("UN", n_items, n_classes), 2 + 15)  # (C-1) + I*C
  expect_equal(count_parameters("MON", n_items, n_classes), 2 + 15)
  expect_equal(count_parameters("IIO", n_items, n_classes), 2 + 15)
  expect_equal(count_parameters("DM", n_items, n_classes), 2 + 15)

  # LCR: (C-1) class probs + C theta (all free) + (I-1) delta (mean = 0)
  # = 2C + I - 2 = 6 + 5 - 2 = 9
  expect_equal(count_parameters("LCR", n_items, n_classes), 9)

  # RM (mirt Rasch): I intercepts + latent variance = I + 1 = 6
  expect_equal(count_parameters("RM", n_items, n_classes), 6)
})

test_that("pava_increasing produces monotonic sequence", {
  x <- c(3, 1, 4, 1, 5)
  result <- pava_increasing(x)
  expect_true(all(diff(result) >= 0))
})

test_that("pava_decreasing produces monotonic sequence", {
  x <- c(3, 1, 4, 1, 5)
  result <- pava_decreasing(x)
  expect_true(all(diff(result) <= 0))
})
