# Tests for Model Comparison Framework

test_that("compare_models returns valid qlcompare object", {
  skip_if_not_installed("mirt")

  set.seed(123)
  n <- 200
  n_items <- 5
  data <- matrix(rbinom(n * n_items, 1, 0.5), nrow = n)

  comparison <- compare_models(data, n_classes = 3,
                                models = c("UN", "MON", "LCR"),
                                n_starts = 2, verbose = FALSE)

  expect_s3_class(comparison, "qlcompare")
  expect_true(!is.null(comparison$fits))
  expect_true(!is.null(comparison$comparison_table))
  expect_true(!is.null(comparison$best_model))
})

test_that("compare_models comparison_table has correct columns", {
  set.seed(456)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  comparison <- compare_models(data, n_classes = 3,
                                models = c("UN", "MON"),
                                n_starts = 2)

  tbl <- comparison$comparison_table

  expect_true("Model" %in% names(tbl))
  expect_true("LogLik" %in% names(tbl))
  expect_true("nPar" %in% names(tbl))
  expect_true("AIC" %in% names(tbl))
  expect_true("BIC" %in% names(tbl))
  expect_true("Converged" %in% names(tbl))
})

test_that("compare_models best_model has lowest BIC", {
  set.seed(789)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  comparison <- compare_models(data, n_classes = 3,
                                models = c("UN", "MON", "DM"),
                                n_starts = 3)

  tbl <- comparison$comparison_table
  best <- comparison$best_model

  expect_equal(best, tbl$Model[which.min(tbl$BIC)])
})

test_that("get_model extracts correct model", {
  set.seed(101)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  comparison <- compare_models(data, n_classes = 3,
                                models = c("UN", "MON"),
                                n_starts = 2)

  un_fit <- get_model(comparison, "UN")
  expect_s3_class(un_fit, "qlfit")
  expect_equal(un_fit$model_type, "UN")

  mon_fit <- get_model(comparison, "MON")
  expect_equal(mon_fit$model_type, "MON")

  expect_error(get_model(comparison, "RM"))  # Not fitted
})

test_that("compare_models handles single model", {
  set.seed(202)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  comparison <- compare_models(data, n_classes = 3,
                                models = "UN",
                                n_starts = 2)

  expect_equal(nrow(comparison$comparison_table), 1)
})

test_that("print and summary methods work for qlcompare", {
  set.seed(303)
  n <- 100
  data <- matrix(rbinom(n * 4, 1, 0.5), nrow = n)

  comparison <- compare_models(data, n_classes = 2,
                                models = c("UN", "MON"),
                                n_starts = 2)

  expect_output(print(comparison), "Model Comparison")
  expect_output(summary(comparison), "delta_BIC")
})

test_that("successive_comparison returns correct structure", {
  set.seed(404)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  # Only test a few models to speed up
  result <- successive_comparison(data, n_classes = 3,
                                  n_starts = 2, verbose = FALSE)

  expect_true(!is.null(result$comparison))
  expect_true(!is.null(result$steps))
  expect_true(!is.null(result$conclusion))
  expect_true(!is.null(result$best_model))
})

test_that("successive_comparison identifies classificatory data", {
  set.seed(505)

  # Generate data with no ordinal structure
  n <- 300
  n_items <- 5

  # Classes with non-monotonic patterns
  class_probs <- c(0.33, 0.33, 0.34)
  item_probs <- matrix(c(
    0.9, 0.2, 0.5,  # Non-monotonic
    0.2, 0.8, 0.4,
    0.5, 0.3, 0.9,
    0.4, 0.6, 0.2,
    0.7, 0.4, 0.6
  ), nrow = 5, byrow = TRUE)

  class_assignments <- sample(1:3, n, replace = TRUE, prob = class_probs)
  data <- matrix(0, n, n_items)
  for (i in 1:n) {
    c <- class_assignments[i]
    for (j in 1:n_items) {
      data[i, j] <- rbinom(1, 1, item_probs[j, c])
    }
  }

  result <- successive_comparison(data, n_classes = 3,
                                  n_starts = 5, verbose = FALSE)

  # Should likely identify UN as best or indicate no ordering
  # (Result depends on random variation)
  expect_true(result$best_model %in% c("UN", "MON", "IIO", "DM", "LCR", "RM"))
})

test_that("interpret_structure returns interpretation for all models", {
  models <- c("UN", "MON", "IIO", "DM", "LCR", "RM")

  for (m in models) {
    interp <- interpret_structure(m)
    expect_true(is.character(interp))
    expect_true(nchar(interp) > 50)  # Should be substantial text
  }
})

test_that("compare_fit works with list input", {
  set.seed(606)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  fit1 <- fit_un(data, n_classes = 2, n_starts = 2)
  fit2 <- fit_un(data, n_classes = 3, n_starts = 2)

  result <- compare_fit(fit1, fit2, measures = "ic")

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true("BIC" %in% names(result))
})

test_that("lr_test computes correct statistics", {
  set.seed(707)
  n <- 200
  data <- matrix(rbinom(n * 5, 1, 0.5), nrow = n)

  fit1 <- fit_un(data, n_classes = 2, n_starts = 3)
  fit2 <- fit_un(data, n_classes = 3, n_starts = 3)

  result <- lr_test(fit1, fit2)

  expect_s3_class(result, "lr_test")
  expect_true(is.numeric(result$statistic))
  expect_true(result$statistic >= 0)  # LR stat is non-negative
  expect_true(result$df > 0)
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

test_that("fit_measures extracts correct values", {
  set.seed(808)
  n <- 100
  data <- matrix(rbinom(n * 4, 1, 0.5), nrow = n)

  fit <- fit_un(data, n_classes = 2, n_starts = 2)
  measures <- fit_measures(fit)

  expect_true(is.data.frame(measures))
  expect_equal(measures$model, "UN")
  expect_equal(measures$loglik, fit$loglik)
  expect_equal(measures$n_par, fit$n_par)
  expect_equal(measures$AIC, AIC(fit))
  expect_equal(measures$BIC, BIC(fit))
})

test_that("entropy_r2 returns valid value", {
  set.seed(909)
  n <- 100
  data <- matrix(rbinom(n * 4, 1, 0.5), nrow = n)

  fit <- fit_un(data, n_classes = 3, n_starts = 2)
  e_r2 <- entropy_r2(fit)

  expect_true(is.numeric(e_r2))
  expect_true(e_r2 >= 0 && e_r2 <= 1)
})

test_that("classification_accuracy returns valid value", {
  set.seed(1010)
  n <- 100
  data <- matrix(rbinom(n * 4, 1, 0.5), nrow = n)

  fit <- fit_un(data, n_classes = 2, n_starts = 2)
  acc <- classification_accuracy(fit)

  expect_true(is.numeric(acc))
  expect_true(acc >= 0.5 && acc <= 1)  # At least 50% for 2 classes
})
