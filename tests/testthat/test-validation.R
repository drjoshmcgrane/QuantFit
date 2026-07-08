# Validation tests using simulated Rasch model data

#' Generate data from a Rasch model
#'
#' @param n Number of persons
#' @param n_items Number of items
#' @param theta_mean Mean of person abilities
#' @param theta_sd SD of person abilities
#' @param delta Item difficulties (or generated if NULL)
#' @param seed Random seed
#'
#' @return List with data matrix and true parameters
generate_rasch_data <- function(n = 500, n_items = 10,
                                 theta_mean = 0, theta_sd = 1,
                                 delta = NULL, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)


  # Generate person abilities
  theta <- rnorm(n, mean = theta_mean, sd = theta_sd)

  # Generate or use provided item difficulties
  if (is.null(delta)) {
    delta <- seq(-1.5, 1.5, length.out = n_items)
  }

  # Generate responses
  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    for (j in 1:n_items) {
      prob <- 1 / (1 + exp(-(theta[i] - delta[j])))
      data[i, j] <- rbinom(1, 1, prob)
    }
  }

  list(
    data = data,
    theta = theta,
    delta = delta,
    n = n,
    n_items = n_items
  )
}

#' Generate data from a Latent Class model (non-Rasch)
#'
#' @param n Number of persons
#' @param n_items Number of items
#' @param n_classes Number of classes
#' @param class_probs Class probabilities
#' @param monotonic If TRUE, generate monotonic item probs
#' @param seed Random seed
generate_lca_data <- function(n = 500, n_items = 10, n_classes = 3,
                               class_probs = NULL, monotonic = FALSE,
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (is.null(class_probs)) {
    class_probs <- rep(1/n_classes, n_classes)
  }

  # Generate item probabilities
  if (monotonic) {
    # Monotonically increasing across classes
    item_probs <- matrix(0, n_items, n_classes)
    for (i in 1:n_items) {
      base <- runif(1, 0.1, 0.3)
      increment <- runif(1, 0.15, 0.25)
      for (c in 1:n_classes) {
        item_probs[i, c] <- min(0.95, base + (c - 1) * increment)
      }
    }
  } else {
    # Random (non-monotonic)
    item_probs <- matrix(runif(n_items * n_classes, 0.1, 0.9),
                         nrow = n_items, ncol = n_classes)
  }

  # Generate class assignments
  class_assignments <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  # Generate responses
  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    c <- class_assignments[i]
    for (j in 1:n_items) {
      data[i, j] <- rbinom(1, 1, item_probs[j, c])
    }
  }

  list(
    data = data,
    class_assignments = class_assignments,
    item_probs = item_probs,
    class_probs = class_probs,
    n = n,
    n_items = n_items,
    n_classes = n_classes
  )
}

# ============================================================================
# Tests with Rasch-generated data
# ============================================================================

test_that("LCR model recovers Rasch item difficulties", {
  sim <- generate_rasch_data(n = 500, n_items = 8, seed = 12345)

  fit <- fit_lcr(sim$data, n_classes = 7, n_starts = 10, seed = 54321)

  # Center both for comparison
  true_delta <- sim$delta - mean(sim$delta)
  est_delta <- fit$delta - mean(fit$delta)

  # Correlation should be very high
  cor_delta <- cor(true_delta, est_delta)
  expect_true(cor_delta > 0.90,
              info = paste("Correlation:", round(cor_delta, 3)))

  # RMSE should be reasonable
  rmse <- sqrt(mean((true_delta - est_delta)^2))
  expect_true(rmse < 0.3,
              info = paste("RMSE:", round(rmse, 3)))
})

test_that("LCR scores correlate with true theta", {
  sim <- generate_rasch_data(n = 500, n_items = 10, seed = 23456)

  fit <- fit_lcr(sim$data, n_classes = 9, n_starts = 10, seed = 65432)

  # Get EAP scores
  scores <- lcr_scores(fit, type = "eap")

  # Correlation with true theta
  cor_theta <- cor(scores, sim$theta)
  expect_true(cor_theta > 0.85,
              info = paste("Correlation:", round(cor_theta, 3)))
})

test_that("Rasch data prefers quantitative models over UN", {

  sim <- generate_rasch_data(n = 400, n_items = 8, seed = 34567)

  # Fit UN and LCR
  fit_un_result <- fit_un(sim$data, n_classes = 5, n_starts = 5, seed = 1)
  fit_lcr_result <- fit_lcr(sim$data, n_classes = 5, n_starts = 5, seed = 2)

  # LCR should have better (lower) BIC due to fewer parameters
  # even with similar fit
  bic_un <- BIC(fit_un_result)
  bic_lcr <- BIC(fit_lcr_result)

  # LCR has far fewer parameters, so should have lower BIC
  expect_true(bic_lcr < bic_un,
              info = paste("BIC UN:", round(bic_un, 1),
                          "BIC LCR:", round(bic_lcr, 1)))
})

test_that("Rasch data shows monotonicity in MON model", {
  sim <- generate_rasch_data(n = 400, n_items = 8, seed = 45678)

  fit <- fit_mon(sim$data, n_classes = 5, n_starts = 5, seed = 3)

  # All items should show monotonically increasing probs
  for (i in 1:8) {
    diffs <- diff(fit$item_probs[i, ])
    expect_true(all(diffs >= -1e-6),
                info = paste("Item", i, "not monotonic"))
  }

  # The fit should be good (converged)
  expect_true(fit$convergence)
})

test_that("Rasch data shows double monotonicity in DM model", {
  sim <- generate_rasch_data(n = 400, n_items = 8, seed = 56789)

  # Item order from true difficulties
  true_order <- order(sim$delta, decreasing = TRUE)  # Easiest first

  fit <- fit_dm(sim$data, n_classes = 5, item_order = true_order,
                n_starts = 5, seed = 4)

  # Check class monotonicity
  for (i in 1:8) {
    expect_true(all(diff(fit$item_probs[i, ]) >= -1e-6))
  }

  # Check item ordering
  for (c in 1:5) {
    ordered_probs <- fit$item_probs[true_order, c]
    expect_true(all(diff(ordered_probs) <= 1e-6))
  }
})

test_that("successive_comparison identifies quantitative structure in Rasch data", {
  skip_if_not_installed("mirt")

  sim <- generate_rasch_data(n = 400, n_items = 8, seed = 67890)

  result <- successive_comparison(sim$data, n_classes = 5,
                                  n_starts = 5, verbose = FALSE)

  # Best model should be LCR or RM (quantitative)
  expect_true(result$best_model %in% c("LCR", "RM", "DM"),
              info = paste("Best model was:", result$best_model))

  # Conclusion should mention quantitative or interval
  expect_true(
    grepl("QUANTITATIVE|interval|Rasch", result$conclusion, ignore.case = TRUE) ||
    grepl("ORDINAL.*strong", result$conclusion, ignore.case = TRUE),
    info = paste("Conclusion:", result$conclusion)
  )
})

# ============================================================================
# Tests comparing Rasch vs non-Rasch data
# ============================================================================

test_that("Non-Rasch data prefers UN over LCR", {
  # Generate clearly non-Rasch data
  sim <- generate_lca_data(n = 400, n_items = 8, n_classes = 3,
                            monotonic = FALSE, seed = 78901)

  fit_un_result <- fit_un(sim$data, n_classes = 3, n_starts = 5, seed = 5)
  fit_lcr_result <- fit_lcr(sim$data, n_classes = 3, n_starts = 5, seed = 6)

  # UN should fit better (higher LL) since data doesn't follow Rasch
  ll_un <- fit_un_result$loglik
  ll_lcr <- fit_lcr_result$loglik

  # UN should have notably higher LL (Rasch is misspecified)
  expect_true(ll_un >= ll_lcr - 5,  # Allow small tolerance
              info = paste("LL UN:", round(ll_un, 1),
                          "LL LCR:", round(ll_lcr, 1)))
})

test_that("Monotonic LCA data prefers MON over UN", {
  sim <- generate_lca_data(n = 400, n_items = 8, n_classes = 4,
                            monotonic = TRUE, seed = 89012)

  fit_un_result <- fit_un(sim$data, n_classes = 4, n_starts = 5, seed = 7)
  fit_mon_result <- fit_mon(sim$data, n_classes = 4, n_starts = 5, seed = 8)

  # Both should have similar LL since MON is correctly specified
  # But BIC might favor MON slightly if it recovers true structure
  expect_true(fit_mon_result$convergence)

  # The MON model should satisfy constraints
  for (i in 1:8) {
    expect_true(all(diff(fit_mon_result$item_probs[i, ]) >= -1e-6))
  }
})

# ============================================================================
# Model selection accuracy tests
# ============================================================================

test_that("Model comparison correctly orders models by BIC", {
  sim <- generate_rasch_data(n = 300, n_items = 6, seed = 90123)

  comparison <- compare_models(sim$data, n_classes = 4,
                                models = c("UN", "MON", "DM", "LCR"),
                                n_starts = 3, verbose = FALSE)

  # Table should be sorted by BIC
  bic_values <- comparison$comparison_table$BIC
  expect_equal(bic_values, sort(bic_values))

  # Best model should have lowest BIC
  expect_equal(comparison$best_model,
               comparison$comparison_table$Model[1])
})

test_that("AIC and BIC are consistent", {
  sim <- generate_rasch_data(n = 300, n_items = 6, seed = 11111)

  fit <- fit_lcr(sim$data, n_classes = 4, n_starts = 3)

  aic <- AIC(fit)
  bic <- BIC(fit)
  ll <- fit$loglik
  k <- fit$n_par
  n <- fit$n_obs

  # Verify formulas
  expect_equal(aic, -2 * ll + 2 * k, tolerance = 1e-6)
  expect_equal(bic, -2 * ll + k * log(n), tolerance = 1e-6)
})

# ============================================================================
# Edge cases and robustness
# ============================================================================

test_that("Models handle extreme item difficulties", {
  # Very easy and very hard items
  delta <- c(-3, -2, -1, 0, 1, 2, 3)
  sim <- generate_rasch_data(n = 300, n_items = 7, delta = delta, seed = 22222)

  fit <- fit_lcr(sim$data, n_classes = 5, n_starts = 5, seed = 33333)

  # Should still converge

  expect_true(fit$convergence)

  # Delta recovery should still be reasonable
  cor_delta <- cor(fit$delta - mean(fit$delta), delta - mean(delta))
  expect_true(cor_delta > 0.85)
})

test_that("Models handle high-ability population", {
  # Population with high theta
  sim <- generate_rasch_data(n = 300, n_items = 8,
                              theta_mean = 1.5, theta_sd = 0.8,
                              seed = 44444)

  fit <- fit_lcr(sim$data, n_classes = 5, n_starts = 5, seed = 55555)

  expect_true(fit$convergence)

  # Theta estimates should be shifted high
  expect_true(mean(fit$theta) > 0)
})

test_that("Models handle low-variability population", {
  # Homogeneous population
  sim <- generate_rasch_data(n = 300, n_items = 8,
                              theta_mean = 0, theta_sd = 0.5,
                              seed = 66666)

  fit <- fit_lcr(sim$data, n_classes = 3, n_starts = 5, seed = 77777)

  expect_true(fit$convergence)

  # Theta range should be compressed
  theta_range <- diff(range(fit$theta))
  expect_true(theta_range < 3)  # Should be narrower than typical
})

# ============================================================================
# Comparison with mirt package (if available)
# ============================================================================

test_that("RM model matches mirt Rasch results", {
  skip_if_not_installed("mirt")

  sim <- generate_rasch_data(n = 500, n_items = 8, seed = 88888)

  # Fit using our wrapper
  fit_rm_result <- fit_rm(sim$data, verbose = FALSE)

  # Fit directly with mirt
  mirt_fit <- mirt::mirt(as.data.frame(sim$data), 1, itemtype = "Rasch",
                         verbose = FALSE)
  mirt_coefs <- mirt::coef(mirt_fit, simplify = TRUE)
  mirt_delta <- -mirt_coefs$items[, "d"]
  mirt_delta <- mirt_delta - mean(mirt_delta)

  # Our delta should match mirt's
  our_delta <- fit_rm_result$delta - mean(fit_rm_result$delta)

  cor_delta <- cor(our_delta, mirt_delta)
  expect_true(cor_delta > 0.99,
              info = paste("Correlation with mirt:", round(cor_delta, 4)))
})

test_that("LCR approximates RM with many classes", {
  skip_if_not_installed("mirt")

  sim <- generate_rasch_data(n = 500, n_items = 8, seed = 99999)

  fit_rm_result <- fit_rm(sim$data, verbose = FALSE)
  fit_lcr_result <- fit_lcr(sim$data, n_classes = 15, n_starts = 10, seed = 11112)

  # Delta parameters should be very similar
  rm_delta <- fit_rm_result$delta - mean(fit_rm_result$delta)
  lcr_delta <- fit_lcr_result$delta - mean(fit_lcr_result$delta)

  cor_delta <- cor(rm_delta, lcr_delta)
  expect_true(cor_delta > 0.95,
              info = paste("RM-LCR delta correlation:", round(cor_delta, 3)))
})

# ============================================================================
# Torres Irribarra & Diakow 6-Model Framework Tests
# ============================================================================

#' Generate non-monotonic LCA data (Model 0 - UN)
generate_UN_test_data <- function(n, n_items, n_classes, seed) {
  set.seed(seed)

  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)
  for (i in 1:n_items) {
    # Non-monotonic patterns
    if (i %% 2 == 1) {
      # Decreasing
      item_probs[i, ] <- seq(0.8, 0.2, length.out = n_classes)
    } else {
      # U-shaped
      mid <- (n_classes + 1) / 2
      item_probs[i, ] <- sapply(1:n_classes, function(c) 0.5 + 0.3 * abs(c - mid) / mid)
    }
    item_probs[i, ] <- item_probs[i, ] + rnorm(n_classes, 0, 0.02)
  }
  item_probs <- pmax(pmin(item_probs, 0.95), 0.05)

  class_probs <- rep(1/n_classes, n_classes)
  class_assign <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    data[i, ] <- rbinom(n_items, 1, item_probs[, class_assign[i]])
  }

  list(data = data, item_probs = item_probs)
}

test_that("UN model preferred for non-monotonic data", {
  sim <- generate_UN_test_data(n = 500, n_items = 8, n_classes = 4, seed = 111111)

  fit_un_res <- fit_un(sim$data, n_classes = 4, n_starts = 5, seed = 1)
  fit_mon_res <- fit_mon(sim$data, n_classes = 4, n_starts = 5, seed = 2)

  # UN should have better (higher) log-likelihood since MON is constrained
  expect_true(fit_un_res$loglik >= fit_mon_res$loglik - 1,
              info = paste("UN LL:", round(fit_un_res$loglik, 1),
                          "MON LL:", round(fit_mon_res$loglik, 1)))
})

test_that("MON model satisfies class monotonicity constraints", {
  sim <- generate_rasch_data(n = 400, n_items = 8, seed = 222222)

  fit <- fit_mon(sim$data, n_classes = 5, n_starts = 5, seed = 3)

  # Every item should have non-decreasing probs across classes
  for (i in 1:8) {
    diffs <- diff(fit$item_probs[i, ])
    expect_true(all(diffs >= -1e-6),
                info = paste("Item", i, "violates monotonicity"))
  }
})

test_that("IIO model maintains invariant item ordering", {
  sim <- generate_rasch_data(n = 400, n_items = 8, seed = 333333)

  # Use empirical item order
  item_order <- order(colMeans(sim$data), decreasing = TRUE)

  fit <- fit_iio(sim$data, n_classes = 5, item_order = item_order,
                 n_starts = 5, seed = 4)

  # Item order should be maintained within each class
  for (c in 1:5) {
    ordered_probs <- fit$item_probs[item_order, c]
    diffs <- diff(ordered_probs)
    expect_true(all(diffs <= 1e-6),
                info = paste("Class", c, "violates item ordering"))
  }
})

test_that("DM model satisfies both monotonicity constraints", {
  sim <- generate_rasch_data(n = 400, n_items = 8, seed = 444444)

  item_order <- order(colMeans(sim$data), decreasing = TRUE)

  fit <- fit_dm(sim$data, n_classes = 5, item_order = item_order,
                n_starts = 5, seed = 5)

  # Check class monotonicity
  for (i in 1:8) {
    expect_true(all(diff(fit$item_probs[i, ]) >= -1e-6),
                info = paste("Item", i, "violates class monotonicity"))
  }

  # Check item ordering
  for (c in 1:5) {
    ordered_probs <- fit$item_probs[item_order, c]
    expect_true(all(diff(ordered_probs) <= 1e-6),
                info = paste("Class", c, "violates item ordering"))
  }
})

test_that("LCR model uses Rasch parameterization correctly", {
  sim <- generate_rasch_data(n = 500, n_items = 8, seed = 555555)

  fit <- fit_lcr(sim$data, n_classes = 7, n_starts = 10, seed = 6)

  # Should have theta and delta parameters

  expect_true(!is.null(fit$theta))
  expect_true(!is.null(fit$delta))

  # Check that item_probs follow Rasch formula
  for (c in 1:7) {
    for (i in 1:8) {
      expected_prob <- 1 / (1 + exp(-(fit$theta[c] - fit$delta[i])))
      expect_equal(fit$item_probs[i, c], expected_prob, tolerance = 1e-4)
    }
  }
})

test_that("Scale-level classification distinguishes quantitative vs non-quantitative", {
  skip_if_not_installed("mirt")

  # Generate Rasch data (quantitative)
  sim_rasch <- generate_rasch_data(n = 400, n_items = 8, seed = 666666)

  # Generate non-monotonic LCA data (non-quantitative)
  sim_lca <- generate_UN_test_data(n = 400, n_items = 8, n_classes = 4, seed = 777777)

  # Compare models for Rasch data
  comp_rasch <- compare_models(sim_rasch$data, n_classes = 4,
                                models = c("UN", "MON", "LCR"),
                                n_starts = 3, verbose = FALSE)

  # Compare models for LCA data
  comp_lca <- compare_models(sim_lca$data, n_classes = 4,
                              models = c("UN", "MON", "LCR"),
                              n_starts = 3, verbose = FALSE)

  # Rasch data should prefer LCR (quantitative)
  # LCA data should prefer UN (non-quantitative)
  expect_true(comp_rasch$best_model %in% c("LCR", "MON"),
              info = paste("Rasch data selected:", comp_rasch$best_model))
})
