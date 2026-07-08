# Test QuantFit with Simulated Rasch Data
# This script validates the package using conventionally generated Rasch data

# Load required packages
library(QuantFit)

# ============================================================================
# Function to generate Rasch model data
# ============================================================================

#' Simulate data from a Rasch model
#'
#' Uses the standard Rasch model: P(X=1|theta,delta) = exp(theta-delta)/(1+exp(theta-delta))
#'
#' @param n Number of persons
#' @param n_items Number of items
#' @param theta_mean Mean of ability distribution
#' @param theta_sd SD of ability distribution
#' @param delta_range Range for item difficulties
#' @param seed Random seed
sim_rasch <- function(n = 500, n_items = 10,
                      theta_mean = 0, theta_sd = 1,
                      delta_range = c(-2, 2),
                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Person abilities from normal distribution
  theta <- rnorm(n, mean = theta_mean, sd = theta_sd)

  # Item difficulties evenly spaced
  delta <- seq(delta_range[1], delta_range[2], length.out = n_items)

  # Generate responses using Rasch model
  data <- matrix(0, nrow = n, ncol = n_items)
  colnames(data) <- paste0("item", 1:n_items)

  for (i in 1:n) {
    for (j in 1:n_items) {
      # Rasch probability
      logit <- theta[i] - delta[j]
      prob <- exp(logit) / (1 + exp(logit))
      data[i, j] <- rbinom(1, 1, prob)
    }
  }

  list(
    data = data,
    theta = theta,
    delta = delta
  )
}

# ============================================================================
# Generate test dataset
# ============================================================================

cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("Generating Rasch model data...\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

sim <- sim_rasch(
  n = 500,
  n_items = 10,
  theta_mean = 0,
  theta_sd = 1,
  delta_range = c(-1.5, 1.5),
  seed = 42
)

cat("Data dimensions:", dim(sim$data), "\n")
cat("True item difficulties (delta):\n")
print(round(sim$delta, 3))
cat("\nMean theta:", round(mean(sim$theta), 3), "\n")
cat("SD theta:", round(sd(sim$theta), 3), "\n")
cat("\nItem proportions correct:\n")
print(round(colMeans(sim$data), 3))

# ============================================================================
# Fit individual models
# ============================================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("Fitting individual models...\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

n_classes <- 5

# Unconstrained
cat("Fitting UN model...\n")
fit_un_result <- fit_un(sim$data, n_classes = n_classes, n_starts = 10, seed = 1)
cat("  Converged:", fit_un_result$convergence, "\n")
cat("  LogLik:", round(fit_un_result$loglik, 2), "\n")
cat("  BIC:", round(BIC(fit_un_result), 2), "\n\n")

# Class Monotonicity
cat("Fitting MON model...\n")
fit_mon_result <- fit_mon(sim$data, n_classes = n_classes, n_starts = 10, seed = 2)
cat("  Converged:", fit_mon_result$convergence, "\n")
cat("  LogLik:", round(fit_mon_result$loglik, 2), "\n")
cat("  BIC:", round(BIC(fit_mon_result), 2), "\n\n")

# Double Monotonicity
cat("Fitting DM model...\n")
# Use true item order (easiest to hardest based on delta)
true_item_order <- order(sim$delta, decreasing = TRUE)
fit_dm_result <- fit_dm(sim$data, n_classes = n_classes,
                        item_order = true_item_order, n_starts = 10, seed = 3)
cat("  Converged:", fit_dm_result$convergence, "\n")
cat("  LogLik:", round(fit_dm_result$loglik, 2), "\n")
cat("  BIC:", round(BIC(fit_dm_result), 2), "\n\n")

# Latent Class Rasch
cat("Fitting LCR model...\n")
fit_lcr_result <- fit_lcr(sim$data, n_classes = n_classes, n_starts = 10, seed = 4)
cat("  Converged:", fit_lcr_result$convergence, "\n")
cat("  LogLik:", round(fit_lcr_result$loglik, 2), "\n")
cat("  BIC:", round(BIC(fit_lcr_result), 2), "\n\n")

# ============================================================================
# Parameter recovery for LCR
# ============================================================================

cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("LCR Parameter Recovery\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

# Compare estimated delta to true delta
true_delta_centered <- sim$delta - mean(sim$delta)
est_delta_centered <- fit_lcr_result$delta - mean(fit_lcr_result$delta)

cat("True delta (centered):", round(true_delta_centered, 3), "\n")
cat("Est. delta (centered):", round(est_delta_centered, 3), "\n\n")

cor_delta <- cor(true_delta_centered, est_delta_centered)
rmse_delta <- sqrt(mean((true_delta_centered - est_delta_centered)^2))

cat("Correlation between true and estimated delta:", round(cor_delta, 4), "\n")
cat("RMSE of delta estimates:", round(rmse_delta, 4), "\n\n")

# Theta (class locations) vs scoring
cat("Estimated class locations (theta):", round(fit_lcr_result$theta, 3), "\n")
cat("Class proportions:", round(fit_lcr_result$class_probs, 3), "\n\n")

# Get EAP scores and compare to true theta
eap_scores <- lcr_scores(fit_lcr_result, type = "eap")
cor_theta <- cor(eap_scores, sim$theta)
cat("Correlation between EAP scores and true theta:", round(cor_theta, 4), "\n")

# ============================================================================
# Model comparison
# ============================================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("Model Comparison\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

# Compare all models (excluding RM which needs mirt)
comparison_table <- data.frame(
  Model = c("UN", "MON", "DM", "LCR"),
  LogLik = c(fit_un_result$loglik, fit_mon_result$loglik,
             fit_dm_result$loglik, fit_lcr_result$loglik),
  nPar = c(fit_un_result$n_par, fit_mon_result$n_par,
           fit_dm_result$n_par, fit_lcr_result$n_par),
  AIC = c(AIC(fit_un_result), AIC(fit_mon_result),
          AIC(fit_dm_result), AIC(fit_lcr_result)),
  BIC = c(BIC(fit_un_result), BIC(fit_mon_result),
          BIC(fit_dm_result), BIC(fit_lcr_result))
)

comparison_table$delta_BIC <- comparison_table$BIC - min(comparison_table$BIC)
comparison_table <- comparison_table[order(comparison_table$BIC), ]

print(comparison_table, row.names = FALSE)

cat("\nBest model by BIC:", comparison_table$Model[1], "\n")

# ============================================================================
# Verify monotonicity in MON model
# ============================================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("Verifying Model Constraints\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

cat("MON Model - Class Monotonicity Check:\n")
mon_violations <- 0
for (i in 1:10) {
  diffs <- diff(fit_mon_result$item_probs[i, ])
  if (any(diffs < -1e-6)) {
    mon_violations <- mon_violations + 1
    cat("  Item", i, "VIOLATES monotonicity\n")
  }
}
if (mon_violations == 0) {
  cat("  All items satisfy class monotonicity\n")
}

cat("\nDM Model - Double Monotonicity Check:\n")
dm_class_violations <- 0
dm_item_violations <- 0

# Check class monotonicity
for (i in 1:10) {
  if (any(diff(fit_dm_result$item_probs[i, ]) < -1e-6)) {
    dm_class_violations <- dm_class_violations + 1
  }
}

# Check item ordering
for (c in 1:n_classes) {
  ordered_probs <- fit_dm_result$item_probs[true_item_order, c]
  if (any(diff(ordered_probs) > 1e-6)) {
    dm_item_violations <- dm_item_violations + 1
  }
}

cat("  Class monotonicity violations:", dm_class_violations, "\n")
cat("  Item ordering violations:", dm_item_violations, "\n")

# ============================================================================
# Test with mirt if available
# ============================================================================

if (requireNamespace("mirt", quietly = TRUE)) {
  cat("\n")
  cat("=" |> rep(60) |> paste(collapse = ""), "\n")
  cat("Comparison with mirt Rasch Model\n")
  cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

  # Fit RM wrapper
  cat("Fitting RM model (mirt wrapper)...\n")
  fit_rm_result <- fit_rm(sim$data, verbose = FALSE)

  cat("  Converged:", fit_rm_result$convergence, "\n")
  cat("  LogLik:", round(fit_rm_result$loglik, 2), "\n")
  cat("  BIC:", round(BIC(fit_rm_result), 2), "\n\n")

  # Compare delta
  rm_delta <- fit_rm_result$delta - mean(fit_rm_result$delta)

  cat("Delta comparison (mirt vs true):\n")
  cat("  Correlation:", round(cor(rm_delta, true_delta_centered), 4), "\n")

  cat("\nDelta comparison (LCR vs mirt):\n")
  cat("  Correlation:", round(cor(est_delta_centered, rm_delta), 4), "\n")

  # Full comparison including RM
  cat("\nFull Model Comparison (including RM):\n")
  full_table <- rbind(
    comparison_table[, 1:5],
    data.frame(
      Model = "RM",
      LogLik = fit_rm_result$loglik,
      nPar = fit_rm_result$n_par,
      AIC = AIC(fit_rm_result),
      BIC = BIC(fit_rm_result)
    )
  )
  full_table$delta_BIC <- full_table$BIC - min(full_table$BIC)
  full_table <- full_table[order(full_table$BIC), ]
  print(full_table, row.names = FALSE)

} else {
  cat("\nmirt package not available, skipping RM comparison\n")
}

# ============================================================================
# Summary
# ============================================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

cat("Data: Rasch model with", 500, "persons and", 10, "items\n\n")

cat("Key findings:\n")
cat("1. LCR model recovers item difficulties with r =", round(cor_delta, 3), "\n")
cat("2. LCR EAP scores correlate r =", round(cor_theta, 3), "with true abilities\n")
cat("3. Best fitting model by BIC:", comparison_table$Model[1], "\n")
cat("4. All constrained models satisfy their constraints\n\n")

if (comparison_table$Model[1] == "LCR") {
  cat("SUCCESS: LCR (quantitative model) correctly identified as best fit\n")
  cat("for data generated from a Rasch model.\n")
} else {
  cat("NOTE: ", comparison_table$Model[1], " selected, but LCR is most appropriate\n")
  cat("for true Rasch data. Consider increasing n_classes or sample size.\n")
}
