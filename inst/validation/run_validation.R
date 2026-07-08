# ============================================================
# QuantFit Package Validation
# Simulates data from each model and tests model selection
# ============================================================

# Load package functions
source("R/utils.R")
source("R/classes.R")
source("R/em_algorithm.R")
source("R/constraints.R")
source("R/fit_measures.R")
source("R/fit_un.R")
source("R/fit_mon.R")
source("R/fit_iio.R")
source("R/fit_dm.R")
source("R/fit_lcr.R")
source("R/fit_rm.R")

# ============================================================
# DATA GENERATION FUNCTIONS
# ============================================================

#' Generate data from Unconstrained LCA (no ordering structure)
#' Creates data that STRONGLY VIOLATES class monotonicity.
#' Half the items have REVERSED monotonicity, others have non-monotonic patterns.
#' This represents truly qualitative/nominal latent structure.
#' @param n Sample size
#' @param n_items Number of items
#' @param n_classes Number of latent classes
#' @param seed Random seed
#' @return List with data matrix and true parameters
generate_UN_data <- function(n, n_items, n_classes, seed) {
  set.seed(seed)

  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)

  # Create STRONGLY unordered structure:
  # - Half the items have REVERSED monotonicity (prob DECREASES with class)
  # - Other items have clearly non-monotonic patterns
  # - This ensures constrained models will fit much worse

  for (i in 1:n_items) {
    if (i <= n_items / 2) {
      # REVERSED monotonicity: clearly decreasing across classes
      # These items STRONGLY violate Rasch/monotonicity assumptions
      base <- 0.85
      decrement <- 0.6 / (n_classes - 1)  # Large decrease
      for (c in 1:n_classes) {
        item_probs[i, c] <- base - (c - 1) * decrement
      }
    } else if (i <= 3 * n_items / 4) {
      # Non-monotonic: peaked in middle class
      for (c in 1:n_classes) {
        mid <- (n_classes + 1) / 2
        dist_from_mid <- abs(c - mid) / mid
        item_probs[i, c] <- 0.8 - 0.5 * dist_from_mid
      }
    } else {
      # Non-monotonic: valley in middle class
      for (c in 1:n_classes) {
        mid <- (n_classes + 1) / 2
        dist_from_mid <- abs(c - mid) / mid
        item_probs[i, c] <- 0.2 + 0.5 * dist_from_mid
      }
    }

    # Add small noise
    item_probs[i, ] <- item_probs[i, ] + rnorm(n_classes, 0, 0.03)
  }

  # Bound probabilities (ensure matrix structure is preserved)
  item_probs[item_probs < 0.10] <- 0.10
  item_probs[item_probs > 0.90] <- 0.90

  class_probs <- rep(1/n_classes, n_classes)
  class_assign <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  # Generate responses
  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    for (j in 1:n_items) {
      data[i, j] <- rbinom(1, 1, item_probs[j, class_assign[i]])
    }
  }

  list(
    data = data,
    item_probs = item_probs,
    class_probs = class_probs,
    class_assign = class_assign,
    model = "UN"
  )
}

#' Generate data from Class Monotonicity model (MON)
#' Items have monotonically increasing probs across classes (ordinal structure)
#' BUT items have EXTREMELY DIFFERENT discriminations - strongly violating
#' Rasch's equal discrimination assumption.
#' Item ordering CHANGES across classes (items cross each other).
#' @param n Sample size
#' @param n_items Number of items
#' @param n_classes Number of latent classes
#' @param seed Random seed
#' @return List with data matrix and true parameters
generate_MON_data <- function(n, n_items, n_classes, seed) {
  set.seed(seed)

  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)

  for (i in 1:n_items) {
    # Each item has VERY DIFFERENT discrimination
    # This strongly violates Rasch's equal discrimination assumption
    # AND causes items to cross (violating invariant item ordering)

    if (i %% 3 == 1) {
      # Very high discrimination: nearly 0 to nearly 1
      item_probs[i, 1] <- 0.08
      item_probs[i, n_classes] <- 0.92
      # Interpolate
      for (c in 2:(n_classes-1)) {
        item_probs[i, c] <- 0.08 + 0.84 * (c - 1) / (n_classes - 1)
      }
    } else if (i %% 3 == 2) {
      # Very low discrimination: stays around 0.5
      for (c in 1:n_classes) {
        item_probs[i, c] <- 0.40 + 0.15 * (c - 1) / (n_classes - 1)
      }
    } else {
      # Medium but offset threshold (starts high, increases little)
      for (c in 1:n_classes) {
        item_probs[i, c] <- 0.55 + 0.30 * (c - 1) / (n_classes - 1)
      }
    }

    # Add noise
    item_probs[i, ] <- item_probs[i, ] + rnorm(n_classes, 0, 0.02)
  }

  # Bound probabilities
  item_probs[item_probs < 0.05] <- 0.05
  item_probs[item_probs > 0.95] <- 0.95

  # Enforce monotonicity
  for (i in 1:n_items) {
    item_probs[i, ] <- sort(item_probs[i, ])
  }

  class_probs <- rep(1/n_classes, n_classes)
  class_assign <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    for (j in 1:n_items) {
      data[i, j] <- rbinom(1, 1, item_probs[j, class_assign[i]])
    }
  }

  list(
    data = data,
    item_probs = item_probs,
    class_probs = class_probs,
    class_assign = class_assign,
    model = "MON"
  )
}

#' Generate data from Invariant Item Ordering model (IIO)
#' Items maintain same relative ordering across ALL classes.
#' BUT classes are NOT monotonically ordered - class ordering is SCRAMBLED.
#' This is ordinal item structure without ordinal class structure.
#' @param n Sample size
#' @param n_items Number of items
#' @param n_classes Number of latent classes
#' @param seed Random seed
#' @return List with data matrix and true parameters
generate_IIO_data <- function(n, n_items, n_classes, seed) {
  set.seed(seed)

  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)

  # Create structure where:
  # - Item ordering is INVARIANT (item 1 always easiest, item n always hardest)
  # - BUT class abilities are SCRAMBLED (not monotonic)

  # Item difficulties (invariant across classes)
  item_difficulty <- seq(0, 1, length.out = n_items)

  # Non-monotonic class abilities: scrambled order
  # For 4 classes, use order like: Class 3 > Class 1 > Class 4 > Class 2
  class_abilities <- c(0.55, 0.25, 0.75, 0.40)  # Scrambled, not monotonic
  if (n_classes != 4) {
    # For other class numbers, create random scrambled order
    class_abilities <- sample(seq(0.2, 0.8, length.out = n_classes))
  }

  for (c in 1:n_classes) {
    for (i in 1:n_items) {
      # Item prob = class ability - item difficulty effect
      # This preserves item ordering within each class
      item_probs[i, c] <- class_abilities[c] + 0.15 - 0.35 * item_difficulty[i]
    }
  }

  # Add noise
  item_probs <- item_probs + matrix(rnorm(n_items * n_classes, 0, 0.03),
                                     nrow = n_items, ncol = n_classes)

  # Bound probabilities
  item_probs[item_probs < 0.08] <- 0.08
  item_probs[item_probs > 0.92] <- 0.92

  # Ensure item ordering within each class (PAVA decreasing)
  for (c in 1:n_classes) {
    item_probs[, c] <- sort(item_probs[, c], decreasing = TRUE)
  }

  class_probs <- rep(1/n_classes, n_classes)
  class_assign <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    for (j in 1:n_items) {
      data[i, j] <- rbinom(1, 1, item_probs[j, class_assign[i]])
    }
  }

  list(
    data = data,
    item_probs = item_probs,
    class_probs = class_probs,
    class_assign = class_assign,
    item_order = 1:n_items,
    model = "IIO"
  )
}

#' Generate data from Double Monotonicity model (DM)
#' Both class monotonicity AND invariant item ordering hold.
#' Uses LINEAR spacing with EQUAL discrimination (Mokken-style).
#' Different from Rasch because the probability function is linear, not logistic.
#' @param n Sample size
#' @param n_items Number of items
#' @param n_classes Number of latent classes
#' @param seed Random seed
#' @return List with data matrix and true parameters
generate_DM_data <- function(n, n_items, n_classes, seed) {
  set.seed(seed)

  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)

  # Create LINEAR (not logistic) probability structure with EQUAL discrimination
  # This satisfies DM and is Mokken-like but differs from Rasch in functional form
  for (i in 1:n_items) {
    # Item difficulty (linear spacing)
    item_difficulty <- (i - 1) / (n_items - 1)  # 0 to 1

    for (c in 1:n_classes) {
      # Class ability (linear spacing, monotonically increasing)
      class_ability <- (c - 1) / (n_classes - 1)  # 0 to 1

      # LINEAR probability function (not logistic)
      # P = base + discrimination * (ability - difficulty)
      # With equal discrimination for all items
      item_probs[i, c] <- 0.10 + 0.80 * class_ability - 0.50 * item_difficulty
    }
  }

  # Add small noise
  item_probs <- item_probs + matrix(rnorm(n_items * n_classes, 0, 0.02),
                                     nrow = n_items, ncol = n_classes)

  # Bound probabilities
  item_probs[item_probs < 0.05] <- 0.05
  item_probs[item_probs > 0.95] <- 0.95

  # Project to ensure DM constraints hold
  for (i in 1:n_items) item_probs[i, ] <- sort(item_probs[i, ])
  for (c in 1:n_classes) item_probs[, c] <- sort(item_probs[, c], decreasing = TRUE)

  class_probs <- rep(1/n_classes, n_classes)
  class_assign <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    for (j in 1:n_items) {
      data[i, j] <- rbinom(1, 1, item_probs[j, class_assign[i]])
    }
  }

  list(
    data = data,
    item_probs = item_probs,
    class_probs = class_probs,
    class_assign = class_assign,
    item_order = 1:n_items,
    model = "DM"
  )
}

#' Generate data from Latent Class Rasch model (LCR)
#' Rasch parameterization with discrete ability classes
#' @param n Sample size
#' @param n_items Number of items
#' @param n_classes Number of latent classes
#' @param seed Random seed
#' @return List with data matrix and true parameters
generate_LCR_data <- function(n, n_items, n_classes, seed) {
  set.seed(seed)

  # Discrete ability levels
  theta <- seq(-1.5, 1.5, length.out = n_classes)

  # Item difficulties
  delta <- seq(-1.5, 1.5, length.out = n_items)

  class_probs <- rep(1/n_classes, n_classes)
  class_assign <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  # Compute item probabilities from Rasch formula
  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)
  for (j in 1:n_items) {
    for (c in 1:n_classes) {
      item_probs[j, c] <- 1 / (1 + exp(-(theta[c] - delta[j])))
    }
  }

  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    for (j in 1:n_items) {
      prob <- 1 / (1 + exp(-(theta[class_assign[i]] - delta[j])))
      data[i, j] <- rbinom(1, 1, prob)
    }
  }

  list(
    data = data,
    theta = theta,
    delta = delta,
    item_probs = item_probs,
    class_probs = class_probs,
    class_assign = class_assign,
    model = "LCR"
  )
}

#' Generate data from continuous Rasch Model (RM)
#' @param n Sample size
#' @param n_items Number of items
#' @param theta_mean Mean of ability distribution
#' @param theta_sd SD of ability distribution
#' @param seed Random seed
#' @return List with data matrix and true parameters
generate_RM_data <- function(n, n_items, theta_mean = 0, theta_sd = 1, seed) {
  set.seed(seed)

  # Continuous ability
  theta <- rnorm(n, mean = theta_mean, sd = theta_sd)

  # Item difficulties
  delta <- seq(-2, 2, length.out = n_items)

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
    model = "RM"
  )
}

# ============================================================
# DIAGNOSTIC FUNCTION
# ============================================================

#' Check the structure of simulated data
#' Verifies that data generation is producing expected patterns
#' @param sim_data Simulated data object from generate_*_data()
#' @return Prints diagnostics
check_data_structure <- function(sim_data) {
  cat("\n--- Data Structure Check for", sim_data$model, "---\n")

  if (!is.null(sim_data$item_probs)) {
    probs <- sim_data$item_probs
    n_items <- nrow(probs)
    n_classes <- ncol(probs)

    cat("Item probability matrix (rows=items, cols=classes):\n")
    print(round(probs, 2))

    # Check class monotonicity (each row should be non-decreasing)
    mon_violations <- 0
    for (i in 1:n_items) {
      if (any(diff(probs[i, ]) < -0.01)) {
        mon_violations <- mon_violations + 1
      }
    }
    cat("\nClass monotonicity violations:", mon_violations, "of", n_items, "items\n")

    # Check if items cross (item ordering changes across classes)
    item_order_class1 <- order(probs[, 1], decreasing = TRUE)
    item_order_classC <- order(probs[, n_classes], decreasing = TRUE)
    orders_match <- all(item_order_class1 == item_order_classC)
    cat("Item ordering invariant:", orders_match, "\n")

    # Discrimination variability
    discriminations <- apply(probs, 1, function(x) max(x) - min(x))
    cat("Discrimination range: [", round(min(discriminations), 2), ",",
        round(max(discriminations), 2), "]\n")
    cat("Discrimination SD:", round(sd(discriminations), 3), "\n")
  }

  cat("\n")
}

# ============================================================
# VALIDATION FUNCTION
# ============================================================

#' Run full model comparison on a dataset
#' Fits all 6 models from Torres Irribarra & Diakow framework
#' @param data Binary response matrix
#' @param n_classes Number of classes for discrete models
#' @param n_starts Number of random starts
#' @return Data frame with fit results
fit_all_models <- function(data, n_classes, n_starts = 5) {

  res_un <- tryCatch(
    suppressWarnings(fit_un(data, n_classes, n_starts = n_starts, seed = 1)),
    error = function(e) NULL
  )

  res_mon <- tryCatch(
    suppressWarnings(fit_mon(data, n_classes, n_starts = n_starts, seed = 2)),
    error = function(e) NULL
  )

  res_iio <- tryCatch(
    suppressWarnings(fit_iio(data, n_classes, n_starts = n_starts, seed = 3)),
    error = function(e) NULL
  )

  res_dm <- tryCatch(
    suppressWarnings(fit_dm(data, n_classes, n_starts = n_starts, seed = 4)),
    error = function(e) NULL
  )

  res_lcr <- tryCatch(
    suppressWarnings(fit_lcr(data, n_classes, n_starts = n_starts, seed = 5)),
    error = function(e) NULL
  )

  res_rm <- tryCatch(
    suppressWarnings(fit_rm(data, verbose = FALSE)),
    error = function(e) NULL
  )

  data.frame(
    LL_UN = if (!is.null(res_un)) res_un$loglik else NA,
    LL_MON = if (!is.null(res_mon)) res_mon$loglik else NA,
    LL_IIO = if (!is.null(res_iio)) res_iio$loglik else NA,
    LL_DM = if (!is.null(res_dm)) res_dm$loglik else NA,
    LL_LCR = if (!is.null(res_lcr)) res_lcr$loglik else NA,
    LL_RM = if (!is.null(res_rm)) res_rm$loglik else NA,
    nPar_UN = if (!is.null(res_un)) res_un$n_par else NA,
    nPar_MON = if (!is.null(res_mon)) res_mon$n_par else NA,
    nPar_IIO = if (!is.null(res_iio)) res_iio$n_par else NA,
    nPar_DM = if (!is.null(res_dm)) res_dm$n_par else NA,
    nPar_LCR = if (!is.null(res_lcr)) res_lcr$n_par else NA,
    nPar_RM = if (!is.null(res_rm)) res_rm$n_par else NA,
    BIC_UN = if (!is.null(res_un)) BIC(res_un) else NA,
    BIC_MON = if (!is.null(res_mon)) BIC(res_mon) else NA,
    BIC_IIO = if (!is.null(res_iio)) BIC(res_iio) else NA,
    BIC_DM = if (!is.null(res_dm)) BIC(res_dm) else NA,
    BIC_LCR = if (!is.null(res_lcr)) BIC(res_lcr) else NA,
    BIC_RM = if (!is.null(res_rm)) BIC(res_rm) else NA
  )
}

# ============================================================
# RUN VALIDATION
# ============================================================

run_validation <- function(n = 400, n_items = 10, n_classes = 4, n_reps = 5,
                           save_results = TRUE, output_dir = "inst/validation") {

  cat("QuantFit Validation Study\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("Settings: n =", n, ", items =", n_items, ", classes =", n_classes,
      ", reps =", n_reps, "\n\n")

  all_results <- data.frame()
  all_simulated_data <- list()

  # All 6 models from Torres Irribarra & Diakow
  models_to_test <- c("UN", "MON", "IIO", "DM", "LCR", "RM")

  for (gen_model in models_to_test) {
    cat("\nGenerating from:", gen_model, "\n")
    cat(rep("-", 40), "\n", sep = "")

    for (rep in 1:n_reps) {
      seed <- match(gen_model, models_to_test) * 1000 + rep * 100

      # Generate data
      sim_data <- switch(gen_model,
        UN = generate_UN_data(n, n_items, n_classes, seed),
        MON = generate_MON_data(n, n_items, n_classes, seed),
        IIO = generate_IIO_data(n, n_items, n_classes, seed),
        DM = generate_DM_data(n, n_items, n_classes, seed),
        LCR = generate_LCR_data(n, n_items, n_classes, seed),
        RM = generate_RM_data(n, n_items, seed = seed)
      )

      # Store simulated data
      data_id <- paste0(gen_model, "_rep", rep)
      all_simulated_data[[data_id]] <- sim_data

      # Fit all models
      fit_results <- fit_all_models(sim_data$data, n_classes)

      # Determine best model
      bics <- c(
        UN = fit_results$BIC_UN,
        MON = fit_results$BIC_MON,
        IIO = fit_results$BIC_IIO,
        DM = fit_results$BIC_DM,
        LCR = fit_results$BIC_LCR,
        RM = fit_results$BIC_RM
      )
      best_model <- names(which.min(bics))

      # Combine results
      result_row <- data.frame(
        gen_model = gen_model,
        rep = rep,
        seed = seed,
        best_model = best_model,
        correct = (best_model == gen_model),
        fit_results
      )
      all_results <- rbind(all_results, result_row)

      # Print progress
      match_marker <- ifelse(best_model == gen_model, " <-- MATCH", "")
      cat(sprintf("  Rep %d: Best=%-3s (UN=%.0f MON=%.0f IIO=%.0f DM=%.0f LCR=%.0f RM=%.0f)%s\n",
          rep, best_model,
          fit_results$BIC_UN, fit_results$BIC_MON, fit_results$BIC_IIO,
          fit_results$BIC_DM, fit_results$BIC_LCR, fit_results$BIC_RM,
          match_marker))
    }
  }

  # Save results
  if (save_results) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    saveRDS(all_results, file.path(output_dir, "validation_results.rds"))
    saveRDS(all_simulated_data, file.path(output_dir, "simulated_data.rds"))

    cat("\n\nResults saved to:", output_dir, "\n")
  }

  # Print summary
  cat("\n\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("VALIDATION SUMMARY\n")
  cat(rep("=", 60), "\n", sep = "")

  cat("\n1. Selection Frequency (rows = generating, cols = selected):\n")
  print(table(all_results$gen_model, all_results$best_model))

  cat("\n2. Accuracy by Generating Model:\n")
  for (gm in models_to_test) {
    subset <- all_results[all_results$gen_model == gm, ]
    n_correct <- sum(subset$best_model == gm)
    cat(sprintf("   %-4s: %d/%d (%.0f%%)\n", gm, n_correct, n_reps, 100 * n_correct / n_reps))
  }

  total_correct <- sum(all_results$correct)
  total_n <- nrow(all_results)
  cat(sprintf("\n   TOTAL: %d/%d (%.0f%%)\n", total_correct, total_n, 100 * total_correct / total_n))

  cat("\n3. Mean BIC by Generating Model:\n")
  cat(sprintf("   %-6s %7s %7s %7s %7s %7s %7s\n", "Gen", "UN", "MON", "IIO", "DM", "LCR", "RM"))
  for (gm in models_to_test) {
    subset <- all_results[all_results$gen_model == gm, ]
    cat(sprintf("   %-6s %7.0f %7.0f %7.0f %7.0f %7.0f %7.0f\n", gm,
        mean(subset$BIC_UN, na.rm = TRUE),
        mean(subset$BIC_MON, na.rm = TRUE),
        mean(subset$BIC_IIO, na.rm = TRUE),
        mean(subset$BIC_DM, na.rm = TRUE),
        mean(subset$BIC_LCR, na.rm = TRUE),
        mean(subset$BIC_RM, na.rm = TRUE)))
  }

  cat("\n4. Mean Log-Likelihood by Generating Model:\n")
  cat(sprintf("   %-6s %7s %7s %7s %7s %7s %7s\n", "Gen", "UN", "MON", "IIO", "DM", "LCR", "RM"))
  for (gm in models_to_test) {
    subset <- all_results[all_results$gen_model == gm, ]
    cat(sprintf("   %-6s %7.0f %7.0f %7.0f %7.0f %7.0f %7.0f\n", gm,
        mean(subset$LL_UN, na.rm = TRUE),
        mean(subset$LL_MON, na.rm = TRUE),
        mean(subset$LL_IIO, na.rm = TRUE),
        mean(subset$LL_DM, na.rm = TRUE),
        mean(subset$LL_LCR, na.rm = TRUE),
        mean(subset$LL_RM, na.rm = TRUE)))
  }

  invisible(list(
    results = all_results,
    simulated_data = all_simulated_data
  ))
}

# ============================================================
# DIAGNOSTIC TEST
# ============================================================

#' Quick diagnostic test of data generation
#' Generates one dataset from each model and checks structure
diagnostic_test <- function(n = 400, n_items = 10, n_classes = 4) {
  cat("=" |> rep(60) |> paste(collapse = ""), "\n")
  cat("DIAGNOSTIC TEST: Data Generation Verification\n")
  cat("=" |> rep(60) |> paste(collapse = ""), "\n")

  cat("\nGenerating one dataset from each model type...\n")

  # All 6 models from Torres Irribarra & Diakow
  models <- c("UN", "MON", "IIO", "DM", "LCR", "RM")

  for (m in models) {
    sim <- switch(m,
      UN = generate_UN_data(n, n_items, n_classes, seed = 12345),
      MON = generate_MON_data(n, n_items, n_classes, seed = 12345),
      IIO = generate_IIO_data(n, n_items, n_classes, seed = 12345),
      DM = generate_DM_data(n, n_items, n_classes, seed = 12345),
      LCR = generate_LCR_data(n, n_items, n_classes, seed = 12345),
      RM = generate_RM_data(n, n_items, seed = 12345)
    )

    check_data_structure(sim)
  }

  cat("\n" |> rep(2) |> paste(collapse = ""))
  cat("KEY EXPECTATIONS:\n")
  cat("-" |> rep(40) |> paste(collapse = ""), "\n")
  cat("UN:  Many monotonicity violations (>50%), items cross\n")
  cat("MON: 0 monotonicity violations, high discrim variability, items cross\n")
  cat("IIO: Monotonicity violations OK, invariant item ordering\n")
  cat("DM:  0 monotonicity violations, invariant item ordering\n")
  cat("LCR: 0 monotonicity violations, Rasch structure (logistic)\n")
  cat("RM:  Continuous theta, Rasch structure (logistic)\n")
}

# Run if executed directly
if (interactive() || !exists("VALIDATION_SOURCED")) {
  VALIDATION_SOURCED <- TRUE
  cat("To run validation, call: run_validation()\n")
  cat("For diagnostics, call: diagnostic_test()\n")
  cat("Or with custom settings: run_validation(n = 500, n_reps = 10)\n")
}
