# ============================================================
# QuantFit Package Validation
# Following Torres Irribarra & Diakow Simulation Design
# ============================================================
#
# This script replicates the simulation study from the position paper:
# - 6 generating models (UN, MON, IIO, DM, LCR, RM)
# - n = 5000 persons (as in original study)
# - 10 items
# - Varying number of classes (2-6)
# - Tests whether QuantFit correctly recovers the generating model
# ============================================================

# ============================================================
# DATA GENERATION FUNCTIONS
# Matching the paper's simulation approach
# ============================================================

#' Generate data from Unconstrained LCA (Model 0)
#' No ordering constraints - purely nominal/classificatory structure
#' Classes differ qualitatively, not quantitatively
generate_model_0_UN <- function(n, n_items, n_classes, seed) {
  set.seed(seed)

  # Create item probabilities with NO systematic ordering
  # Each class has a distinct profile that doesn't follow monotonicity
  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)

  for (i in 1:n_items) {
    # Random baseline for each item
    base <- runif(1, 0.3, 0.7)

    # Generate NON-MONOTONIC class-specific deviations
    # Use different patterns for different items to ensure no systematic ordering
    if (i %% 4 == 1) {
      # Decreasing pattern
      deviations <- seq(0.25, -0.25, length.out = n_classes)
    } else if (i %% 4 == 2) {
      # U-shaped pattern
      mid <- (n_classes + 1) / 2
      deviations <- sapply(1:n_classes, function(c) 0.2 * (abs(c - mid) / mid - 0.5))
    } else if (i %% 4 == 3) {
      # Inverted U-shaped pattern
      mid <- (n_classes + 1) / 2
      deviations <- sapply(1:n_classes, function(c) -0.2 * (abs(c - mid) / mid - 0.5))
    } else {
      # Random permutation
      deviations <- sample(seq(-0.2, 0.2, length.out = n_classes))
    }

    item_probs[i, ] <- base + deviations + rnorm(n_classes, 0, 0.03)
  }

  # Bound probabilities
  item_probs <- pmax(pmin(item_probs, 0.95), 0.05)

  # Generate class memberships and responses
  class_probs <- rep(1/n_classes, n_classes)
  class_assign <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    data[i, ] <- rbinom(n_items, 1, item_probs[, class_assign[i]])
  }

  list(
    data = data,
    item_probs = item_probs,
    class_probs = class_probs,
    true_class = class_assign,
    model_code = 0,
    model_name = "UN"
  )
}

#' Generate data from Class Monotonicity model (Model 1)
#' Item probs increase monotonically across classes
#' But items have DIFFERENT discriminations (violates Rasch)
#' And item ordering is NOT invariant (items cross)
generate_model_1_MON <- function(n, n_items, n_classes, seed) {
  set.seed(seed)

  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)

  for (i in 1:n_items) {
    # Each item has different discrimination (slope)
    # This ensures items cross and IIO is violated
    discrimination <- runif(1, 0.3, 0.9)
    threshold <- runif(1, 0.2, 0.8)

    for (c in 1:n_classes) {
      class_ability <- (c - 1) / (n_classes - 1)
      item_probs[i, c] <- threshold + discrimination * (class_ability - 0.5)
    }

    # Ensure monotonicity by sorting
    item_probs[i, ] <- sort(item_probs[i, ])
  }

  # Add noise and bound
  item_probs <- item_probs + matrix(rnorm(n_items * n_classes, 0, 0.02),
                                     nrow = n_items)
  item_probs <- pmax(pmin(item_probs, 0.95), 0.05)

  # Re-enforce monotonicity after noise
  for (i in 1:n_items) {
    item_probs[i, ] <- sort(item_probs[i, ])
  }

  # Generate data
  class_probs <- rep(1/n_classes, n_classes)
  class_assign <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    data[i, ] <- rbinom(n_items, 1, item_probs[, class_assign[i]])
  }

  list(
    data = data,
    item_probs = item_probs,
    class_probs = class_probs,
    true_class = class_assign,
    model_code = 1,
    model_name = "MON"
  )
}

#' Generate data from Invariant Item Ordering model (Model 2)
#' Items maintain same relative difficulty ordering across ALL classes
#' But class ordering is NOT monotonic (scrambled class abilities)
#' AND items have DIFFERENT discriminations (violates Rasch equal discrim)
generate_model_2_IIO <- function(n, n_items, n_classes, seed) {
  set.seed(seed)

  # Item difficulties (fixed ordering from easiest to hardest)
  item_difficulties <- seq(0, 1, length.out = n_items)

  # Non-monotonic class abilities (STRONGLY scrambled)
  # Use a specific scrambled pattern that clearly violates monotonicity
  if (n_classes == 4) {
    class_abilities <- c(0.6, 0.3, 0.8, 0.45)  # Class 3 > Class 1 > Class 4 > Class 2
  } else {
    # Create scrambled order
    ordered_abilities <- seq(0.25, 0.75, length.out = n_classes)
    scramble_idx <- c(seq(2, n_classes, by = 2), seq(1, n_classes, by = 2))
    class_abilities <- ordered_abilities[scramble_idx[1:n_classes]]
  }

  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)

  # Items have DIFFERENT discriminations (violates Rasch)
  # But item ORDER is preserved within each class
  item_discriminations <- seq(0.3, 0.7, length.out = n_items)

  for (c in 1:n_classes) {
    for (i in 1:n_items) {
      # Different discriminations per item - violates Rasch equal discrimination
      item_probs[i, c] <- class_abilities[c] - item_discriminations[i] * item_difficulties[i] + 0.2
    }
  }

  # Add noise
  item_probs <- item_probs + matrix(rnorm(n_items * n_classes, 0, 0.015),
                                     nrow = n_items)
  item_probs <- pmax(pmin(item_probs, 0.92), 0.08)

  # Enforce IIO within each class (items sorted by difficulty)
  for (c in 1:n_classes) {
    item_probs[, c] <- sort(item_probs[, c], decreasing = TRUE)
  }

  # Generate data
  class_probs <- rep(1/n_classes, n_classes)
  class_assign <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    data[i, ] <- rbinom(n_items, 1, item_probs[, class_assign[i]])
  }

  list(
    data = data,
    item_probs = item_probs,
    class_probs = class_probs,
    true_class = class_assign,
    item_order = 1:n_items,
    model_code = 2,
    model_name = "IIO"
  )
}

#' Generate data from Double Monotonicity model (Model 3)
#' Both class monotonicity AND invariant item ordering hold
#' Mokken-style: UNEQUAL discriminations that still preserve double monotonicity
#' This differs from Rasch (equal discrim) but satisfies DM constraints
generate_model_3_DM <- function(n, n_items, n_classes, seed) {
  set.seed(seed)

  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)

  # Item difficulties (ordered from easy to hard)
  item_difficulties <- seq(0, 1, length.out = n_items)

  # Class abilities (monotonically ordered)
  class_abilities <- seq(0.2, 0.8, length.out = n_classes)

  # DIFFERENT discriminations per item (violates Rasch equal discrimination)
  # But structured so DM constraints are still satisfied
  item_discriminations <- seq(0.4, 0.8, length.out = n_items)

  for (i in 1:n_items) {
    for (c in 1:n_classes) {
      # Each item has different discrimination
      # P = base + discrim_i * (ability_c - difficulty_i)
      base_prob <- 0.5
      item_probs[i, c] <- base_prob + item_discriminations[i] * (class_abilities[c] - item_difficulties[i])
    }
  }

  # Add small noise
  item_probs <- item_probs + matrix(rnorm(n_items * n_classes, 0, 0.015),
                                     nrow = n_items)
  item_probs <- pmax(pmin(item_probs, 0.92), 0.08)

  # Enforce double monotonicity (project onto constraint space)
  # First ensure class monotonicity (each row non-decreasing)
  for (i in 1:n_items) {
    item_probs[i, ] <- sort(item_probs[i, ])
  }
  # Then ensure item ordering (each column non-increasing by item difficulty)
  for (c in 1:n_classes) {
    item_probs[, c] <- sort(item_probs[, c], decreasing = TRUE)
  }

  # Generate data
  class_probs <- rep(1/n_classes, n_classes)
  class_assign <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    data[i, ] <- rbinom(n_items, 1, item_probs[, class_assign[i]])
  }

  list(
    data = data,
    item_probs = item_probs,
    class_probs = class_probs,
    true_class = class_assign,
    item_order = 1:n_items,
    model_code = 3,
    model_name = "DM"
  )
}

#' Generate data from Latent Class Rasch model (Model 4)
#' Rasch parameterization: P(X=1) = exp(theta_c - delta_i) / (1 + exp(theta_c - delta_i))
#' Discrete ability classes with Rasch item response function
generate_model_4_LCR <- function(n, n_items, n_classes, seed) {
  set.seed(seed)

  # Discrete class ability levels (theta)
  theta <- seq(-2, 2, length.out = n_classes)

  # Item difficulties (delta)
  delta <- seq(-1.5, 1.5, length.out = n_items)

  # Compute item probabilities using Rasch formula
  item_probs <- matrix(0, nrow = n_items, ncol = n_classes)
  for (i in 1:n_items) {
    for (c in 1:n_classes) {
      item_probs[i, c] <- 1 / (1 + exp(-(theta[c] - delta[i])))
    }
  }

  # Generate data
  class_probs <- rep(1/n_classes, n_classes)
  class_assign <- sample(1:n_classes, n, replace = TRUE, prob = class_probs)

  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    prob_vec <- 1 / (1 + exp(-(theta[class_assign[i]] - delta)))
    data[i, ] <- rbinom(n_items, 1, prob_vec)
  }

  list(
    data = data,
    item_probs = item_probs,
    class_probs = class_probs,
    true_class = class_assign,
    theta = theta,
    delta = delta,
    model_code = 4,
    model_name = "LCR"
  )
}

#' Generate data from continuous Rasch Model (Model 5)
#' Standard Rasch with continuous theta ~ N(0, 1)
generate_model_5_RM <- function(n, n_items, seed) {
  set.seed(seed)

  # Continuous ability
  theta <- rnorm(n, mean = 0, sd = 1)

  # Item difficulties
  delta <- seq(-2, 2, length.out = n_items)

  # Generate data using Rasch model
  data <- matrix(0, nrow = n, ncol = n_items)
  for (i in 1:n) {
    prob_vec <- 1 / (1 + exp(-(theta[i] - delta)))
    data[i, ] <- rbinom(n_items, 1, prob_vec)
  }

  list(
    data = data,
    theta = theta,
    delta = delta,
    model_code = 5,
    model_name = "RM"
  )
}

# ============================================================
# MAIN VALIDATION FUNCTION
# ============================================================

#' Run validation study following paper's simulation design
#'
#' @param n_datasets Number of datasets to generate (paper used 90)
#' @param n_persons Number of persons per dataset (paper used 5000)
#' @param n_items Number of items (paper used 10)
#' @param n_classes_range Range of classes to use (paper used 2-6)
#' @param n_starts Number of random starts for estimation
#' @param verbose Print progress
#' @return List with results and generated datasets
run_paper_validation <- function(n_datasets = 30,
                                  n_persons = 1000,
                                  n_items = 10,
                                  n_classes_range = 3:5,
                                  n_starts = 5,
                                  verbose = TRUE) {

  if (verbose) {
    cat("=" |> rep(70) |> paste(collapse = ""), "\n")
    cat("QuantFit Validation Study\n")
    cat("Following Torres Irribarra & Diakow Simulation Design\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")
    cat("Settings:\n")
    cat("  Datasets:", n_datasets, "\n")
    cat("  Persons per dataset:", n_persons, "\n")
    cat("  Items:", n_items, "\n")
    cat("  Classes range:", paste(range(n_classes_range), collapse = "-"), "\n")
    cat("  Random starts:", n_starts, "\n\n")
  }

  # Model codes matching the paper
  model_codes <- 0:5
  model_names <- c("UN", "MON", "IIO", "DM", "LCR", "RM")

  # Allocate datasets across models (roughly equal)
  datasets_per_model <- ceiling(n_datasets / 6)

  all_results <- data.frame()
  all_datasets <- list()

  dataset_id <- 0

  for (model_code in model_codes) {
    model_name <- model_names[model_code + 1]

    if (verbose) {
      cat("\n--- Generating from Model", model_code, "(", model_name, ") ---\n")
    }

    for (rep in 1:datasets_per_model) {
      dataset_id <- dataset_id + 1
      if (dataset_id > n_datasets) break

      # Random number of classes (or use n_persons for RM)
      n_classes <- if (model_code == 5) {
        n_persons  # RM uses continuous theta
      } else {
        sample(n_classes_range, 1)
      }

      # Generate data
      seed <- dataset_id * 100 + model_code

      sim_data <- switch(as.character(model_code),
        "0" = generate_model_0_UN(n_persons, n_items, n_classes, seed),
        "1" = generate_model_1_MON(n_persons, n_items, n_classes, seed),
        "2" = generate_model_2_IIO(n_persons, n_items, n_classes, seed),
        "3" = generate_model_3_DM(n_persons, n_items, n_classes, seed),
        "4" = generate_model_4_LCR(n_persons, n_items, n_classes, seed),
        "5" = generate_model_5_RM(n_persons, n_items, seed)
      )

      # Store dataset
      all_datasets[[paste0("D", dataset_id)]] <- sim_data

      # Fit all models using QuantFit
      fit_n_classes <- if (model_code == 5) 4 else n_classes

      fit_results <- tryCatch({
        fit_all_quantfit_models(sim_data$data, fit_n_classes, n_starts)
      }, error = function(e) {
        if (verbose) cat("  Error fitting dataset", dataset_id, ":", e$message, "\n")
        NULL
      })

      if (is.null(fit_results)) next

      # Determine selected model
      bics <- c(
        UN = fit_results$BIC_UN,
        MON = fit_results$BIC_MON,
        IIO = fit_results$BIC_IIO,
        DM = fit_results$BIC_DM,
        LCR = fit_results$BIC_LCR,
        RM = fit_results$BIC_RM
      )
      selected_model <- names(which.min(bics))
      selected_code <- match(selected_model, model_names) - 1

      # Record result
      result <- data.frame(
        dataset_id = dataset_id,
        gen_model_code = model_code,
        gen_model_name = model_name,
        n_classes = if (model_code == 5) NA else n_classes,
        sel_model_code = selected_code,
        sel_model_name = selected_model,
        correct = (selected_model == model_name),
        fit_results,
        stringsAsFactors = FALSE
      )
      all_results <- rbind(all_results, result)

      if (verbose) {
        marker <- if (selected_model == model_name) "CORRECT" else ""
        cat(sprintf("  D%02d: Gen=%s, C=%s, Sel=%s %s\n",
            dataset_id, model_name,
            if (model_code == 5) "cont" else as.character(n_classes),
            selected_model, marker))
      }
    }
    if (dataset_id > n_datasets) break
  }

  # Print summary
  if (verbose) {
    cat("\n\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n")
    cat("VALIDATION RESULTS SUMMARY\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

    # Confusion matrix
    cat("Selection Matrix (rows=generating, cols=selected):\n\n")
    conf_matrix <- table(
      factor(all_results$gen_model_name, levels = model_names),
      factor(all_results$sel_model_name, levels = model_names)
    )
    print(conf_matrix)

    # Accuracy by model
    cat("\n\nAccuracy by Generating Model:\n")
    cat("-" |> rep(40) |> paste(collapse = ""), "\n")
    for (mn in model_names) {
      subset <- all_results[all_results$gen_model_name == mn, ]
      if (nrow(subset) > 0) {
        n_correct <- sum(subset$correct)
        n_total <- nrow(subset)
        cat(sprintf("  %s: %d/%d (%.1f%%)\n", mn, n_correct, n_total,
            100 * n_correct / n_total))
      }
    }

    total_correct <- sum(all_results$correct)
    total_n <- nrow(all_results)
    cat(sprintf("\n  OVERALL: %d/%d (%.1f%%)\n", total_correct, total_n,
        100 * total_correct / total_n))

    # Scale-level accuracy (Quantitative vs Non-Quantitative)
    cat("\n\nScale-Level Classification:\n")
    cat("-" |> rep(40) |> paste(collapse = ""), "\n")
    all_results$gen_scale <- ifelse(all_results$gen_model_code >= 4, "Quantitative", "Non-Quantitative")
    all_results$sel_scale <- ifelse(all_results$sel_model_code >= 4, "Quantitative", "Non-Quantitative")
    scale_matrix <- table(all_results$gen_scale, all_results$sel_scale)
    print(scale_matrix)
    scale_acc <- sum(diag(scale_matrix)) / sum(scale_matrix)
    cat(sprintf("\nScale-level accuracy: %.1f%%\n", 100 * scale_acc))
  }

  invisible(list(
    results = all_results,
    datasets = all_datasets,
    confusion_matrix = conf_matrix
  ))
}

#' Fit all QuantFit models to a dataset
#' @keywords internal
fit_all_quantfit_models <- function(data, n_classes, n_starts) {

  # Fit each model
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
    BIC_RM = if (!is.null(res_rm)) BIC(res_rm) else NA,
    Conv_UN = if (!is.null(res_un)) res_un$convergence else NA,
    Conv_MON = if (!is.null(res_mon)) res_mon$convergence else NA,
    Conv_IIO = if (!is.null(res_iio)) res_iio$convergence else NA,
    Conv_DM = if (!is.null(res_dm)) res_dm$convergence else NA,
    Conv_LCR = if (!is.null(res_lcr)) res_lcr$convergence else NA,
    Conv_RM = if (!is.null(res_rm)) res_rm$convergence else NA
  )
}

# ============================================================
# QUICK TEST FUNCTION
# ============================================================

#' Quick validation test with small sample
quick_validation_test <- function() {
  cat("Running quick validation test...\n\n")

  run_paper_validation(
    n_datasets = 12,      # 2 per model
    n_persons = 500,      # Smaller for speed
    n_items = 10,
    n_classes_range = 4:4, # Fixed 4 classes
    n_starts = 3,
    verbose = TRUE
  )
}

# ============================================================
# USAGE
# ============================================================

cat("QuantFit Paper Validation Script Loaded\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")
cat("Usage:\n")
cat("  quick_validation_test()  - Fast test with 12 datasets\n")
cat("  run_paper_validation()   - Full validation (30 datasets)\n")
cat("  run_paper_validation(n_datasets=90, n_persons=5000)  - Paper settings\n")
cat("\n")
