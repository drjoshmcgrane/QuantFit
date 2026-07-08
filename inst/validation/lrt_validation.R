# ============================================================
# QuantFit Validation with Likelihood Ratio Tests
# Using Paper's Exact Simulation Functions
# ============================================================

# Load QuantFit package
if (!requireNamespace("QuantFit", quietly = TRUE)) {
  if (file.exists("NAMESPACE")) {
    pkgload::load_all(".")
  }
} else {
  library(QuantFit)
}

# ============================================================
# PAPER'S SIMULATION FUNCTIONS (from functions/simFunc3.r)
# ============================================================

#' Generate latent parameters - Paper's exact method
#' @param nC Number of classes (or persons for model 5)
#' @param nI Number of items
#' @param model Model type: 0=UN, 1=MON, 2=IIO, 3=DM, 4=LCR, 5=RM
#' @param slope Discrimination multiplier
sim.LParam <- function(nC = 4, nI = 20, model = 0, slope = 1,
                       logit = NULL, lC = NULL, lI = NULL) {

  if (model > 5) stop("Invalid model")

  output <- list()
  output$model <- model
  output$nC <- nC
  output$nI <- nI

  # Model 5: Continuous Rasch
  if (model == 5) {
    nR <- nC
    output$nR <- nR
    output$nC <- 1

    if (is.null(lC)) lC <- rnorm(nR)
    if (is.null(lI)) lI <- runif(nI, min = -3, max = 3)

    logit <- matrix(rep(sort(-lI), nR), ncol = nI, byrow = TRUE)
    logit <- sweep(logit, 1, lC, "+")

    output$lC <- lC
    output$lI <- lI
  }

  # Model 4: Latent Class Rasch
  if (model == 4) {
    if (is.null(lC)) {
      # Ensure minimum separation between class locations
      sep <- TRUE
      while (sep) {
        lC <- runif(nC, min = -4, max = 4)
        diffs <- abs(outer(lC, lC, "-"))
        sep <- any(diffs > 0 & diffs < 0.5)
      }
    }
    if (is.null(lI)) lI <- runif(nI, min = -4, max = 4)

    logit <- matrix(rep(sort(-lI), nC), ncol = nI, byrow = TRUE)
    logit <- sweep(logit, 1, lC, "+")

    output$lC <- lC
    output$lI <- lI
  }

  # Models 0-3: Start with random logits from Uniform(-4, 4)
  if (model < 4 && is.null(logit)) {
    logit <- matrix(runif(nC * nI, min = -4, max = 4), ncol = nI)
  }

  # Apply constraints by sorting
  if (model == 3) {  # DM: Sort both
    logit <- apply(logit, 2, sort)
    logit <- t(apply(logit, 1, sort))
  }
  if (model == 2) {  # IIO: Sort rows
    logit <- t(apply(logit, 1, sort))
  }
  if (model == 1) {  # MON: Sort columns
    logit <- apply(logit, 2, sort)
  }

  # Convert to probabilities
  prob <- 1 / (1 + exp(-slope * logit))

  colnames(logit) <- paste0("i", 1:ncol(logit))
  colnames(prob) <- paste0("i", 1:ncol(prob))
  rownames(logit) <- paste0("c", 1:nrow(logit))
  rownames(prob) <- paste0("c", 1:nrow(prob))

  output$logit <- slope * logit
  output$prob <- prob

  output
}

#' Generate observed responses - Paper's exact method
sim.OResp <- function(LPar, prC = NULL, nR = NULL) {
  output <- list()

  if (LPar$model < 5) {
    if (is.null(prC) || is.null(nR)) stop("Need prC and nR for models 0-4")
    if (LPar$nC != length(prC)) stop("prC length mismatch")

    # Assign respondents to classes
    cR <- rep(1:length(prC), rmultinom(1, nR, prC))
    cR <- sample(cR)

    prRI <- LPar$prob[cR, , drop = FALSE]
    genD <- cbind(cR, prRI)
  } else {
    nR <- LPar$nR
    prRI <- LPar$prob
    genD <- prRI
  }

  # Generate binary responses
  respRI <- matrix(rbinom(length(prRI), 1, prRI), ncol = ncol(prRI))
  colnames(respRI) <- paste0("i", 1:ncol(respRI))

  output$obsData <- respRI
  output$genData <- genD
  output$model <- LPar$model
  output$nC <- LPar$nC
  output$nI <- LPar$nI
  output$nR <- nR
  output$logit <- LPar$logit
  output$prob <- LPar$prob

  if (LPar$model < 5) {
    output$cR <- table(cR)
    output$prC <- prC
  }
  if (LPar$model >= 4) {
    output$lC <- LPar$lC
    output$lI <- LPar$lI
  }

  output
}

# ============================================================
# CONSTRAINT CHECKING FUNCTIONS
# ============================================================

#' Check if probability matrix satisfies class monotonicity (MON)
#' MON: For each item, P(X=1|c) increases with c
#' Note: We first order classes by their average probability
#' @param prob Matrix of probabilities (classes x items)
#' @param tolerance Proportion of violations allowed
#' @return TRUE if constraint is approximately satisfied
check_mon_constraint <- function(prob, tolerance = 0.15) {
  # Order classes by their average probability (theta proxy)
  class_means <- rowMeans(prob)
  class_order <- order(class_means)
  prob_ordered <- prob[class_order, , drop = FALSE]

  n_classes <- nrow(prob_ordered)
  n_items <- ncol(prob_ordered)

  violations <- 0
  comparisons <- 0

  for (j in 1:n_items) {
    for (c in 1:(n_classes - 1)) {
      comparisons <- comparisons + 1
      if (prob_ordered[c + 1, j] < prob_ordered[c, j]) {
        violations <- violations + 1
      }
    }
  }

  violation_rate <- violations / comparisons
  violation_rate <= tolerance
}

#' Check if probability matrix satisfies invariant item ordering (IIO)
#' IIO: Items have same difficulty ordering across all classes
#' @param prob Matrix of probabilities (classes x items)
#' @param tolerance Proportion of violations allowed
#' @return TRUE if constraint is approximately satisfied
check_iio_constraint <- function(prob, tolerance = 0.15) {
  n_classes <- nrow(prob)
  n_items <- ncol(prob)

  # Count ordering disagreements across class pairs
  violations <- 0
  comparisons <- 0

  for (c1 in 1:(n_classes - 1)) {
    for (c2 in (c1 + 1):n_classes) {
      # Compare item orderings between classes c1 and c2
      for (i in 1:(n_items - 1)) {
        for (j in (i + 1):n_items) {
          comparisons <- comparisons + 1
          # Check if items i,j have same relative ordering in both classes
          diff_c1 <- prob[c1, i] - prob[c1, j]
          diff_c2 <- prob[c2, i] - prob[c2, j]
          # Violation if signs differ (one positive, one negative)
          if (sign(diff_c1) != sign(diff_c2) && abs(diff_c1) > 0.05 && abs(diff_c2) > 0.05) {
            violations <- violations + 1
          }
        }
      }
    }
  }

  if (comparisons == 0) return(TRUE)
  violation_rate <- violations / comparisons
  violation_rate <= tolerance
}

#' Check if probability matrix satisfies double monotonicity (DM)
#' DM: Both MON and IIO constraints hold
check_dm_constraint <- function(prob, tolerance = 0.10) {
  check_mon_constraint(prob, tolerance) && check_iio_constraint(prob, tolerance)
}

#' Check if probability matrix is consistent with Rasch structure
#' Rasch: prob[c,j] = f(theta_c - beta_j) for some theta and beta
#' Approximate check: look for additive structure in logits
check_rasch_constraint <- function(prob, tolerance = 0.25) {
  # Clip probabilities to avoid infinite logits
  prob_clipped <- pmax(pmin(prob, 0.999), 0.001)
  logit <- log(prob_clipped / (1 - prob_clipped))

  # Order classes by mean logit (theta proxy)
  class_means <- rowMeans(logit)
  class_order <- order(class_means)
  logit_ordered <- logit[class_order, , drop = FALSE]

  n_classes <- nrow(logit_ordered)
  n_items <- ncol(logit_ordered)

  # If Rasch, logit[c,j] = theta_c - beta_j
  # So for each pair of classes (c1, c2):
  # logit[c2,j] - logit[c1,j] = theta_c2 - theta_c1 (constant across items)

  # Compute coefficient of variation of inter-class logit differences
  cv_values <- c()
  for (c in 1:(n_classes - 1)) {
    row_diffs <- logit_ordered[c + 1, ] - logit_ordered[c, ]
    if (sd(row_diffs) > 0 && abs(mean(row_diffs)) > 0.1) {
      cv <- sd(row_diffs) / abs(mean(row_diffs))
      cv_values <- c(cv_values, cv)
    }
  }

  if (length(cv_values) == 0) return(FALSE)

  mean(cv_values) < tolerance
}

# ============================================================
# LIKELIHOOD RATIO TEST FUNCTIONS
# ============================================================

#' Compute likelihood ratio statistic
#' @param ll_restricted Log-likelihood of restricted (constrained) model
#' @param ll_full Log-likelihood of full (unconstrained) model
#' @return LR statistic (positive if full model fits better)
lr_statistic <- function(ll_restricted, ll_full) {
  -2 * (ll_restricted - ll_full)
}

#' Test if constrained model fits adequately vs unconstrained
#' For models with same parameters, use tolerance-based approach
#' @param ll_constrained Log-likelihood of constrained model
#' @param ll_unconstrained Log-likelihood of unconstrained model
#' @param n_obs Number of observations
#' @param alpha Significance level for the test
#' @return TRUE if constrained model is adequate (prefer simpler interpretation)
constraint_adequate <- function(ll_constrained, ll_unconstrained, n_obs, alpha = 0.05) {
  # Handle NA values
  if (is.na(ll_constrained) || is.na(ll_unconstrained)) {
    return(FALSE)
  }

  # LR statistic
  lr <- -2 * (ll_constrained - ll_unconstrained)

  # For same-parameter models, use a practical threshold

  # Based on AIC logic: prefer constrained if LL difference < 2
  # Or based on per-observation threshold
  threshold <- 2  # AIC-style threshold

  # Constrained is adequate if it doesn't fit much worse
  lr < threshold
}

#' Vuong test for non-nested model comparison
#' @param ll1 Vector of observation-level log-likelihoods for model 1
#' @param ll2 Vector of observation-level log-likelihoods for model 2
#' @return List with z-statistic and p-value
vuong_test <- function(ll1, ll2) {
  n <- length(ll1)
  diff <- ll1 - ll2
  mean_diff <- mean(diff)
  var_diff <- var(diff)

  if (var_diff == 0) {
    return(list(z = 0, p_value = 1, preferred = "neither"))
  }

  z <- sqrt(n) * mean_diff / sqrt(var_diff)
  p_value <- 2 * pnorm(-abs(z))

  preferred <- if (z > 1.96) "model1" else if (z < -1.96) "model2" else "neither"

  list(z = z, p_value = p_value, preferred = preferred)
}

#' Model selection using constraint checking on UN estimates
#' @param data Binary response matrix
#' @param n_classes Number of classes for discrete models
#' @param n_starts Random starts
#' @param tolerance Tolerance for constraint violations
#' @return List with selected model and all fits
select_model_constraint <- function(data, n_classes, n_starts = 5, tolerance = 0.15) {

  n_obs <- nrow(data)

  # Fit UN model to get unconstrained estimates
  un_fit <- tryCatch(
    suppressWarnings(fit_un(data, n_classes, n_starts = n_starts, seed = 1)),
    error = function(e) NULL)

  if (is.null(un_fit)) {
    return(list(selected = NA, interpretation = "FIT_FAILED",
                constraints = list(mon = NA, iio = NA, dm = NA, rasch = NA),
                bics = c(UN = NA, RM = NA)))
  }

  # Get estimated probabilities from UN fit
  # item_probs is items x classes, we need classes x items
  prob <- t(un_fit$item_probs)

  # Check constraints on estimated probabilities
  mon_satisfied <- check_mon_constraint(prob, tolerance)
  iio_satisfied <- check_iio_constraint(prob, tolerance)
  dm_satisfied <- mon_satisfied && iio_satisfied
  rasch_satisfied <- check_rasch_constraint(prob, tolerance)

  # Also fit RM for comparison
  rm_fit <- tryCatch(
    suppressWarnings(fit_rm(data, verbose = FALSE)),
    error = function(e) NULL)

  # Get BICs
  un_bic <- BIC(un_fit)
  rm_bic <- if (!is.null(rm_fit)) BIC(rm_fit) else Inf

  # Also fit LCR for BIC comparison
  lcr_fit <- tryCatch(
    suppressWarnings(fit_lcr(data, n_classes, n_starts = n_starts, seed = 5)),
    error = function(e) NULL)
  lcr_bic <- if (!is.null(lcr_fit)) BIC(lcr_fit) else Inf

  # Selection logic based on constraints + BIC
  # 1. If both DM and Rasch structure satisfied, prefer RM if BIC supports it
  # 2. If only DM satisfied, select DM (or LCR/RM if BIC strongly supports)
  # 3. If only MON or IIO, select that
  # 4. Otherwise UN

  best_quant_bic <- min(rm_bic, lcr_bic)

  if (dm_satisfied && rasch_satisfied) {
    # Strong Rasch structure - check if quantitative model has better BIC
    if (best_quant_bic < un_bic) {
      selected <- if (rm_bic < lcr_bic) "RM" else "LCR"
      interpretation <- "QUANTITATIVE (DM + Rasch structure)"
    } else {
      selected <- "DM"
      interpretation <- "ORDINAL (double monotonicity with Rasch-like structure)"
    }
  } else if (dm_satisfied) {
    # DM but not Rasch - only select quantitative if BIC strongly better
    if (best_quant_bic < un_bic - 20) {
      selected <- if (rm_bic < lcr_bic) "RM" else "LCR"
      interpretation <- "QUANTITATIVE (DM + BIC preference)"
    } else {
      selected <- "DM"
      interpretation <- "ORDINAL (double monotonicity)"
    }
  } else if (rasch_satisfied && best_quant_bic < un_bic) {
    # Rasch structure without full DM
    selected <- if (rm_bic < lcr_bic) "RM" else "LCR"
    interpretation <- "QUANTITATIVE (Rasch structure detected)"
  } else if (mon_satisfied && iio_satisfied) {
    selected <- "DM"
    interpretation <- "ORDINAL (double monotonicity)"
  } else if (iio_satisfied) {
    selected <- "IIO"
    interpretation <- "ORDINAL (invariant item ordering)"
  } else if (mon_satisfied) {
    selected <- "MON"
    interpretation <- "ORDINAL (class monotonicity)"
  } else {
    selected <- "UN"
    interpretation <- "CLASSIFICATORY (no ordinal structure)"
  }

  list(
    selected = selected,
    interpretation = interpretation,
    constraints = list(
      mon = mon_satisfied,
      iio = iio_satisfied,
      dm = dm_satisfied,
      rasch = rasch_satisfied
    ),
    bics = c(UN = un_bic, LCR = lcr_bic, RM = rm_bic)
  )
}

#' Model selection using LRT-based successive testing
#' @param data Binary response matrix
#' @param n_classes Number of classes for discrete models
#' @param n_starts Random starts
#' @param alpha Significance level
#' @return List with selected model and all fits
select_model_lrt <- function(data, n_classes, n_starts = 5, alpha = 0.05) {

  n_obs <- nrow(data)

  # Fit all models
  fits <- list()

  fits$UN <- tryCatch(
    suppressWarnings(fit_un(data, n_classes, n_starts = n_starts, seed = 1)),
    error = function(e) NULL)

  fits$MON <- tryCatch(
    suppressWarnings(fit_mon(data, n_classes, n_starts = n_starts, seed = 2)),
    error = function(e) NULL)

  fits$IIO <- tryCatch(
    suppressWarnings(fit_iio(data, n_classes, n_starts = n_starts, seed = 3)),
    error = function(e) NULL)

  fits$DM <- tryCatch(
    suppressWarnings(fit_dm(data, n_classes, n_starts = n_starts, seed = 4)),
    error = function(e) NULL)

  fits$LCR <- tryCatch(
    suppressWarnings(fit_lcr(data, n_classes, n_starts = n_starts, seed = 5)),
    error = function(e) NULL)

  fits$RM <- tryCatch(
    suppressWarnings(fit_rm(data, verbose = FALSE)),
    error = function(e) NULL)

  # Get log-likelihoods
  get_ll <- function(m) if (!is.null(fits[[m]])) fits[[m]]$loglik else NA
  lls <- sapply(c("UN", "MON", "IIO", "DM", "LCR", "RM"), get_ll)

  # Get BICs
  get_bic <- function(m) if (!is.null(fits[[m]])) BIC(fits[[m]]) else NA
  bics <- sapply(c("UN", "MON", "IIO", "DM", "LCR", "RM"), get_bic)

  # LRT-based selection strategy
  # Step 1: Test for class monotonicity (UN vs MON)
  mon_adequate <- constraint_adequate(lls["MON"], lls["UN"], n_obs, alpha)

  # Step 2: Test for item ordering (UN vs IIO)
  iio_adequate <- constraint_adequate(lls["IIO"], lls["UN"], n_obs, alpha)

  # Step 3: If both MON and IIO adequate, test DM
  if (isTRUE(mon_adequate) && isTRUE(iio_adequate)) {
    # Test DM against the better of MON/IIO
    best_ordinal <- max(lls["MON"], lls["IIO"], na.rm = TRUE)
    if (!is.finite(best_ordinal)) {
      dm_adequate <- FALSE
    } else {
      dm_adequate <- constraint_adequate(lls["DM"], best_ordinal, n_obs, alpha)
    }
  } else {
    dm_adequate <- FALSE
  }

  # Step 4: Determine best non-quantitative model
  # Use overall best BIC among UN, MON, IIO, DM
  non_quant_bics <- bics[c("UN", "MON", "IIO", "DM")]
  best_non_quant <- names(which.min(non_quant_bics))
  best_non_quant_bic <- min(non_quant_bics, na.rm = TRUE)

  # Step 5: Compare best non-quantitative to best quantitative
  # Require substantial improvement (BIC difference > 10) to prefer quantitative
  quant_bics <- bics[c("LCR", "RM")]
  best_quant <- names(which.min(quant_bics))
  best_quant_bic <- min(quant_bics, na.rm = TRUE)

  # Quantitative only wins if substantially better (not just marginally)
  quant_threshold <- 10  # Require BIC improvement of at least 10
  quant_preferred <- is.finite(best_quant_bic) &&
                     is.finite(best_non_quant_bic) &&
                     (best_non_quant_bic - best_quant_bic) > quant_threshold

  # Which quantitative model?
  rm_preferred_over_lcr <- !is.na(bics["RM"]) && !is.na(bics["LCR"]) && bics["RM"] < bics["LCR"]

  # Determine selected model based on LRT logic
  if (quant_preferred) {
    # Quantitative model substantially better
    if (rm_preferred_over_lcr) {
      selected <- "RM"
      interpretation <- "QUANTITATIVE (continuous)"
    } else {
      selected <- "LCR"
      interpretation <- "QUANTITATIVE (discrete)"
    }
  } else if (isTRUE(dm_adequate)) {
    selected <- "DM"
    interpretation <- "ORDINAL (double monotonicity)"
  } else if (isTRUE(mon_adequate) && isTRUE(iio_adequate)) {
    # Both orderings but not DM - multidimensional?
    selected <- if (!is.na(lls["MON"]) && !is.na(lls["IIO"]) && lls["MON"] > lls["IIO"]) "MON" else "IIO"
    interpretation <- "ORDINAL (partial)"
  } else if (isTRUE(mon_adequate)) {
    selected <- "MON"
    interpretation <- "ORDINAL (class monotonicity)"
  } else if (isTRUE(iio_adequate)) {
    selected <- "IIO"
    interpretation <- "ORDINAL (item ordering)"
  } else {
    selected <- "UN"
    interpretation <- "CLASSIFICATORY"
  }

  # Also compute BIC-based selection for comparison
  bic_selected <- names(which.min(bics))

  list(
    lrt_selected = selected,
    bic_selected = bic_selected,
    interpretation = interpretation,
    log_likelihoods = lls,
    bics = bics,
    fits = fits,
    tests = list(
      mon_adequate = mon_adequate,
      iio_adequate = iio_adequate,
      dm_adequate = dm_adequate,
      quant_preferred = quant_preferred,
      rm_preferred = rm_preferred_over_lcr
    )
  )
}

# ============================================================
# VALIDATION FUNCTION
# ============================================================

#' Run validation comparing BIC and LRT selection
#' @param n_datasets Number of datasets per model
#' @param n_persons Number of persons
#' @param n_items Number of items
#' @param n_classes Number of classes
#' @param n_starts Random starts
run_lrt_validation <- function(n_datasets = 4,
                                n_persons = 1000,
                                n_items = 10,
                                n_classes = 4,
                                n_starts = 5,
                                verbose = TRUE) {

  model_names <- c("UN", "MON", "IIO", "DM", "LCR", "RM")

  if (verbose) {
    cat("=" |> rep(70) |> paste(collapse = ""), "\n")
    cat("QuantFit Validation: BIC vs LRT Selection\n")
    cat("Using Paper's Exact Simulation Functions\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")
    cat("Settings: n=", n_persons, ", items=", n_items, ", classes=", n_classes, "\n\n")
  }

  results <- data.frame()

  for (model_code in 0:5) {
    model_name <- model_names[model_code + 1]

    if (verbose) cat("--- Generating from", model_name, "---\n")

    for (rep in 1:n_datasets) {
      seed <- model_code * 1000 + rep * 100

      # Generate data using paper's method
      set.seed(seed)
      if (model_code == 5) {
        LPar <- sim.LParam(nC = n_persons, nI = n_items, model = model_code)
        LOResp <- sim.OResp(LPar)
      } else {
        LPar <- sim.LParam(nC = n_classes, nI = n_items, model = model_code)
        prC <- rep(1/n_classes, n_classes)
        LOResp <- sim.OResp(LPar, prC = prC, nR = n_persons)
      }

      data <- LOResp$obsData

      # Select model using both methods
      selection <- select_model_lrt(data, n_classes, n_starts)

      lrt_correct <- selection$lrt_selected == model_name
      bic_correct <- selection$bic_selected == model_name

      result <- data.frame(
        gen_model = model_name,
        rep = rep,
        lrt_selected = selection$lrt_selected,
        bic_selected = selection$bic_selected,
        lrt_correct = lrt_correct,
        bic_correct = bic_correct,
        stringsAsFactors = FALSE
      )
      results <- rbind(results, result)

      if (verbose) {
        lrt_mark <- if (lrt_correct) "v" else " "
        bic_mark <- if (bic_correct) "v" else " "
        cat(sprintf("  Rep %d: LRT=%s[%s] BIC=%s[%s]\n",
            rep, selection$lrt_selected, lrt_mark,
            selection$bic_selected, bic_mark))
      }
    }
  }

  # Summary
  if (verbose) {
    cat("\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n")
    cat("RESULTS SUMMARY\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

    cat("LRT Selection Confusion Matrix:\n")
    print(table(factor(results$gen_model, levels = model_names),
                factor(results$lrt_selected, levels = model_names)))

    cat("\nBIC Selection Confusion Matrix:\n")
    print(table(factor(results$gen_model, levels = model_names),
                factor(results$bic_selected, levels = model_names)))

    cat("\n\nAccuracy Comparison:\n")
    cat("-" |> rep(50) |> paste(collapse = ""), "\n")
    cat(sprintf("%-8s %8s %8s\n", "Model", "LRT", "BIC"))
    cat("-" |> rep(50) |> paste(collapse = ""), "\n")

    for (m in model_names) {
      sub <- results[results$gen_model == m, ]
      lrt_acc <- mean(sub$lrt_correct) * 100
      bic_acc <- mean(sub$bic_correct) * 100
      cat(sprintf("%-8s %7.1f%% %7.1f%%\n", m, lrt_acc, bic_acc))
    }

    cat("-" |> rep(50) |> paste(collapse = ""), "\n")
    cat(sprintf("%-8s %7.1f%% %7.1f%%\n", "OVERALL",
        mean(results$lrt_correct) * 100,
        mean(results$bic_correct) * 100))

    # Scale-level accuracy
    cat("\n\nScale-Level Accuracy:\n")
    results$gen_scale <- ifelse(results$gen_model %in% c("LCR", "RM"), "Quant", "NonQuant")
    results$lrt_scale <- ifelse(results$lrt_selected %in% c("LCR", "RM"), "Quant", "NonQuant")
    results$bic_scale <- ifelse(results$bic_selected %in% c("LCR", "RM"), "Quant", "NonQuant")

    lrt_scale_acc <- mean(results$gen_scale == results$lrt_scale) * 100
    bic_scale_acc <- mean(results$gen_scale == results$bic_scale) * 100

    cat(sprintf("  LRT: %.1f%%\n", lrt_scale_acc))
    cat(sprintf("  BIC: %.1f%%\n", bic_scale_acc))
  }

  invisible(results)
}

# Quick test
quick_lrt_test <- function() {
  run_lrt_validation(n_datasets = 2, n_persons = 500, n_starts = 3)
}

#' Run validation using constraint-based selection
run_constraint_validation <- function(n_datasets = 4,
                                       n_persons = 1000,
                                       n_items = 10,
                                       n_classes = 4,
                                       n_starts = 5,
                                       verbose = TRUE) {

  model_names <- c("UN", "MON", "IIO", "DM", "LCR", "RM")

  if (verbose) {
    cat("=" |> rep(70) |> paste(collapse = ""), "\n")
    cat("QuantFit Validation: Constraint-Based Selection\n")
    cat("Using Paper's Exact Simulation Functions\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")
    cat("Settings: n=", n_persons, ", items=", n_items, ", classes=", n_classes, "\n\n")
  }

  results <- data.frame()

  for (model_code in 0:5) {
    model_name <- model_names[model_code + 1]

    if (verbose) cat("--- Generating from", model_name, "---\n")

    for (rep in 1:n_datasets) {
      seed <- model_code * 1000 + rep * 100

      # Generate data using paper's method
      set.seed(seed)
      if (model_code == 5) {
        LPar <- sim.LParam(nC = n_persons, nI = n_items, model = model_code)
        LOResp <- sim.OResp(LPar)
      } else {
        LPar <- sim.LParam(nC = n_classes, nI = n_items, model = model_code)
        prC <- rep(1/n_classes, n_classes)
        LOResp <- sim.OResp(LPar, prC = prC, nR = n_persons)
      }

      data <- LOResp$obsData

      # Select model using constraint checking
      selection <- select_model_constraint(data, n_classes, n_starts)

      constraint_correct <- selection$selected == model_name

      result <- data.frame(
        gen_model = model_name,
        rep = rep,
        constraint_selected = selection$selected,
        constraint_correct = constraint_correct,
        mon = selection$constraints$mon,
        iio = selection$constraints$iio,
        dm = selection$constraints$dm,
        rasch = selection$constraints$rasch,
        stringsAsFactors = FALSE
      )
      results <- rbind(results, result)

      if (verbose) {
        mark <- if (constraint_correct) "v" else " "
        cat(sprintf("  Rep %d: Sel=%s[%s] (MON=%s IIO=%s DM=%s Rasch=%s)\n",
            rep, selection$selected, mark,
            if (selection$constraints$mon) "Y" else "N",
            if (selection$constraints$iio) "Y" else "N",
            if (selection$constraints$dm) "Y" else "N",
            if (selection$constraints$rasch) "Y" else "N"))
      }
    }
  }

  # Summary
  if (verbose) {
    cat("\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n")
    cat("RESULTS SUMMARY\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

    cat("Constraint Selection Confusion Matrix:\n")
    print(table(factor(results$gen_model, levels = model_names),
                factor(results$constraint_selected, levels = model_names)))

    cat("\n\nAccuracy by Model:\n")
    cat("-" |> rep(50) |> paste(collapse = ""), "\n")

    for (m in model_names) {
      sub <- results[results$gen_model == m, ]
      acc <- mean(sub$constraint_correct) * 100
      cat(sprintf("  %s: %.1f%%\n", m, acc))
    }

    cat("-" |> rep(50) |> paste(collapse = ""), "\n")
    cat(sprintf("  OVERALL: %.1f%%\n", mean(results$constraint_correct) * 100))

    # Scale-level accuracy
    cat("\n\nScale-Level Accuracy:\n")
    results$gen_scale <- ifelse(results$gen_model %in% c("LCR", "RM"), "Quant", "NonQuant")
    results$sel_scale <- ifelse(results$constraint_selected %in% c("LCR", "RM"), "Quant", "NonQuant")

    scale_acc <- mean(results$gen_scale == results$sel_scale) * 100
    cat(sprintf("  Constraint: %.1f%%\n", scale_acc))
  }

  invisible(results)
}

cat("\nLRT Validation Script Loaded\n")
cat("Usage: run_lrt_validation() or run_constraint_validation()\n\n")
