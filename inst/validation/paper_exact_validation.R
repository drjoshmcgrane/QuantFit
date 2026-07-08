# ============================================================
# QuantFit Package Validation
# Using EXACT Paper Simulation Functions
# From Torres Irribarra & Diakow
# ============================================================
#
# This script uses the original simulation functions from the paper
# to generate data and validate that QuantFit correctly recovers
# the generating model structure.
# ============================================================

# ============================================================
# ORIGINAL PAPER SIMULATION FUNCTIONS
# Copied from functions/simFunc3.r
# ============================================================

#' Generate latent parameters for simulation
#'
#' @param nC Number of classes (or number of persons for model 5)
#' @param nI Number of items
#' @param model Model type: 0=UN, 1=MON, 2=IIO, 3=DM, 4=LCR, 5=RM
#' @param slope Discrimination parameter (default 1)
#' @param logit Optional pre-specified logit matrix
#' @param lC Optional class locations (for models 4-5)
#' @param lI Optional item difficulties (for models 4-5)
#' @return List with model parameters
sim.LParam <- function(nC = 4, nI = 20, model = 0, slope = 1,
                       logit = NULL, lC = NULL, lI = NULL) {

  if (model > 5) { return(print('Error: Invalid Model.')) }

  if (!is.null(lC)) { nC <- length(lC) }
  if (!is.null(lI)) { nI <- length(lI) }

  output <- list()

  output["model"] <- model
  output["nC"]    <- nC
  output["nI"]    <- nI

  if (model == 5) {
    # Rasch Model: nC is number of persons, continuous theta
    nR <- nC
    output["nR"] <- nR
    output["nC"] <- 1

    if (is.null(lC)) { lC <- rnorm(nC) }
    if (is.null(lI)) { lI <- runif(nI, min = -3, max = 3) }

    logit <- matrix(rep(sort(-lI), nR), ncol = nI, byrow = TRUE)
    logit <- apply(logit, 2, '+', lC)

    output[["lC"]] <- lC
    output[["lI"]] <- lI
  }

  if (model == 4) {
    # Latent Class Rasch: discrete class locations
    if (is.null(lC)) {
      sep <- TRUE
      while (sep) {
        lC <- runif(nC, min = -4, max = 4)
        # Ensure minimum separation of 0.5 between classes
        sep <- any(0 < abs(kronecker(lC, lC, FUN = '-')) &
                     abs(kronecker(lC, lC, FUN = '-')) < 0.5)
      }
    }

    if (is.null(lI)) { lI <- runif(nI, min = -4, max = 4) }

    logit <- matrix(rep(sort(-lI), nC), ncol = nI, byrow = TRUE)
    logit <- apply(logit, 2, '+', lC)

    output[["lC"]] <- lC
    output[["lI"]] <- lI
  }

  # Models 0-3: Start with random logits
  if (model < 4 & is.null(logit)) {
    logit <- matrix(runif(nC * nI, min = -4, max = 4), ncol = nI)
  }

  # Model 3 (DM): Sort both ways for double monotonicity
  if (model == 3) {
    logit <- apply(logit, 2, sort)      # Sort columns (class monotonicity)
    logit <- t(apply(logit, 1, sort))   # Sort rows (item ordering)
  }

  # Model 2 (IIO): Sort rows for invariant item ordering
  if (model == 2) {
    logit <- t(apply(logit, 1, sort))
  }

  # Model 1 (MON): Sort columns for class monotonicity
  if (model == 1) {
    logit <- apply(logit, 2, sort)
  }

  # Convert logits to probabilities
  prob <- 1 / (1 + exp(-slope * logit))

  colnames(logit) <- paste("i", seq(1:ncol(logit)), sep = '')
  colnames(prob)  <- paste("i", seq(1:ncol(logit)), sep = '')
  rownames(logit) <- paste("c", seq(1:nrow(logit)), sep = '')
  rownames(prob)  <- paste("c", seq(1:nrow(logit)), sep = '')

  output[["logit"]] <- slope * logit
  output[["prob"]]  <- prob

  return(output)
}

#' Generate observed responses from latent parameters
#'
#' @param LPar Output from sim.LParam
#' @param prC Class proportions (vector summing to 1)
#' @param nR Number of respondents
#' @return List with observed data and generating parameters
sim.OResp <- function(LPar, prC = NULL, nR = NULL) {

  output <- list()

  if (LPar["model"] < 5) {
    if (is.null(prC) | is.null(nR)) {
      return(print("Error: Class proportions or number of cases missing."))
    }

    if (LPar["nC"] != length(prC)) {
      return(print("Error: Class proportions do not match nClass."))
    }

    # Assign respondents to classes
    cR <- rep(seq(1:length(prC)), rmultinom(1, nR, prC))
    cR <- sample(cR, size = nR)

    # Get probabilities for each respondent based on class
    prRI <- LPar[["prob"]][cR, ]
    genD <- cbind(cR, prRI)
  }

  if (LPar["model"] == 5) {
    nR <- LPar[["nR"]]
    prRI <- LPar[["prob"]]
    genD <- prRI
  }

  # Generate binary responses
  respRI <- matrix(sapply(c(prRI), rbinom, n = 1, size = 1), ncol = ncol(prRI))
  obsD <- respRI

  colnames(obsD) <- paste("i", seq(1:ncol(obsD)), sep = "")

  output[["obsData"]] <- obsD
  output[["genData"]] <- genD
  output["model"] <- LPar["model"]
  output["nC"] <- LPar["nC"]
  output["nI"] <- LPar["nI"]
  output["nR"] <- nR

  if (LPar["model"] < 5) {
    output[["cR"]] <- table(cR)
    output[["prC"]] <- prC
  }

  output[["logit"]] <- LPar[["logit"]]
  output[["prob"]] <- LPar[["prob"]]

  if (LPar["model"] >= 4) {
    output[["lC"]] <- LPar[["lC"]]
    output[["lI"]] <- LPar[["lI"]]
  }

  return(output)
}

# ============================================================
# VALIDATION FUNCTIONS
# ============================================================

#' Generate a single dataset using paper's method
#'
#' @param model Model code (0-5)
#' @param nR Number of respondents
#' @param nI Number of items
#' @param nC Number of classes
#' @param seed Random seed
#' @return List with data and parameters
generate_paper_data <- function(model, nR, nI, nC, seed) {
  set.seed(seed)

  model_names <- c("UN", "MON", "IIO", "DM", "LCR", "RM")

  if (model == 5) {
    # Rasch model: nC parameter becomes nR
    LPar <- sim.LParam(nC = nR, nI = nI, model = model)
    LOResp <- sim.OResp(LPar)
  } else {
    LPar <- sim.LParam(nC = nC, nI = nI, model = model)
    prC <- rep(1/nC, nC)
    LOResp <- sim.OResp(LPar, prC = prC, nR = nR)
  }

  list(
    data = LOResp$obsData,
    model_code = model,
    model_name = model_names[model + 1],
    n_classes = nC,
    prob = LOResp$prob,
    logit = LOResp$logit
  )
}

#' Fit all QuantFit models and return BIC values
#'
#' @param data Binary response matrix
#' @param n_classes Number of classes for discrete models
#' @param n_starts Number of random starts
#' @return Data frame with fit results
fit_all_models_paper <- function(data, n_classes, n_starts = 5) {

  results <- list()

  # UN
  res <- tryCatch(
    suppressWarnings(fit_un(data, n_classes, n_starts = n_starts, seed = 1)),
    error = function(e) NULL
  )
  results$UN <- if (!is.null(res)) list(ll = res$loglik, npar = res$n_par,
                                         bic = BIC(res), conv = res$convergence) else NULL

  # MON
  res <- tryCatch(
    suppressWarnings(fit_mon(data, n_classes, n_starts = n_starts, seed = 2)),
    error = function(e) NULL
  )
  results$MON <- if (!is.null(res)) list(ll = res$loglik, npar = res$n_par,
                                          bic = BIC(res), conv = res$convergence) else NULL

  # IIO
  res <- tryCatch(
    suppressWarnings(fit_iio(data, n_classes, n_starts = n_starts, seed = 3)),
    error = function(e) NULL
  )
  results$IIO <- if (!is.null(res)) list(ll = res$loglik, npar = res$n_par,
                                          bic = BIC(res), conv = res$convergence) else NULL

  # DM
  res <- tryCatch(
    suppressWarnings(fit_dm(data, n_classes, n_starts = n_starts, seed = 4)),
    error = function(e) NULL
  )
  results$DM <- if (!is.null(res)) list(ll = res$loglik, npar = res$n_par,
                                         bic = BIC(res), conv = res$convergence) else NULL

  # LCR
  res <- tryCatch(
    suppressWarnings(fit_lcr(data, n_classes, n_starts = n_starts, seed = 5)),
    error = function(e) NULL
  )
  results$LCR <- if (!is.null(res)) list(ll = res$loglik, npar = res$n_par,
                                          bic = BIC(res), conv = res$convergence) else NULL

  # RM
  res <- tryCatch(
    suppressWarnings(fit_rm(data, verbose = FALSE)),
    error = function(e) NULL
  )
  results$RM <- if (!is.null(res)) list(ll = res$loglik, npar = res$n_par,
                                         bic = BIC(res), conv = res$convergence) else NULL

  results
}

#' Run validation study using paper's exact simulation
#'
#' @param n_datasets Total number of datasets
#' @param n_persons Number of persons per dataset
#' @param n_items Number of items
#' @param n_classes Number of classes for discrete models
#' @param n_starts Random starts for model fitting
#' @param verbose Print progress
#' @return List with results
run_exact_paper_validation <- function(n_datasets = 30,
                                        n_persons = 1000,
                                        n_items = 10,
                                        n_classes = 4,
                                        n_starts = 5,
                                        verbose = TRUE) {

  model_names <- c("UN", "MON", "IIO", "DM", "LCR", "RM")

  if (verbose) {
    cat("=" |> rep(70) |> paste(collapse = ""), "\n")
    cat("QuantFit Validation - Using Paper's Exact Simulation Method\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")
    cat("Settings:\n")
    cat("  Datasets per model:", ceiling(n_datasets/6), "\n")
    cat("  Persons:", n_persons, "\n")
    cat("  Items:", n_items, "\n")
    cat("  Classes:", n_classes, "\n\n")
  }

  all_results <- data.frame()
  datasets_per_model <- ceiling(n_datasets / 6)

  for (model_code in 0:5) {
    model_name <- model_names[model_code + 1]

    if (verbose) {
      cat("\n--- Model", model_code, "(", model_name, ") ---\n")
    }

    for (rep in 1:datasets_per_model) {
      seed <- model_code * 1000 + rep * 100

      # Generate data using paper's exact method
      sim <- generate_paper_data(model_code, n_persons, n_items, n_classes, seed)

      # Fit all models
      fits <- fit_all_models_paper(sim$data, n_classes, n_starts)

      # Get BIC values
      bics <- sapply(model_names, function(m) {
        if (!is.null(fits[[m]])) fits[[m]]$bic else NA
      })

      # Select best model
      selected <- names(which.min(bics))
      correct <- (selected == model_name)

      # Store result
      result <- data.frame(
        gen_model = model_name,
        gen_code = model_code,
        rep = rep,
        seed = seed,
        sel_model = selected,
        correct = correct,
        BIC_UN = bics["UN"],
        BIC_MON = bics["MON"],
        BIC_IIO = bics["IIO"],
        BIC_DM = bics["DM"],
        BIC_LCR = bics["LCR"],
        BIC_RM = bics["RM"],
        stringsAsFactors = FALSE
      )
      all_results <- rbind(all_results, result)

      if (verbose) {
        marker <- if (correct) "CORRECT" else ""
        cat(sprintf("  Rep %d: Sel=%s %s\n", rep, selected, marker))
      }
    }
  }

  # Summary
  if (verbose) {
    cat("\n\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n")
    cat("RESULTS SUMMARY\n")
    cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

    # Confusion matrix
    cat("Selection Matrix (rows=generating, cols=selected):\n\n")
    conf <- table(
      factor(all_results$gen_model, levels = model_names),
      factor(all_results$sel_model, levels = model_names)
    )
    print(conf)

    # Accuracy
    cat("\n\nAccuracy by Model:\n")
    cat("-" |> rep(40) |> paste(collapse = ""), "\n")
    for (m in model_names) {
      sub <- all_results[all_results$gen_model == m, ]
      n_correct <- sum(sub$correct)
      n_total <- nrow(sub)
      if (n_total > 0) {
        cat(sprintf("  %s: %d/%d (%.1f%%)\n", m, n_correct, n_total,
            100 * n_correct / n_total))
      }
    }

    total <- sum(all_results$correct)
    n <- nrow(all_results)
    cat(sprintf("\n  OVERALL: %d/%d (%.1f%%)\n", total, n, 100 * total / n))

    # Scale-level (Quantitative = models 4,5; Non-Quant = models 0-3)
    cat("\n\nScale-Level Classification:\n")
    cat("-" |> rep(40) |> paste(collapse = ""), "\n")
    all_results$gen_scale <- ifelse(all_results$gen_code >= 4, "Quant", "Non-Quant")
    all_results$sel_scale <- ifelse(all_results$sel_model %in% c("LCR", "RM"), "Quant", "Non-Quant")
    scale_conf <- table(all_results$gen_scale, all_results$sel_scale)
    print(scale_conf)
    scale_acc <- sum(diag(scale_conf)) / sum(scale_conf)
    cat(sprintf("\nScale accuracy: %.1f%%\n", 100 * scale_acc))
  }

  invisible(list(
    results = all_results,
    confusion = conf
  ))
}

#' Quick test with paper's method
quick_paper_test <- function() {
  run_exact_paper_validation(
    n_datasets = 12,
    n_persons = 1000,
    n_items = 10,
    n_classes = 4,
    n_starts = 5,
    verbose = TRUE
  )
}

# ============================================================
# DIAGNOSTIC: Check data generation matches paper
# ============================================================

#' Diagnostic to verify data generation
diagnose_paper_generation <- function() {
  cat("=" |> rep(60) |> paste(collapse = ""), "\n")
  cat("Diagnostic: Paper Data Generation Verification\n")
  cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

  model_names <- c("UN", "MON", "IIO", "DM", "LCR", "RM")

  for (model in 0:5) {
    cat("\n--- Model", model, "(", model_names[model + 1], ") ---\n")

    set.seed(12345)

    if (model == 5) {
      LPar <- sim.LParam(nC = 100, nI = 10, model = model)
    } else {
      LPar <- sim.LParam(nC = 4, nI = 10, model = model)
    }

    prob <- LPar$prob

    if (model < 5) {
      cat("Probability matrix (classes x items):\n")
      print(round(prob, 2))

      # Check class monotonicity (each column should be non-decreasing)
      class_mon <- all(apply(prob, 2, function(x) all(diff(x) >= -0.001)))
      cat("Class monotonicity satisfied:", class_mon, "\n")

      # Check item ordering (each row should be non-decreasing)
      item_ord <- all(apply(prob, 1, function(x) all(diff(x) >= -0.001)))
      cat("Item ordering satisfied:", item_ord, "\n")
    } else {
      cat("Rasch model - continuous theta\n")
      cat("Theta range:", round(range(LPar$lC), 2), "\n")
      cat("Delta range:", round(range(LPar$lI), 2), "\n")
    }
  }
}

# ============================================================
# USAGE
# ============================================================

cat("\n")
cat("Paper-Exact Validation Script Loaded\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n")
cat("Functions:\n")
cat("  quick_paper_test()            - Quick test (12 datasets)\n")
cat("  run_exact_paper_validation()  - Full validation\n")
cat("  diagnose_paper_generation()   - Check data generation\n")
cat("\n")
