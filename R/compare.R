#' Model Comparison Framework for QuantFit
#'
#' @name compare
#' @description Functions for comparing latent structure models following the
#'   Torres Irribarra & Diakow framework.
NULL

#' Compare Multiple Latent Structure Models
#'
#' @param data A matrix or data frame of binary responses (0/1)
#' @param n_classes Number of latent classes for discrete models (UN, MON, IIO, DM, LCR)
#' @param models Character vector of models to fit. Default is all six:
#'   c("UN", "MON", "IIO", "DM", "LCR", "RM")
#' @param item_order Item order for IIO and DM models. If NULL, estimated from data.
#' @param n_starts Number of random starts for each model (default 10)
#' @param verbose Print progress messages (default FALSE)
#' @param ... Additional arguments passed to individual fitting functions
#'
#' @return A qlcompare object containing:
#' \describe{
#'   \item{fits}{List of fitted model objects}
#'   \item{comparison_table}{Data frame with fit statistics for all models}
#'   \item{best_model}{Name of best model by BIC}
#' }
#'
#' @details
#' This function fits all specified models and returns a comparison table with
#' fit statistics including log-likelihood, number of parameters, AIC, BIC, and SABIC.
#'
#' The six models represent different levels of structure:
#' \itemize{
#'   \item \strong{UN}: No ordering (classificatory)
#'   \item \strong{MON}: Classes are ordered
#'   \item \strong{IIO}: Items maintain same relative difficulty across classes
#'   \item \strong{DM}: Both MON and IIO (strong ordinal)
#'   \item \strong{LCR}: Rasch parameterization with discrete classes (quantitative)
#'   \item \strong{RM}: Continuous Rasch (fully quantitative)
#' }
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 500
#' data <- matrix(rbinom(n * 10, 1, 0.5), nrow = n)
#'
#' # Compare all models
#' comparison <- compare_models(data, n_classes = 3)
#' print(comparison)
#' summary(comparison)
#' plot(comparison)
#'
#' # Compare only discrete models
#' comparison2 <- compare_models(data, n_classes = 3,
#'                                models = c("UN", "MON", "DM", "LCR"))
#' }
#'
#' @seealso \code{\link{successive_comparison}} for stepwise comparison
#'
#' @export
compare_models <- function(data, n_classes,
                           models = c("UN", "MON", "IIO", "DM", "LCR", "RM"),
                           item_order = NULL,
                           n_starts = 10,
                           verbose = FALSE,
                           ...) {

  # Capture call
  call <- match.call()

  # Validate inputs
  data <- validate_data(data)

  valid_models <- c("UN", "MON", "IIO", "DM", "LCR", "RM")
  models <- match.arg(models, valid_models, several.ok = TRUE)

  n_obs <- nrow(data)
  n_items <- ncol(data)

  # Estimate item order if needed and not provided
  if ((any(c("IIO", "DM") %in% models)) && is.null(item_order)) {
    item_order <- estimate_item_order(data)
    if (verbose) {
      cat("Estimated item order (easiest to hardest):", item_order, "\n\n")
    }
  }

  # Fit each model
  fits <- list()
  results <- list()

  for (model in models) {
    if (verbose) cat("Fitting", model, "model...\n")

    fit <- tryCatch({
      switch(model,
        UN = fit_un(data, n_classes, n_starts = n_starts, verbose = verbose > 1, ...),
        MON = fit_mon(data, n_classes, n_starts = n_starts, verbose = verbose > 1, ...),
        IIO = fit_iio(data, n_classes, item_order = item_order,
                      n_starts = n_starts, verbose = verbose > 1, ...),
        DM = fit_dm(data, n_classes, item_order = item_order,
                    n_starts = n_starts, verbose = verbose > 1, ...),
        LCR = fit_lcr(data, n_classes, n_starts = n_starts, verbose = verbose > 1, ...),
        RM = fit_rm(data, verbose = verbose > 1, ...)
      )
    }, error = function(e) {
      warning("Failed to fit ", model, " model: ", e$message)
      NULL
    })

    if (!is.null(fit)) {
      fits[[model]] <- fit

      results[[model]] <- data.frame(
        Model = model,
        LogLik = fit$loglik,
        nPar = fit$n_par,
        AIC = AIC(fit),
        BIC = BIC(fit),
        SABIC = sabic(fit),
        Converged = fit$convergence,
        stringsAsFactors = FALSE
      )

      if (verbose) {
        cat("  LogLik:", round(fit$loglik, 2),
            " AIC:", round(AIC(fit), 2),
            " BIC:", round(BIC(fit), 2), "\n\n")
      }
    }
  }

  if (length(fits) == 0) {
    stop("All model fits failed")
  }

  # Combine results
  comparison_table <- do.call(rbind, results)
  rownames(comparison_table) <- NULL

  # Sort by BIC
  comparison_table <- comparison_table[order(comparison_table$BIC), ]

  # Identify best model
  best_model <- comparison_table$Model[1]

  # Create qlcompare object
  result <- new_qlcompare(
    fits = fits,
    comparison_table = comparison_table,
    best_model = best_model,
    call = call
  )

  result
}

#' Successive Model Comparison Strategy
#'
#' Implements the stepwise comparison strategy from Torres Irribarra & Diakow,
#' progressively testing whether stronger structure is supported by the data.
#'
#' @param data A matrix or data frame of binary responses (0/1)
#' @param n_classes Number of latent classes for discrete models
#' @param item_order Item order for IIO and DM models. If NULL, estimated from data.
#' @param n_starts Number of random starts for each model (default 10)
#' @param bic_threshold BIC difference threshold for preferring simpler model (default 2)
#' @param verbose Print progress and interpretation (default TRUE)
#' @param ... Additional arguments passed to fitting functions
#'
#' @return A list containing:
#' \describe{
#'   \item{comparison}{A qlcompare object with all fitted models}
#'   \item{steps}{Data frame with step-by-step comparison results}
#'   \item{conclusion}{Character string describing the supported structure}
#'   \item{best_model}{Name of best model overall}
#' }
#'
#' @details
#' The successive comparison follows these steps:
#' \enumerate{
#'   \item \strong{Step 1 - Ordering}: Compare UN vs MON and UN vs IIO.
#'     If constrained models fit as well or better, there is evidence for ordering.
#'   \item \strong{Step 2 - Unidimensionality}: Compare MON+IIO vs DM.
#'     If DM fits as well, there is evidence for double monotonicity/unidimensionality.
#'   \item \strong{Step 3 - Interval Scale}: Compare DM vs LCR.
#'     If LCR fits as well, there is evidence for interval-level measurement.
#'   \item \strong{Step 4 - Continuity}: Compare LCR vs RM.
#'     Tests whether discrete or continuous latent trait is preferred.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' data <- matrix(rbinom(500 * 10, 1, 0.5), nrow = 500)
#'
#' result <- successive_comparison(data, n_classes = 4)
#' print(result$conclusion)
#' }
#'
#' @export
successive_comparison <- function(data, n_classes,
                                  item_order = NULL,
                                  n_starts = 10,
                                  bic_threshold = 2,
                                  verbose = TRUE,
                                  ...) {

  # Validate inputs
  data <- validate_data(data)

  n_obs <- nrow(data)
  n_items <- ncol(data)

  # Estimate item order if needed
  if (is.null(item_order)) {
    item_order <- estimate_item_order(data)
    if (verbose) {
      cat("Estimated item order:", item_order, "\n\n")
    }
  }

  # Fit all models
  if (verbose) cat("Fitting all models...\n\n")

  comparison <- compare_models(
    data, n_classes,
    models = c("UN", "MON", "IIO", "DM", "LCR", "RM"),
    item_order = item_order,
    n_starts = n_starts,
    verbose = verbose,
    ...
  )

  # Extract BIC values
  fits <- comparison$fits
  get_bic <- function(m) if (!is.null(fits[[m]])) BIC(fits[[m]]) else NA

  bic_un <- get_bic("UN")
  bic_mon <- get_bic("MON")
  bic_iio <- get_bic("IIO")
  bic_dm <- get_bic("DM")
  bic_lcr <- get_bic("LCR")
  bic_rm <- get_bic("RM")

  # Successive comparisons
  steps <- list()

  # Step 1: Is there ordering?
  if (verbose) cat("\n=== Step 1: Test for Ordering ===\n")

  # UN vs MON
  delta_un_mon <- bic_un - bic_mon
  mon_preferred <- !is.na(delta_un_mon) && delta_un_mon > -bic_threshold

  # UN vs IIO
  delta_un_iio <- bic_un - bic_iio
  iio_preferred <- !is.na(delta_un_iio) && delta_un_iio > -bic_threshold

  ordering_supported <- mon_preferred || iio_preferred

  steps$step1 <- data.frame(
    Step = 1,
    Question = "Is there ordering?",
    Comparison = "UN vs MON, UN vs IIO",
    delta_BIC_MON = round(delta_un_mon, 2),
    delta_BIC_IIO = round(delta_un_iio, 2),
    Result = if (ordering_supported) "Yes - ordering supported" else "No - keep UN"
  )

  if (verbose) {
    cat("  UN BIC:", round(bic_un, 2), "\n")
    cat("  MON BIC:", round(bic_mon, 2), " (delta:", round(delta_un_mon, 2), ")\n")
    cat("  IIO BIC:", round(bic_iio, 2), " (delta:", round(delta_un_iio, 2), ")\n")
    cat("  --> ", steps$step1$Result, "\n")
  }

  # Step 2: Is there unidimensionality (double monotonicity)?
  if (verbose) cat("\n=== Step 2: Test for Unidimensionality ===\n")

  if (ordering_supported) {
    # Compare best of MON/IIO to DM
    best_ordered_bic <- min(bic_mon, bic_iio, na.rm = TRUE)
    delta_ordered_dm <- best_ordered_bic - bic_dm
    dm_preferred <- !is.na(delta_ordered_dm) && delta_ordered_dm > -bic_threshold

    steps$step2 <- data.frame(
      Step = 2,
      Question = "Is there unidimensionality?",
      Comparison = "Best(MON,IIO) vs DM",
      delta_BIC = round(delta_ordered_dm, 2),
      Result = if (dm_preferred) "Yes - double monotonicity supported" else "No - keep separate MON/IIO"
    )

    if (verbose) {
      cat("  Best ordered BIC:", round(best_ordered_bic, 2), "\n")
      cat("  DM BIC:", round(bic_dm, 2), " (delta:", round(delta_ordered_dm, 2), ")\n")
      cat("  --> ", steps$step2$Result, "\n")
    }
  } else {
    dm_preferred <- FALSE
    steps$step2 <- data.frame(
      Step = 2,
      Question = "Is there unidimensionality?",
      Comparison = "Skipped (no ordering)",
      delta_BIC = NA,
      Result = "Skipped"
    )
    if (verbose) cat("  Skipped (no ordering evidence)\n")
  }

  # Step 3: Is there interval-level measurement?
  if (verbose) cat("\n=== Step 3: Test for Interval Scale ===\n")

  # ALWAYS check LCR against the best non-quantitative model
  # This ensures we don't miss quantitative structure when ordinal steps fail
  if (dm_preferred) {
    compare_bic <- bic_dm
    compare_label <- "DM"
  } else if (ordering_supported) {
    compare_bic <- min(bic_mon, bic_iio, na.rm = TRUE)
    compare_label <- "Best(MON,IIO)"
  } else {
    compare_bic <- bic_un
    compare_label <- "UN"
  }

  delta_to_lcr <- compare_bic - bic_lcr
  lcr_preferred <- !is.na(delta_to_lcr) && delta_to_lcr > bic_threshold

  steps$step3 <- data.frame(
    Step = 3,
    Question = "Is there interval-level measurement?",
    Comparison = paste(compare_label, "vs LCR"),
    delta_BIC = round(delta_to_lcr, 2),
    Result = if (lcr_preferred) "Yes - Rasch structure supported" else "No - keep ordinal/nominal"
  )

  if (verbose) {
    cat("  ", compare_label, " BIC:", round(compare_bic, 2), "\n")
    cat("  LCR BIC:", round(bic_lcr, 2), " (delta:", round(delta_to_lcr, 2), ")\n")
    cat("  --> ", steps$step3$Result, "\n")
  }

  # Step 4: Is continuous theta preferred?
  if (verbose) cat("\n=== Step 4: Test for Continuous Theta ===\n")

  # Always compare LCR vs RM when quantitative structure is possible
  delta_lcr_rm <- bic_lcr - bic_rm
  rm_preferred <- !is.na(delta_lcr_rm) && delta_lcr_rm > bic_threshold

  # Also check if RM beats the best non-quantitative model directly
  rm_beats_ordinal <- !is.na(bic_rm) && bic_rm < compare_bic - bic_threshold

  steps$step4 <- data.frame(
    Step = 4,
    Question = "Is continuous theta preferred?",
    Comparison = "LCR vs RM",
    delta_BIC = round(delta_lcr_rm, 2),
    Result = if (rm_preferred) "Yes - continuous model preferred" else "No - discrete classes sufficient"
  )

  if (verbose) {
    cat("  LCR BIC:", round(bic_lcr, 2), "\n")
    cat("  RM BIC:", round(bic_rm, 2), " (delta:", round(delta_lcr_rm, 2), ")\n")
    cat("  --> ", steps$step4$Result, "\n")
  }

  # Determine conclusion based on the overall best model by BIC
  # This ensures we pick the right model even when hierarchical steps give mixed signals
  all_bics <- c(UN = bic_un, MON = bic_mon, IIO = bic_iio,
                DM = bic_dm, LCR = bic_lcr, RM = bic_rm)
  overall_best <- names(which.min(all_bics))

  # Generate conclusion based on overall best model
  if (overall_best == "UN") {
    conclusion <- "CLASSIFICATORY: Data support a qualitative/nominal latent structure (UN model preferred)"
    best_model <- "UN"
  } else if (overall_best == "MON") {
    conclusion <- "ORDINAL (classes ordered): MON model preferred - classes are ordered but items do not maintain ordering"
    best_model <- "MON"
  } else if (overall_best == "IIO") {
    conclusion <- "ORDINAL (items ordered): IIO model preferred - items maintain difficulty ordering"
    best_model <- "IIO"
  } else if (overall_best == "DM") {
    conclusion <- "ORDINAL (strong): Double monotonicity supported - unidimensional ordinal scale"
    best_model <- "DM"
  } else if (overall_best == "LCR") {
    conclusion <- "QUANTITATIVE (discrete): Latent Class Rasch model preferred - interval scale with discrete levels"
    best_model <- "LCR"
  } else {
    conclusion <- "QUANTITATIVE (continuous): Full Rasch model preferred - interval scale with continuous theta"
    best_model <- "RM"
  }

  # Compile steps summary (simplified - just step/question/result)
  steps_df <- data.frame(
    Step = sapply(steps, function(s) s$Step),
    Question = sapply(steps, function(s) s$Question),
    Result = sapply(steps, function(s) s$Result),
    stringsAsFactors = FALSE
  )
  rownames(steps_df) <- NULL

  if (verbose) {
    cat("\n========================================\n")
    cat("CONCLUSION:", conclusion, "\n")
    cat("Best model:", best_model, "\n")
    cat("========================================\n")
  }

  list(
    comparison = comparison,
    steps = steps_df,
    conclusion = conclusion,
    best_model = best_model,
    interpretation = interpret_structure(best_model)
  )
}

#' Interpret latent structure based on best model
#'
#' @param model Model code (UN, MON, IIO, DM, LCR, RM)
#'
#' @return Character string with interpretation
#' @keywords internal
interpret_structure <- function(model) {
  interpretations <- list(
    UN = paste(
      "The Unconstrained model suggests a CLASSIFICATORY interpretation.",
      "The latent variable represents distinct categories without inherent ordering.",
      "Classes differ qualitatively rather than quantitatively.",
      "Examples: diagnostic categories, learning styles, personality types."
    ),
    MON = paste(
      "The Class Monotonicity model suggests an ORDINAL interpretation.",
      "Latent classes can be ordered from low to high.",
      "Higher classes have uniformly higher item response probabilities.",
      "However, items may not maintain consistent relative difficulty across classes."
    ),
    IIO = paste(
      "The Invariant Item Ordering model suggests an ORDINAL interpretation.",
      "Items maintain the same difficulty ordering across all latent classes.",
      "This is consistent with a single underlying dimension.",
      "However, class ordering is not constrained."
    ),
    DM = paste(
      "The Double Monotonicity model suggests a strong ORDINAL interpretation.",
      "Both class ordering and item ordering are preserved.",
      "This is consistent with a unidimensional construct.",
      "The scale is ordinal - distances between classes are not necessarily equal."
    ),
    LCR = paste(
      "The Latent Class Rasch model suggests a QUANTITATIVE interpretation.",
      "The Rasch parameterization implies interval-level measurement.",
      "Differences between class locations (theta) are meaningful.",
      "The discrete classes may reflect a naturally categorical construct",
      "or may be a sufficient approximation to a continuous trait."
    ),
    RM = paste(
      "The Rasch Model suggests a fully QUANTITATIVE interpretation.",
      "The latent trait is continuous and interval-scaled.",
      "Person ability differences are directly comparable.",
      "This supports the strongest measurement claims."
    )
  )

  interpretations[[model]]
}

#' Plot comparison of model fit indices
#'
#' @param x A qlcompare object
#' @param criterion Which criterion to plot: "BIC" (default), "AIC", or "both"
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the plot data
#' @export
plot_comparison <- function(x, criterion = c("BIC", "AIC", "both"), ...) {
  if (!inherits(x, "qlcompare")) {
    stop("x must be a qlcompare object")
  }

  criterion <- match.arg(criterion)

  tbl <- x$comparison_table

  # Base R plotting
  if (criterion == "both") {
    par(mfrow = c(1, 2))

    # AIC plot
    barplot(tbl$AIC, names.arg = tbl$Model,
            main = "Model Comparison: AIC",
            ylab = "AIC", col = "steelblue",
            las = 2)
    abline(h = min(tbl$AIC), lty = 2, col = "red")

    # BIC plot
    barplot(tbl$BIC, names.arg = tbl$Model,
            main = "Model Comparison: BIC",
            ylab = "BIC", col = "steelblue",
            las = 2)
    abline(h = min(tbl$BIC), lty = 2, col = "red")

    par(mfrow = c(1, 1))

  } else {
    values <- if (criterion == "BIC") tbl$BIC else tbl$AIC

    barplot(values, names.arg = tbl$Model,
            main = paste("Model Comparison:", criterion),
            ylab = criterion,
            col = ifelse(values == min(values), "darkgreen", "steelblue"),
            las = 2)
    abline(h = min(values), lty = 2, col = "red")

    # Add delta values as text
    delta_values <- values - min(values)
    text(seq(0.7, by = 1.2, length.out = nrow(tbl)),
         values + diff(range(values)) * 0.05,
         labels = paste0("Δ=", round(delta_values, 1)),
         cex = 0.8)
  }

  invisible(tbl)
}
