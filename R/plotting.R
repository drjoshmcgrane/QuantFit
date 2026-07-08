#' Visualization Functions for QuantFit
#'
#' @name plotting
#' @description Functions for visualizing latent structure models and comparisons
NULL

#' Plot Item Response Functions (IRFs)
#'
#' Creates a plot showing item response probabilities across latent classes,
#' similar to Figure 3.3 in Torres Irribarra & Diakow.
#'
#' @param object A qlfit object
#' @param items Which items to plot. Default is all items. Can be numeric indices
#'   or item names.
#' @param log_odds If TRUE, plot log-odds instead of probabilities (default FALSE)
#' @param legend_pos Position of legend: "right" (default), "bottom", "none"
#' @param colors Optional vector of colors for items
#' @param main Optional title for the plot
#' @param ... Additional arguments passed to plot()
#'
#' @return Invisibly returns the plot data
#'
#' @examples
#' \dontrun{
#' fit <- fit_mon(data, n_classes = 4)
#' plot_irfs(fit)
#' plot_irfs(fit, log_odds = TRUE)
#' plot_irfs(fit, items = 1:5)
#' }
#'
#' @export
plot_irfs <- function(object, items = NULL, log_odds = FALSE,
                      legend_pos = c("right", "bottom", "none"),
                      colors = NULL, main = NULL, ...) {

  if (!inherits(object, "qlfit")) {
    stop("object must be a qlfit object")
  }

  legend_pos <- match.arg(legend_pos)

  # Handle RM model differently (continuous theta)
  if (object$model_type == "RM") {
    return(plot_irfs_rm(object, items, legend_pos, colors, main, ...))
  }

  # Extract item probabilities
  item_probs <- object$item_probs
  n_items <- nrow(item_probs)
  n_classes <- ncol(item_probs)

  # Select items to plot
  if (is.null(items)) {
    items <- 1:n_items
  } else if (is.character(items)) {
    items <- match(items, rownames(item_probs))
    items <- items[!is.na(items)]
  }

  if (length(items) == 0) {
    stop("No valid items to plot")
  }

  # Subset item probabilities
  plot_probs <- item_probs[items, , drop = FALSE]
  n_plot_items <- nrow(plot_probs)

  # Convert to log-odds if requested
  if (log_odds) {
    plot_values <- log(plot_probs / (1 - plot_probs))
    ylab <- "Log-odds"
  } else {
    plot_values <- plot_probs
    ylab <- "P(X = 1)"
  }

  # Set up colors
  if (is.null(colors)) {
    colors <- rainbow(n_plot_items, s = 0.7, v = 0.8)
  }

  # Set up plot title
  if (is.null(main)) {
    model_names <- c(
      UN = "Unconstrained",
      MON = "Class Monotonicity",
      IIO = "Item Ordering",
      DM = "Double Monotonicity",
      LCR = "Latent Class Rasch"
    )
    main <- paste("Item Response Functions:", model_names[object$model_type])
  }

  # Create plot
  x_vals <- 1:n_classes

  # Set up margins for legend
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (legend_pos == "right") {
    par(mar = c(5, 4, 4, 8), xpd = TRUE)
  } else if (legend_pos == "bottom") {
    par(mar = c(8, 4, 4, 2))
  }

  # Initialize plot
  y_range <- range(plot_values, na.rm = TRUE)
  if (!log_odds) y_range <- c(0, 1)

  plot(x_vals, plot_values[1, ], type = "n",
       xlim = c(0.5, n_classes + 0.5),
       ylim = y_range,
       xlab = "Latent Class", ylab = ylab,
       main = main,
       xaxt = "n", ...)

  axis(1, at = x_vals, labels = paste0("C", x_vals))

  # Add grid
  abline(h = if (log_odds) 0 else 0.5, lty = 2, col = "gray70")
  grid(nx = NA, ny = NULL, col = "gray90")

  # Plot each item
  for (i in 1:n_plot_items) {
    lines(x_vals, plot_values[i, ], col = colors[i], lwd = 2)
    points(x_vals, plot_values[i, ], col = colors[i], pch = 19, cex = 1.2)
  }

  # Add legend
  item_labels <- if (!is.null(rownames(plot_probs))) {
    rownames(plot_probs)
  } else {
    paste0("Item ", items)
  }

  if (legend_pos == "right") {
    legend("topright", inset = c(-0.25, 0),
           legend = item_labels, col = colors,
           lwd = 2, pch = 19, cex = 0.8, bty = "n")
  } else if (legend_pos == "bottom") {
    legend("bottom", inset = c(0, -0.35),
           legend = item_labels, col = colors,
           lwd = 2, pch = 19, cex = 0.7,
           horiz = TRUE, bty = "n", xpd = TRUE)
  }

  invisible(data.frame(
    item = rep(items, each = n_classes),
    class = rep(1:n_classes, n_plot_items),
    value = as.vector(t(plot_values))
  ))
}

#' Plot IRFs for Rasch Model (continuous theta)
#'
#' @keywords internal
plot_irfs_rm <- function(object, items = NULL, legend_pos = "right",
                         colors = NULL, main = NULL, ...) {

  delta <- object$delta
  n_items <- length(delta)

  # Select items
  if (is.null(items)) items <- 1:n_items

  # Theta range
  theta <- seq(-4, 4, length.out = 100)

  # Compute probabilities
  plot_probs <- matrix(0, nrow = length(items), ncol = length(theta))
  for (i in seq_along(items)) {
    plot_probs[i, ] <- inv_logit(theta - delta[items[i]])
  }

  # Colors
  if (is.null(colors)) {
    colors <- rainbow(length(items), s = 0.7, v = 0.8)
  }

  # Title
  if (is.null(main)) main <- "Item Characteristic Curves: Rasch Model"

  # Plot
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (legend_pos == "right") {
    par(mar = c(5, 4, 4, 8), xpd = TRUE)
  }

  plot(theta, plot_probs[1, ], type = "n",
       xlim = c(-4, 4), ylim = c(0, 1),
       xlab = expression(theta), ylab = "P(X = 1)",
       main = main, ...)

  abline(h = 0.5, lty = 2, col = "gray70")
  grid(col = "gray90")

  for (i in seq_along(items)) {
    lines(theta, plot_probs[i, ], col = colors[i], lwd = 2)
    # Mark item difficulty
    points(delta[items[i]], 0.5, col = colors[i], pch = "|", cex = 1.5)
  }

  # Legend
  item_labels <- if (!is.null(names(delta))) {
    names(delta)[items]
  } else {
    paste0("Item ", items)
  }

  if (legend_pos == "right") {
    legend("topright", inset = c(-0.25, 0),
           legend = item_labels, col = colors,
           lwd = 2, cex = 0.8, bty = "n")
  }

  invisible(plot_probs)
}

#' Plot Class Profiles
#'
#' Shows the probability profile for each latent class across items.
#'
#' @param object A qlfit object
#' @param classes Which classes to plot (default all)
#' @param colors Optional vector of colors for classes
#' @param main Optional title
#' @param ... Additional arguments passed to plot()
#'
#' @return Invisibly returns the plot data
#' @export
plot_class_profiles <- function(object, classes = NULL, colors = NULL,
                                main = NULL, ...) {

  if (!inherits(object, "qlfit")) {
    stop("object must be a qlfit object")
  }

  if (object$model_type == "RM") {
    stop("Class profiles not available for continuous Rasch model")
  }

  item_probs <- object$item_probs
  n_items <- nrow(item_probs)
  n_classes <- ncol(item_probs)

  # Select classes
  if (is.null(classes)) classes <- 1:n_classes

  # Colors
  if (is.null(colors)) {
    colors <- rainbow(length(classes), s = 0.7, v = 0.8)
  }

  # Title
  if (is.null(main)) {
    main <- paste("Class Response Profiles:", object$model_type)
  }

  # Plot setup
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 4, 4, 8), xpd = TRUE)

  x_vals <- 1:n_items

  plot(x_vals, item_probs[, classes[1]], type = "n",
       xlim = c(0.5, n_items + 0.5), ylim = c(0, 1),
       xlab = "Item", ylab = "P(X = 1)",
       main = main, xaxt = "n", ...)

  axis(1, at = x_vals)
  grid(col = "gray90")
  abline(h = 0.5, lty = 2, col = "gray70")

  for (i in seq_along(classes)) {
    c <- classes[i]
    lines(x_vals, item_probs[, c], col = colors[i], lwd = 2)
    points(x_vals, item_probs[, c], col = colors[i], pch = 19)
  }

  legend("topright", inset = c(-0.2, 0),
         legend = paste0("Class ", classes, " (",
                        round(object$class_probs[classes] * 100), "%)"),
         col = colors, lwd = 2, pch = 19, cex = 0.8, bty = "n")

  invisible(item_probs[, classes])
}

#' Plot Posterior Class Membership
#'
#' Visualizes the distribution of posterior class memberships.
#'
#' @param object A qlfit object
#' @param type Type of plot: "histogram" (default), "stacked", or "heatmap"
#' @param ... Additional arguments
#'
#' @return Invisibly returns posterior data
#' @export
plot_posteriors <- function(object, type = c("histogram", "stacked", "heatmap"), ...) {
  if (!inherits(object, "qlfit")) {
    stop("object must be a qlfit object")
  }

  if (is.null(object$posteriors)) {
    stop("Posteriors not available for this model")
  }

  type <- match.arg(type)
  posteriors <- object$posteriors
  n_classes <- ncol(posteriors)

  if (type == "histogram") {
    # Histogram of modal class assignments
    assignments <- apply(posteriors, 1, which.max)

    barplot(table(factor(assignments, levels = 1:n_classes)),
            main = "Modal Class Assignments",
            xlab = "Latent Class",
            ylab = "Frequency",
            col = rainbow(n_classes, s = 0.7, v = 0.8),
            names.arg = paste0("Class ", 1:n_classes))

  } else if (type == "stacked") {
    # Stacked bar for each observation (useful for small n)
    n_obs <- nrow(posteriors)
    if (n_obs > 100) {
      # Sample for visibility
      idx <- sample(n_obs, 100)
      posteriors <- posteriors[idx, ]
    }

    barplot(t(posteriors),
            main = "Posterior Class Memberships",
            xlab = "Observation",
            ylab = "Probability",
            col = rainbow(n_classes, s = 0.7, v = 0.8),
            border = NA,
            space = 0)

    legend("topright",
           legend = paste0("Class ", 1:n_classes),
           fill = rainbow(n_classes, s = 0.7, v = 0.8),
           cex = 0.8, bty = "n")

  } else if (type == "heatmap") {
    # Sort observations by modal class
    assignments <- apply(posteriors, 1, which.max)
    ord <- order(assignments, -apply(posteriors, 1, max))
    posteriors_sorted <- posteriors[ord, ]

    image(t(posteriors_sorted),
          col = colorRampPalette(c("white", "darkblue"))(100),
          main = "Posterior Class Membership Heatmap",
          xlab = "Latent Class",
          ylab = "Observation (sorted)",
          xaxt = "n", yaxt = "n")

    axis(1, at = seq(0, 1, length.out = n_classes),
         labels = paste0("C", 1:n_classes))
  }

  invisible(posteriors)
}

#' Plot Model Comparison
#'
#' Visualizes comparison of multiple models' fit statistics.
#'
#' @param object A qlcompare object or list of qlfit objects
#' @param criterion Which criterion to highlight: "BIC" (default), "AIC", or "SABIC"
#' @param ... Additional arguments
#'
#' @return Invisibly returns comparison data
#' @export
plot_model_comparison <- function(object, criterion = c("BIC", "AIC", "SABIC"), ...) {
  criterion <- match.arg(criterion)

  if (inherits(object, "qlcompare")) {
    tbl <- object$comparison_table
  } else {
    stop("object must be a qlcompare object")
  }

  n_models <- nrow(tbl)

  # Set up multi-panel plot
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

  # Panel 1: Log-likelihood and parameters
  barplot(
    rbind(tbl$LogLik / max(abs(tbl$LogLik)), tbl$nPar / max(tbl$nPar)),
    beside = TRUE,
    names.arg = tbl$Model,
    col = c("steelblue", "coral"),
    main = "Log-Likelihood and Parameters",
    ylab = "Scaled Value",
    las = 2
  )
  legend("topright",
         legend = c("Log-Lik (scaled)", "nPar (scaled)"),
         fill = c("steelblue", "coral"),
         cex = 0.8, bty = "n")

  # Panel 2: Information criteria
  ic_values <- switch(criterion,
    BIC = tbl$BIC,
    AIC = tbl$AIC,
    SABIC = tbl$SABIC
  )

  colors <- rep("steelblue", n_models)
  colors[which.min(ic_values)] <- "darkgreen"

  barplot(ic_values,
          names.arg = tbl$Model,
          col = colors,
          main = paste("Model Comparison:", criterion),
          ylab = criterion,
          las = 2)

  abline(h = min(ic_values), lty = 2, col = "red")

  # Add delta labels
  delta <- ic_values - min(ic_values)
  mtext(paste0("Δ=", round(delta, 1)),
        side = 1, line = 3, at = seq(0.7, by = 1.2, length.out = n_models),
        cex = 0.7)

  invisible(tbl)
}

#' Plot Latent Structure Interpretation
#'
#' Creates a visual representation of the latent structure hierarchy
#' and where the data falls.
#'
#' @param result A result from successive_comparison()
#' @param ... Additional arguments
#'
#' @return Invisibly returns the result
#' @export
plot_structure_hierarchy <- function(result, ...) {
  if (!is.list(result) || is.null(result$best_model)) {
    stop("result must be output from successive_comparison()")
  }

  # Define hierarchy
  models <- c("UN", "MON", "IIO", "DM", "LCR", "RM")
  labels <- c("Unconstrained\n(Classificatory)",
              "Class\nMonotonicity",
              "Item\nOrdering",
              "Double\nMonotonicity",
              "Latent Class\nRasch",
              "Rasch\n(Quantitative)")
  structure_type <- c("Nominal", "Ordinal", "Ordinal",
                      "Ordinal", "Interval", "Interval")

  # Positions (x-coords representing constraint level)
  x_pos <- c(1, 2, 2, 3, 4, 5)
  y_pos <- c(1, 1.5, 0.5, 1, 1, 1)

  # Set up plot
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(2, 2, 3, 2))

  plot(x_pos, y_pos, type = "n",
       xlim = c(0.5, 5.5), ylim = c(0, 2),
       xlab = "", ylab = "",
       main = "Latent Structure Hierarchy",
       axes = FALSE, ...)

  # Draw arrows showing nesting
  arrows(1.2, 1, 1.8, 1.4, length = 0.1, lwd = 1.5, col = "gray50")  # UN -> MON
  arrows(1.2, 1, 1.8, 0.6, length = 0.1, lwd = 1.5, col = "gray50")  # UN -> IIO
  arrows(2.2, 1.4, 2.8, 1.1, length = 0.1, lwd = 1.5, col = "gray50")  # MON -> DM
  arrows(2.2, 0.6, 2.8, 0.9, length = 0.1, lwd = 1.5, col = "gray50")  # IIO -> DM
  arrows(3.2, 1, 3.8, 1, length = 0.1, lwd = 1.5, col = "gray50")  # DM -> LCR
  arrows(4.2, 1, 4.8, 1, length = 0.1, lwd = 1.5, col = "gray50")  # LCR -> RM

  # Draw boxes for each model
  for (i in 1:6) {
    is_best <- models[i] == result$best_model
    box_col <- if (is_best) "lightgreen" else "lightgray"
    border_col <- if (is_best) "darkgreen" else "gray50"
    lwd <- if (is_best) 3 else 1

    rect(x_pos[i] - 0.35, y_pos[i] - 0.25,
         x_pos[i] + 0.35, y_pos[i] + 0.25,
         col = box_col, border = border_col, lwd = lwd)

    text(x_pos[i], y_pos[i], labels[i], cex = 0.7)
  }

  # Add structure type labels at bottom
  text(1, -0.1, "Nominal", cex = 0.8, font = 2)
  text(2.5, -0.1, "Ordinal", cex = 0.8, font = 2)
  text(4.5, -0.1, "Interval", cex = 0.8, font = 2)

  # Add arrow showing increasing constraint
  arrows(0.7, 1.8, 5.3, 1.8, length = 0.15, lwd = 2, col = "darkblue")
  text(3, 1.95, "Increasing structural constraints", cex = 0.8, col = "darkblue")

  # Add best model label
  text(3, -0.4, paste("Best model:", result$best_model), cex = 1, font = 2)

  invisible(result)
}
