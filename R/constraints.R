#' Constraint Specification Functions for QuantFit
#'
#' @name constraints
#' @description Functions to build and manage inequality constraints for ordered latent class models
NULL

#' Specify constraints for model fitting
#'
#' @param class_monotonicity Logical, enforce class monotonicity (β_ic ≤ β_ic' for c < c')
#' @param item_ordering Logical, enforce invariant item ordering (β_ic ≤ β_i'c for i < i')
#' @param item_order Optional vector specifying item order (indices from easiest to hardest).
#'   Required if item_ordering = TRUE.
#'
#' @return A constraints specification object
#'
#' @examples
#' # Class monotonicity only
#' spec <- specify_constraints(class_monotonicity = TRUE)
#'
#' # Both constraints (Double Monotonicity)
#' spec <- specify_constraints(
#'   class_monotonicity = TRUE,
#'   item_ordering = TRUE,
#'   item_order = c(1, 3, 2, 5, 4)
#' )
#'
#' @export
specify_constraints <- function(class_monotonicity = FALSE,
                                item_ordering = FALSE,
                                item_order = NULL) {

  if (item_ordering && is.null(item_order)) {
    message("Note: item_order not specified. Will be estimated from data.")
  }

  structure(
    list(
      class_monotonicity = class_monotonicity,
      item_ordering = item_ordering,
      item_order = item_order
    ),
    class = "ql_constraints"
  )
}

#' Build class monotonicity constraint function
#'
#' Constraint: β_ic ≤ β_ic' for all items i and classes c < c'
#' This means item probabilities increase across classes.
#'
#' @param n_items Number of items
#' @param n_classes Number of latent classes
#'
#' @return Function that takes item_probs matrix and returns constraint values (>= 0 if satisfied)
#' @keywords internal
build_class_monotonicity_fn <- function(n_items, n_classes) {
  function(item_probs) {
    # For each item and each pair of adjacent classes,
    # constraint is: p_{i,c+1} - p_{i,c} >= 0

    constraints <- numeric(0)

    for (i in 1:n_items) {
      for (c in 1:(n_classes - 1)) {
        # p_{i,c+1} >= p_{i,c}
        constraints <- c(constraints, item_probs[i, c + 1] - item_probs[i, c])
      }
    }

    constraints
  }
}

#' Build item ordering constraint function
#'
#' Constraint: β_ic ≤ β_i'c for items ordered i < i' (easier items have higher probs)
#' This means within each class, items maintain the same difficulty ordering.
#'
#' @param n_items Number of items
#' @param n_classes Number of latent classes
#' @param item_order Vector of item indices from easiest to hardest
#'
#' @return Function that takes item_probs matrix and returns constraint values (>= 0 if satisfied)
#' @keywords internal
build_item_ordering_fn <- function(n_items, n_classes, item_order) {
  function(item_probs) {
    # For each class and each pair of adjacent items in order,
    # constraint is: p_{easier,c} - p_{harder,c} >= 0

    constraints <- numeric(0)

    for (c in 1:n_classes) {
      for (k in 1:(length(item_order) - 1)) {
        easier_item <- item_order[k]
        harder_item <- item_order[k + 1]
        # p_{easier,c} >= p_{harder,c}
        constraints <- c(constraints,
                        item_probs[easier_item, c] - item_probs[harder_item, c])
      }
    }

    constraints
  }
}

#' Build double monotonicity constraint function
#'
#' Combines class monotonicity and item ordering constraints.
#'
#' @param n_items Number of items
#' @param n_classes Number of latent classes
#' @param item_order Vector of item indices from easiest to hardest
#'
#' @return Function that takes item_probs matrix and returns constraint values (>= 0 if satisfied)
#' @keywords internal
build_double_monotonicity_fn <- function(n_items, n_classes, item_order) {
  class_mon_fn <- build_class_monotonicity_fn(n_items, n_classes)
  item_ord_fn <- build_item_ordering_fn(n_items, n_classes, item_order)

  function(item_probs) {
    c(class_mon_fn(item_probs), item_ord_fn(item_probs))
  }
}

#' Build constraint function from specification
#'
#' @param constraints A ql_constraints object
#' @param n_items Number of items
#' @param n_classes Number of latent classes
#' @param item_order Item order (may be estimated if not provided)
#'
#' @return Constraint function
#' @keywords internal
build_constraint_fn <- function(constraints, n_items, n_classes, item_order = NULL) {
  if (constraints$class_monotonicity && constraints$item_ordering) {
    # Double monotonicity
    if (is.null(item_order)) {
      stop("item_order must be specified for item ordering constraints")
    }
    return(build_double_monotonicity_fn(n_items, n_classes, item_order))

  } else if (constraints$class_monotonicity) {
    # Class monotonicity only
    return(build_class_monotonicity_fn(n_items, n_classes))

  } else if (constraints$item_ordering) {
    # Item ordering only
    if (is.null(item_order)) {
      stop("item_order must be specified for item ordering constraints")
    }
    return(build_item_ordering_fn(n_items, n_classes, item_order))

  } else {
    # No constraints - return function that always satisfies
    return(function(item_probs) numeric(0))
  }
}

#' Count number of constraints
#'
#' @param constraints A ql_constraints object
#' @param n_items Number of items
#' @param n_classes Number of latent classes
#'
#' @return Number of inequality constraints
#' @keywords internal
count_constraints <- function(constraints, n_items, n_classes) {
  n_constraints <- 0

  if (constraints$class_monotonicity) {
    # I items x (C-1) adjacent class pairs
    n_constraints <- n_constraints + n_items * (n_classes - 1)
  }

  if (constraints$item_ordering) {
    # C classes x (I-1) adjacent item pairs
    n_constraints <- n_constraints + n_classes * (n_items - 1)
  }

  n_constraints
}

#' Check if constraints are satisfied
#'
#' @param item_probs Item probability matrix (I x C)
#' @param constraints A ql_constraints object
#' @param item_order Item order vector (if needed)
#' @param tol Tolerance for constraint satisfaction
#'
#' @return List with 'satisfied' logical and 'violations' details
#' @keywords internal
check_constraints <- function(item_probs, constraints, item_order = NULL, tol = 1e-6) {
  n_items <- nrow(item_probs)
  n_classes <- ncol(item_probs)

  violations <- list()

  # Check class monotonicity
  if (constraints$class_monotonicity) {
    for (i in 1:n_items) {
      for (c in 1:(n_classes - 1)) {
        if (item_probs[i, c + 1] < item_probs[i, c] - tol) {
          violations <- c(violations, list(list(
            type = "class_monotonicity",
            item = i,
            classes = c(c, c + 1),
            probs = c(item_probs[i, c], item_probs[i, c + 1]),
            violation = item_probs[i, c] - item_probs[i, c + 1]
          )))
        }
      }
    }
  }

  # Check item ordering
  if (constraints$item_ordering && !is.null(item_order)) {
    for (c in 1:n_classes) {
      for (k in 1:(length(item_order) - 1)) {
        easier <- item_order[k]
        harder <- item_order[k + 1]
        if (item_probs[easier, c] < item_probs[harder, c] - tol) {
          violations <- c(violations, list(list(
            type = "item_ordering",
            class = c,
            items = c(easier, harder),
            probs = c(item_probs[easier, c], item_probs[harder, c]),
            violation = item_probs[harder, c] - item_probs[easier, c]
          )))
        }
      }
    }
  }

  list(
    satisfied = length(violations) == 0,
    n_violations = length(violations),
    violations = violations
  )
}

#' Project item probabilities onto constraint space
#'
#' Uses isotonic regression to enforce monotonicity constraints.
#'
#' @param item_probs Item probability matrix (I x C)
#' @param constraints A ql_constraints object
#' @param item_order Item order vector (if needed)
#'
#' @return Projected item probability matrix
#' @keywords internal
project_constraints <- function(item_probs, constraints, item_order = NULL) {
  n_items <- nrow(item_probs)
  n_classes <- ncol(item_probs)

  proj_probs <- item_probs

  # Project class monotonicity using PAVA (pool adjacent violators)
  if (constraints$class_monotonicity) {
    for (i in 1:n_items) {
      proj_probs[i, ] <- pava_increasing(proj_probs[i, ])
    }
  }

  # Project item ordering
  if (constraints$item_ordering && !is.null(item_order)) {
    for (c in 1:n_classes) {
      # Get probs in item order
      ordered_probs <- proj_probs[item_order, c]
      # Apply PAVA (decreasing because easier items have higher probs)
      ordered_probs <- pava_decreasing(ordered_probs)
      # Put back
      proj_probs[item_order, c] <- ordered_probs
    }
  }

  # Bound probabilities
  proj_probs <- bound_probs(proj_probs)

  proj_probs
}

#' Pool Adjacent Violators Algorithm (increasing)
#'
#' @param x Numeric vector
#'
#' @return Isotonically increasing vector
#' @keywords internal
pava_increasing <- function(x) {
  n <- length(x)
  if (n <= 1) return(x)

  # Initialize blocks
  result <- x
  weights <- rep(1, n)

  repeat {
    # Find first violation
    violation_found <- FALSE
    for (i in 1:(n - 1)) {
      if (result[i] > result[i + 1]) {
        # Pool adjacent blocks
        pooled_value <- (result[i] * weights[i] + result[i + 1] * weights[i + 1]) /
                        (weights[i] + weights[i + 1])
        result[i] <- pooled_value
        result[i + 1] <- pooled_value
        weights[i] <- weights[i] + weights[i + 1]
        weights[i + 1] <- weights[i]
        violation_found <- TRUE
        break
      }
    }

    if (!violation_found) break
  }

  result
}

#' Pool Adjacent Violators Algorithm (decreasing)
#'
#' @param x Numeric vector
#'
#' @return Isotonically decreasing vector
#' @keywords internal
pava_decreasing <- function(x) {
  # Reverse, apply increasing PAVA, reverse back
  rev(pava_increasing(rev(x)))
}

#' Print method for constraint specifications
#'
#' @param x A ql_constraints object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns x
#' @export
print.ql_constraints <- function(x, ...) {
  cat("QuantFit Constraint Specification\n")
  cat("------------------------------------\n")
  cat("Class monotonicity:", x$class_monotonicity, "\n")
  cat("Item ordering:", x$item_ordering, "\n")

  if (!is.null(x$item_order)) {
    cat("Item order:", paste(x$item_order, collapse = " -> "), "\n")
  }

  invisible(x)
}
