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
#' Uses (weighted) isotonic regression to enforce monotonicity constraints.
#' This is a convenience wrapper around
#' \code{\link{project_constraints_weighted}} with unit weights.
#'
#' @param item_probs Item probability matrix (I x C)
#' @param constraints A ql_constraints object
#' @param item_order Item order vector (if needed)
#'
#' @return Projected item probability matrix
#' @keywords internal
project_constraints <- function(item_probs, constraints, item_order = NULL) {
  project_constraints_weighted(item_probs, constraints, item_order,
                               class_weights = NULL)
}

#' Weighted L2 projection of item probabilities onto constraint space
#'
#' Computes the exact projection of the item probability matrix onto the
#' constraint set under the weighted L2 norm with cell weights
#' \eqn{w_{ic} = N_c} (the expected class counts). By the exponential-family
#' isotonic MLE theorem (Robertson, Wright & Dykstra), this projection of the
#' unconstrained M-step means equals the exact constrained M-step maximizer
#' of the expected complete-data log-likelihood.
#'
#' \itemize{
#'   \item Class monotonicity only: weighted PAVA per item (row) with
#'     weights \eqn{N_c}. Exact.
#'   \item Item ordering only: PAVA per class (column); weights are constant
#'     within a column so plain PAVA is exact.
#'   \item Both (double monotonicity): 2D weighted isotonic regression via
#'     Dykstra's alternating-projections algorithm
#'     (\code{\link{dykstra_dm_projection}}). Exact.
#' }
#'
#' @param item_probs Item probability matrix (I x C), typically the
#'   unconstrained M-step weighted means
#' @param constraints A ql_constraints object
#' @param item_order Item order vector (if needed)
#' @param class_weights Vector of class weights (length C), typically the
#'   expected class counts N_c from the E-step. NULL means unit weights.
#'
#' @return Projected item probability matrix
#' @keywords internal
project_constraints_weighted <- function(item_probs, constraints,
                                         item_order = NULL,
                                         class_weights = NULL) {
  n_items <- nrow(item_probs)
  n_classes <- ncol(item_probs)

  if (is.null(class_weights)) {
    class_weights <- rep(1, n_classes)
  }
  # Guard against fully collapsed classes (zero weight breaks weighted means)
  class_weights <- pmax(class_weights, 1e-10)

  if (constraints$class_monotonicity && constraints$item_ordering) {
    # Double monotonicity: exact 2D weighted isotonic regression via Dykstra
    if (is.null(item_order)) {
      stop("item_order must be specified for item ordering constraints")
    }
    proj_probs <- dykstra_dm_projection(item_probs, item_order, class_weights)
    return(bound_probs(proj_probs))
  }

  proj_probs <- item_probs

  # Class monotonicity: weighted PAVA per item across ordered classes
  if (constraints$class_monotonicity) {
    for (i in 1:n_items) {
      proj_probs[i, ] <- pava_increasing(proj_probs[i, ], class_weights)
    }
  }

  # Item ordering: PAVA per class across ordered items. Cell weights are
  # constant within a class (all equal N_c), so plain PAVA is the exact
  # weighted projection.
  if (constraints$item_ordering && !is.null(item_order)) {
    for (c in 1:n_classes) {
      ordered_probs <- proj_probs[item_order, c]
      ordered_probs <- pava_decreasing(ordered_probs)
      proj_probs[item_order, c] <- ordered_probs
    }
  }

  bound_probs(proj_probs)
}

#' 2D isotonic regression via Dykstra's alternating projections
#'
#' Computes the exact weighted L2 projection of a matrix onto the set of
#' matrices that are (a) non-decreasing along each row (over ordered classes)
#' and (b) non-increasing along each column when rows are taken in
#' \code{item_order} (easier items keep higher probabilities). Cell weights
#' are \eqn{w_{ic} = class\_weights[c]}, constant across items within a
#' class, so both partial projections are one-dimensional PAVA problems:
#' weighted PAVA for rows, plain PAVA for columns.
#'
#' Dykstra's algorithm cycles the two projections with increment (correction)
#' matrices; unlike naive alternating projections it converges to the exact
#' projection onto the intersection.
#'
#' @param item_probs Item probability matrix (I x C)
#' @param item_order Vector of item indices from easiest to hardest
#' @param class_weights Vector of class weights (length C); NULL = unit
#' @param tol Convergence tolerance on the iterate change (default 1e-10)
#' @param max_cycles Maximum number of projection cycles (default 500)
#'
#' @return Projected matrix satisfying both constraint sets
#' @keywords internal
dykstra_dm_projection <- function(item_probs, item_order,
                                  class_weights = NULL,
                                  tol = 1e-10, max_cycles = 500) {
  n_items <- nrow(item_probs)
  n_classes <- ncol(item_probs)

  if (is.null(class_weights)) {
    class_weights <- rep(1, n_classes)
  }
  class_weights <- pmax(class_weights, 1e-10)

  x <- item_probs
  incr_row <- matrix(0, n_items, n_classes)  # correction for row constraint set
  incr_col <- matrix(0, n_items, n_classes)  # correction for column constraint set

  for (cycle in 1:max_cycles) {
    x_old <- x

    # Project onto row-monotone set (weighted PAVA per row)
    z <- x + incr_row
    y <- z
    for (i in 1:n_items) {
      y[i, ] <- pava_increasing(z[i, ], class_weights)
    }
    incr_row <- z - y

    # Project onto column-ordered set (plain PAVA per column in item_order;
    # weights constant within a column under the weighted norm)
    z <- y + incr_col
    x <- z
    for (c in 1:n_classes) {
      x[item_order, c] <- pava_decreasing(z[item_order, c])
    }
    incr_col <- z - x

    if (max(abs(x - x_old)) < tol) break
  }

  x
}

#' Pool Adjacent Violators Algorithm (increasing)
#'
#' Weighted isotonic regression: minimizes \eqn{\sum_k w_k (y_k - x_k)^2}
#' subject to \eqn{y_1 \le y_2 \le \dots \le y_n}. Uses the classical
#' block-merging PAVA: maintains a stack of blocks (weighted mean, summed
#' weight, count) and merges while the last two blocks violate monotonicity.
#' Matches \code{stats::isoreg} for unit weights.
#'
#' @param x Numeric vector
#' @param w Optional non-negative weights (same length as x). NULL = unit.
#'
#' @return Isotonically increasing vector (weighted L2 projection of x)
#' @keywords internal
pava_increasing <- function(x, w = NULL) {
  n <- length(x)
  if (n <= 1) return(x)

  if (is.null(w)) {
    w <- rep(1, n)
  } else {
    if (length(w) != n) stop("weights must have same length as x")
    w <- pmax(w, 0)
  }

  # Stack of blocks: value = weighted mean, weight = summed weight, count
  vals <- numeric(n)
  wts <- numeric(n)
  cnts <- integer(n)
  nb <- 0L

  for (i in 1:n) {
    nb <- nb + 1L
    vals[nb] <- x[i]
    wts[nb] <- w[i]
    cnts[nb] <- 1L

    # Merge while the last two blocks violate monotonicity
    while (nb > 1L && vals[nb - 1L] > vals[nb]) {
      w_sum <- wts[nb - 1L] + wts[nb]
      if (w_sum > 0) {
        vals[nb - 1L] <- (vals[nb - 1L] * wts[nb - 1L] + vals[nb] * wts[nb]) / w_sum
      } else {
        # Both blocks have zero weight: use simple average
        vals[nb - 1L] <- (vals[nb - 1L] + vals[nb]) / 2
      }
      wts[nb - 1L] <- w_sum
      cnts[nb - 1L] <- cnts[nb - 1L] + cnts[nb]
      nb <- nb - 1L
    }
  }

  # Expand blocks back to the full vector
  rep(vals[seq_len(nb)], cnts[seq_len(nb)])
}

#' Pool Adjacent Violators Algorithm (decreasing)
#'
#' Antitonic (non-increasing) weighted isotonic regression, implemented by
#' reversing, applying \code{\link{pava_increasing}}, and reversing back.
#'
#' @param x Numeric vector
#' @param w Optional non-negative weights (same length as x). NULL = unit.
#'
#' @return Isotonically decreasing vector
#' @keywords internal
pava_decreasing <- function(x, w = NULL) {
  if (is.null(w)) {
    rev(pava_increasing(rev(x)))
  } else {
    rev(pava_increasing(rev(x), rev(w)))
  }
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
