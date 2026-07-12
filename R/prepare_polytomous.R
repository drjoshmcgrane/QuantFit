#' Recode polytomous items into adjacent-category dichotomous sub-items
#'
#' Splits each polytomous item into its adjacent-category ("partial credit")
#' dichotomisations. An item with categories \eqn{0, 1, \ldots, m} becomes
#' \eqn{m} binary sub-items; sub-item \eqn{k} (for \eqn{k = 1, \ldots, m})
#' contrasts category \eqn{k} against category \eqn{k-1} and is defined only
#' for respondents who scored \eqn{k-1} or \eqn{k}:
#' \deqn{Y_{jk} = 1 \textrm{ if } X_j = k, \quad 0 \textrm{ if } X_j = k-1,
#'       \quad \textrm{NA otherwise.}}
#'
#' This is the polytomous Rasch (partial credit) decomposition: the adjacent
#' category logits of a partial credit model are exactly the Rasch item steps.
#' It is the recoding used to bring polytomous data into the conjoint checks and
#' the latent-structure models. Cumulative ("item-step" \eqn{X \ge k})
#' dichotomisation is a *different* construction and is deliberately not used
#' here.
#'
#' @param resp Integer matrix or data frame of polytomous responses, persons in
#'   rows and items in columns. Each item must be scored with consecutive
#'   integer categories starting at 0. Complete data only (no `NA`).
#'
#' @return A numeric matrix with one column per (item, step) sub-item, holding
#'   `1`, `0`, or `NA`. The integer attributes `"item"` and `"step"` give the
#'   originating item index and category step for each column.
#'
#' @seealso [prepare_polytomous()], [PrepareChecks()]
#' @keywords internal
#' @export
recode_adjacent <- function(resp) {
  resp <- as.matrix(resp)
  if (any(is.na(resp))) {
    stop("recode_adjacent() requires complete data (no NA in 'resp').")
  }
  if (any(resp != round(resp)) || any(resp < 0)) {
    stop("Polytomous responses must be non-negative integer category scores ",
         "starting at 0.")
  }
  J <- ncol(resp)
  cols <- vector("list", 0L)
  item <- step <- integer(0)
  labs <- character(0)
  cn <- colnames(resp)
  if (is.null(cn)) cn <- paste0("I", seq_len(J))
  for (j in seq_len(J)) {
    m <- max(resp[, j])
    if (m < 1L) {
      stop("Item ", cn[j], " has a single category; every item must use at ",
           "least two categories.")
    }
    x <- resp[, j]
    for (k in seq_len(m)) {
      cols[[length(cols) + 1L]] <-
        ifelse(x == k, 1L, ifelse(x == k - 1L, 0L, NA_integer_))
      item <- c(item, j)
      step <- c(step, k)
      labs <- c(labs, sprintf("%s.s%d", cn[j], k))
    }
  }
  out <- do.call(cbind, cols)
  colnames(out) <- labs
  attr(out, "item") <- item
  attr(out, "step") <- step
  out
}

#' Prepare polytomous data for the conjoint checks
#'
#' Builds the `N` and `n` count matrices required by [ConjointChecks()],
#' [HiConjointChecks()], and [KaraChecks()] from *polytomous* response data,
#' via the adjacent-category recoding of [recode_adjacent()]. Each polytomous
#' item is split into its adjacent-category dichotomous sub-items, respondents
#' are grouped by their total (sum) score across the original polytomous items
#' -- the sufficient statistic for the partial credit model -- and the
#' score-by-sub-item count matrix is formed exactly as [PrepareChecks()] does
#' for dichotomous data.
#'
#' Adjacent-category coding leaves some score-group-by-sub-item cells
#' structurally empty (a low-scoring group contains nobody at a high category
#' boundary, and vice versa), whereas [ConjointChecks()] requires every cell to
#' be populated. The score groups (rows) and sub-items (columns) are therefore
#' pruned to the largest dense block in which every cell has at least
#' `cell.lower` in-play respondents, iteratively dropping the row or column
#' carrying the most under-populated cells. This extends the sum-score-group
#' filtering that [PrepareChecks()] already applies through `ss.lower` to both
#' dimensions.
#'
#' @param resp Integer matrix or data frame of polytomous responses (persons in
#'   rows, items in columns), each item scored with consecutive integer
#'   categories starting at 0. Complete data only.
#' @param ss.lower Only total-score groups with at least this many respondents
#'   are used (as in [PrepareChecks()]).
#' @param cell.lower Minimum number of in-play respondents required in every
#'   retained cell. Cells with fewer are pruned away; larger values give a
#'   smaller but more stably estimated matrix.
#' @param order.columns If `TRUE` (default) the retained sub-items are ordered
#'   by ascending overall proportion (`sum(n)/sum(N)`), mirroring the
#'   difficulty ordering [PrepareChecks()] applies to items.
#'
#' @return A list with elements `N` and `n` -- the number of in-play responses
#'   and the number scoring the higher category, for each score-group (row) by
#'   sub-item (column) cell -- ready to pass to [ConjointChecks()] or
#'   [KaraChecks()]. The attributes `"dropped"` (a list of the removed score
#'   groups and sub-items) and `"sub_items"` (the surviving item/step map) record
#'   what the dense-core pruning kept.
#'
#' @seealso [recode_adjacent()], [PrepareChecks()], [ConjointChecks()],
#'   [KaraChecks()]
#'
#' @examples
#' # simulate a small partial-credit data set (4 categories, 6 items)
#' set.seed(1)
#' nP <- 800; J <- 6; M <- 3
#' theta <- rnorm(nP)
#' delta <- matrix(seq(-1, 1, length.out = J), J, M) +
#'   matrix(seq(-1, 1, length.out = M), J, M, byrow = TRUE)
#' resp <- matrix(0L, nP, J)
#' for (j in seq_len(J)) {
#'   num <- cbind(0, t(apply(outer(theta, delta[j, ], "-"), 1, cumsum)))
#'   P <- exp(num) / rowSums(exp(num))
#'   resp[, j] <- rowSums(runif(nP) > t(apply(P, 1, cumsum)))
#' }
#' pc <- prepare_polytomous(resp, ss.lower = 10, cell.lower = 3)
#' dim(pc$N)
#'
#' @export
prepare_polytomous <- function(resp, ss.lower = 10, cell.lower = 5,
                               order.columns = TRUE) {
  resp <- as.matrix(resp)
  if (any(is.na(resp))) {
    stop("Checks will only work with complete data. Suggestion: remove ",
         "respondents with missing responses.")
  }
  if (ss.lower < 2) {
    message("ss.lower must be greater than 1, setting to 2.")
    ss.lower <- 2
  }
  if (cell.lower < 1) cell.lower <- 1

  score <- rowSums(resp)
  rec <- recode_adjacent(resp)                     # persons x sub-items (w/ NA)

  tab <- table(score)
  lev <- as.numeric(names(tab))[tab >= ss.lower]
  if (length(lev) < 3) {
    stop("Fewer than three total-score groups have at least ss.lower = ",
         ss.lower, " respondents; cannot form a conjoint matrix.")
  }
  Nl <- nl <- vector("list", length(lev))
  names(Nl) <- names(nl) <- as.character(lev)
  for (i in seq_along(lev)) {
    sub <- rec[score == lev[i], , drop = FALSE]
    Nl[[i]] <- colSums(!is.na(sub))
    nl[[i]] <- colSums(sub == 1, na.rm = TRUE)
  }
  N <- do.call(rbind, Nl)
  n <- do.call(rbind, nl)
  item <- attr(rec, "item")
  step <- attr(rec, "step")

  # prune to the largest dense block: every retained cell has N >= cell.lower
  dropped_rows <- character(0)
  dropped_cols <- character(0)
  repeat {
    bad <- N < cell.lower
    if (!any(bad)) break
    rowbad <- rowSums(bad)
    colbad <- colSums(bad)
    if (max(rowbad) >= max(colbad)) {
      d <- which.max(rowbad)
      dropped_rows <- c(dropped_rows, rownames(N)[d])
      N <- N[-d, , drop = FALSE]; n <- n[-d, , drop = FALSE]
    } else {
      d <- which.max(colbad)
      dropped_cols <- c(dropped_cols, colnames(N)[d])
      N <- N[, -d, drop = FALSE]; n <- n[, -d, drop = FALSE]
      item <- item[-d]; step <- step[-d]
    }
    if (nrow(N) < 3 || ncol(N) < 3) {
      stop("Dense-core pruning left fewer than a 3x3 matrix (got ", nrow(N),
           "x", ncol(N), "). Try lowering cell.lower or ss.lower, or supply ",
           "more items/respondents.")
    }
  }

  if (order.columns) {
    ord <- order(colSums(n) / colSums(N))
    N <- N[, ord, drop = FALSE]; n <- n[, ord, drop = FALSE]
    item <- item[ord]; step <- step[ord]
  }

  out <- list(N = N, n = n)
  attr(out, "dropped") <- list(score_groups = dropped_rows,
                               sub_items = dropped_cols)
  attr(out, "sub_items") <- data.frame(sub_item = colnames(N),
                                       item = item, step = step,
                                       stringsAsFactors = FALSE)
  out
}
