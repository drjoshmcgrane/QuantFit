#' Simulate item responses from a specified latent-structure model
#'
#' Generates dichotomous or polytomous item response data from any of the six
#' models in the Torres Irribarra & Diakow hierarchy, for validation and power
#' studies. The dichotomous case reproduces the generators of the original paper
#' (random logits on \eqn{U(-4, 4)}, sorted across classes for class
#' monotonicity, across items for invariant item ordering, or both for double
#' monotonicity; a partial-credit/Rasch structure for the quantitative models),
#' and the polytomous case generalises them through a cumulative-threshold
#' (for UN/MON/IIO/DM) or partial-credit (for LCR/RM) formulation.
#'
#' @details
#' For the non-parametric models each item-by-class combination is given a
#' location \eqn{L_{cj}} and the category structure is built from shared,
#' increasing threshold offsets, so that
#' \eqn{P(X \ge k \mid c, j) = \mathrm{logit}^{-1}(L_{cj} - \tau_k)}:
#' \describe{
#'   \item{UN}{\eqn{L} left unsorted (class and item profiles may cross).}
#'   \item{MON}{each column of \eqn{L} sorted, so category distributions are
#'     stochastically ordered across classes.}
#'   \item{IIO}{each row of \eqn{L} sorted, so expected item scores share one
#'     ordering across classes.}
#'   \item{DM}{both, so the data are doubly monotone but not necessarily
#'     partial-credit.}
#' }
#' The quantitative models use the partial credit model
#' \eqn{P(X = x \mid \theta, j) \propto \exp\{\sum_{l \le x}(\theta - \delta_{jl})\}}
#' with discrete, separated class locations (LCR) or a continuous
#' \eqn{\theta \sim N(0, 1)} (RM). When `n_cat = 2` every model reduces to the
#' dichotomous generator.
#'
#' @section The partial-order model (`"PO"`):
#' TI&D's generators are one \eqn{C \times J} matrix of \eqn{U(-4,4)} logits
#' with each column (item) fully sorted for MON and left unsorted for UN.
#' Sorting a column imposes the CHAIN (total class order); not sorting is the
#' ANTICHAIN. `"PO"` is the natural interpolation within the same idiom: each
#' column's values are sorted descending and assigned to classes along a
#' **random linear extension of a specified partial order**, drawn independently
#' per item. Whenever \eqn{a \succeq b} in the poset, \eqn{a} precedes \eqn{b}
#' in *every* linear extension, so \eqn{L_{aj} \ge L_{bj}} for all items -
#' \eqn{a} dominates \eqn{b} exactly. Incomparable pairs receive random relative
#' order per item and therefore cross (a rejection loop - the same device TI&D
#' use for LCR class separation - guarantees each incomparable pair crosses in
#' both directions by at least `po_margin` on the probability scale). The chain
#' poset reproduces `"MON"` exactly and the antichain reproduces `"UN"`, so
#' `"PO"` nests both endpoints. The item side is left unconstrained (as in UN
#' and MON), so the class structure is the only signal: the resulting data
#' violate class monotonicity (a crossing pair exists) yet carry real dominance
#' structure that the six-model set cannot express - the ground truth for
#' testing the partial-order refinement of [select_model_hybrid()].
#'
#' @param model One of `"UN"`, `"MON"`, `"IIO"`, `"DM"`, `"LCR"`, `"RM"`,
#'   `"PO"` (partial class order; see the dedicated section), or `"PO_ITEMS"`
#'   (the mirror: classes totally ordered - MON holds exactly - while ITEMS are
#'   partially ordered; `poset` is then over items, keywords `"layers2"` /
#'   `"single"` or a J x J dominance matrix; routes to the MON cell, where the
#'   item-poset refinement must catch what the six-model verdict cannot say).
#' @param n_persons Number of respondents.
#' @param n_items Number of items.
#' @param n_classes Number of latent classes (ignored for `"RM"`).
#' @param n_cat Number of ordered response categories (2 = dichotomous).
#' @param class_probs Optional class mixing proportions (length `n_classes`);
#'   defaults to equal.
#' @param poset For `model = "PO"` only: the class partial order. Either a
#'   keyword - `"V"` (one maximal class dominates all others, which are mutually
#'   incomparable; the default), `"Lambda"` (all classes dominate one minimal
#'   class), `"single"` (class 1 dominates class 2; everything else
#'   incomparable) - or a `n_classes x n_classes` logical matrix `D` with
#'   `D[a, b] = TRUE` meaning class `a` dominates class `b` (must be
#'   irreflexive and acyclic; transitive closure is taken). Must be neither a
#'   chain (use `"MON"`) nor an antichain (use `"UN"`).
#' @param po_margin For `"PO"`/`"PO_ITEMS"`: minimum crossing depth, on the
#'   probability scale, required of every incomparable pair in both directions
#'   (rejection-sampled; default 0.05).
#' @param item_order For `model = "PO"` only: `"free"` (default) leaves the item
#'   side unconstrained, routing to the UN cell; `"invariant"` builds class
#'   lines over one shared sorted item spine (`L = a_c s_j + b_c`, `a_c > 0`),
#'   so invariant item ordering HOLDS while class dominance follows the poset -
#'   routing to the IIO cell, where only the class-poset refinement can express
#'   the truth.
#' @param seed Optional random seed.
#'
#' @return An integer matrix of responses (`n_persons` rows, `n_items` columns),
#'   scored `0..n_cat-1`, with attributes `"model"` and `"params"` (the
#'   generating parameters and, for the mixture models, the class memberships).
#'
#' @examples
#' # doubly monotone dichotomous data, 3 classes
#' d <- simulate_responses("DM", n_persons = 300, n_items = 8, n_classes = 3,
#'                         seed = 1)
#' # partial-credit (Rasch) data with 4 categories
#' p <- simulate_responses("RM", n_persons = 300, n_items = 6, n_cat = 4,
#'                         seed = 1)
#'
#' @seealso [select_model_ll()], [quant_fit()]
#' @export
simulate_responses <- function(model = c("UN", "MON", "IIO", "DM", "LCR", "RM",
                                         "PO", "PO_ITEMS"),
                               n_persons = 500, n_items = 10, n_classes = 3,
                               n_cat = 2L, class_probs = NULL,
                               poset = "V", po_margin = 0.05,
                               item_order = c("free", "invariant"),
                               seed = NULL) {
  model <- match.arg(model); item_order <- match.arg(item_order)
  if (!is.null(seed)) set.seed(seed)
  n_cat <- as.integer(n_cat)
  if (n_cat < 2L) stop("n_cat must be >= 2")
  m <- n_cat - 1L
  tau <- if (m == 1L) 0 else seq(-1.2, 1.2, length.out = m)  # increasing offsets

  # cumulative P(X>=k) (length m) -> category probs (length m+1)
  cum_to_cat <- function(cum) c(1, cum) - c(cum, 0)

  params <- list()

  if (model %in% c("UN", "MON", "IIO", "DM", "PO", "PO_ITEMS")) {
    if (is.null(class_probs)) class_probs <- rep(1 / n_classes, n_classes)
    L <- matrix(stats::runif(n_classes * n_items, -4, 4), n_classes, n_items)
    if (model %in% c("MON", "DM")) L <- apply(L, 2L, sort)          # classes ordered
    if (model %in% c("IIO", "DM")) L <- t(apply(L, 1L, sort))       # items ordered
    if (model == "PO") {
      D <- .po_resolve_poset(poset, n_classes)
      # Rejection loop (TI&D's own device, cf. their LCR separation loop):
      # regenerate until every incomparable pair crosses in BOTH directions by
      # at least po_margin on the probability scale, so no incomparable pair can
      # be read as dominated by accident.
      if (item_order == "free") {
        repeat {
          L <- .po_assign(matrix(stats::runif(n_classes * n_items, -4, 4),
                                 n_classes, n_items), D)
          if (.po_crossing_ok(stats::plogis(L), D, po_margin)) break
        }
      } else {
        # item_order = "invariant": class lines over one shared sorted item
        # spine, L[c, j] = a_c * s_j + b_c with a_c > 0, so each row is sorted
        # (invariant item ordering HOLDS) while class dominance follows the
        # poset: comparable pairs must not cross within the spine's range,
        # incomparable pairs must cross. This is the PO analogue of the IIO
        # row-sort, and it routes to the IIO cell - where the class-poset
        # refinement must catch what the six-model verdict cannot express.
        # class crossings must clear the MON axis's own detection floor
        # (resample-q05 vs eps ~ 0.01 treats mean-violations < ~0.1 as noise),
        # otherwise the data route to DM instead of the IIO cell.
        # CONSTRUCTIVE search, not blind rejection (blind draws satisfy the
        # joint exact-dominance + balanced-crossing geometry with near-zero
        # probability - minutes per dataset): slopes get guaranteed pairwise
        # separation (shuffled equal-spaced grid), intercepts are placed in
        # topological order so every comparable pair is dominated EXACTLY by
        # construction, and each incomparable pair's crossing point is placed
        # deliberately inside the spine. Capped best-of loop; achieved margin
        # recorded, warning if it falls short (mirrors PO_ITEMS).
        m_eff <- max(po_margin, 0.12)
        n_dom <- rowSums(D)                     # how many classes each dominates
        best_L <- NULL; best_m <- -Inf
        inc_pairs <- which(upper.tri(D) & !(D | t(D)), arr.ind = TRUE)
        for (try_i in seq_len(400L)) {
          s_j <- sort(stats::runif(n_items, -2.5, 2.5))
          a <- sample(seq(0.4, 1.6, length.out = n_classes))
          if (nrow(inc_pairs)) {                # extreme slope separation for
            a[inc_pairs[1L, ]] <- sample(c(0.4, 1.6))   # the crossing pair
            a[-inc_pairs[1L, ]] <- sample(seq(0.7, 1.3,
              length.out = n_classes - 2L))
          }
          b <- numeric(n_classes)
          for (cc in order(n_dom)) {            # bottoms first
            below <- which(D[cc, ])            # classes cc must dominate
            placed <- below[match(below, order(n_dom)) < match(cc, order(n_dom))]
            inc <- setdiff(which(!(D[cc, ] | D[, cc]) & seq_len(n_classes) != cc),
                           cc)
            inc <- inc[match(inc, order(n_dom)) < match(cc, order(n_dom))]
            b[cc] <- stats::runif(1, -2, 2)
            if (length(inc)) {                  # engineer a crossing at the
              # item-mass median (balanced masses on both sides)
              u <- stats::median(s_j) + stats::runif(1, -0.4, 0.4)
              x <- inc[1L]
              b[cc] <- a[x] * u + b[x] - a[cc] * u
            }
            if (length(placed))                 # exact dominance (hard floor)
              b[cc] <- max(b[cc], max(vapply(placed, function(y)
                max(a[y] * s_j + b[y] - a[cc] * s_j), numeric(1))) + 0.05)
          }
          L <- outer(a, s_j) + b
          # recentre on the INCOMPARABLE classes (constant shift keeps
          # dominance/crossings; puts the crossing pair at mid-probability
          # where the margin lives - dominant classes may saturate, they only
          # need dominance)
          ctr <- if (nrow(inc_pairs)) unique(as.vector(inc_pairs))
                 else seq_len(n_classes)
          L <- L - mean(L[ctr, ])
          P <- stats::plogis(L)
          pr <- inc_pairs
          m_got <- if (nrow(pr)) min(vapply(seq_len(nrow(pr)), function(k) {
            x <- pr[k, 1]; y <- pr[k, 2]
            min(mean(pmax(0, P[x, ] - P[y, ])), mean(pmax(0, P[y, ] - P[x, ])))
          }, numeric(1))) else Inf
          if (m_got > best_m) { best_m <- m_got; best_L <- L }
          if (best_m >= m_eff) break
        }
        L <- best_L
        if (best_m < m_eff)
          warning(sprintf(
            "PO invariant: achieved crossing margin %.4f < target %.2f (best of 400 draws)",
            best_m, m_eff))
      }
    }
    if (model == "PO_ITEMS") {
      # The MIRROR: classes totally ordered (MON holds), items PARTIALLY
      # ordered. theta_c is a well-separated chain; per class, item effects are
      # assigned along a random linear extension of the ITEM poset with spread
      # strictly smaller than the class spacing, so every column stays sorted
      # (exact class monotonicity) while item pairs dominate iff comparable in
      # the poset. Routes to the MON cell, where the item-poset refinement must
      # catch the structure. `poset` here is over ITEMS (keywords: "layers2" =
      # top half dominates bottom half, halves internally incomparable;
      # "single"; or a J x J dominance matrix).
      Di <- .po_resolve_item_poset(poset, n_items)
      delta <- 2.2                                   # class spacing
      theta <- (seq_len(n_classes) - (n_classes + 1) / 2) * delta
      # capped best-of search (joint crossing over many pairs is a rare event
      # for a single draw); the achieved margin is recorded in params
      best <- NULL; best_m <- -Inf
      for (it in seq_len(2000L)) {
        E <- t(.po_assign(t(matrix(stats::runif(n_classes * n_items,
               -delta / 2 + 0.08, delta / 2 - 0.08), n_classes, n_items)), Di))
        Lc <- matrix(theta, n_classes, n_items) + E
        mc <- .po_item_min_crossing(stats::plogis(Lc), Di)
        if (mc > best_m) { best_m <- mc; best <- Lc }
        if (mc >= max(po_margin, 0.10)) break
      }
      L <- best
      if (best_m < max(po_margin, 0.10))
        warning("PO_ITEMS: achieved crossing margin ", signif(best_m, 3),
                " < po_margin ", po_margin, " (best of 2000 draws)")
    }
    cls <- sample.int(n_classes, n_persons, replace = TRUE, prob = class_probs)
    resp <- matrix(0L, n_persons, n_items)
    for (j in seq_len(n_items)) {
      cumj <- plogis(outer(L[, j], tau, "-"))                       # C x m: P(X>=k|c)
      P <- cbind(1, cumj) - cbind(cumj, 0)                          # C x (m+1)
      for (c in seq_len(n_classes)) {
        who <- which(cls == c)
        if (length(who))
          resp[who, j] <- sample.int(m + 1L, length(who), replace = TRUE,
                                     prob = P[c, ]) - 1L
      }
    }
    params <- list(L = L, tau = tau, class = cls, class_probs = class_probs)
    if (model == "PO") params$poset <- D

  } else {  # LCR, RM: partial credit model
    b <- stats::runif(n_items, -2, 2)                              # item locations
    delta_list <- lapply(b, function(bj) bj + tau)                 # item step params
    if (model == "LCR") {
      if (is.null(class_probs)) class_probs <- rep(1 / n_classes, n_classes)
      repeat {                                                     # separated classes
        a <- sort(stats::runif(n_classes, -3, 3))
        if (all(diff(a) > 0.5)) break
      }
      cls <- sample.int(n_classes, n_persons, replace = TRUE, prob = class_probs)
      theta <- a[cls]
      params <- list(theta_class = a, item = b, tau = tau, class = cls)
    } else {
      theta <- stats::rnorm(n_persons)
      params <- list(item = b, tau = tau)
    }
    ip <- compute_pcm_probs(theta, delta_list, use_cpp = TRUE)      # list of n x (m+1)
    resp <- matrix(0L, n_persons, n_items)
    for (j in seq_len(n_items)) {
      cdf <- t(apply(ip[[j]], 1L, cumsum))
      resp[, j] <- rowSums(stats::runif(n_persons) > cdf)
    }
  }

  colnames(resp) <- paste0("Item", seq_len(n_items))
  attr(resp, "model") <- model
  attr(resp, "params") <- params
  resp
}

# --- helpers for the "PO" (partial class order) generator -------------------
# Resolve a poset spec to a strict-dominance logical matrix D (D[a,b]: a >= b),
# transitively closed, validated as irreflexive/acyclic and neither chain nor
# antichain.
.po_resolve_poset <- function(poset, C) {
  if (is.character(poset) && length(poset) == 1L) {
    D <- matrix(FALSE, C, C)
    if (identical(poset, "V")) {
      if (C < 3L) stop("poset \"V\" needs n_classes >= 3")
      D[1L, 2:C] <- TRUE                       # class 1 dominates all others
    } else if (identical(poset, "Lambda")) {
      if (C < 3L) stop("poset \"Lambda\" needs n_classes >= 3")
      D[2:C, 1L] <- TRUE                       # all others dominate class 1
    } else if (identical(poset, "single")) {
      if (C < 3L) stop("poset \"single\" needs n_classes >= 3")
      D[1L, 2L] <- TRUE                        # only 1 >= 2; rest incomparable
    } else stop("unknown poset keyword: ", poset)
  } else {
    D <- matrix(as.logical(poset), C, C)
    if (any(is.na(D))) stop("poset matrix must be logical")
    if (any(diag(D))) stop("poset must be irreflexive")
    # transitive closure (Warshall)
    for (k in seq_len(C)) for (a in seq_len(C)) if (D[a, k])
      D[a, ] <- D[a, ] | D[k, ]
    if (any(D & t(D))) stop("poset contains a cycle")
  }
  n_comp <- sum(D[upper.tri(D)] | t(D)[upper.tri(D)])
  if (n_comp == 0L) stop("poset is an antichain - use model = \"UN\"")
  if (n_comp == C * (C - 1L) / 2L) stop("poset is a chain - use model = \"MON\"")
  D
}

# Assign each column's values (sorted descending) to classes along a random
# linear extension of D, drawn independently per item: comparable pairs then
# dominate identically on every item; incomparable pairs get random relative
# order per item. Chain -> exactly the MON column sort; antichain -> UN.
.po_assign <- function(L, D) {
  C <- nrow(L)
  for (j in seq_len(ncol(L))) {
    vals <- sort(L[, j], decreasing = TRUE)
    remaining <- seq_len(C); ext <- integer(0)
    while (length(remaining)) {                 # random topological order
      maximal <- remaining[vapply(remaining, function(b)
        !any(D[remaining, b]), logical(1))]
      pick <- if (length(maximal) == 1L) maximal else sample(maximal, 1L)
      ext <- c(ext, pick); remaining <- setdiff(remaining, pick)
    }
    L[ext, j] <- vals                           # largest value -> first placed
  }
  L
}

# Every incomparable pair must cross in BOTH directions by >= margin on the
# probability scale (mean over items of the one-sided exceedance).
.po_crossing_ok <- function(P, D, margin) {
  C <- nrow(P)
  for (a in seq_len(C - 1L)) for (b in (a + 1L):C) {
    if (D[a, b] || D[b, a]) next
    if (mean(pmax(0, P[a, ] - P[b, ])) < margin ||
        mean(pmax(0, P[b, ] - P[a, ])) < margin) return(FALSE)
  }
  TRUE
}

# Item-side poset resolution for "PO_ITEMS": keywords over ITEMS, or a J x J
# dominance matrix (validated as for classes). "layers2": the first half of the
# items dominates the second half; halves internally incomparable.
.po_resolve_item_poset <- function(poset, J) {
  if (is.character(poset) && length(poset) == 1L) {
    D <- matrix(FALSE, J, J)
    if (identical(poset, "layers2") || identical(poset, "V")) {
      if (J < 4L) stop("item poset \"layers2\" needs n_items >= 4")
      top <- seq_len(J %/% 2L); bot <- setdiff(seq_len(J), top)
      D[top, bot] <- TRUE
    } else if (identical(poset, "single")) {
      if (J < 3L) stop("item poset \"single\" needs n_items >= 3")
      D[1L, 2L] <- TRUE
    } else stop("unknown item poset keyword: ", poset,
                " (use \"layers2\", \"single\", or a J x J matrix)")
    D
  } else .po_resolve_poset(poset, J)
}

# Crossing depth for an incomparable ITEM pair: the MAX single-class reversal
# in each direction (with only C classes, a mean-over-classes requirement is
# jointly near-unsatisfiable across many pairs; a crossing exists with the
# depth of its largest reversal). Returns the minimum depth over all
# incomparable pairs and both directions.
.po_item_min_crossing <- function(P, D) {
  J <- ncol(P); m <- Inf
  for (j in seq_len(J - 1L)) for (k in (j + 1L):J) {
    if (D[j, k] || D[k, j]) next
    m <- min(m, max(P[, j] - P[, k]), max(P[, k] - P[, j]))
  }
  m
}
