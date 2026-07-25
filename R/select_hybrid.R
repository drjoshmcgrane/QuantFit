# HYBRID selector: property-based 2x2 ordinal layer + LR quantitative edge.
#
# The ordinal/nominal layer is decided by testing the two defining constraints
# directly; the DM -> quantitative decision is delegated to select_model_ll()'s
# likelihood-ratio edge. Both layers are latent-model based, deliberately: TI&D's
# hierarchy is defined on the latent class x item table and class monotonicity
# has no faithful manifest proxy (manifest_coherence_finding.md). An earlier
# variant closed DM -> quant with manifest double cancellation; that mixed an
# observable-conjoint axiom into a latent hierarchy and was removed. Double
# cancellation remains its own route (cc_bootstrap_null / cc_bootstrap_hierarchy).
#
# The ordinal / nominal layer (UN, MON, IIO, DM) is decided by testing the two
# defining PROPERTIES directly against the data, rather than comparing two
# near-equivalent constrained-model fits (the LR-edge approach, whose null
# degenerates when models coincide - notably DM vs IIO on doubly-monotone
# data, giving near-zero power). This mirrors Torres Irribarra & Diakow's own
# constraint-presence graphical logic, automated and calibrated:
#
#   IIO axis  invariant item ordering. Do item response functions cross when
#             persons are ordered by rest-score? A model-free crossing
#             magnitude, calibrated by a parametric null in which item ordering
#             is imposed (simulate under the fitted DM). Small p => IIO
#             violated.
#   MON axis  class / person monotonicity. MON holds iff SOME class ordering
#             makes every item's class-probability monotone, so the statistic
#             is the MINIMUM total downward movement over all class orderings -
#             the exact read of the "exists an order" definition, rather than
#             the mean-ordering heuristic. (Empirically the two coincide on
#             property-HOLDS data and differ only occasionally on violated data,
#             where perm-min is the slightly more conservative choice; it is
#             kept for definitional correctness, NOT because it fixed the
#             low-separation size problem - the per-item normalisation and the
#             eps recalibration did that.) Normalised per item so the scale is
#             J-invariant.
#             Calibrated by DATA RESAMPLING (the statistic's own sampling
#             distribution): a genuine crossing is structural and survives
#             resampling (q05 stays high); ordering noise on near-flat classes
#             shakes out (q05 collapses toward 0). Property holds if the q05
#             lower bound is at or below a per-item tolerance. A parametric
#             monotone null instead reintroduces unconstrained-fit noise and,
#             worse, fit_mon absorbs the IIO crossing - collapsing power.
#
#   IIO holds & MON holds -> DM ;  IIO holds & MON violated -> IIO
#   IIO violated & MON holds -> MON ;  both violated -> UN
#
# When DM is reached the quantitative sequence DM -> LCR -> RM is entered using
# the same calibrated machinery as [select_model_ll()].

# --- IIO axis --------------------------------------------------------------
.manifest_iio_stat <- function(data, min_group = 20L) {
  data <- as.matrix(data); J <- ncol(data); N <- nrow(data)
  p <- colMeans(data, na.rm = TRUE); ord <- order(p, decreasing = TRUE)
  total <- rowSums(data, na.rm = TRUE); mag <- 0
  for (ai in seq_len(J - 1L)) for (bi in (ai + 1L):J) {
    a <- ord[ai]; b <- ord[bi]; rest <- total - data[, a] - data[, b]
    for (r in sort(unique(rest))) {
      g <- rest == r; if (sum(g) < min_group) next
      d <- mean(data[g, b]) - mean(data[g, a])   # overall-harder item now easier
      if (d > 0) mag <- mag + d * sum(g)
    }
  }
  mag / N
}

.manifest_iio_holds <- function(data, C, B, n_starts, use_cpp, seed) {
  obs <- .manifest_iio_stat(data)
  dm <- tryCatch(refit_model_type("DM", as.matrix(data), C, n_starts, use_cpp),
                 error = function(e) NULL)
  if (is.null(dm)) return(list(p = NA_real_, holds = NA))
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(data)
  null <- vapply(seq_len(B), function(b)
    .manifest_iio_stat(.impose_mask(simulate_from_qlfit(dm, n), data)), numeric(1))
  p <- (1 + sum(null >= obs)) / (B + 1)
  list(stat = obs, p = p, holds = p > 0.05)
}

# --- POSET refinement of the UN cell ----------------------------------------
# MON asks whether a TOTAL order on classes exists. Its negation conflates two
# very different things: an ANTICHAIN (no class dominates any other -> genuinely
# nominal) and a PARTIAL ORDER (some classes dominate others, just not a chain
# -> strictly MORE structure than nominal). TI&D's model set has no rung for the
# latter, so a partial order is currently reported as "nominal", which
# understates what the data support.
#
# Dominance is defined consistently with the MON statistic - class a dominates b
# iff the per-item mean downward violation is negligible - so this reuses the
# calibrated mon_eps instead of adding another threshold.
#
# CALIBRATION. Class labels are not identified and permute across bootstrap
# refits, so the dominance MATRIX cannot be compared across resamples. The
# statistic is therefore the label-invariant COUNT of comparable pairs, and the
# reported structure is its resampling lower bound (q05), mirroring the MON
# axis: a dominance relation counts only if it survives resampling.
#
# On real TI&D data this cleanly separated correctly- from wrongly-nominal
# verdicts: of the datasets BOTH selectors called nominal, all 19 antichains
# were true UN, while all 5 misclassified datasets (4 IIO, 1 MON) came out as
# partial orders. See partial_order_finding.md.
.class_dominance <- function(P, eps) {           # P = items x classes
  C <- ncol(P); D <- matrix(FALSE, C, C)
  for (a in seq_len(C)) for (b in seq_len(C)) if (a != b)
    D[a, b] <- mean(pmax(0, P[, b] - P[, a])) <= eps
  n <- 0L
  for (a in seq_len(C - 1L)) for (b in (a + 1L):C)
    if (xor(D[a, b], D[b, a])) n <- n + 1L       # strictly comparable pair
  n
}
.poset_shape <- function(data, C, B, n_starts, use_cpp, eps, seed) {
  data <- as.matrix(data); n <- nrow(data)
  tot <- C * (C - 1L) / 2L
  obs <- .class_dominance(refit_model_type("UN", data, C, n_starts,
                                           use_cpp)$item_probs, eps)
  if (!is.null(seed)) set.seed(seed)
  boot <- vapply(seq_len(B), function(b) {
    f <- tryCatch(refit_model_type("UN",
           data[sample.int(n, n, replace = TRUE), , drop = FALSE], C,
           n_starts, use_cpp), error = function(e) NULL)
    if (is.null(f)) NA_real_ else .class_dominance(f$item_probs, eps)
  }, numeric(1))
  boot <- boot[!is.na(boot)]
  lo <- if (length(boot)) stats::quantile(boot, 0.05, names = FALSE) else 0
  shape <- if (lo >= tot) "chain" else if (lo >= 1) "partial" else "antichain"
  list(comparable = obs, total = tot, lo = lo, shape = shape)
}

# --- MON axis --------------------------------------------------------------
# All orderings of C classes (rows), for the exact "exists an order" read of
# MON. C is small in the manifest layer (2-4); guard larger C with the
# mean-ordering fallback so cost stays bounded.
.class_orderings <- function(C) {
  if (C <= 1L) return(matrix(1L, 1L, 1L))
  do.call(rbind, lapply(seq_len(C), function(i)
    cbind(i, matrix(setdiff(seq_len(C), i)[t(.class_orderings(C - 1L))],
                    ncol = C - 1L, byrow = TRUE))))
}
.manifest_mon_stat <- function(data, C, n_starts, use_cpp) {
  un <- refit_model_type("UN", as.matrix(data), C, n_starts, use_cpp)
  P <- un$item_probs                              # items x classes
  J <- nrow(P)
  down <- function(o) sum(pmax(0, -t(apply(P[, o, drop = FALSE], 1, diff))))
  ords <- if (C <= 5L) .class_orderings(C) else matrix(order(colMeans(P)), 1L)
  min(apply(ords, 1, down)) / J                   # min over orders, per item
}

.manifest_mon_holds <- function(data, C, B, n_starts, use_cpp, eps, seed) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(data)
  boot <- vapply(seq_len(B), function(b)
    .manifest_mon_stat(data[sample.int(n, n, replace = TRUE), , drop = FALSE],
                       C, n_starts, use_cpp), numeric(1))
  lo <- stats::quantile(boot, 0.05, names = FALSE)  # lower bound of the statistic
  list(stat = .manifest_mon_stat(data, C, n_starts, use_cpp), lo = lo,
       holds = lo <= eps)
}

# --- Quantitative edge: DM -> LCR -> RM -------------------------------------
# Runs ONLY the quantitative succession, which is all the 2x2 needs once it has
# established DM. Deliberately does NOT call select_model_ll(): that would
# re-select the class count by BIC, refit all six models and re-run the entire
# ordinal lattice with bootstraps, only to have its ordinal verdict discarded.
# That was pure waste, and the discarded verdict was not a usable second opinion
# anyway - the lattice's ordinal layer runs at a BIC-selected C while the 2x2
# runs at a fixed C, so any disagreement confounds class count with structure.
#
# The DM-vs-LCR comparison is made at the Lindsay bridge grain
# ceiling((score_max + 1) / 2) - the support-point count at which a latent-class
# Rasch model can represent the continuous Rasch model - so grain selection is
# removed from the scale test. Only if the bridge supports quantity is LCR
# profiled over grain and compared with the continuous RM.
.hybrid_quant_edge <- function(data, grain_range, B, n_starts, boot_n_starts,
                               alpha, use_cpp, mc.cores, seed, min_effect = 1,
                               verbose = FALSE) {
  data <- as.matrix(data)
  score_max <- if (!.is_polytomous(data)) ncol(data) else
    sum(vapply(seq_len(ncol(data)), function(j) {
      x <- data[, j]; if (all(is.na(x))) 0L else as.integer(max(x, na.rm = TRUE))
    }, integer(1)))
  bridge_C <- max(2L, as.integer(ceiling((score_max + 1) / 2)))
  grain_grid <- sort(unique(c(as.integer(grain_range[grain_range >= 2L]),
                              bridge_C)))
  sd_ <- function(off) if (is.null(seed)) NULL else seed + off
  ok <- function(lbl, expr) tryCatch(suppressWarnings(expr), error = function(e) {
    if (verbose) cat("Fit of", lbl, "failed:", conditionMessage(e), "\n"); NULL })

  dm_bridge  <- ok("DM bridge",  fit_dm(data, bridge_C, n_starts = n_starts,
                                        use_cpp = use_cpp, seed = sd_(5000L)))
  lcr_bridge <- ok("LCR bridge", fit_lcr(data, bridge_C, n_starts = n_starts,
                                         use_cpp = use_cpp, seed = sd_(6000L)))
  out <- list(supports_quant = NA, selected = NA_character_,
              bridge_C = bridge_C, grain_grid = grain_grid,
              lcr_vs_dm = NULL, rm_vs_lcr = NULL)
  if (is.null(dm_bridge) || is.null(lcr_bridge)) return(out)

  # LCR vs DM at the bridge grain, with the same degenerate-null retention the
  # lattice uses: when even the null's 95th percentile is a trivially small LR,
  # the bootstrap itself certifies the models indistinguishable, so retain the
  # constrained model by parsimony rather than reject on a negligible effect.
  if (verbose) cat("Quantitative edge: LCR vs DM at C =", bridge_C, "...\n")
  t <- tryCatch(ll_equivalence_test(data, lcr_bridge, dm_bridge, B = B,
           n_starts = boot_n_starts, seed = sd_(4000L), use_cpp = use_cpp,
           mc.cores = mc.cores), error = function(e) NULL)
  if (is.null(t)) return(out)
  adequate <- t$p_value > alpha
  q95 <- if (length(t$null_distribution))
    stats::quantile(t$null_distribution, 0.95, names = FALSE) else NA_real_
  if (!adequate && is.finite(q95) && q95 < min_effect) adequate <- TRUE
  out$lcr_vs_dm <- t
  out$supports_quant <- adequate
  if (!adequate) { out$selected <- "DM"; return(out) }

  # Quantitative: profile LCR over grain by BIC, then located/discrete LCR
  # versus continuous RM.
  lcr_profile <- lapply(grain_grid, function(C)
    if (C == bridge_C) lcr_bridge else
      ok(paste("LCR", C), fit_lcr(data, C, n_starts = n_starts,
                                  use_cpp = use_cpp, seed = sd_(7000L + C))))
  bic <- vapply(lcr_profile, function(f) if (is.null(f)) NA_real_ else BIC(f),
                numeric(1))
  usable <- which(is.finite(bic))
  pick <- if (length(usable)) usable[which.min(bic[usable])] else
    match(bridge_C, grain_grid)
  profile_C <- grain_grid[pick]; lcr_best <- lcr_profile[[pick]]
  rm_fit <- ok("RM", fit_rm(data, verbose = FALSE))
  out$selected <- "LCR"; out$C <- profile_C; out$LCR <- lcr_best
  if (is.null(rm_fit) || is.null(lcr_best)) return(out)
  if (verbose) cat("Quantitative edge: RM vs LCR ...\n")
  rl <- tryCatch(rm_vs_lcr_test(data, rm_fit, lcr_best, profile_C, B = B,
           C_range = grain_grid, alpha = alpha, n_starts = boot_n_starts,
           use_cpp = use_cpp, observed_fits = lcr_profile, mc.cores = mc.cores,
           seed = sd_(8000L)), error = function(e) NULL)
  out$rm_vs_lcr <- rl
  if (!is.null(rl)) {
    if (isTRUE(rl$available) && !is.na(rl$profiled_C)) out$C <- rl$profiled_C
    rm_pref <- if (!isTRUE(rl$available))
      (BIC(rm_fit) <= BIC(lcr_best)) else !rl$select_lcr
    if (isTRUE(rm_pref)) { out$selected <- "RM"; out$RM <- rm_fit }
  }
  out
}

#' Hybrid latent-structure selector (property-based 2x2 + LR quantitative edge)
#'
#' Decides the ordinal / nominal layer (UN, MON, IIO, DM) by testing the
#' invariant-item-ordering and class-monotonicity PROPERTIES directly against the
#' data (a 2x2 on the two defining constraints, mirroring Torres Irribarra &
#' Diakow's own constraint-presence logic), then delegates the DM -> quantitative
#' decision to the likelihood-ratio edge of [select_model_ll()].
#'
#' The 2x2 roughly doubles IIO recovery relative to the LR-edge lattice, which
#' loses power precisely where DM and IIO coincide (paired K=20 audit: IIO 12/12
#' vs 6/12, McNemar p = 0.031); the LR edge in turn holds the DM/quantitative
#' boundary that a double-cancellation gate leaks across (DM 10/10 vs 7/10).
#'
#' Both layers are latent-model based, which is deliberate. TI&D's six-model
#' hierarchy is defined by constraints on the latent class x item table, and
#' class monotonicity has NO faithful manifest proxy (Mokken monotone
#' homogeneity is a strictly weaker property and is near-blind to these
#' violations - see `manifest_coherence_finding.md`). An earlier variant closed
#' DM -> quant with a manifest double-cancellation test instead; that mixed an
#' observable-conjoint-array axiom into a latent hierarchy and was both
#' conceptually and empirically worse. It has been removed. Double cancellation
#' remains available as its own route via [cc_bootstrap_null()].
#'
#' @param data Binary response matrix (persons x items).
#' @param n_classes Number of latent classes (or a range; the median is used
#'   for the constraint fits). Default 3.
#' @param B Bootstrap replicates per axis and per quantitative edge (default 49).
#' @param n_starts Random starts for the constraint fits (default 5).
#' @param mon_eps Per-item tolerance on the perm-min downward-movement q05
#'   below which class monotonicity is treated as holding (default 0.01,
#'   anchored at N = 1500 and N-scaled internally). Set deliberately tight: on
#'   the standard well-separated generator a real IIO class-monotonicity
#'   violation gives q05 ~0.01-0.06, so a larger eps silently classifies IIO
#'   data as DM (0.04 collapsed IIO recovery to ~20\%). A tight eps keeps IIO
#'   detection; the cost is only at artificially low class separation
#'   (< ~1 logit), where weak IIO and near-collapsed DM are genuinely
#'   indistinguishable - an identifiability limit, not a tuning target.
#' @param lr_boot_n_starts Multistart count for the LR-edge bootstrap null refits
#'   (default 2, matching the validated lattice run; keeps a single quant-edge
#'   dataset tractable).
#' @param lr_n_classes Latent-class range passed to the LR delegation
#'   (default \code{2:6}, matching the standalone
#'   \code{select_model_ll} audit). The lattice selects the class count by BIC
#'   over this range; forcing a single C (e.g. the 2x2's C = 3) misspecifies the
#'   quant grain and biases the LCR-vs-DM edge toward retaining DM on
#'   short/finely-graded data (it cost ~3/4 of the J=8 quant misses). The LCR
#'   bridge grain ceiling((J+1)/2) is always added internally, so the quant edge
#'   is never capped below the Lindsay bound regardless of this range.
#' @param use_cpp Use the compiled EM engine.
#' @param mc.cores Cores for the bootstraps.
#' @param seed Optional integer seed.
#' @param verbose Print progress.
#' @return A list with `selected`, `interpretation`, `scale`, `n_classes`, the
#'   two axis results (`iio`, `mon`) and, when DM was reached, `quant` holding
#'   the LR-edge verdict (`supports_quant`) and the full [select_model_ll()]
#'   object (`lr`).
#' @seealso [select_model_ll()] for the LR-edge lattice on its own;
#'   [cc_bootstrap_null()] for the conjoint-cancellation route.
#' @export
select_model_hybrid <- function(data, n_classes = 3L, B = 49L, n_starts = 5L,
                                mon_eps = 0.01, alpha = 0.05,
                                lr_boot_n_starts = 2L, lr_n_classes = 2:6,
                                use_cpp = TRUE,
                                mc.cores = 1L, seed = NULL, verbose = FALSE) {
  data <- validate_data_any(data, allow_na = FALSE)
  C <- as.integer(stats::median(n_classes)); if (C < 2L) C <- 2L
  J <- ncol(data)
  s <- function(off) if (is.null(seed)) NULL else seed + off

  if (verbose) cat("Manifest ordinal layer: IIO axis...\n")
  iio <- .manifest_iio_holds(data, C, B, n_starts, use_cpp, s(0L))
  if (verbose) cat("Manifest ordinal layer: MON axis...\n")
  # mon_eps is a PER-ITEM tolerance on the perm-min downward-movement q05.
  # N-scale it: the resampling noise floor of the class-probability decrease
  # scales ~ 1/sqrt(N), so anchor mon_eps at N = 1500 and widen it at smaller
  # N (keeps size ~alpha across N; the calibration study showed the holds/
  # violated q05 gap sits at ~0.02 vs ~0.06 per item at N = 1500).
  mon_eps_N <- mon_eps * sqrt(1500 / nrow(data))
  mon <- .manifest_mon_holds(data, C, max(B %/% 2L, 20L), n_starts, use_cpp,
                             mon_eps_N, s(1000L))

  ih <- isTRUE(iio$holds); mh <- isTRUE(mon$holds)
  ordinal <- if (ih && mh) "DM" else if (ih && !mh) "IIO" else
             if (!ih && mh) "MON" else "UN"
  interp <- c(UN = "CLASSIFICATORY (no ordinal structure)",
              MON = "ORDINAL (class monotonicity)",
              IIO = "ORDINAL (invariant item ordering)",
              DM = "ORDINAL (double monotonicity)")[ordinal]
  selected <- ordinal; interpretation <- interp

  # Refine the UN cell: distinguish a genuine antichain (nominal) from a partial
  # order (real, incomplete ordinal structure that the six-model set cannot
  # express). This does NOT change `selected` or `scale` - the six models are
  # TI&D's and a partial order is not one of them - it is reported as additional
  # structure in `poset` and in the interpretation.
  poset <- NULL
  if (identical(ordinal, "UN")) {
    if (verbose) cat("UN cell: class-dominance poset ...\n")
    poset <- tryCatch(.poset_shape(data, C, max(B %/% 2L, 20L), n_starts,
               use_cpp, mon_eps_N, s(4000L)), error = function(e) NULL)
    if (!is.null(poset) && identical(poset$shape, "partial"))
      interpretation <- paste0("PARTIAL ORDER (", poset$lo, " of ",
        poset$total, " class pairs comparable) - more than nominal, but no ",
        "total order; not representable in the six-model set")
  }

  # Quantitative sequence only from DM (order before quantity).
  # The DM -> quant decision is delegated to the LR-edge lattice
  # (select_model_ll). We take the lattice verdict ONLY when it is quantitative;
  # if the lattice finds no quant support the 2x2's DM stands (its ordinal-layer
  # call is authoritative).
  if (ordinal == "DM") {
    if (verbose) cat("DM vs quant: LR edge ...\n")
    add <- tryCatch(.hybrid_quant_edge(data, grain_range = lr_n_classes,
             B = B, n_starts = n_starts, boot_n_starts = lr_boot_n_starts,
             alpha = alpha, use_cpp = use_cpp, mc.cores = mc.cores,
             seed = if (!is.null(seed)) seed + 2000L else NULL,
             verbose = verbose), error = function(e) NULL)
    if (isTRUE(add$supports_quant) && !is.na(add$selected)) {
      selected <- add$selected
      interpretation <- if (identical(selected, "LCR"))
        "QUANTITATIVE (discrete: latent class Rasch)" else
        "QUANTITATIVE (continuous: Rasch model)"
      C <- add$C
    }                                   # else stays DM (2x2 ordinal call)
  }

  scale <- c(UN = "nominal", MON = "ordinal", IIO = "ordinal", DM = "ordinal",
             LCR = "quant", RM = "quant")[selected]
  structure(list(selected = selected, interpretation = interpretation,
                 scale = unname(scale), n_classes = C,
                 iio = iio, mon = mon,
                 quant = if (exists("add", inherits = FALSE)) add else NULL,
                 poset = poset,
                 method = "hybrid-2x2+lr"),
            class = "qlselect_hybrid")
}

#' @export
print.qlselect_hybrid <- function(x, ...) {
  cat("\nHybrid latent-structure selection (2x2 ordinal layer + LR quant edge)\n")
  cat("---------------------------------------------------------------------\n")
  cat(sprintf("Selected        : %s  [%s]\n", x$selected, x$scale))
  cat(sprintf("Interpretation  : %s\n", x$interpretation))
  cat(sprintf("IIO axis        : stat %.4f, p %.3f -> %s\n",
              x$iio$stat, x$iio$p, if (isTRUE(x$iio$holds)) "holds" else "violated"))
  cat(sprintf("MON axis        : stat %.4f, lo %.4f -> %s\n",
              x$mon$stat, x$mon$lo, if (isTRUE(x$mon$holds)) "holds" else "violated"))
  if (!is.null(x$poset))
    cat(sprintf("Class poset     : %d/%d pairs comparable (q05 %g) -> %s\n",
                x$poset$comparable, x$poset$total, x$poset$lo, x$poset$shape))
  if (!is.null(x$quant))
    cat(sprintf("DM vs quant (LR): supports_quant = %s -> %s\n",
                isTRUE(x$quant$supports_quant),
                if (isTRUE(x$quant$supports_quant)) "quantitative" else "stays DM"))
  invisible(x)
}
