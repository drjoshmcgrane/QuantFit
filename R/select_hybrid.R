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
  selected <- ordinal; interpretation <- interp; rl <- NULL

  # Quantitative sequence only from DM (order before quantity).
  # The DM -> quant decision is delegated to the LR-edge lattice
  # (select_model_ll). We take the lattice verdict ONLY when it is quantitative;
  # if the lattice finds no quant support the 2x2's DM stands (its ordinal-layer
  # call is authoritative).
  if (ordinal == "DM") {
    if (verbose) cat("DM vs quant: LR edge (lattice) ...\n")
    ll <- tryCatch(select_model_ll(data, n_classes = lr_n_classes, alpha = alpha,
             alpha_quant = alpha, B = B, n_starts = n_starts,
             boot_n_starts = lr_boot_n_starts, use_cpp = use_cpp,
             mc.cores = mc.cores,
             seed = if (!is.null(seed)) seed + 2000L else NULL,
             verbose = FALSE), error = function(e) NULL)
    add <- list(supports_quant = if (is.null(ll)) NA else ll$selected %in% c("LCR", "RM"),
                lr = ll)
    if (isTRUE(add$supports_quant)) {
      selected <- ll$selected; interpretation <- ll$interpretation
    }                                   # else stays DM (2x2 ordinal call)
  }

  scale <- c(UN = "nominal", MON = "ordinal", IIO = "ordinal", DM = "ordinal",
             LCR = "quant", RM = "quant")[selected]
  structure(list(selected = selected, interpretation = interpretation,
                 scale = unname(scale), n_classes = C,
                 iio = iio, mon = mon,
                 quant = if (exists("add", inherits = FALSE)) add else NULL,
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
  if (!is.null(x$quant))
    cat(sprintf("DM vs quant (LR): supports_quant = %s -> %s\n",
                isTRUE(x$quant$supports_quant),
                if (isTRUE(x$quant$supports_quant)) "quantitative" else "stays DM"))
  invisible(x)
}
