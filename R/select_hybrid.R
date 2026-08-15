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
#             is imposed (simulate under the fitted IIO model - IIO alone, the
#             hypothesis under test). Small p => IIO violated.
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
    a <- ord[ai]; b <- ord[bi]
    okp <- !is.na(data[, a]) & !is.na(data[, b])  # pairwise complete (masked
    rest <- total[okp] - data[okp, a] - data[okp, b]  # data; identical to the
    da <- data[okp, a]; db <- data[okp, b]        # old path when nothing is NA)
    for (r in sort(unique(rest))) {
      g <- rest == r; if (sum(g) < min_group) next
      d <- mean(db[g]) - mean(da[g])              # overall-harder item now easier
      if (d > 0) mag <- mag + d * sum(g)
    }
  }
  mag / N
}

.manifest_iio_holds <- function(data, C, B, n_starts, use_cpp, seed,
                                null_C_range = 2:6, B_refine = 999L,
                                refine_band = NULL) {
  obs <- .manifest_iio_stat(data)
  # Null simulated under the FITTED IIO model - the hypothesis is IIO alone.
  # (Simulating from DM, as before the 2026-07-28 audit, additionally imposed
  # class monotonicity: a misspecified, overly restrictive composite null for
  # IIO-but-non-MON data.)
  #
  # NULL CLASS COUNT (2026-08-14): the null's class count is selected by BIC
  # over `null_C_range`, floored at the routing C - it is NOT the 2x2's fixed
  # C. Evidence: on 4-class DM data at J = 48 the C = 3 null is too tight
  # (the fitted 3-class IIO model cannot reproduce the data's crossing level,
  # so simulated statistics sit systematically below the observed one) and
  # the axis falsely rejects: TA173 p = 0.020 at C = 3 vs 0.260 at C = 4,
  # TA220 p = 0.020 vs 0.220, with the OBSERVED statistic identical (it is
  # model-free). The bias exists at every J but only becomes decisive when
  # the statistic aggregates over many item pairs (1128 at J = 48 vs 28 at
  # J = 8), which is why fixed C = 3 passed validation at J <= 24 (DM
  # recovery 97/100/93%) and collapses at J = 48. Routing still uses the
  # fixed C; only the reference distribution adapts, so the stability
  # rationale for a small fixed C in the property tests is preserved.
  fit_at <- function(cc) tryCatch(
    refit_model_type("IIO", as.matrix(data), cc, n_starts, use_cpp),
    error = function(e) NULL)
  cand <- sort(unique(c(C, as.integer(null_C_range))))
  cand <- cand[cand >= C]
  fits <- lapply(cand, fit_at)
  bics <- vapply(fits, function(f) if (is.null(f)) NA_real_ else BIC(f),
                 numeric(1))
  usable <- which(is.finite(bics))
  if (!length(usable)) return(list(p = NA_real_, holds = NA, null_C = NA))
  pick <- usable[which.min(bics[usable])]
  iio_fit <- fits[[pick]]; null_C <- cand[pick]
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(data)
  draw <- function(m) vapply(seq_len(m), function(b)
    .manifest_iio_stat(.impose_mask(simulate_from_qlfit(iio_fit, n), data)),
    numeric(1))
  null <- draw(B)
  p <- (1 + sum(null >= obs)) / (length(null) + 1)
  # ADAPTIVE PRECISION (2026-08-15): a bootstrap p-value near the decision
  # threshold has SE ~ sqrt(p(1-p)/B) = 0.031 at B = 49, so borderline
  # decisions are not reproducible across seeds (measured: only 2 of 9
  # borderline TI&D datasets gave the same verdict over 5 seeds at B = 49;
  # 6 of 9 at B = 199). When p lands in the indeterminate band around alpha,
  # extend the same null to B_refine draws (SE ~ 0.007 at 999) and decide on
  # the pooled distribution. Decisive datasets pay nothing.
  refined <- FALSE
  # ANALYTIC refinement band - no tuned constants. A bootstrap p-value has
  # SE = sqrt(p(1-p)/B), so decisions are unreliable when p lies within a
  # few SEs of the threshold. Lower edge: the Monte-Carlo floor guard
  # (1.5/(B+1)) - a p-value at the floor is the most decisive rejection
  # available and cannot benefit from more draws. Upper edge:
  # alpha + 3*SE(alpha, B). Both derive from alpha and B alone.
  alpha_axis <- 0.05
  if (is.null(refine_band))
    refine_band <- c(1.5 / (B + 1),
                     alpha_axis + 3 * sqrt(alpha_axis * (1 - alpha_axis) / B))
  lo <- max(refine_band[1], 1.5 / (B + 1))
  if (!is.null(B_refine) && B_refine > B &&
      p >= lo && p <= refine_band[2]) {
    null <- c(null, draw(B_refine - B))
    p <- (1 + sum(null >= obs)) / (length(null) + 1)
    refined <- TRUE
  }
  list(stat = obs, p = p, holds = p > 0.05, null_C = null_C,
       B_eff = length(null), refined = refined)
}

# item_probs comes back as a J x C matrix from the dichotomous engine but as a
# per-item list of class x category matrices from the masked/polytomous engine
# (any NA in the data routes there). Normalise to the J x C expected-score
# matrix on [0, 1] (expected score / (K - 1)) - identical to P(X = 1) in the
# dichotomous case, and keeping every dominance tolerance (eps) on one scale
# for polytomous items.
.hyb_item_probs <- function(fit) {
  P <- fit$item_probs
  if (is.matrix(P)) return(P)
  do.call(rbind, lapply(P, function(m)
    as.numeric(m %*% (seq_len(ncol(m)) - 1L)) / max(1L, ncol(m) - 1L)))
}

# --- POSET refinements (class side AND item side) ---------------------------
# The six-model set asserts TOTAL orders: MON on classes, IIO on items. Each
# negation conflates an ANTICHAIN (genuinely unordered) with a PARTIAL ORDER
# (real, incomplete order structure with no rung in the framework). The
# refinement is therefore SYMMETRIC:
#
#   UN cell  (neither order)  -> class poset AND item poset
#   IIO cell (items ordered)  -> class poset (is the class side truly unordered?)
#   MON cell (classes ordered)-> item poset  (is the item side truly unordered?)
#   DM cell  (both ordered)   -> nothing to refine
#
# Dominance mirrors the calibrated statistics: class a dominates b iff the
# per-item mean downward violation is <= eps; item j dominates k iff the
# per-class mean downward violation is <= eps. The CLASS statistic must be
# label-invariant (class labels permute across bootstrap refits), so it is the
# count of strictly comparable pairs; the ITEM statistic is invariant for free
# (item identity is observed, and quantifying over all classes cancels class
# relabelling). Both report the resampling lower bound (q05) of the count, and
# both are CAPPED at "partial": the refinement only runs in cells where the
# corresponding total order was already rejected, so claiming a full chain
# would contradict the axis that routed us there.
#
# For C = 3 the class poset also reports its isomorphism TYPE from the
# label-invariant degree profile of the point-estimate dominance digraph:
# "single" (one comparable pair), "V" (one class dominates the two others),
# "Lambda" (two classes dominate one). Point-estimate based - treat as
# descriptive, not calibrated.
.class_dominance_matrix <- function(P, eps) {    # P = items x classes
  C <- ncol(P); D <- matrix(FALSE, C, C)
  for (a in seq_len(C)) for (b in seq_len(C)) if (a != b)
    D[a, b] <- mean(pmax(0, P[, b] - P[, a])) <= eps
  D & !t(D)                                       # strict dominance only
}
.item_dominance_count <- function(P, eps) {      # P = items x classes
  J <- nrow(P); n <- 0L
  for (j in seq_len(J - 1L)) for (k in (j + 1L):J) {
    jk <- mean(pmax(0, P[k, ] - P[j, ])) <= eps   # j dominates k
    kj <- mean(pmax(0, P[j, ] - P[k, ])) <= eps
    if (xor(jk, kj)) n <- n + 1L
  }
  n
}
.class_poset_type <- function(D) {               # D strict dominance, C = 3
  if (nrow(D) != 3L) return(NA_character_)
  npairs <- sum(D)
  if (npairs == 1L) return("single")
  if (npairs == 2L) {
    if (any(rowSums(D) == 2L)) return("V")
    if (any(colSums(D) == 2L)) return("Lambda")
  }
  NA_character_
}
# DOMINANCE-DEMONSTRATION POSET TESTS (redesigned 2026-07-29, external review).
# The earlier "dominance-asymmetry" statistic algebraically collapsed to the
# absolute difference of side-average probabilities (|mean(d+)| - |mean(d-)|
# = |mean(d)|) and carried NO crossing information; its permutation null
# tested profile exchangeability, not the antichain hypothesis. Both are gone.
# The redesign tests the ESTIMAND directly:
#
#   For a directed pair (x, y), the violation mass V(x >= y) - the mean over
#   the opposing side of pmax(0, P_y - P_x) - is ZERO iff x weakly dominates y.
#   ONE parametric bootstrap from the fitted UN table (B refits, class labels
#   aligned back to the observed fit by best-permutation profile matching)
#   yields the sampling distribution of every directed V. A pair is declared
#   COMPARABLE iff, under SIMULTANEOUS max-statistic bootstrap bounds over all
#   directed pairs of its side (one calibrated (1 - alpha/2) quantile of the
#   max deviation per direction - resolvable at modest B, unlike per-pair
#   Bonferroni tails), the upper bound of one direction's violation is <= eps
#   (dominance demonstrated to the axis tolerance) while the lower bound of
#   the reverse is > eps (a real gap - two identical profiles are
#   incomparable, not mutually dominant).
#
# "partial" = at least one demonstrated pair; "antichain" = none demonstrated
# (a failure to demonstrate, not a certified absence - stated as such in the
# docs). There is no permutation null: the verdict is a simultaneous
# demonstration of specific dominances with familywise control, conservative
# under ANY antichain - crossing profiles, unequal averages, or identical
# profiles alike. eps is the axis-calibrated tolerance (mon_eps, N-scaled),
# entering exactly as in the MON axis. Refits that fail are dropped and
# reported via b_eff; below max(20, B/2) the refinement refuses to answer.
# The demonstrated relation is PAIRWISE eps-dominance, which need not be
# transitive; this check verifies it assembles into a strict partial order
# (acyclic; transitive closure adds no edge beyond the demonstrated set).
# Package-level so tests exercise the PRODUCTION logic, not a copy.
.poset_relation_ok <- function(pr, nn) {
  if (!nrow(pr)) return(TRUE)
  D <- matrix(FALSE, nn, nn); D[cbind(pr$dominant, pr$dominated)] <- TRUE
  R <- D
  for (i in seq_len(nn)) R <- R | ((R %*% R) > 0)    # transitive closure
  if (any(diag(R))) return(FALSE)                    # cycle
  all((R & !D) == FALSE)                             # closure adds no edge
}

.poset_align <- function(Pb, P0) {
  C <- ncol(P0)
  perms <- .class_orderings(C)
  best <- perms[1L, ]; bestd <- Inf
  for (i in seq_len(nrow(perms))) {
    d <- sum(abs(Pb[, perms[i, ], drop = FALSE] - P0))
    if (d < bestd) { bestd <- d; best <- perms[i, ] }
  }
  Pb[, best, drop = FALSE]
}

.poset_refine <- function(data, C, B, n_starts, use_cpp, eps, alpha, seed,
                          sides = c("class", "item"), poset_B = NULL) {
  data <- as.matrix(data); n <- nrow(data)
  # class side needs C >= 3: at C = 2 there is nothing between antichain and
  # chain, so "partial" does not exist and reporting it would be misleading
  if (C < 3L) sides <- setdiff(sides, "class")
  if (!length(sides)) return(NULL)
  fit0 <- refit_model_type("UN", data, C, n_starts, use_cpp)
  P0 <- .hyb_item_probs(fit0)                      # items x classes
  pc <- pmax(fit0$class_probs, 0); pc <- pc / sum(pc)
  J <- nrow(P0)
  if (!is.null(seed)) set.seed(seed)
  # Inferential resolution: the simultaneous construction below calibrates ONE
  # (1 - alpha/2) quantile of a max-deviation statistic, so B need not scale
  # with the pair count - but the top order statistic of B = 49 is thin, so
  # the refinement enforces its own floor of 99 replicates.
  if (!is.null(poset_B) && (length(poset_B) != 1L || !is.finite(poset_B) ||
                            poset_B < 1))
    stop("poset_B must be a finite positive scalar")
  B_pos <- if (is.null(poset_B)) max(B, 99L) else max(as.integer(poset_B), 99L)
  poly <- !is.matrix(fit0$item_probs)              # masked / polytomous engine
  sim_one <- function() {
    cls <- sample.int(C, n, replace = TRUE, prob = pc)
    if (poly) {
      # categorical sampling from the fitted per-category probabilities -
      # expected scores are NOT probabilities and must never feed rbinom
      sim <- matrix(0L, n, J)
      for (j in seq_len(J)) {
        Pm <- fit0$item_probs[[j]]                 # C x K categories
        for (cc in seq_len(C)) {
          who <- which(cls == cc)
          if (length(who))
            sim[who, j] <- sample.int(ncol(Pm), length(who), replace = TRUE,
                                      prob = Pm[cc, ]) - 1L
        }
      }
    } else {
      sim <- matrix(stats::rbinom(n * J, 1, t(P0)[cls, , drop = FALSE]), n, J)
      storage.mode(sim) <- "integer"
    }
    .impose_mask(sim, data)                        # share observed missingness
  }
  boots <- vector("list", B_pos)
  for (b in seq_len(B_pos)) {
    f <- tryCatch(refit_model_type("UN", sim_one(), C, n_starts, use_cpp),
                  error = function(e) NULL)
    boots[[b]] <- if (is.null(f)) NULL else .poset_align(.hyb_item_probs(f), P0)
  }
  boots <- Filter(Negate(is.null), boots)
  b_eff <- length(boots)
  if (b_eff < max(20L, B_pos %/% 2L))
    stop("poset refinement: only ", b_eff, "/", B_pos, " bootstrap UN refits ",
         "succeeded - refusing to report a poset verdict from a broken ",
         "reference distribution")
  # SIMULTANEOUS STUDENTIZED MAX-STATISTIC BOUNDS (2026-07-29 round 5;
  # replaces per-pair Bonferroni tails, which finite B cannot resolve once the
  # pair count is large): with Vhat the observed directed violation masses,
  # V^(b) their aligned bootstrap replicates and s_k each direction's
  # bootstrap sd,
  #   u* = q_{1-alpha/2} of max_k (V^(b)_k - Vhat_k) / s_k
  #   l* = q_{1-alpha/2} of max_k (Vhat_k - V^(b)_k) / s_k
  # give simultaneous bounds Vhat_k + u* s_k and Vhat_k - l* s_k for the whole
  # family at once. x-dominates-y is DEMONSTRATED iff
  #   Vhat_xy + u* s_xy <= eps  (violation of "x >= y" excluded to tolerance)
  #   Vhat_yx - l* s_yx >  eps  (the reverse direction carries a real gap).
  # One calibrated quantile per family, resolvable at B ~ 99 regardless of m;
  # studentization stops loose (crossing) directions inflating the bound for
  # tight (dominated) ones.
  decide_side <- function(idx, viol) {
    m <- nrow(idx)
    # collect observed and bootstrap V for all 2m directed pairs
    vhat <- numeric(2L * m)
    for (k in seq_len(m)) {
      vhat[2L * k - 1L] <- viol(P0, idx[k, 1], idx[k, 2])
      vhat[2L * k]      <- viol(P0, idx[k, 2], idx[k, 1])
    }
    Vb <- matrix(0, b_eff, 2L * m)
    for (b in seq_len(b_eff)) {
      Pb <- boots[[b]]
      for (k in seq_len(m)) {
        Vb[b, 2L * k - 1L] <- viol(Pb, idx[k, 1], idx[k, 2])
        Vb[b, 2L * k]      <- viol(Pb, idx[k, 2], idx[k, 1])
      }
    }
    # STUDENTIZED max-statistic: per-direction bootstrap scale, one
    # simultaneous quantile. Without studentization the loose (crossing)
    # directions inflate the max and destroy power for the tight (dominated)
    # ones; with it each direction gets a bound proportional to its own
    # sampling noise. Zero-variance directions (violation identically 0
    # across refits) get a floor so the bound collapses to vhat itself.
    sk <- pmax(apply(Vb, 2L, stats::sd), 1e-6)
    dev <- sweep(sweep(Vb, 2L, vhat, "-"), 2L, sk, "/")
    u_q <- stats::quantile(apply(dev,  1L, max), 1 - alpha / 2, names = FALSE)
    l_q <- stats::quantile(apply(-dev, 1L, max), 1 - alpha / 2, names = FALSE)
    up <- vhat + u_q * sk                          # simultaneous upper bounds
    lo <- vhat - l_q * sk                          # simultaneous lower bounds
    pr <- data.frame(dominant = integer(0), dominated = integer(0))
    for (k in seq_len(m)) {
      x <- idx[k, 1]; y <- idx[k, 2]
      x_dom <- up[2L * k - 1L] <= eps && lo[2L * k] > eps
      y_dom <- up[2L * k] <= eps && lo[2L * k - 1L] > eps
      if (xor(x_dom, y_dom))
        pr <- rbind(pr, data.frame(dominant = if (x_dom) x else y,
                                   dominated = if (x_dom) y else x))
    }
    pr
  }
  out <- list()
  if ("class" %in% sides) {
    idx <- t(utils::combn(C, 2L))
    viol <- function(P, x, y) mean(pmax(0, P[, y] - P[, x]))   # V(x >= y)
    pr <- decide_side(idx, viol)
    D <- matrix(FALSE, C, C)
    if (nrow(pr)) D[cbind(pr$dominant, pr$dominated)] <- TRUE
    tr_ok <- .poset_relation_ok(pr, C)
    # a non-transitive demonstrated set is a collection of pairwise
    # eps-dominances, NOT a partial order - labelled as its own shape
    shape <- if (nrow(pr) == 0L) "antichain" else
             if (tr_ok) "partial" else "nontransitive_dominance"
    out$class <- list(comparable = nrow(pr), total = nrow(idx), pairs = pr,
      b_eff = b_eff, eps = eps, shape = shape, transitive = tr_ok,
      type = if (identical(shape, "partial"))
        .class_poset_type(D) else NA_character_)
  }
  if ("item" %in% sides) {
    idx <- t(utils::combn(J, 2L))
    viol <- function(P, x, y) mean(pmax(0, P[y, ] - P[x, ]))   # V(x >= y)
    pr <- decide_side(idx, viol)
    tr_ok <- .poset_relation_ok(pr, J)
    out$item <- list(comparable = nrow(pr), total = nrow(idx), pairs = pr,
      b_eff = b_eff, eps = eps, transitive = tr_ok,
      shape = if (nrow(pr) == 0L) "antichain" else
              if (tr_ok) "partial" else "nontransitive_dominance")
  }
  out
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
  P <- .hyb_item_probs(un)                        # items x classes
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
              lcr_vs_dm = NULL, rm_vs_lcr = NULL, rm_stage_failed = FALSE)
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
  if (!adequate && is.finite(q95) && q95 < min_effect &&
      t$statistic < min_effect) adequate <- TRUE   # both scales negligible
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
  if (is.null(rm_fit) || is.null(lcr_best)) {
    out$rm_stage_failed <- TRUE      # LCR verdict reflects the bridge test
    return(out)                      # only; RM continuity was NOT assessed
  }
  if (verbose) cat("Quantitative edge: RM vs LCR ...\n")
  rl <- tryCatch(rm_vs_lcr_test(data, rm_fit, lcr_best, profile_C, B = B,
           C_range = grain_grid, alpha = alpha, n_starts = boot_n_starts,
           use_cpp = use_cpp, observed_fits = lcr_profile, mc.cores = mc.cores,
           seed = sd_(8000L)), error = function(e) NULL)
  out$rm_vs_lcr <- rl
  # the calibrated comparison is unavailable BOTH when the call errors (NULL)
  # and when it returns available = FALSE (e.g. every bootstrap replicate
  # failed) - in either case the discrete/continuous verdict rests on the raw
  # BIC fallback, not the calibrated test, and must be flagged
  if (is.null(rl) || !isTRUE(rl$available)) out$rm_stage_failed <- TRUE
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
#' Diakow's own constraint-presence logic), then decides DM -> quantitative with
#' the same likelihood-ratio machinery as [select_model_ll()]'s quantitative
#' edge (`ll_equivalence_test` at the Lindsay bridge grain, then
#' `rm_vs_lcr_test`), run directly - without refitting the ordinal lattice.
#'
#' The 2x2 doubles IIO recovery relative to the LR-edge lattice, which loses
#' power precisely where DM and IIO coincide (complete paired K=20 grid, J=8:
#' IIO 18/20 vs 9/20, McNemar 9:0 discordant, p = 0.0039; on real TI&D data at
#' full N = 5000, 13/18 vs 7/18, p = 0.070); the LR edge in turn holds the
#' DM/quantitative boundary that a double-cancellation gate leaks across
#' (DM 10/10 vs 7/10).
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
#' @param data Binary response matrix (persons x items). `NA` entries are
#'   permitted and treated as missing at random: all fits use the masked
#'   likelihood, the manifest IIO statistic is computed on pairwise-complete
#'   responses, and every bootstrap / permutation-null replicate has the
#'   observed missingness mask re-imposed (rank-matched, as in
#'   [select_model_ll()]), so the calibrations see the same masking the
#'   observed statistics do.
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
#' @param alpha Significance level for the quantitative-edge decisions
#'   (LCR vs DM and RM vs LCR; default 0.05). The 2x2 axes use their own
#'   calibrations (`mon_eps` for MON; the IIO axis's bootstrap p at 0.05).
#' @param iio_B_refine Bootstrap draws used to RE-DECIDE an IIO-axis
#'   p-value that lands in `iio_refine_band` around the threshold (default
#'   999; `NULL` disables). A bootstrap p-value near 0.05 has standard error
#'   ~0.031 at B = 49, so borderline decisions are not reproducible across
#'   seeds - measured on TI&D data, only 2 of 9 borderline datasets gave the
#'   same verdict over 5 seeds at B = 49 (6 of 9 at B = 199). Decisive
#'   datasets never pay this cost.
#' @param iio_refine_band Indeterminate p-value band triggering the refined
#'   draw. Default `NULL` = ANALYTIC: lower edge `1.5/(B+1)` (the
#'   Monte-Carlo floor guard; a floor p-value is the most decisive rejection
#'   possible), upper edge `alpha + 3*sqrt(alpha(1-alpha)/B)` (three
#'   binomial standard errors above the threshold). Both derive from alpha
#'   and B alone - no data-tuned constants.
#' @param iio_null_C_range Class counts considered for the IIO axis's NULL
#'   model, selected among them by BIC and floored at `n_classes` (default
#'   `2:6`). The axis statistic is model-free; only its reference
#'   distribution uses this fit. A null class count below the data's true
#'   heterogeneity is too tight and inflates false rejections - decisively so
#'   on long tests, where the statistic aggregates over many item pairs (see
#'   `tid_ni48_finding.md`). Set to `n_classes` to restore the pre-2026-08-14
#'   fixed-C behaviour.
#' @param poset_B Bootstrap replicates for the poset refinement's aligned
#'   UN bootstrap (default `NULL` = `max(B, 99)`). The 99-replicate floor is
#'   ENFORCED - explicit smaller values are raised to 99. Raise to 499+ for a
#'   final analysis - the simultaneous quantile is a top-order statistic and
#'   deepens with B.
#' @param min_effect Practical-equivalence threshold (log-likelihood-ratio
#'   units) for the degenerate-null retention on the LCR-vs-DM edge: a
#'   significant rejection is overridden ONLY when both the observed LR and
#'   the null's 95th percentile fall below it - i.e. the effect and its
#'   reference scale are both negligible. Default 1 (matching
#'   [select_model_ll()]); set 0 to disable the retention entirely. A
#'   judgment call, deliberately public.
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
#'   two axis results (`iio`, `mon`), `poset` (partial-order refinements of the
#'   ordinal verdict: `$class` and/or `$item`, each with `shape` ("partial"
#'   iff at least one dominance pair is DEMONSTRATED under simultaneous
#'   max-statistic bootstrap bounds from an aligned parametric bootstrap of
#'   the fitted UN table (internal floor of 99 replicates); "antichain"
#'   records a failure to demonstrate any pair, not a certified absence), `comparable` (demonstrated pair count), `total`,
#'   `pairs` (the demonstrated dominances), `b_eff` (successful bootstrap
#'   refits), `eps` (the axis-calibrated tolerance used), and for the C = 3
#'   class side a descriptive isomorphism `type` of "single"/"V"/"Lambda";
#'   class+item in the UN cell, class only in the IIO cell, item only in the
#'   MON cell), and, when DM was
#'   reached, `quant` holding the quantitative-edge evidence (`supports_quant`,
#'   `bridge_C`, `lcr_vs_dm`, `rm_vs_lcr`, and the selected `LCR`/`RM` fit).
#' @seealso [select_model_ll()] for the LR-edge lattice on its own;
#'   [cc_bootstrap_null()] for the conjoint-cancellation route.
#' @export
select_model_hybrid <- function(data, n_classes = 3L, B = 49L, n_starts = 5L,
                                mon_eps = 0.01, alpha = 0.05, min_effect = 1,
                                iio_null_C_range = 2:6,
                                iio_B_refine = 999L,
                                iio_refine_band = NULL,
                                poset_B = NULL,
                                lr_boot_n_starts = 2L, lr_n_classes = 2:6,
                                use_cpp = TRUE,
                                mc.cores = 1L, seed = NULL, verbose = FALSE) {
  data <- validate_data_any(data, allow_na = TRUE)
  if (anyNA(data)) {
    message("Data contain NA responses: fits use the masked likelihood and ",
            "every bootstrap/null replicate preserves the observed ",
            "missingness mask (MAR assumed; see the missing-data design).")
  }
  C <- as.integer(stats::median(n_classes)); if (C < 2L) C <- 2L
  J <- ncol(data)
  s <- function(off) if (is.null(seed)) NULL else seed + off

  if (verbose) cat("Manifest ordinal layer: IIO axis...\n")
  iio <- .manifest_iio_holds(data, C, B, n_starts, use_cpp, s(0L),
                             null_C_range = iio_null_C_range,
                             B_refine = iio_B_refine,
                             refine_band = iio_refine_band)
  if (verbose) cat("Manifest ordinal layer: MON axis...\n")
  # mon_eps is a PER-ITEM tolerance on the perm-min downward-movement q05.
  # N-scale it: the resampling noise floor of the class-probability decrease
  # scales ~ 1/sqrt(N), so anchor mon_eps at N = 1500 and widen it at smaller
  # N (keeps size ~alpha across N; the calibration study showed the holds/
  # violated q05 gap sits at ~0.02 vs ~0.06 per item at N = 1500).
  mon_eps_N <- mon_eps * sqrt(1500 / nrow(data))
  mon <- .manifest_mon_holds(data, C, max(B %/% 2L, 20L), n_starts, use_cpp,
                             mon_eps_N, s(1000L))

  if (is.na(iio$holds))
    stop("IIO axis calibration failed (constrained IIO fit did not converge); ",
         "try more n_starts - refusing to treat numerical failure as an ",
         "IIO violation")
  ih <- isTRUE(iio$holds); mh <- isTRUE(mon$holds)
  ordinal <- if (ih && mh) "DM" else if (ih && !mh) "IIO" else
             if (!ih && mh) "MON" else "UN"
  interp <- c(UN = "CLASSIFICATORY (no ordinal structure)",
              MON = "ORDINAL (class monotonicity)",
              IIO = "ORDINAL (invariant item ordering)",
              DM = "ORDINAL (double monotonicity)")[ordinal]
  selected <- ordinal; interpretation <- interp

  # Poset refinements: distinguish genuine antichains from PARTIAL orders on
  # whichever side(s) the 2x2 left unordered. Never changes `selected`/`scale` -
  # a partial order is not one of TI&D's six models - it is reported as
  # additional structure in `poset` and in the interpretation.
  poset <- NULL; poset_refusal <- NULL
  poset_sides <- switch(ordinal, UN = c("class", "item"),
                        IIO = "class", MON = "item", NULL)
  if (!is.null(poset_sides) && identical(poset_sides, "class") && C < 3L) {
    poset_refusal <- sprintf(
      "class poset not applicable at C = %d (requires C >= 3)", C)
    poset_sides <- NULL
  }
  if (!is.null(poset_sides)) {
    if (verbose) cat(ordinal, "cell: poset refinement (",
                     paste(poset_sides, collapse = "+"), ") ...\n")
    poset <- tryCatch(.poset_refine(data, C, B, n_starts,
               use_cpp, mon_eps_N, alpha, s(4000L), sides = poset_sides,
               poset_B = poset_B),
               error = function(e) {
                 warning("poset refinement failed (", conditionMessage(e),
                         ") - no partial-order verdict is reported")
                 poset_refusal <<- conditionMessage(e)
                 NULL })
    notes <- character(0)
    nt <- character(0)
    if (!is.null(poset$class) &&
        identical(poset$class$shape, "nontransitive_dominance"))
      nt <- c(nt, sprintf("class side: %d pairwise dominances, NOT transitive",
                          poset$class$comparable))
    if (!is.null(poset$item) &&
        identical(poset$item$shape, "nontransitive_dominance"))
      nt <- c(nt, sprintf("item side: %d pairwise dominances, NOT transitive",
                          poset$item$comparable))
    if (!is.null(poset$class) && identical(poset$class$shape, "partial"))
      notes <- c(notes, paste0(poset$class$comparable, " of ", poset$class$total,
        " class pairs demonstrated (simultaneous ", signif(alpha, 2), ", B_eff ",
        poset$class$b_eff, ")",
        if (!is.na(poset$class$type)) paste0(" (", poset$class$type, ")") else ""))
    if (!is.null(poset$item) && identical(poset$item$shape, "partial"))
      notes <- c(notes, paste0(poset$item$comparable, " of ", poset$item$total,
        " item pairs demonstrated (simultaneous ", signif(alpha, 2), ", B_eff ",
        poset$item$b_eff, ")"))
    if (length(notes))
      interpretation <- paste0(interpretation, " + PARTIAL ORDER [",
        paste(notes, collapse = "; "),
        "] - real order structure the six-model set cannot express")
    if (length(nt))
      interpretation <- paste0(interpretation, " + [",
        paste(nt, collapse = "; "),
        "] - pairwise eps-dominance only, no partial-order claim")
  }

  # Quantitative sequence only from DM (order before quantity).
  # The DM -> quant decision is delegated to the LR-edge lattice
  # (select_model_ll). We take the lattice verdict ONLY when it is quantitative;
  # if the lattice finds no quant support the 2x2's DM stands (its ordinal-layer
  # call is authoritative).
  quant_edge_failed <- FALSE
  if (ordinal == "DM") {
    if (verbose) cat("DM vs quant: LR edge ...\n")
    add <- tryCatch(.hybrid_quant_edge(data, grain_range = lr_n_classes,
             B = B, n_starts = n_starts, boot_n_starts = lr_boot_n_starts,
             alpha = alpha, use_cpp = use_cpp, mc.cores = mc.cores,
             seed = if (!is.null(seed)) seed + 2000L else NULL,
             min_effect = min_effect, verbose = verbose),
             error = function(e) NULL)
    if (is.null(add) || is.na(add$supports_quant)) {
      # numerical failure, NOT evidence against quantity - flag it loudly
      # rather than letting DM stand as if the edge had tested and retained it
      warning("quantitative edge did not complete (bridge fits or LR ",
              "bootstrap failed); the DM verdict is the 2x2 ordinal call ",
              "only - quantitative structure was NOT assessed")
      quant_edge_failed <- TRUE
    }
    if (isTRUE(add$rm_stage_failed)) {
      warning("RM-vs-LCR stage did not complete (fit_rm or the comparison ",
              "failed); a quantitative verdict reflects the LCR-vs-DM bridge ",
              "test only - discrete-vs-continuous was NOT assessed")
      quant_edge_failed <- TRUE
    }
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
                 poset_refusal = poset_refusal,
                 quant_edge_failed = quant_edge_failed,
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
  if (!is.null(x$poset$class))
    cat(sprintf("Class poset     : %d/%d pairs demonstrated (B_eff %d, transitive %s) -> %s%s\n",
                x$poset$class$comparable, x$poset$class$total,
                x$poset$class$b_eff, x$poset$class$transitive,
                x$poset$class$shape,
                if (!is.na(x$poset$class$type))
                  paste0(" [", x$poset$class$type, "]") else ""))
  if (!is.null(x$poset$item))
    cat(sprintf("Item poset      : %d/%d pairs demonstrated (B_eff %d, transitive %s) -> %s\n",
                x$poset$item$comparable, x$poset$item$total,
                x$poset$item$b_eff, x$poset$item$transitive,
                x$poset$item$shape))
  if (!is.null(x$quant))
    cat(sprintf("DM vs quant (LR): supports_quant = %s -> %s\n",
                isTRUE(x$quant$supports_quant),
                if (isTRUE(x$quant$supports_quant)) "quantitative" else "stays DM"))
  invisible(x)
}
