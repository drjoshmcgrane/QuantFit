# Manifest-2x2 ordinal selector (separate from the LR-edge lattice).
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
#   MON axis  class / person monotonicity. Order the fitted unconstrained
#             latent classes by mean success; how much does any item's
#             class-probability DECREASE across the ordered classes? Calibrated
#             by DATA RESAMPLING (the statistic's own sampling distribution),
#             NOT a parametric null - a parametric simulate-and-refit
#             reintroduces the very unconstrained-fit noise the statistic
#             avoids, collapsing power.
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
.manifest_mon_stat <- function(data, C, n_starts, use_cpp) {
  un <- refit_model_type("UN", as.matrix(data), C, n_starts, use_cpp)
  P <- un$item_probs
  P <- P[, order(colMeans(P)), drop = FALSE]     # classes easy -> hard
  sum(pmax(0, -t(apply(P, 1, diff))))            # total downward class movement
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

# --- ADD axis (ordinal DM vs quantitative): additive separability -----------
.manifest_add_stat <- function(un_fit) {
  P <- pmin(pmax(un_fit$item_probs, 1e-3), 1 - 1e-3)   # items x classes
  L <- stats::qlogis(P)
  mu <- mean(L); ri <- rowMeans(L) - mu; cc <- colMeans(L) - mu
  resid <- L - (outer(ri, cc, "+") + mu)               # remove additive main effects
  sqrt(mean(resid^2)) / sqrt(mean((L - mu)^2))         # scale-free interaction size
}
# parametric null under ADDITIVITY (Rasch form imposed on the class x item
# table): simulate class responses from the additive-projected probabilities,
# refit UN, recompute the interaction. Operates on the class structure only -
# no fit_rm. Observed interaction above the null 95th pct => non-additive => DM.
.manifest_add_holds <- function(data, C, B, n_starts, use_cpp, seed) {
  un <- refit_model_type("UN", as.matrix(data), C, n_starts, use_cpp)
  obs <- .manifest_add_stat(un)
  P <- pmin(pmax(un$item_probs, 1e-3), 1 - 1e-3)
  L <- stats::qlogis(P); mu <- mean(L)
  ri <- rowMeans(L) - mu; cc <- colMeans(L) - mu
  Padd <- stats::plogis(outer(ri, cc, "+") + mu)       # additive-projected probs
  pc <- un$class_probs; pc <- pmax(pc, 0); pc <- pc / sum(pc)
  n <- nrow(data); J <- ncol(data)
  if (!is.null(seed)) set.seed(seed)
  null <- vapply(seq_len(B), function(b) {
    cl <- sample.int(C, n, replace = TRUE, prob = pc)
    d2 <- matrix(stats::rbinom(n * J, 1, t(Padd)[cl, , drop = FALSE]), n, J)
    u2 <- tryCatch(refit_model_type("UN", d2, C, n_starts, use_cpp),
                   error = function(e) NULL)
    if (is.null(u2)) NA_real_ else .manifest_add_stat(u2)
  }, numeric(1))
  null <- null[!is.na(null)]
  p <- if (!length(null)) NA_real_ else (1 + sum(null >= obs)) / (length(null) + 1)
  list(stat = obs, p = p, additive = isTRUE(p > 0.05))
}

# --- DIP axis (LCR discrete vs RM continuous): latent-distribution shape -----
.manifest_dip_stat <- function(data, jseed) {
  s <- rowSums(data)
  if (!is.null(jseed)) set.seed(jseed)
  unname(diptest::dip.test(s + stats::runif(length(s), -0.4, 0.4))$statistic)
}
# calibrate the dip against a CONTINUOUS-theta (RM) null at the fitted item
# difficulties: observed dip above the null 95th pct => the score distribution
# is more multimodal than any continuous theta explains => LCR (discrete).
.manifest_dip_lcr <- function(data, rm_fit, B, seed) {
  obs <- .manifest_dip_stat(data, seed)
  beta <- rm_fit$delta; sigma <- sqrt(rasch_latent_var(rm_fit))
  n <- nrow(data); J <- ncol(data)
  if (!is.null(seed)) set.seed(seed)
  null <- vapply(seq_len(B), function(b) {
    th <- stats::rnorm(n, 0, sigma)
    d2 <- matrix(stats::rbinom(n * J, 1, stats::plogis(outer(th, beta, "-"))), n, J)
    .manifest_dip_stat(d2, b)
  }, numeric(1))
  p <- (1 + sum(null >= obs)) / (B + 1)
  list(stat = obs, p = p, discrete = isTRUE(p <= 0.05))   # small p => LCR
}

#' Manifest-2x2 latent-structure selector (property-based ordinal layer)
#'
#' An alternative to [select_model_ll()] that decides the ordinal / nominal
#' layer (UN, MON, IIO, DM) by testing the invariant-item-ordering and
#' class-monotonicity PROPERTIES directly against the data, then enters the
#' quantitative sequence DM -> LCR -> RM. On simulated data this roughly
#' doubles IIO recovery relative to the LR-edge lattice (which loses power when
#' DM and IIO coincide) and is far cheaper. Kept as a SEPARATE selector; it
#' does not change [select_model_ll()].
#'
#' @param data Binary response matrix (persons x items).
#' @param n_classes Number of latent classes (or a range; the median is used
#'   for the constraint fits). Default 3.
#' @param B Bootstrap replicates per axis and per quantitative edge (default 49).
#' @param n_starts Random starts for the constraint fits (default 5).
#' @param mon_eps Negligible class-probability decrease below which monotonicity
#'   is treated as holding (default 0.03).
#' @param alpha Significance level for the quantitative edges (default 0.05).
#' @param use_cpp Use the compiled EM engine.
#' @param mc.cores Cores for the bootstraps.
#' @param seed Optional integer seed.
#' @param verbose Print progress.
#' @return A list with `selected`, `interpretation`, `scale`, and the axis /
#'   edge evidence (`iio`, `mon`, and, when reached, `rm_vs_lcr`).
#' @seealso [select_model_ll()] for the LR-edge lattice.
#' @export
select_model_manifest <- function(data, n_classes = 3L, B = 49L, n_starts = 5L,
                                   mon_eps = 0.03, alpha = 0.05, use_cpp = TRUE,
                                   mc.cores = 1L, seed = NULL, verbose = FALSE) {
  data <- validate_data_any(data, allow_na = FALSE)
  C <- as.integer(stats::median(n_classes)); if (C < 2L) C <- 2L
  J <- ncol(data)
  s <- function(off) if (is.null(seed)) NULL else seed + off

  if (verbose) cat("Manifest ordinal layer: IIO axis...\n")
  iio <- .manifest_iio_holds(data, C, B, n_starts, use_cpp, s(0L))
  if (verbose) cat("Manifest ordinal layer: MON axis...\n")
  mon <- .manifest_mon_holds(data, C, max(B %/% 2L, 20L), n_starts, use_cpp,
                             mon_eps, s(1000L))

  ih <- isTRUE(iio$holds); mh <- isTRUE(mon$holds)
  ordinal <- if (ih && mh) "DM" else if (ih && !mh) "IIO" else
             if (!ih && mh) "MON" else "UN"
  interp <- c(UN = "CLASSIFICATORY (no ordinal structure)",
              MON = "ORDINAL (class monotonicity)",
              IIO = "ORDINAL (invariant item ordering)",
              DM = "ORDINAL (double monotonicity)")[ordinal]
  selected <- ordinal; interpretation <- interp; rl <- NULL

  # Quantitative sequence only from DM (order before quantity)
  if (ordinal == "DM") {
    if (verbose) cat("Quantitative axis: additivity (DM vs quant)...\n")
    # ADD axis: is the class x item logit table additively separable
    # (theta_c - beta_i)? Additive -> quantitative; not -> ordinal (DM).
    add <- .manifest_add_holds(data, C, B, n_starts, use_cpp, s(2000L))
    if (isTRUE(add$additive)) {
      # quantitative. DIP axis decides discrete (LCR) vs continuous (RM) from
      # the shape of the score distribution (sufficient for theta).
      if (verbose) cat("Quantitative axis: latent shape (LCR vs RM)...\n")
      rm_fit <- tryCatch(suppressWarnings(fit_rm(data, verbose = FALSE)),
                         error = function(e) NULL)
      if (!is.null(rm_fit)) {
        dip <- .manifest_dip_lcr(data, rm_fit, B, s(3000L))
        if (isTRUE(dip$discrete)) {
          selected <- "LCR"; interpretation <- "QUANTITATIVE (discrete: latent class Rasch)"
        } else {
          selected <- "RM"; interpretation <- "QUANTITATIVE (continuous: Rasch model)"
        }
      } else {
        selected <- "LCR"; interpretation <- "QUANTITATIVE (discrete: latent class Rasch)"
        dip <- NULL
      }
    } else add <- add   # not additive -> stays DM
  }

  scale <- c(UN = "nominal", MON = "ordinal", IIO = "ordinal", DM = "ordinal",
             LCR = "quant", RM = "quant")[selected]
  structure(list(selected = selected, interpretation = interpretation,
                 scale = unname(scale), n_classes = C,
                 iio = iio, mon = mon,
                 add = if (exists("add", inherits = FALSE)) add else NULL,
                 dip = if (exists("dip", inherits = FALSE)) dip else NULL,
                 method = "manifest-4axis"),
            class = "qlselect_manifest")
}

#' @export
print.qlselect_manifest <- function(x, ...) {
  cat("\nManifest-2x2 latent-structure selection\n")
  cat("---------------------------------------\n")
  cat(sprintf("Selected        : %s  [%s]\n", x$selected, x$scale))
  cat(sprintf("Interpretation  : %s\n", x$interpretation))
  cat(sprintf("IIO axis        : stat %.4f, p %.3f -> %s\n",
              x$iio$stat, x$iio$p, if (isTRUE(x$iio$holds)) "holds" else "violated"))
  cat(sprintf("MON axis        : stat %.4f, lo %.4f -> %s\n",
              x$mon$stat, x$mon$lo, if (isTRUE(x$mon$holds)) "holds" else "violated"))
  if (!is.null(x$add))
    cat(sprintf("ADD axis        : stat %.4f, p %.3f -> %s\n",
                x$add$stat, x$add$p, if (isTRUE(x$add$additive)) "additive (quant)" else "non-additive (DM)"))
  if (!is.null(x$dip))
    cat(sprintf("DIP axis        : stat %.4f, p %.3f -> %s\n",
                x$dip$stat, x$dip$p, if (isTRUE(x$dip$discrete)) "discrete (LCR)" else "continuous (RM)"))
  if (FALSE)
    cat(sprintf("RM vs LCR       : stat %.2f, p %.3f\n",
                x$rm_vs_lcr$statistic, x$rm_vs_lcr$p_value))
  invisible(x)
}
