# Validation of the DOMINANCE-DEMONSTRATION partial-order tests (v4:
# studentized simultaneous max-statistic bounds; polytomous arm; XANTI at
# K=100 per length for familywise-size resolution; NEARANTI boundary class;
# demonstrated pair IDENTITIES recorded so subset claims are verifiable).
#
# History: v2 (collapsed statistic + exchangeability null) withdrawn; v3
# (per-pair Bonferroni tails) superseded - finite B cannot resolve alpha/(2m)
# quantiles at large m. v4 uses STUDENTIZED SIMULTANEOUS MAX-STATISTIC bounds:
# one calibrated (1-alpha/2) quantile of the studentized max deviation per
# family, per-direction scale, internal floor of 99 bootstrap refits.
#
# Truths and what each measures:
#   PO free (m x nI)      class-side POWER + type recovery (UN cell)
#   PO_INV (0.12 x nI)    class-side POWER under invariant items (IIO cell)
#   PO_ITEMS C=4 (nI 8,12) item-side POWER at feasible margins (MON cell)
#   PO_ITEMS C=4 nI=24    item-side BOUNDARY: achievable margin ~0.008 < eps
#                         0.01 - expected NON-detection, recorded as such
#   PO_ITEMS C=3 nI=8     item-side marginal case (margin ~0.012 ~ eps)
#   UN / IIO controls     class-side SIZE on genuine antichains (scored
#                         against POPULATION pairs from params$L)
#   MON control           item side vs population truth
#   CROSSING-ANTICHAIN    the reviewer counterexample class: profiles cross
#                         deeply, averages differ, zero dominance - SIZE for
#                         exactly the case that broke v2
suppressMessages(library(QuantFit))
`%||%` <- function(a, b) if (is.null(a)) b else a
K     <- as.integer(Sys.getenv("PO_K", "8"))
NR    <- as.integer(Sys.getenv("PO_N", "1500"))
B     <- as.integer(Sys.getenv("PO_B", "49"))
cores <- as.integer(Sys.getenv("PO_CORES", "8"))
out <- Sys.getenv("PO_OUT", file.path("..", "qf_evidence", "po_validation4_out"))
if (startsWith(paste0(normalizePath(out, mustWork = FALSE), "/"),
               paste0(normalizePath(getwd()), "/")))
  stop("PO_OUT must live OUTSIDE the git checkout (dirty-tree provenance)")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
# PROVENANCE CANARIES (review round 6): every row records the package git SHA,
# the inference method marker, and the poset bootstrap depth; the skip-if-
# exists resume check REFUSES to reuse a file whose method/SHA disagree with
# the current build, so superseded results can never be silently recycled.
METHOD <- "max-T-studentized-v4"
source("inst/validation/validation_shared.R")
prov <- qf_provenance_check(METHOD)
SHA <- prov$sha
HEAD_SHA <- prov$head


KX <- as.integer(Sys.getenv("PO_KX", "100"))             # antichain-size reps
CONFIG <- sprintf("B=%d;N=%d;K=%d;KX=%d;alpha=0.05;posetfloor=99", B, NR, K, KX)
grid <- rbind(
  expand.grid(truth = "PO", C = 3L, margin = c(0.05, 0.10, 0.15),
              nI = c(6L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "PO_INV", C = 3L, margin = NA,     # accept-best (recorded)
              nI = c(6L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "PO_ITEMS", C = 4L, margin = NA,   # accept-best (recorded)
              nI = c(8L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "PO_ITEMS", C = 3L, margin = NA,
              nI = 8L, rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "PO_POLY", C = 3L, margin = 0.05,  # 4-category PO free
              nI = 8L, rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = c("UN", "MON", "IIO"), C = 3L, margin = NA,
              nI = c(6L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "XANTI", C = 3L, margin = NA,      # reviewer counterexample
              nI = c(8L, 24L), rep = seq_len(KX), stringsAsFactors = FALSE),
  expand.grid(truth = "NEARANTI", C = 3L, margin = NA,   # crossings just above
              nI = 8L, rep = seq_len(KX %/% 2L),         # the tolerance
              stringsAsFactors = FALSE))

gen <- function(tr, C, m, nI, rep) {
  sd <- 51000L + 977L * rep + 13L * nI + round(100 * ifelse(is.na(m), 0, m)) +
        7L * C + match(tr, c("PO","PO_INV","PO_ITEMS","UN","MON","IIO",
                             "XANTI","NEARANTI","PO_POLY"))
  if (tr %in% c("XANTI", "NEARANTI"))
    return(qf_gen_anti(tr, nI, NR, rep, C))
  if (tr == "PO_POLY")
    return(simulate_responses("PO", n_persons = NR, n_items = nI,
             n_classes = C, n_cat = 4, poset = "V", po_margin = m, seed = sd))
  switch(tr,
    PO       = simulate_responses("PO", n_persons = NR, n_items = nI,
                 n_classes = C, poset = "V", po_margin = m, seed = sd),
    PO_INV   = simulate_responses("PO", n_persons = NR, n_items = nI,
                 n_classes = C, poset = "V", item_order = "invariant",
                 po_margin = NULL, seed = sd),
    PO_ITEMS = simulate_responses("PO_ITEMS", n_persons = NR, n_items = nI,
                 n_classes = C, poset = "layers2", po_margin = NULL,
                 seed = sd),
    simulate_responses(tr, n_persons = NR, n_items = nI, n_classes = C,
                       seed = sd))
}
pop_pairs <- function(L, eps = 0.01, tau = 0) {
  # population table on the SAME estimand the test uses: normalised expected
  # scores sum_k plogis(L - tau_k) / (K - 1); dichotomous (tau = 0) reduces
  # to plogis(L)
  m <- length(tau)
  P <- Reduce(`+`, lapply(tau, function(tk) stats::plogis(L - tk))) / max(1L, m)
  C <- nrow(P); J <- ncol(P); cl <- 0L; it <- 0L
  for (a in seq_len(C - 1L)) for (b in (a + 1L):C) {
    ab <- mean(pmax(0, P[b, ] - P[a, ])) <= eps
    ba <- mean(pmax(0, P[a, ] - P[b, ])) <= eps
    if (xor(ab, ba)) cl <- cl + 1L }
  for (j in seq_len(J - 1L)) for (k in (j + 1L):J) {
    jk <- mean(pmax(0, P[, k] - P[, j])) <= eps
    kj <- mean(pmax(0, P[, j] - P[, k])) <= eps
    if (xor(jk, kj)) it <- it + 1L }
  c(class = cl, item = it)
}

run <- function(k) {
  cs <- grid[k, ]
  f <- file.path(out, sprintf("v_%s_C%d_%s_%02d_%02d.csv", cs$truth, cs$C,
                              ifelse(is.na(cs$margin), "na", cs$margin),
                              cs$nI, cs$rep))
  if (file.exists(f)) {
    prev <- tryCatch(read.csv(f, nrows = 1), error = function(e) NULL)
    if (qf_canary_match(prev, SHA, METHOD, CONFIG)) return(invisible())
    unlink(f)                       # stale schema/method/build: regenerate
  }
  d <- gen(cs$truth, cs$C, cs$margin, cs$nI, cs$rep)
  pp <- pop_pairs(attr(d, "params")$L,
                  tau = attr(d, "params")$tau %||% 0)
  ach <- attr(d, "params")$po_achieved
  t0 <- proc.time()[3]
  r <- withCallingHandlers(
         select_model_hybrid(d, n_classes = cs$C, B = B, mc.cores = 1L,
                             seed = 1, verbose = FALSE),
         warning = function(w) {
           if (grepl("poset refinement failed|quantitative edge",
                     conditionMessage(w)))
             stop("selector warning treated as failure: ",
                  conditionMessage(w))
           invokeRestart("muffleWarning")
         })
  if (is.null(r) || is.na(r$selected))
    stop("selector returned no verdict")
  if (r$selected %in% c("UN", "IIO", "MON") && is.null(r$poset))
    stop("poset refinement missing for an ordinal-cell verdict")
  g <- function(x, fld) if (is.null(x) || is.null(x[[fld]])) NA else x[[fld]]
  pstr <- function(x) if (is.null(x) || !nrow(x$pairs)) "" else
    paste(sprintf("%d>%d", x$pairs$dominant, x$pairs$dominated), collapse = ";")
  write.csv(data.frame(method = METHOD, sha = SHA, head = HEAD_SHA,
    config = CONFIG,
    poset_b = if (!is.null(r$poset$class)) g(r$poset$class, "b_eff") else
              g(r$poset$item, "b_eff"),
    truth = cs$truth, C = cs$C, margin = cs$margin,
    nI = cs$nI, rep = cs$rep,
    selected = if (is.null(r)) NA else r$selected,
    class_shape = g(r$poset$class, "shape"),
    class_comp = g(r$poset$class, "comparable"),
    class_type = g(r$poset$class, "type"),
    class_pairs = if (is.null(r)) "" else pstr(r$poset$class),
    class_transitive = g(r$poset$class, "transitive"),
    item_shape = g(r$poset$item, "shape"),
    item_comp = g(r$poset$item, "comparable"),
    item_pairs = if (is.null(r)) "" else pstr(r$poset$item),
    item_transitive = g(r$poset$item, "transitive"),
    b_eff = if (!is.null(r$poset$class)) g(r$poset$class, "b_eff") else
            g(r$poset$item, "b_eff"),
    pop_class_pairs = pp["class"], pop_item_pairs = pp["item"],
    po_achieved = if (is.null(ach)) NA else ach,
    secs = round(proc.time()[3] - t0, 1)), f, row.names = FALSE)
}
cat("PO validation v4:", nrow(grid), "datasets (B =", B, ", poset floor 99)\n")
fails <- parallel::mclapply(seq_len(nrow(grid)), function(k)
  tryCatch({ run(k); NULL }, error = function(e)
    data.frame(row = k, truth = grid$truth[k], nI = grid$nI[k],
               rep = grid$rep[k], error = conditionMessage(e))),
  mc.cores = cores, mc.preschedule = FALSE)
fails <- do.call(rbind, Filter(Negate(is.null), fails))
if (!is.null(fails) && nrow(fails)) {
  write.csv(fails, file.path(out, "FAILURES.csv"), row.names = FALSE)
  cat("PO V4 FAILED:", nrow(fails), "cells errored - see FAILURES.csv\n")
  quit(save = "no", status = 1L)
}
unlink(file.path(out, "FAILURES.csv"))
cat("PO V4 DONE\n")
