# Validation of the DOMINANCE-DEMONSTRATION partial-order tests (v3).
#
# Supersedes v2 entirely (external review 2026-07-29): the v2 "asymmetry"
# statistic algebraically collapsed to |mean difference| and its permutation
# null tested exchangeability, so v2's size/power tables validated the wrong
# thing (raw files remain in po_validation2_results/ as history). v3 validates
# the redesign: per-directed-pair violation masses from ONE aligned parametric
# bootstrap of the fitted UN table, pairs declared comparable at
# Bonferroni-corrected one-sided levels, "partial" = >= 1 demonstrated pair.
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
K     <- as.integer(Sys.getenv("PO_K", "8"))
NR    <- as.integer(Sys.getenv("PO_N", "1500"))
B     <- as.integer(Sys.getenv("PO_B", "49"))
cores <- as.integer(Sys.getenv("PO_CORES", "8"))
out   <- Sys.getenv("PO_OUT", "po_validation3_out"); dir.create(out, showWarnings = FALSE)

grid <- rbind(
  expand.grid(truth = "PO", C = 3L, margin = c(0.05, 0.10, 0.15),
              nI = c(6L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "PO_INV", C = 3L, margin = 0.12,
              nI = c(6L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "PO_ITEMS", C = 4L, margin = NA,   # feasible per nI below
              nI = c(8L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "PO_ITEMS", C = 3L, margin = NA,
              nI = 8L, rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = c("UN", "MON", "IIO"), C = 3L, margin = NA,
              nI = c(6L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "XANTI", C = 3L, margin = NA,      # reviewer counterexample
              nI = c(8L, 24L), rep = seq_len(K), stringsAsFactors = FALSE))
po_items_margin <- function(C, nI)                       # just under ceiling
  if (C == 4L) c(`8` = 0.030, `12` = 0.017, `24` = 0.007)[as.character(nI)] else
  c(`8` = 0.011)[as.character(nI)]

gen <- function(tr, C, m, nI, rep) {
  sd <- 51000L + 977L * rep + 13L * nI + round(100 * ifelse(is.na(m), 0, m)) +
        7L * C + match(tr, c("PO","PO_INV","PO_ITEMS","UN","MON","IIO","XANTI"))
  if (tr == "XANTI") {
    # crossing antichain with unequal averages: increasing / shifted
    # decreasing / alternating logit profiles - zero population dominance
    set.seed(sd)
    base <- seq(-1.6, 1.6, length.out = nI)
    P <- rbind(plogis(base), plogis(rev(base) + 0.35),
               plogis(base * rep_len(c(-1, 1), nI) + 0.15))
    cls <- sample.int(3L, NR, replace = TRUE)
    d <- matrix(rbinom(NR * nI, 1, P[cls, ]), NR, nI)
    storage.mode(d) <- "integer"
    attr(d, "params") <- list(L = qlogis(P))
    return(d)
  }
  switch(tr,
    PO       = simulate_responses("PO", n_persons = NR, n_items = nI,
                 n_classes = C, poset = "V", po_margin = m, seed = sd),
    PO_INV   = simulate_responses("PO", n_persons = NR, n_items = nI,
                 n_classes = C, poset = "V", item_order = "invariant",
                 po_margin = m, seed = sd),
    PO_ITEMS = simulate_responses("PO_ITEMS", n_persons = NR, n_items = nI,
                 n_classes = C, poset = "layers2",
                 po_margin = po_items_margin(C, nI), seed = sd),
    simulate_responses(tr, n_persons = NR, n_items = nI, n_classes = C,
                       seed = sd))
}
pop_pairs <- function(L, eps = 0.01) {
  P <- stats::plogis(L); C <- nrow(P); J <- ncol(P); cl <- 0L; it <- 0L
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
  if (file.exists(f)) return(invisible())
  d <- gen(cs$truth, cs$C, cs$margin, cs$nI, cs$rep)
  pp <- pop_pairs(attr(d, "params")$L)
  ach <- attr(d, "params")$po_achieved
  t0 <- proc.time()[3]
  r <- tryCatch(suppressWarnings(select_model_hybrid(d, n_classes = cs$C,
         B = B, mc.cores = 1L, seed = 1, verbose = FALSE)),
       error = function(e) NULL)
  g <- function(x, fld) if (is.null(x) || is.null(x[[fld]])) NA else x[[fld]]
  write.csv(data.frame(truth = cs$truth, C = cs$C, margin = cs$margin,
    nI = cs$nI, rep = cs$rep,
    selected = if (is.null(r)) NA else r$selected,
    class_shape = g(r$poset$class, "shape"),
    class_comp = g(r$poset$class, "comparable"),
    class_type = g(r$poset$class, "type"),
    item_shape = g(r$poset$item, "shape"),
    item_comp = g(r$poset$item, "comparable"),
    b_eff = if (!is.null(r$poset$class)) g(r$poset$class, "b_eff") else
            g(r$poset$item, "b_eff"),
    pop_class_pairs = pp["class"], pop_item_pairs = pp["item"],
    po_achieved = if (is.null(ach)) NA else ach,
    secs = round(proc.time()[3] - t0, 1)), f, row.names = FALSE)
}
cat("PO validation v3:", nrow(grid), "datasets (B =", B, ")\n")
invisible(parallel::mclapply(seq_len(nrow(grid)),
  function(k) tryCatch(run(k), error = function(e) NULL),
  mc.cores = cores, mc.preschedule = FALSE))
cat("PO V3 DONE\n")
