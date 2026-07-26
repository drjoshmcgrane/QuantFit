# Validation of the CALIBRATED partial-order tests (v2; supersedes the v1
# count/threshold study whose raw results remain in po_validation_results/).
#
# Design mirrors how the 2x2 axes were validated: SIZE and POWER for each
# antichain test, swept across the dimensions that govern them. Key structural
# fact (established in the smoke tests, quantified here): each side's power is
# governed by the size of the OTHER side - the class test has J observations
# per class pair, the item test has C per item pair - so class-side power is
# swept over nI and PO crossing margin, and every case reports both tests where
# its cell runs them.
#
# Truths and their routing:
#   PO free (margin m)   -> UN cell: class test (power) + item test (vs pop)
#   PO invariant         -> IIO cell: class test (power)
#   PO_ITEMS layers2     -> MON cell: item test (power)
#   UN control           -> UN cell: class test SIZE + item test (vs pop)
#   MON control          -> MON cell: item test vs POPULATION truth
#   IIO control          -> IIO cell: class test SIZE
#
# Item-side outcomes are scored against the POPULATION table (params$L): TI&D
# style draws genuinely contain margin-driven item dominance, and the test
# detecting it is correct - the population comparable-pair count (eps = 0.01 on
# the true probabilities) is recorded per dataset for exactly that comparison.
suppressMessages(library(QuantFit))
K     <- as.integer(Sys.getenv("PO_K", "8"))
NR    <- as.integer(Sys.getenv("PO_N", "1500"))
B     <- as.integer(Sys.getenv("PO_B", "49"))
cores <- as.integer(Sys.getenv("PO_CORES", "8"))
out   <- Sys.getenv("PO_OUT", "po_validation2_out"); dir.create(out, showWarnings = FALSE)

grid <- rbind(
  expand.grid(truth = "PO",       margin = c(0.05, 0.10, 0.15),
              nI = c(6L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "PO_INV",   margin = 0.12,
              nI = c(6L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "PO_ITEMS", margin = 0.10,
              nI = c(8L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = c("UN", "MON", "IIO"), margin = NA,
              nI = c(6L, 12L, 24L), rep = seq_len(K), stringsAsFactors = FALSE))

gen <- function(tr, m, nI, rep) {
  sd <- 41000L + 977L * rep + 13L * nI + round(100 * ifelse(is.na(m), 0, m)) +
        match(tr, c("PO","PO_INV","PO_ITEMS","UN","MON","IIO"))
  switch(tr,
    PO       = simulate_responses("PO", n_persons = NR, n_items = nI,
                 n_classes = 3, poset = "V", po_margin = m, seed = sd),
    PO_INV   = simulate_responses("PO", n_persons = NR, n_items = nI,
                 n_classes = 3, poset = "V", item_order = "invariant",
                 po_margin = m, seed = sd),
    PO_ITEMS = suppressWarnings(simulate_responses("PO_ITEMS", n_persons = NR,
                 n_items = nI, n_classes = 3, poset = "layers2",
                 po_margin = m, seed = sd)),
    simulate_responses(tr, n_persons = NR, n_items = nI, n_classes = 3,
                       seed = sd))
}
pop_pairs <- function(L, eps = 0.01) {           # population dominance counts
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
  f <- file.path(out, sprintf("v_%s_%s_%02d_%02d.csv", cs$truth,
                              ifelse(is.na(cs$margin), "na", cs$margin),
                              cs$nI, cs$rep))
  if (file.exists(f)) return(invisible())
  d <- gen(cs$truth, cs$margin, cs$nI, cs$rep)
  pp <- pop_pairs(attr(d, "params")$L); storage.mode(d) <- "integer"
  t0 <- proc.time()[3]
  r <- tryCatch(select_model_hybrid(d, n_classes = 3L, B = B, mc.cores = 1L,
         seed = 1, verbose = FALSE), error = function(e) NULL)
  g <- function(x, fld) if (is.null(x)) NA else x[[fld]]
  write.csv(data.frame(truth = cs$truth, margin = cs$margin, nI = cs$nI,
    rep = cs$rep, selected = if (is.null(r)) NA else r$selected,
    class_shape = g(r$poset$class, "shape"), class_p = g(r$poset$class, "p"),
    class_type = g(r$poset$class, "type"),
    item_shape = g(r$poset$item, "shape"), item_p = g(r$poset$item, "p"),
    pop_class_pairs = pp["class"], pop_item_pairs = pp["item"],
    secs = round(proc.time()[3] - t0, 1)), f, row.names = FALSE)
}
cat("PO validation v2:", nrow(grid), "datasets (B =", B, ")\n")
invisible(parallel::mclapply(seq_len(nrow(grid)),
  function(k) tryCatch(run(k), error = function(e) NULL),
  mc.cores = cores, mc.preschedule = FALSE))
cat("PO V2 DONE\n")
