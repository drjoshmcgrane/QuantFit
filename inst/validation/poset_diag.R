# PARTIAL-ORDER DIAGNOSTIC
#
# The 2x2 asks whether a TOTAL order on classes exists (MON). Its negation lumps
# together an ANTICHAIN (no class dominates another -> genuinely nominal) with a
# PARTIAL ORDER (some classes dominate others, just not a chain -> strictly more
# structure than nominal). This measures which of the two the "UN" verdicts
# actually are.
#
# Dominance is defined consistently with the MON statistic: class a dominates
# class b iff the per-item mean downward violation is negligible,
#     mean_j max(0, P[j,b] - P[j,a]) <= eps
# so it reuses the already-calibrated mon_eps rather than inventing a threshold.
#
# Poset shape for C classes (C(C-1)/2 pairs):
#   all pairs comparable & acyclic -> CHAIN     (a total order exists = MON)
#   0 pairs comparable             -> ANTICHAIN (genuinely nominal)
#   otherwise                      -> PARTIAL   (currently misfiled as nominal)
suppressMessages(library(QuantFit))
env <- asNamespace("QuantFit"); rmt <- get("refit_model_type", env)

poset <- function(P, eps) {                 # P = items x classes
  C <- ncol(P); D <- matrix(FALSE, C, C)
  for (a in seq_len(C)) for (b in seq_len(C)) if (a != b)
    D[a, b] <- mean(pmax(0, P[, b] - P[, a])) <= eps   # a dominates b
  pairs <- 0L; cyc <- 0L
  for (a in 1:(C-1)) for (b in (a+1):C) {
    ab <- D[a,b]; ba <- D[b,a]
    if (xor(ab, ba)) pairs <- pairs + 1L
    else if (ab && ba) cyc <- cyc + 1L      # mutual "dominance" = tied/flat
  }
  tot <- C*(C-1)/2
  shape <- if (pairs == tot) "chain" else if (pairs == 0L) "antichain" else "partial"
  list(shape = shape, comparable = pairs, total = tot, ties = cyc)
}

EPS <- as.numeric(strsplit(Sys.getenv("PO_EPS","0.01,0.03,0.05"), ",")[[1]])
od <- Sys.getenv("PO_OUT","po_out"); dir.create(od, showWarnings = FALSE)

# ---- source A: simulated data from each generating model -------------------
simA <- expand.grid(m = c("UN","MON","IIO","DM"), rep = 1:8, stringsAsFactors = FALSE)
# ---- source B: REAL TI&D datasets the lattice classified UN ----------------
pri <- do.call(rbind, lapply(list.files("tid_results", full.names = TRUE), read.csv))
pri <- pri[!is.na(pri$lc_selected) & pri$lc_selected == "UN" & pri$nI %in% c(6,12,24), ]

run_sim <- function(k) {
  cs <- simA[k, ]; f <- file.path(od, sprintf("s_%s_%02d.csv", cs$m, cs$rep))
  if (file.exists(f)) return(invisible())
  d <- simulate_responses(cs$m, n_persons = 1500, n_items = 12, n_classes = 3,
                          seed = 900 * cs$rep + match(cs$m, c("UN","MON","IIO","DM")))
  d <- if (is.list(d)) d$data else d; storage.mode(d) <- "integer"
  P <- rmt("UN", d, 3L, 5L, TRUE)$item_probs
  rows <- do.call(rbind, lapply(EPS, function(e) { p <- poset(P, e)
    data.frame(src = "sim", truth = cs$m, id = cs$rep, eps = e,
               shape = p$shape, comparable = p$comparable, total = p$total) }))
  write.csv(rows, f, row.names = FALSE)
}
run_tid <- function(k) {
  cs <- pri[k, ]; f <- file.path(od, sprintf("t_%d.csv", cs$id))
  if (file.exists(f)) return(invisible())
  e <- new.env(); load(sprintf("tid_data/TA%d.Rdata", cs$id), envir = e)
  d <- get(ls(e)[1], e)$obsData; storage.mode(d) <- "integer"
  P <- rmt("UN", d, 3L, 5L, TRUE)$item_probs
  rows <- do.call(rbind, lapply(EPS, function(ee) { p <- poset(P, ee)
    data.frame(src = "tid", truth = cs$genM, id = cs$id, eps = ee,
               shape = p$shape, comparable = p$comparable, total = p$total) }))
  write.csv(rows, f, row.names = FALSE)
}
nc <- as.integer(Sys.getenv("PO_CORES","3"))
invisible(parallel::mclapply(seq_len(nrow(simA)), function(k) tryCatch(run_sim(k), error=function(e) NULL), mc.cores = nc))
cat("sim done; TI&D UN-classified datasets:", nrow(pri), "\n")
invisible(parallel::mclapply(seq_len(nrow(pri)), function(k) tryCatch(run_tid(k), error=function(e) NULL), mc.cores = nc))
cat("PO DONE\n")
