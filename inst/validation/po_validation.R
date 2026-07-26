# Validation of the "PO" (partial class order) generator at TI&D spec.
#
# simulate_responses("PO") extends TI&D's own generation idiom (one C x J
# U(-4,4) logit matrix; column-sort = MON chain, unsorted = UN antichain) with
# per-item random linear extensions of a specified partial order - the natural
# interpolation, nesting both endpoints exactly (chain-D reproduces the MON
# column sort identically; verified). Ground truth here is a class structure
# the six-model set cannot express.
#
# Questions:
#   1. What does the hybrid SELECT on true-PO data? (Expected UN - or IIO when
#      the unconstrained item side happens to pass the IIO axis, in which case
#      the poset rung never runs; track that leak.)
#   2. Rung SENSITIVITY: how often does the poset rung flag "partial" on
#      true-PO data reaching the UN cell? (The TI&D-real-data estimate was 2/5
#      on misspecification-induced structure; this is purpose-built structure
#      with guaranteed crossing margins.)
#   3. Rung SPECIFICITY at the same spec: UN controls should stay "antichain".
#   4. (Small lattice arm at nI=12) where does the lattice file PO data?
#
# TI&D spec: N = 5000, C = 3, equal class probs. Grid: poset type
# {V, Lambda, single} x nI {6,12,24} x K reps, + UN controls. Resumable.
suppressMessages(library(QuantFit))
K     <- as.integer(Sys.getenv("PO_K", "10"))
NR    <- as.integer(Sys.getenv("PO_N", "5000"))
cores <- as.integer(Sys.getenv("PO_CORES", "8"))
ARM   <- Sys.getenv("PO_ARM", "hybrid")          # "hybrid" or "lattice"
out   <- Sys.getenv("PO_OUT", "po_validation_out"); dir.create(out, showWarnings = FALSE)

grid <- rbind(
  expand.grid(truth = c("V", "Lambda", "single"), nI = c(6L, 12L, 24L),
              rep = seq_len(K), stringsAsFactors = FALSE),
  expand.grid(truth = "UN", nI = c(6L, 12L, 24L), rep = seq_len(K),
              stringsAsFactors = FALSE))
if (ARM == "lattice")                             # small arm: nI=12 only
  grid <- grid[grid$nI == 12L & grid$rep <= min(K, 6L), ]

gen <- function(truth, nI, rep) {
  sd <- 31000L + 611L * rep + 7L * nI + match(truth, c("V","Lambda","single","UN"))
  if (truth == "UN")
    simulate_responses("UN", n_persons = NR, n_items = nI, n_classes = 3,
                       seed = sd)
  else
    simulate_responses("PO", n_persons = NR, n_items = nI, n_classes = 3,
                       poset = truth, seed = sd)
}

run <- function(k) {
  cs <- grid[k, ]
  f <- file.path(out, sprintf("%s_%s_%02d_%02d.csv", ARM, cs$truth, cs$nI, cs$rep))
  if (file.exists(f)) return(invisible())
  d <- gen(cs$truth, cs$nI, cs$rep); storage.mode(d) <- "integer"
  t0 <- proc.time()[3]
  if (ARM == "hybrid") {
    r <- tryCatch(select_model_hybrid(d, n_classes = 3L, B = 49L,
             mc.cores = 1L, seed = 1, verbose = FALSE), error = function(e) NULL)
    write.csv(data.frame(arm = ARM, truth = cs$truth, nI = cs$nI, rep = cs$rep,
      selected = if (is.null(r)) NA else r$selected,
      shape    = if (is.null(r) || is.null(r$poset)) NA else r$poset$shape,
      comparable = if (is.null(r) || is.null(r$poset)) NA else r$poset$comparable,
      lo       = if (is.null(r) || is.null(r$poset)) NA else r$poset$lo,
      secs = round(proc.time()[3] - t0, 1)), f, row.names = FALSE)
  } else {
    r <- tryCatch(select_model_ll(d, n_classes = 2:6, B = 49L, boot_n_starts = 2L,
             mc.cores = 1L, seed = 1, verbose = FALSE), error = function(e) NULL)
    write.csv(data.frame(arm = ARM, truth = cs$truth, nI = cs$nI, rep = cs$rep,
      selected = if (is.null(r)) NA else r$selected,
      shape = NA, comparable = NA, lo = NA,
      secs = round(proc.time()[3] - t0, 1)), f, row.names = FALSE)
  }
}
cat("PO validation [", ARM, "]:", nrow(grid), "datasets (N =", NR, ")\n")
invisible(parallel::mclapply(seq_len(nrow(grid)),
  function(k) tryCatch(run(k), error = function(e) NULL),
  mc.cores = cores, mc.preschedule = FALSE))
cat("PO VALIDATION DONE\n")
