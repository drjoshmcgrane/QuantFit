# mon_eps sensitivity sweep for select_model_hybrid's ordinal layer.
#
# mon_eps enters routing ONLY as a threshold on the MON axis's resampling q05
# (mon$lo <= eps_N => MON holds), so the 2x2 verdict under ANY eps is derivable
# from one run per dataset: record (iio$p, mon$lo), then evaluate the whole eps
# grid analytically. The quantitative edge is eps-independent given a DM
# verdict, so ordinal-layer routing is the object of study; `selected` at the
# shipped default (0.01) is recorded as a check.
#
# Truths: all six. Correct ordinal verdict: UN->UN, MON->MON, IIO->IIO,
# DM->DM, LCR/RM->DM (quant truths must route to DM to reach the edge).
suppressMessages(library(QuantFit))
K     <- as.integer(Sys.getenv("MES_K", "10"))
B     <- as.integer(Sys.getenv("MES_B", "49"))
cores <- as.integer(Sys.getenv("MES_CORES", "8"))
out   <- Sys.getenv("MES_OUT", "mon_eps_sweep_out"); dir.create(out, showWarnings = FALSE)
models <- c("UN", "MON", "IIO", "DM", "LCR", "RM")
grid <- expand.grid(model = models, rep = seq_len(K), stringsAsFactors = FALSE)

run <- function(k) {
  cs <- grid[k, ]
  f <- file.path(out, sprintf("e_%s_%02d.csv", cs$model, cs$rep))
  if (file.exists(f)) return(invisible())
  d <- simulate_responses(cs$model, n_persons = 1500, n_items = 8,
                          n_classes = 3, seed = 46000L + 977L * cs$rep +
                            match(cs$model, models))
  d <- if (is.list(d)) d$data else d; storage.mode(d) <- "integer"
  r <- tryCatch(select_model_hybrid(d, n_classes = 3L, B = B, mc.cores = 1L,
         seed = 1, verbose = FALSE), error = function(e) NULL)
  if (is.null(r)) return(invisible())
  write.csv(data.frame(truth = cs$model, rep = cs$rep, selected = r$selected,
    iio_p = r$iio$p, iio_stat = r$iio$stat,
    mon_lo = r$mon$lo, mon_stat = r$mon$stat), f, row.names = FALSE)
}
cat("mon_eps sweep:", nrow(grid), "datasets (one run each, B =", B, ")\n")
invisible(parallel::mclapply(seq_len(nrow(grid)),
  function(k) tryCatch(run(k), error = function(e) NULL),
  mc.cores = cores, mc.preschedule = FALSE))

# --- post-hoc eps grid -------------------------------------------------------
d <- do.call(rbind, lapply(list.files(out, "^e_", full.names = TRUE), read.csv))
eps_grid <- c(0.0025, 0.005, 0.01, 0.02, 0.04, 0.08)
target <- c(UN = "UN", MON = "MON", IIO = "IIO", DM = "DM", LCR = "DM", RM = "DM")
cat(sprintf("\n%d/%d datasets | ordinal-layer routing accuracy by mon_eps\n",
            nrow(d), nrow(grid)))
cat("(N = 1500 so eps_N = eps; correct verdict: quant truths -> DM)\n\n")
hdr <- sprintf("%-6s", "eps"); for (m in models) hdr <- paste0(hdr, sprintf("%6s", m))
cat(hdr, "  overall\n")
for (eps in eps_grid) {
  ih <- d$iio_p > 0.05; mh <- d$mon_lo <= eps
  v <- ifelse(ih & mh, "DM", ifelse(ih, "IIO", ifelse(mh, "MON", "UN")))
  row <- sprintf("%-6.4g", eps); tot <- 0L
  for (m in models) {
    s <- d$truth == m; okn <- sum(v[s] == target[m]); tot <- tot + okn
    row <- paste0(row, sprintf("%4d/%d", okn, sum(s)))
  }
  cat(row, sprintf("  %d/%d\n", tot, nrow(d)))
}
cat("\nq05 (mon_lo) distributions by truth:\n")
for (m in models) {
  s <- d[d$truth == m, ]
  cat(sprintf("%-3s: min %.4f  median %.4f  max %.4f\n", m,
      min(s$mon_lo), median(s$mon_lo), max(s$mon_lo)))
}
cat("MES DONE\n")
