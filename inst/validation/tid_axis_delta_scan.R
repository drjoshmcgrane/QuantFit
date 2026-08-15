# AXIS DELTA SCAN: which TI&D datasets can change verdict under the v0.3.4
# IIO-axis null fix? The fix touches ONLY that axis, so a hybrid verdict can
# differ only where the axis decision flips. Evaluating the axis alone (both
# configs, identical data) is ~1-2 orders cheaper than re-running the full
# pipeline, and yields the exact regeneration set.
suppressMessages(library(QuantFit))
DD    <- Sys.getenv("TAX_DATA", "../tid_data")
out   <- Sys.getenv("TAX_OUT", "../qf_evidence/axis_delta_scan")
B     <- as.integer(Sys.getenv("TAX_B", "49"))
cores <- as.integer(Sys.getenv("TAX_CORES", "6"))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
bi <- quantfit_build_info()
lv <- c("UN","MON","IIO","DM","LCR","RM")
gm <- read.csv(file.path(DD, "generatingModels.csv")); gm$genM <- lv[gm$genM + 1L]
gc_ <- read.csv(file.path(DD, "generatingConditions.csv")); names(gc_)[1] <- "id"
sel <- merge(gm[, c("id","genM")], gc_[, c("id","nI")], by = "id")
sel <- sel[sel$nI <= 24 & !(sel$id %in% c(218L, 376L)), ]   # crashers excluded
run <- function(k) {
  cs <- sel[k, ]; f <- file.path(out, sprintf("s_%d.csv", cs$id))
  if (file.exists(f)) return(invisible())
  e <- new.env(); load(file.path(DD, sprintf("TA%d.Rdata", cs$id)), envir = e)
  d <- get(ls(e)[1], e)$obsData; storage.mode(d) <- "integer"
  # OLD = the configuration that produced the committed tables (fixed null
  # C, no adaptive precision). NEW = current shipped behaviour (BIC null C +
  # refinement of borderline p-values), so the flip set is reproducible.
  old <- QuantFit:::.manifest_iio_holds(d, 3L, B, 5L, TRUE, 1L,
           null_C_range = 3L, B_refine = NULL)
  new <- QuantFit:::.manifest_iio_holds(d, 3L, B, 5L, TRUE, 1L,
           null_C_range = 2:6)
  write.csv(data.frame(sha = as.character(bi$sha), id = cs$id,
    truth = cs$genM, nI = cs$nI, stat = old$stat, p_old = old$p,
    holds_old = old$holds, p_new = new$p, holds_new = new$holds,
    null_C = new$null_C, B_eff_new = new$B_eff, refined = new$refined,
    flip = !identical(old$holds, new$holds)),
    f, row.names = FALSE)
}
cat("axis delta scan:", nrow(sel), "TI&D datasets\n")
invisible(parallel::mclapply(seq_len(nrow(sel)),
  function(k) tryCatch(run(k), error = function(e) NULL),
  mc.cores = cores, mc.preschedule = FALSE))
cat("SCAN DONE\n")
