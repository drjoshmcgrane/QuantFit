# IIO-AXIS NULL CLASS-COUNT STUDY (2026-08-14)
#
# Measures size and power of the IIO axis under the OLD (null at the fixed
# routing C) and NEW (BIC-selected null C, floored at routing C) reference
# distributions, across test length and true class count. The axis statistic
# is model-free and identical under both; only the null changes.
#
#   size  : truths where IIO HOLDS (IIO, DM) - rejection = false rejection
#   power : truths where IIO is VIOLATED (UN, MON) - rejection = correct
#
# Motivation: on 4-class DM data at J = 48 the fixed C = 3 null falsely
# rejected (TA173 p = .020 vs .280 at C = 4). The risk of the fix is the
# opposite error: an over-rich null costing power on UN/MON.
suppressMessages(library(QuantFit))
K     <- as.integer(Sys.getenv("AX_K", "20"))
NR    <- as.integer(Sys.getenv("AX_N", "1500"))
B     <- as.integer(Sys.getenv("AX_B", "49"))
cores <- as.integer(Sys.getenv("AX_CORES", "6"))
out   <- Sys.getenv("AX_OUT", "../qf_evidence/axis_null_study")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
bi <- quantfit_build_info()
grid <- expand.grid(truth = c("IIO","DM","MON","UN"), J = c(6L,12L,24L,48L),
                    trueC = c(3L,4L), rep = seq_len(K),
                    stringsAsFactors = FALSE)
run <- function(k) {
  cs <- grid[k, ]
  f <- file.path(out, sprintf("a_%s_J%02d_C%d_%02d.csv", cs$truth, cs$J,
                              cs$trueC, cs$rep))
  if (file.exists(f)) return(invisible())
  d <- simulate_responses(cs$truth, n_persons = NR, n_items = cs$J,
         n_classes = cs$trueC, seed = 91000L + 733L*cs$rep + 17L*cs$J +
                                       3L*cs$trueC + match(cs$truth,
                                       c("IIO","DM","MON","UN")))
  d <- if (is.list(d)) d$data else d; storage.mode(d) <- "integer"
  old <- QuantFit:::.manifest_iio_holds(d, 3L, B, 5L, TRUE, 1L,
                                        null_C_range = 3L)
  new <- QuantFit:::.manifest_iio_holds(d, 3L, B, 5L, TRUE, 1L,
                                        null_C_range = 2:6)
  write.csv(data.frame(sha = as.character(bi$sha), truth = cs$truth,
    J = cs$J, trueC = cs$trueC, rep = cs$rep, stat = old$stat,
    p_old = old$p, holds_old = old$holds,
    p_new = new$p, holds_new = new$holds, null_C_new = new$null_C),
    f, row.names = FALSE)
}
cat("axis-null study:", nrow(grid), "datasets\n")
fails <- parallel::mclapply(seq_len(nrow(grid)), function(k)
  tryCatch({ run(k); NULL }, error = function(e)
    data.frame(row = k, error = conditionMessage(e))),
  mc.cores = cores, mc.preschedule = FALSE)
fails <- do.call(rbind, Filter(Negate(is.null), fails))
if (!is.null(fails) && nrow(fails))
  write.csv(fails, file.path(out, "FAILURES.csv"), row.names = FALSE)
cat("AXIS STUDY DONE\n")
