# MON-AXIS TRUE-CLASS-COUNT STUDY (2026-08-14)
#
# Companion to iio_axis_null_study.R. The MON axis computes its statistic
# from a UN fit at the FIXED routing C and calibrates by data resampling
# (not model simulation), so it cannot fail the way the IIO null did - but
# the statistic rests on that fit, and a C below the data's true class count
# merges classes, which could MASK a monotonicity violation living between
# the merged pair (under-rejection -> MON-violating data called DM/IIO) or
# manufacture one. Never tested with true C > routing C: exactly the blind
# spot that hid the IIO-axis defect.
#
#   size  : truths where MON HOLDS (MON, DM) - rejection = false rejection
#   power : truths where MON is VIOLATED (IIO, UN) - rejection = correct
# Routing C fixed at 3 throughout; true C varied 3:5.
suppressMessages(library(QuantFit))
K     <- as.integer(Sys.getenv("MX_K", "20"))
NR    <- as.integer(Sys.getenv("MX_N", "1500"))
B     <- as.integer(Sys.getenv("MX_B", "25"))
cores <- as.integer(Sys.getenv("MX_CORES", "6"))
out   <- Sys.getenv("MX_OUT", "../qf_evidence/mon_axis_study")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
bi <- quantfit_build_info()
grid <- expand.grid(truth = c("MON","DM","IIO","UN"), J = c(8L,24L,48L),
                    trueC = c(3L,4L,5L), rep = seq_len(K),
                    stringsAsFactors = FALSE)
run <- function(k) {
  cs <- grid[k, ]
  f <- file.path(out, sprintf("m_%s_J%02d_C%d_%02d.csv", cs$truth, cs$J,
                              cs$trueC, cs$rep))
  if (file.exists(f)) return(invisible())
  d <- simulate_responses(cs$truth, n_persons = NR, n_items = cs$J,
         n_classes = cs$trueC, seed = 77000L + 811L*cs$rep + 19L*cs$J +
                                       5L*cs$trueC + match(cs$truth,
                                       c("MON","DM","IIO","UN")))
  d <- if (is.list(d)) d$data else d; storage.mode(d) <- "integer"
  eps <- 0.01 * sqrt(1500 / nrow(d))
  r3 <- QuantFit:::.manifest_mon_holds(d, 3L, B, 5L, TRUE, eps, 1L)
  r4 <- QuantFit:::.manifest_mon_holds(d, cs$trueC, B, 5L, TRUE, eps, 1L)
  write.csv(data.frame(sha = as.character(bi$sha), truth = cs$truth,
    J = cs$J, trueC = cs$trueC, rep = cs$rep,
    stat_C3 = r3$stat, lo_C3 = r3$lo, holds_C3 = r3$holds,
    stat_Ct = r4$stat, lo_Ct = r4$lo, holds_Ct = r4$holds),
    f, row.names = FALSE)
}
cat("MON-axis true-C study:", nrow(grid), "datasets\n")
fails <- parallel::mclapply(seq_len(nrow(grid)), function(k)
  tryCatch({ run(k); NULL }, error = function(e)
    data.frame(row = k, error = conditionMessage(e))),
  mc.cores = cores, mc.preschedule = FALSE)
fails <- do.call(rbind, Filter(Negate(is.null), fails))
if (!is.null(fails) && nrow(fails))
  write.csv(fails, file.path(out, "FAILURES.csv"), row.names = FALSE)
cat("MON STUDY DONE\n")
