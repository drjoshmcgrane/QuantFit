# Snapshot of the TI&D validation as results accumulate.
res_dir <- "/Users/josh/Documents/Claude_Code/QuantModelFitting/tid_results"
fs <- list.files(res_dir, full.names = TRUE, pattern = "\\.csv$")
if (!length(fs)) { cat("no results yet\n"); quit(save = "no") }
r <- do.call(rbind, lapply(fs, read.csv))
models <- c("UN", "MON", "IIO", "DM", "LCR", "RM")
r$genM <- factor(r$genM, models); r$lc_selected <- factor(r$lc_selected, models)
scale_of <- c(UN="nominal",MON="ordinal",IIO="ordinal",DM="ordinal",LCR="quant",RM="quant")
cat(sprintf("=== %d / 1080 datasets done (clean: %d, contaminated: %d) ===\n",
            nrow(r), sum(r$clean), sum(!r$clean)))
cl <- r[r$clean & !is.na(r$lc_selected), ]
if (nrow(cl)) {
  cat("\n-- LC route, CLEAN unidimensional data --\n")
  print(table(true = cl$genM, selected = cl$lc_selected))
  sc <- table(true = scale_of[as.character(cl$genM)],
              sel = scale_of[as.character(cl$lc_selected)])
  cat(sprintf("exact %.1f%% | scale-type %.1f%%\n",
      100*mean(as.character(cl$genM)==as.character(cl$lc_selected)),
      100*mean(scale_of[as.character(cl$genM)]==scale_of[as.character(cl$lc_selected)])))
  print(sc)
  cat("\n-- CC route (reject rate by true model, clean) --\n")
  print(round(tapply(cl$cc_reject, cl$genM, mean, na.rm = TRUE), 2))
  cat("-- Kara route (reject rate by true model, clean) --\n")
  print(round(tapply(cl$kara_reject, cl$genM, mean, na.rm = TRUE), 2))
}
ct <- r[!r$clean & !is.na(r$lc_selected), ]
if (nrow(ct) > 5) {
  cat("\n-- LC route, CONTAMINATED (multidim) data --\n")
  print(table(true = ct$genM, selected = ct$lc_selected))
}
cat(sprintf("\nmedian secs/dataset: LC %.0f  CC %.0f  Kara %.0f | errors: %d\n",
    median(r$secs_lc, na.rm=TRUE), median(r$secs_cc, na.rm=TRUE),
    median(r$secs_kara, na.rm=TRUE), sum(r$err != "", na.rm=TRUE)))
