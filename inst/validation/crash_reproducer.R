# DETERMINISTIC SEGFAULT REPRODUCER - known-crasher TA datasets (2026-08-12)
#
# TA218 (LCR, nI 24) and TA376 (IIO, nI 24) crash the selector with SIGSEGV
# ('invalid permissions', wild address - heap corruption) at the LCR-vs-DM
# bridge stage of the quantitative edge. Reproduced forked and non-forked,
# use_cpp TRUE and FALSE, multiple seeds, on two machines - but INTERMITTENT,
# not deterministic: TA552, which crashed the same way many times, later
# completed normally (see known_crashers.md). Implicates compiled code shared
# by both engines (constrained-fit optimizers are the prime suspect), with
# manifestation depending on allocation state. Run this script REPEATEDLY;
# a single clean pass does not mean the bug is gone.
#
# Run (expect the R process to DIE, so run in a throwaway process):
#   Rscript inst/validation/crash_reproducer.R /path/to/tid_data 552
args <- commandArgs(TRUE)
DD <- args[1]; id <- as.integer(args[2] %||% "552")
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
suppressMessages(library(QuantFit))
e <- new.env(); load(file.path(DD, sprintf("TA%d.Rdata", id)), envir = e)
d <- get(ls(e)[1], e)$obsData; storage.mode(d) <- "integer"
cat("running TA", id, " (", nrow(d), "x", ncol(d), ") - expect SIGSEGV at the",
    "LCR-vs-DM bridge...\n")
r <- select_model_hybrid(d, n_classes = 3L, B = 49, lr_boot_n_starts = 2L,
                         mc.cores = 1L, seed = 1, verbose = TRUE)
cat("NO CRASH - selected:", r$selected, "(bug may be fixed; update",
    "known_crashers.md)\n")
