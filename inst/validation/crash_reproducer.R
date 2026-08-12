# DETERMINISTIC SEGFAULT REPRODUCER - known-crasher TA datasets (2026-08-12)
#
# TA376 and TA552 (both IIO-truth, nI = 12/24, TI&D archive) crash the
# selector with SIGSEGV ('invalid permissions', wild address - heap
# corruption) at the LCR-vs-DM bridge stage of the quantitative edge.
# Reproduced: forked and non-forked, use_cpp TRUE and FALSE, seeds 1 and 2 -
# deterministic in the data, engine-independent, which implicates compiled
# code shared by both paths (constrained-fit optimizers are the prime
# suspect). Signature matches the long-standing fork-isolated Kara GC
# segfault (see kara-gc-segfault memory / finding).
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
