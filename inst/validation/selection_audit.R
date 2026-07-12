# Selection-accuracy audit for select_model_ll().
#
# Simulates data from each of the six latent-structure models with
# simulate_responses(), runs the calibrated likelihood-ratio-bootstrap selector,
# and tabulates (a) the exact-model confusion matrix and (b) the coarser
# scale-type accuracy (nominal = UN, ordinal = MON/IIO/DM, quantitative =
# LCR/RM). Works for dichotomous (n_cat = 2) and polytomous (n_cat > 2) data.
#
# Usage:  Rscript selection_audit.R <n_cat> <reps> <B> [mc.cores]

suppressMessages(library(QuantFit))

a <- commandArgs(trailingOnly = TRUE)
n_cat  <- if (length(a) >= 1) as.integer(a[1]) else 2L
K      <- if (length(a) >= 2) as.integer(a[2]) else 30L
B      <- if (length(a) >= 3) as.integer(a[3]) else 99L
cores  <- if (length(a) >= 4) as.integer(a[4]) else max(1L, parallel::detectCores() - 2L)

n_items <- if (n_cat == 2L) 8L else 6L
models  <- c("UN", "MON", "IIO", "DM", "LCR", "RM")
scale_of <- c(UN = "nominal", MON = "ordinal", IIO = "ordinal", DM = "ordinal",
              LCR = "quant", RM = "quant")

select_one <- function(true, rep) {
  d <- simulate_responses(true, n_persons = 1500, n_items = n_items,
                          n_classes = 3, n_cat = n_cat,
                          seed = 7000 * rep + match(true, models))
  sel <- tryCatch(
    select_model_ll(d, n_classes = 3, B = B, n_starts = 3, boot_n_starts = 2,
                    seed = 1, mc.cores = 1),
    error = function(e) NULL)
  if (is.null(sel)) NA_character_ else sel$selected
}

conf <- matrix(0L, 6, 6, dimnames = list(true = models, selected = models))
for (true in models) {
  sels <- unlist(parallel::mclapply(seq_len(K), function(r) select_one(true, r),
                                    mc.cores = cores))
  for (s in sels[!is.na(sels)]) conf[true, s] <- conf[true, s] + 1L
}

cat(sprintf("Selection audit: %s data, K = %d reps, B = %d\n\n",
            if (n_cat == 2) "dichotomous" else paste0(n_cat, "-category"), K, B))
print(conf)

sc_true <- scale_of[rep(models, times = rowSums(conf))]
sc_sel  <- scale_of[unlist(lapply(models, function(t) rep(colnames(conf), conf[t, ])))]
cat(sprintf("\nExact-model recovery: %.1f%%   Scale-type accuracy: %.1f%%\n",
            100 * sum(diag(conf)) / sum(conf), 100 * mean(sc_true == sc_sel)))
cat("\nScale-type confusion:\n")
print(table(true = sc_true, selected = sc_sel))
