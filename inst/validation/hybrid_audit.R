# Selection-accuracy audit for the HYBRID selector
# (select_model_manifest(dm_quant="lr"): manifest 2x2 ordinal layer + LR-edge
# DM->quant), on the standard large audit grid used for select_model_ll:
#   6 models x K reps, N=1500, J=8 dichotomous, n_classes=3, B=49
# (J=8 here vs J=12 in the head-to-head -> a generalisation check.) One CSV per
# (model,rep) => RESUMABLE across the machine's background-job kills. Config
# matches the validated head-to-head hybrid (B=49, lr_boot_n_starts=2).
suppressMessages(library(QuantFit))
K       <- as.integer(Sys.getenv("AUDIT_K", "30"))
B       <- 49L
cores   <- as.integer(Sys.getenv("AUDIT_CORES", "4"))
n_items <- 8L
models  <- c("UN", "MON", "IIO", "DM", "LCR", "RM")
scale_of <- c(UN="nominal", MON="ordinal", IIO="ordinal", DM="ordinal",
              LCR="quant", RM="quant")
out <- Sys.getenv("HAUDIT_OUT", "hybrid_audit_out"); dir.create(out, showWarnings=FALSE)
grid <- expand.grid(model=models, rep=seq_len(K), stringsAsFactors=FALSE)

run <- function(k) {
  cs <- grid[k,]; f <- file.path(out, sprintf("a_%s_%03d.csv", cs$model, cs$rep))
  if (file.exists(f)) return(invisible())
  d <- simulate_responses(cs$model, n_persons=1500, n_items=n_items,
                          n_classes=3, seed=7000*cs$rep + match(cs$model, models))
  d <- if (is.list(d)) d$data else d; storage.mode(d) <- "integer"
  sel <- tryCatch(select_model_manifest(d, n_classes=3L, B=B, dm_quant="lr",
             lr_boot_n_starts=2L, mc.cores=1L, seed=1, verbose=FALSE),
           error=function(e) NULL)
  sm <- if (is.null(sel)) NA_character_ else sel$selected
  write.csv(data.frame(truth=cs$model, truth_scale=scale_of[cs$model], rep=cs$rep,
    selected=sm, selected_scale=if(is.na(sm)) NA else scale_of[sm]),
    f, row.names=FALSE)
}
cat("Hybrid audit:", nrow(grid), "datasets (6 models x K=", K, ", J=8, N=1500)\n")
invisible(parallel::mclapply(seq_len(nrow(grid)),
  function(k) tryCatch(run(k), error=function(e) NULL), mc.cores=cores))
cat("HYBRID AUDIT DONE\n")
