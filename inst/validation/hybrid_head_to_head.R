# Hybrid selector pass for the three-way head-to-head: manifest 2x2 ordinal
# layer + LR-edge DM->quant (select_model_manifest(dm_quant="lr")), on the SAME
# shared-seed datasets as the manifest and lattice passes. Resumable (one CSV per
# dataset). Join by id with the manifest_/lattice_ results.
suppressMessages(library(QuantFit))
gen <- function(m, J, N, seed) { set.seed(seed)
  if (m == "RM") { b <- sort(runif(J,-2,2)); th <- rnorm(N)
    d <- matrix(rbinom(N*J,1,plogis(outer(th,b,"-"))),N,J) }
  else { d <- simulate_responses(m, n_persons=N, n_items=J, n_classes=3L, seed=seed)
    d <- if (is.list(d)) d$data else d }
  storage.mode(d) <- "integer"; d }
scale_of <- function(m) c(UN="nominal", MON="ordinal", IIO="ordinal",
  DM="ordinal", LCR="quant", RM="quant")[m]
# identical grid/seeds to head_to_head.R
cases <- expand.grid(model=c("UN","MON","IIO","DM","LCR","RM"), J=12L,
                     N=c(1500L,3000L), rep=1:5, stringsAsFactors=FALSE)
cases$id <- seq_len(nrow(cases)); set.seed(80808)
cases$seed <- sample.int(.Machine$integer.max, nrow(cases))
out <- Sys.getenv("H2H_OUT","h2h_out"); dir.create(out, showWarnings=FALSE)
NCORE <- as.integer(Sys.getenv("H2H_CORES", "6"))
run <- function(k) { cs <- cases[k,]; f <- file.path(out, sprintf("hybrid_h%03d.csv", k))
  if (file.exists(f)) return(invisible())
  d <- gen(cs$model, cs$J, cs$N, cs$seed)
  sel <- tryCatch(select_model_manifest(d, n_classes=3L, B=49L, dm_quant="lr",
             mc.cores=1L, seed=cs$id, verbose=FALSE), error=function(e) NULL)
  sm <- if (is.null(sel)) NA_character_ else sel$selected
  write.csv(data.frame(id=cs$id, truth=cs$model, truth_scale=scale_of(cs$model),
    J=cs$J, N=cs$N, rep=cs$rep, selector="hybrid",
    selected=sm, selected_scale=if(is.na(sm)) NA else scale_of(sm)),
    f, row.names=FALSE) }
cat("Hybrid pass:", nrow(cases), "datasets\n")
invisible(parallel::mclapply(seq_len(nrow(cases)),
  function(k) tryCatch(run(k), error=function(e) NULL), mc.cores=NCORE))
cat("H2H DONE\n")
