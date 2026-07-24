# Six-model head-to-head: MANIFEST 2x2 selector (select_model_manifest) vs
# LR-edge LATTICE selector (select_model_ll), on the SAME datasets (shared
# seeds), across the full generating hierarchy. Production settings (manifest DC
# n.mat=500; lattice B=99).
#
# The two selectors are run as SEPARATE passes (env SELECTOR=lattice|manifest),
# because the manifest DC step (~12 min/dataset) is ~50x slower than the lattice
# selector; running them together makes each expensive dataset ~30 min and the
# job never survives machine-sleep. Split => the lattice table finishes in ~20
# min while the manifest pass grinds resumably. One CSV per (selector,dataset)
# => fully RESUMABLE: re-run to continue after a kill; done datasets are skipped.
# Join the two passes by dataset id for the comparison.
#
#   generating models : UN, MON, IIO, DM (ordinal via sorted random logits,
#                       DM doubly-monotone & non-additive), LCR, RM (quant/additive)
#   truth scale-type  : UN->nominal ; MON/IIO/DM->ordinal ; LCR/RM->quant
suppressMessages(library(QuantFit))
SELECTOR <- Sys.getenv("SELECTOR", "lattice")   # "lattice" or "manifest"

gen <- function(m, J, N, seed) { set.seed(seed)
  if (m == "RM") { b <- sort(runif(J,-2,2)); th <- rnorm(N)
    d <- matrix(rbinom(N*J,1,plogis(outer(th,b,"-"))),N,J) }
  else { d <- simulate_responses(m, n_persons=N, n_items=J, n_classes=3L, seed=seed)
    d <- if (is.list(d)) d$data else d }
  storage.mode(d) <- "integer"; d }

scale_of <- function(m) c(UN="nominal", MON="ordinal", IIO="ordinal",
  DM="ordinal", LCR="quant", RM="quant")[m]

cases <- expand.grid(model=c("UN","MON","IIO","DM","LCR","RM"), J=12L,
                     N=c(1500L,3000L), rep=1:5, stringsAsFactors=FALSE)
cases$id <- seq_len(nrow(cases)); set.seed(80808)
cases$seed <- sample.int(.Machine$integer.max, nrow(cases))
out <- Sys.getenv("H2H_OUT","h2h_out"); dir.create(out, showWarnings=FALSE)

run <- function(k) {
  cs <- cases[k,]; f <- file.path(out, sprintf("%s_h%03d.csv", SELECTOR, k))
  if (file.exists(f)) return(invisible())
  d <- gen(cs$model, cs$J, cs$N, cs$seed)
  sel <- if (SELECTOR == "manifest")
      tryCatch(select_model_manifest(d, n_classes=3L, B=49L, cc_B=49L,
                 cc_n_mat=500L, mc.cores=1L, seed=cs$id, verbose=FALSE),
             error=function(e) NULL)
    else
      tryCatch(select_model_ll(d, n_classes=3L, B=99L, mc.cores=1L,
                 seed=cs$id, verbose=FALSE), error=function(e) NULL)
  sm <- if (is.null(sel)) NA_character_ else sel$selected
  write.csv(data.frame(id=cs$id, truth=cs$model, truth_scale=scale_of(cs$model),
    J=cs$J, N=cs$N, rep=cs$rep, selector=SELECTOR,
    selected=sm, selected_scale=if(is.na(sm)) NA else scale_of(sm)),
    f, row.names=FALSE) }

NCORE <- as.integer(Sys.getenv("H2H_CORES", "4"))
cat("Head-to-head [", SELECTOR, "]:", nrow(cases), "datasets (production n.mat=500,",
    NCORE, "cores)\n")
invisible(parallel::mclapply(seq_len(nrow(cases)),
  function(k) tryCatch(run(k), error=function(e) NULL), mc.cores=NCORE))
cat("H2H DONE\n")
