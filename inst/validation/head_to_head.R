# LATTICE arm of the six-model head-to-head (select_model_ll). The hybrid arm is
# hybrid_head_to_head.R, which uses the SAME grid and seeds (shared
# seeds), so the two are paired. Historical note: a third "manifest" arm
# (2x2 + double cancellation) existed here; that selector has been removed.
#
# One CSV per dataset => fully RESUMABLE: re-run to continue; done datasets are
# skipped. Join with the hybrid pass by dataset id for the paired comparison.
#
#   generating models : UN, MON, IIO, DM (ordinal via sorted random logits,
#                       DM doubly-monotone & non-additive), LCR, RM (quant/additive)
#   truth scale-type  : UN->nominal ; MON/IIO/DM->ordinal ; LCR/RM->quant
suppressMessages(library(QuantFit))
SELECTOR <- "lattice"

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
  # B=49, boot_n_starts=2 to match the hybrid arm's LR delegation exactly, so the
  # paired comparison is at matched settings. (An earlier note here claimed
  # production B=99 was "infeasible on this machine"; that was based on a
  # since-disproven belief about job kill limits - see tid_hybrid.R.)
  sel <- tryCatch(select_model_ll(d, n_classes=3L, B=49L, boot_n_starts=2L,
                 mc.cores=1L, seed=cs$id, verbose=FALSE), error=function(e) NULL)
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
