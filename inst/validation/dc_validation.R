# Double-cancellation (DM->quant) validation at scale, exactly as the manifest
# selector invokes it: cc_bootstrap_null(check="double", defaults).
#   POWER: DM  (simulate_responses sorted random logits = doubly monotone,
#              non-additive) -> want reject=TRUE (non-additive -> stays DM)
#   SIZE : LCR, RM (additive Rasch) -> want reject=FALSE (~alpha) -> supports_quant
suppressMessages(library(QuantFit))
cc_B <- 49L; cc_n_mat <- 500L; alpha <- 0.05

gen <- function(m, J, N, seed) { set.seed(seed)
  if (m == "RM") { b <- sort(runif(J,-2,2)); th <- rnorm(N)
    d <- matrix(rbinom(N*J,1,plogis(outer(th,b,"-"))),N,J) }
  else { d <- simulate_responses(m, n_persons=N, n_items=J, n_classes=3L, seed=seed)
    d <- if (is.list(d)) d$data else d }
  storage.mode(d) <- "integer"; d }

cases <- expand.grid(model=c("DM","LCR","RM"), J=c(6L,12L), N=c(1500L,3000L),
                     rep=1:8, stringsAsFactors=FALSE)
cases$id <- seq_len(nrow(cases)); set.seed(3131)
cases$seed <- sample.int(.Machine$integer.max, nrow(cases))
out <- Sys.getenv("DC_OUT","dc_out"); dir.create(out, showWarnings=FALSE)
run <- function(k){ cs <- cases[k,]; f <- file.path(out, sprintf("d%04d.csv",k))
  if (file.exists(f)) return(invisible())
  d <- gen(cs$model, cs$J, cs$N, cs$seed)
  cc <- tryCatch(cc_bootstrap_null(d, check="double", B=cc_B, n.mat=cc_n_mat,
           alpha=alpha, mc.cores=1L, seed=cs$id, verbose=FALSE),
         error=function(e) NULL)
  if (is.null(cc)) return(invisible())
  write.csv(data.frame(model=cs$model, J=cs$J, N=cs$N, rep=cs$rep,
    reject=cc$reject, p=round(cc$p_value,4), obs=round(cc$observed,4)),
    f, row.names=FALSE) }
cat("Double-cancellation validation:", nrow(cases), "datasets\n")
invisible(parallel::mclapply(seq_len(nrow(cases)),
  function(k) tryCatch(run(k), error=function(e) NULL), mc.cores=4))
cat("DC VALIDATION DONE\n")
