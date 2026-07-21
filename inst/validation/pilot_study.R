# Small fresh PILOT (pre-locked-study shakedown): 6 models x J=12 x 3 reps.
suppressMessages(library(QuantFit))
local({ set.seed(41); th <- rnorm(400); b <- runif(8,-2,2)+0.7
  d <- matrix(rbinom(400*8,1,plogis(outer(th,b,"-"))),400,8); storage.mode(d)<-"integer"
  dl <- QuantFit:::.cml_fit_general(d)
  stopifnot(cor(vapply(dl,sum,1), b) > 0.95, max(abs(unlist(dl))) < 20) })
out <- "pilot_results"; dir.create(out, showWarnings = FALSE)
grid <- expand.grid(model = c("UN","MON","IIO","DM","LCR","RM"), rep = 1:3,
                    stringsAsFactors = FALSE)
grid$id <- seq_len(nrow(grid))
set.seed(31415); grid$seed <- sample.int(.Machine$integer.max, nrow(grid))
run1 <- function(i) {
  g <- grid[grid$id == i, ]; f <- file.path(out, sprintf("P%02d.csv", i))
  if (file.exists(f)) return(invisible())
  set.seed(g$seed)
  d <- simulate_responses(g$model, n_persons = 1500, n_items = 12, n_classes = 3)
  d <- if (is.list(d)) d$data else d; storage.mode(d) <- "integer"
  ms <- g$id * 7919L
  row <- data.frame(id=i, model=g$model, rep=g$rep, lc=NA,
                    gate_p=NA, gate_override=NA, rmlcr_p=NA, prof_C=NA,
                    cc_rej=NA, cc_p=NA, cc_B_eff=NA,
                    omni_rej=NA, omni_p=NA, omni_B_eff=NA,
                    ess_obs=NA, ess_null_med=NA, ess_null_low=NA,
                    secs=NA, err="")
  t0 <- proc.time()[3]
  lco <- tryCatch(suppressWarnings(select_model_ll(d, n_classes=1:5, B=99,
        n_starts=5, boot_n_starts=3, seed=ms, mc.cores=1)),
        error=function(e) NULL)
  if (!is.null(lco)) {
    row$lc <- lco$selected
    if (!is.null(lco$quant_gate)) {
      row$gate_p <- round(lco$quant_gate$p_value, 4)
      row$gate_override <- isTRUE(lco$quant_gate$severity_override)
    }
    if (!is.null(lco$rm_vs_lcr)) {
      row$rmlcr_p <- round(lco$rm_vs_lcr$p_value, 4)
      row$prof_C <- lco$rm_vs_lcr$profiled_C
    }
  } else row$err <- "lc;"
  cc <- tryCatch(suppressWarnings(cc_bootstrap_hierarchy(d, B=99, n.mat=500,
        seed=ms, mc.cores=1, verbose=FALSE)), error=function(e) NULL)
  if (!is.null(cc)) { row$cc_rej <- !cc$supports_quant
    row$cc_p <- round(min(cc$p_adjusted),4)
    row$cc_B_eff <- min(vapply(cc$levels, function(l) l$B, 1))
  } else row$err <- paste0(row$err,"cc;")
  kj <- parallel::mcparallel(tryCatch(suppressWarnings(omni_bootstrap_null(d,
        B=99, S=3000, N_synth=25, seed=ms, mc.cores=1, verbose=FALSE)),
        error=function(e) NULL), silent=TRUE)
  ka <- tryCatch(parallel::mccollect(kj, wait=TRUE)[[1]], error=function(e) NULL)
  if (inherits(ka,"try-error")) ka <- NULL
  if (!is.null(ka) && is.list(ka)) { row$omni_rej <- ka$reject
    row$omni_p <- round(ka$p_value,4); row$omni_B_eff <- ka$B
    row$ess_obs <- round(ka$ess_min,1)
    row$ess_null_med <- round(ka$null_ess_median,1)
    row$ess_null_low <- ka$null_ess_low_n
  } else row$err <- paste0(row$err,"omni;")
  row$secs <- round(proc.time()[3]-t0,1)
  write.csv(row, f, row.names=FALSE)
}
invisible(parallel::mclapply(grid$id, function(i)
  tryCatch(run1(i), error=function(e) NULL), mc.cores = 5L))
cat("PILOT", sum(file.exists(file.path(out, sprintf("P%02d.csv", grid$id)))),
    "/", nrow(grid), "DONE\n")
