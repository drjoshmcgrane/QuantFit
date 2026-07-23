# Cross-parameter robustness / calibration study for the four manifest axes.
# For each axis, on data where the property HOLDS we want the false-violation
# rate ~ alpha across N/J; on data where it is VIOLATED we want power. Records
# the raw statistic and calibrated p (or resampling lower bound) per dataset so
# cutoffs can be checked/tuned rather than assumed. One CSV per dataset.
#
# Property map:
#   IIO holds : IIO, DM, LCR, RM   (violated: UN, MON)
#   MON holds : MON, DM, LCR, RM   (violated: UN, IIO)
#   ADD holds : LCR, RM            (violated: DM)      [only tested when doubly monotone]
#   DIP       : LCR discrete       (RM continuous)
suppressMessages({ library(QuantFit); library(diptest) })
env <- asNamespace("QuantFit")
iio_h <- get(".manifest_iio_holds", env)
mon_h <- get(".manifest_mon_holds", env)
add_h <- get(".manifest_add_holds", env)
rmt   <- get("refit_model_type", env)
dip_l <- get(".manifest_dip_lcr", env)

gen <- function(m, J, N, seed) { set.seed(seed)
  if (m == "RM") { b <- runif(J,-2,2); th <- rnorm(N)
    d <- matrix(rbinom(N*J,1,plogis(outer(th,b,"-"))),N,J) }
  else { d <- simulate_responses(m, n_persons=N, n_items=J, n_classes=3, seed=seed)
    d <- if(is.list(d)) d$data else d }
  storage.mode(d) <- "integer"; d }

B <- 49L; C <- 3L; ns <- 5L
cases <- expand.grid(model=c("UN","MON","IIO","DM","LCR","RM"),
                     J=c(6L,12L,24L), N=c(750L,1500L,3000L), rep=1:8,
                     stringsAsFactors=FALSE)
cases$id <- seq_len(nrow(cases)); set.seed(72301)
cases$seed <- sample.int(.Machine$integer.max, nrow(cases))
out <- Sys.getenv("ROBUST_OUT", "manifest_robust_out"); dir.create(out, showWarnings=FALSE)

run <- function(k) {
  cs <- cases[k,]; f <- file.path(out, sprintf("r%04d.csv", k))
  if (file.exists(f)) return(invisible())
  d <- gen(cs$model, cs$J, cs$N, cs$seed); ms <- cs$id
  iio <- tryCatch(iio_h(d, C, B, ns, TRUE, ms),            error=function(e) list(p=NA,holds=NA))
  mon <- tryCatch(mon_h(d, C, max(B%/%2L,20L), ns, TRUE, 0.03, ms+7L), error=function(e) list(stat=NA,lo=NA,holds=NA))
  add <- tryCatch(add_h(d, C, B, ns, TRUE, ms+11L),        error=function(e) list(stat=NA,p=NA,additive=NA))
  rmf <- tryCatch(suppressWarnings(fit_rm(d, verbose=FALSE)), error=function(e) NULL)
  dip <- if (is.null(rmf)) list(stat=NA,p=NA,discrete=NA) else
         tryCatch(dip_l(d, rmf, B, ms+13L), error=function(e) list(stat=NA,p=NA,discrete=NA))
  write.csv(data.frame(model=cs$model, J=cs$J, N=cs$N, rep=cs$rep,
    iio_p=round(iio$p,4), iio_holds=iio$holds,
    mon_stat=round(mon$stat,4), mon_lo=round(mon$lo,4), mon_holds=mon$holds,
    add_stat=round(add$stat,4), add_p=round(add$p,4), add_additive=add$additive,
    dip_stat=round(dip$stat,5), dip_p=round(dip$p,4), dip_discrete=dip$discrete),
    f, row.names=FALSE)
}
cat("manifest robustness:", nrow(cases), "datasets\n")
invisible(parallel::mclapply(seq_len(nrow(cases)),
  function(k) tryCatch(run(k), error=function(e) NULL), mc.cores=8))
cat("MANIFEST ROBUST DONE\n")
