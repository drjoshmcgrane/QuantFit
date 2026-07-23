# Small fresh PILOT (pre-locked-study shakedown): 6 models x J=12 x 3 reps.
# Set PILOT_LC_ONLY=true to validate the changed TI&D route without rerunning
# the unchanged CC/Omni routes. PILOT_B, PILOT_CORES (datasets), and
# PILOT_INNER_CORES (bootstrap replicates within a dataset) are configurable.
suppressMessages(library(QuantFit))
local({ set.seed(41); th <- rnorm(400); b <- runif(8,-2,2)+0.7
  d <- matrix(rbinom(400*8,1,plogis(outer(th,b,"-"))),400,8); storage.mode(d)<-"integer"
  dl <- QuantFit:::.cml_fit_general(d)
  stopifnot(cor(vapply(dl,sum,1), b) > 0.95, max(abs(unlist(dl))) < 20) })
lc_only <- tolower(Sys.getenv("PILOT_LC_ONLY", "false")) %in%
  c("1", "true", "yes")
pilot_B <- as.integer(Sys.getenv("PILOT_B", "99"))
pilot_cores <- as.integer(Sys.getenv("PILOT_CORES", "5"))
pilot_inner_cores <- as.integer(Sys.getenv("PILOT_INNER_CORES", "1"))
pilot_reps <- as.integer(Sys.getenv("PILOT_REPS", "3"))
pilot_tag <- Sys.getenv("PILOT_TAG", "")
stopifnot(is.finite(pilot_B), pilot_B >= 19L,
          is.finite(pilot_cores), pilot_cores >= 1L,
          is.finite(pilot_inner_cores), pilot_inner_cores >= 1L,
          is.finite(pilot_reps), pilot_reps >= 1L)
out <- if (lc_only) "pilot_results_tid032_lc" else "pilot_results_tid032"
if (nzchar(pilot_tag)) out <- paste0(out, "_", pilot_tag)
dir.create(out, showWarnings = FALSE)
grid <- expand.grid(model = c("UN","MON","IIO","DM","LCR","RM"),
                    rep = seq_len(pilot_reps),
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
                    mon_p=NA, iio_p=NA, dmiio_p=NA, dmmon_p=NA,
                    mon_raw=NA, iio_raw=NA, dmiio_raw=NA, dmmon_raw=NA,
                    mon_pre=NA, iio_pre=NA, dmiio_pre=NA, dmmon_pre=NA,
                    mon_warm=NA, iio_warm=NA, dmiio_warm=NA, dmmon_warm=NA,
                    lcrdm_p=NA, lcrdm_B_eff=NA, bridge_C=NA,
                    lcrdm_raw=NA, lcrdm_pre=NA, lcrdm_warm=NA,
                    lcrdm_boot_warm=NA,
                    rmlcr_p=NA, rmlcr_B_eff=NA, prof_C=NA,
                    cc_rej=NA, cc_p=NA, cc_B_eff=NA,
                    omni_rej=NA, omni_p=NA, omni_B_eff=NA,
                    ess_obs=NA, ess_null_med=NA, ess_null_low=NA,
                    secs=NA, err="")
  t0 <- proc.time()[3]
  lco <- tryCatch(suppressWarnings(select_model_ll(d, n_classes=1:5, B=pilot_B,
        n_starts=5, boot_n_starts=5, method="lattice", severity=FALSE,
        seed=ms, mc.cores=pilot_inner_cores)),
        error=function(e) NULL)
  if (!is.null(lco)) {
    row$lc <- lco$selected
    edge_p <- function(label) {
      z <- lco$tests$p_value[lco$tests$comparison == label]
      if (length(z)) round(z[[1L]], 4) else NA_real_
    }
    edge_raw <- function(label) {
      z <- lco$tests$raw_statistic[lco$tests$comparison == label]
      if (length(z)) round(z[[1L]], 6) else NA_real_
    }
    edge_warm <- function(label) {
      z <- lco$tests$general_warm_started[lco$tests$comparison == label]
      if (length(z)) z[[1L]] else NA
    }
    edge_pre <- function(label) {
      z <- lco$tests$pre_refinement_statistic[
        lco$tests$comparison == label]
      if (length(z)) round(z[[1L]], 6) else NA_real_
    }
    row$mon_p <- edge_p("MON vs UN")
    row$iio_p <- edge_p("IIO vs UN")
    row$dmiio_p <- edge_p("DM vs IIO")
    row$dmmon_p <- edge_p("DM vs MON")
    row$mon_raw <- edge_raw("MON vs UN")
    row$iio_raw <- edge_raw("IIO vs UN")
    row$dmiio_raw <- edge_raw("DM vs IIO")
    row$dmmon_raw <- edge_raw("DM vs MON")
    row$mon_pre <- edge_pre("MON vs UN")
    row$iio_pre <- edge_pre("IIO vs UN")
    row$dmiio_pre <- edge_pre("DM vs IIO")
    row$dmmon_pre <- edge_pre("DM vs MON")
    row$mon_warm <- edge_warm("MON vs UN")
    row$iio_warm <- edge_warm("IIO vs UN")
    row$dmiio_warm <- edge_warm("DM vs IIO")
    row$dmmon_warm <- edge_warm("DM vs MON")
    if (!is.null(lco$lcr_vs_dm)) {
      row$lcrdm_p <- round(lco$lcr_vs_dm$p_value, 4)
      row$lcrdm_B_eff <- lco$lcr_vs_dm$B_effective
      row$lcrdm_raw <- lco$lcr_vs_dm$raw_statistic
      row$lcrdm_pre <- lco$lcr_vs_dm$pre_refinement_statistic
      row$lcrdm_warm <- isTRUE(lco$lcr_vs_dm$general_warm_started)
      row$lcrdm_boot_warm <- lco$lcr_vs_dm$bootstrap_general_warm_started
    }
    if (!is.null(lco$quant_fits)) row$bridge_C <- lco$quant_fits$bridge_C
    if (!is.null(lco$rm_vs_lcr)) {
      row$rmlcr_p <- round(lco$rm_vs_lcr$p_value, 4)
      row$rmlcr_B_eff <- lco$rm_vs_lcr$B_effective
      row$prof_C <- lco$rm_vs_lcr$profiled_C
    }
  } else row$err <- "lc;"
  if (!lc_only) {
    cc <- tryCatch(suppressWarnings(cc_bootstrap_hierarchy(d, B=99, n.mat=500,
          seed=ms, mc.cores=1, verbose=FALSE)), error=function(e) NULL)
    if (!is.null(cc)) { row$cc_rej <- !cc$supports_quant
      row$cc_p <- round(min(cc$p_adjusted),4)
      row$cc_B_eff <- min(vapply(cc$levels, function(l) l$B, 1))
    } else row$err <- paste0(row$err,"cc;")
    ka <- tryCatch(suppressWarnings(omni_bootstrap_null(d,
          B=99, S=12000, N_synth=25, seed=ms, mc.cores=1, verbose=FALSE)),
          error=function(e) NULL)
    if (!is.null(ka) && is.list(ka)) { row$omni_rej <- ka$reject
      row$omni_p <- round(ka$p_value,4); row$omni_B_eff <- ka$B
      row$ess_obs <- round(ka$ess_min,1)
      row$ess_null_med <- round(ka$null_ess_median,1)
      row$ess_null_low <- ka$null_ess_low_n
    } else row$err <- paste0(row$err,"omni;")
  }
  row$secs <- round(proc.time()[3]-t0,1)
  write.csv(row, f, row.names=FALSE)
}
errlog <- file.path(out, "errors.log")
invisible(parallel::mclapply(grid$id, function(i)
  tryCatch(run1(i), error=function(e) {
    cat(sprintf("[%s] id %d: %s\n", format(Sys.time()), i,
        conditionMessage(e)), file = errlog, append = TRUE)
    NULL
  }), mc.cores = pilot_cores,
  mc.preschedule = FALSE))
cat("PILOT", sum(file.exists(file.path(out, sprintf("P%02d.csv", grid$id)))),
    "/", nrow(grid), "DONE\n")
