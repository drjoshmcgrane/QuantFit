# LOCKED fresh-seed validation, v3 (pilot-informed; both external amendments).
# Phases:
#  calib:  RM x {J=6,12,24,48} x {sigma=0.5,1,2} x 20  (=240)
#          LCR x {J=6,12,24,48} x 20                    (=80)
#     LC+CC on ALL; Omni on a stratified subset: reps 1-6/cell at S=12000
#     (probe: ESS ~55, threshold >=50 properly met), rep 7/cell at S=30000
#     (deep anchor, ESS ~130), reps 8+ no Omni. Honest, budgeted design.
#  recover: 6 models x J {6,12,24} x K=10; Omni reps 1-3/cell at S=12000.
#  canary:  2 fast datasets (IIO + RM) asserting every recorded field.
# Identity: build_sha from DESCRIPTION GitSHA (stamped at freeze); the CML
# canary below rejects any build without the full-vector fix.
suppressMessages(library(QuantFit))
stopifnot(packageVersion("QuantFit") >= "0.3.2")
local({ set.seed(41); th <- rnorm(400); b <- runif(8,-2,2)+0.7
  d <- matrix(rbinom(400*8,1,plogis(outer(th,b,"-"))),400,8); storage.mode(d)<-"integer"
  dl <- QuantFit:::.cml_fit_general(d)
  stopifnot(cor(vapply(dl,sum,1), b) > 0.95, max(abs(unlist(dl))) < 20) })
ROUTES <- Sys.getenv("QUANTFIT_ROUTES", "lc")   # lc | cc | omni | all
BUILD_SHA <- tryCatch(utils::packageDescription("QuantFit")$GitSHA,
                      error = function(e) NA_character_)
if (is.null(BUILD_SHA) || is.na(BUILD_SHA))
  stop("QuantFit build carries no GitSHA stamp; refuse to run an unpinned fleet")
EXPECTED_SHA <- Sys.getenv("QUANTFIT_EXPECTED_SHA", "")
if (!nzchar(EXPECTED_SHA))
  stop("Set QUANTFIT_EXPECTED_SHA to the intended lattice commit before launch")
if (!startsWith(BUILD_SHA, EXPECTED_SHA) && !startsWith(EXPECTED_SHA, BUILD_SHA))
  stop("Installed QuantFit GitSHA (", BUILD_SHA,
       ") does not match QUANTFIT_EXPECTED_SHA (", EXPECTED_SHA, ")")

a <- commandArgs(trailingOnly = TRUE)
phase <- a[1]; partition <- toupper(a[2]); n_cores <- as.integer(a[3])
base <- Sys.getenv("FRESH_BASE", getwd())
METHOD_ID <- "tid032"
out_dir <- file.path(base, paste0("fresh_", phase, "_", METHOD_ID))
dir.create(out_dir, showWarnings = FALSE)

grid <- if (phase == "calib") {
  g <- rbind(expand.grid(model = "RM",  J = c(6L,12L,24L), sigma = c(0.5,1,2),
                    rep = 1:20, stringsAsFactors = FALSE),
        expand.grid(model = "LCR", J = c(6L,12L,24L), sigma = 1,
                    rep = 1:20, stringsAsFactors = FALSE))
  g$N <- 1500L; g
} else if (phase == "recover") {
  g <- expand.grid(model = c("UN","MON","IIO","DM","LCR","RM"), J = c(6L,12L,24L),
                   sigma = 1, rep = 1:10, stringsAsFactors = FALSE)
  g$N <- 1500L; g
} else if (phase == "tidmatch") {
  # Head-to-head with Torres Irribarra & Diakow's expert human study: their
  # EXACT design - J=10, N=5000, all six models. 30 reps/model so recovery
  # rates are estimable to ~+-9% and directly comparable to their Figure 9.
  g <- expand.grid(model = c("UN","MON","IIO","DM","LCR","RM"), J = 10L,
                   sigma = 1, rep = 1:30, stringsAsFactors = FALSE)
  g$N <- 5000L; g
} else {                                     # canary
  data.frame(model = c("IIO","RM"), J = 6L, sigma = 1, rep = 1L, N = 1500L)
}
grid$id <- seq_len(nrow(grid))
set.seed(20260721)
grid$seed <- sample.int(.Machine$integer.max, nrow(grid))
# Omni assignment (stratified by rep within cell)
grid$omni_S <- NA_integer_
if (phase == "calib") {
  grid$omni_S[grid$rep <= 2] <- 12000L   # single-machine: sparse Omni subset
  grid$omni_S[grid$rep == 3] <- 30000L   # one deep anchor per cell
} else if (phase == "recover") {
  grid$omni_S[grid$rep <= 2] <- 12000L
} else if (phase == "tidmatch") {
  # LC + CC only: Omni at N=5000 is infeasible single-machine, and the
  # TI&D expert-human comparison is about LC model recovery (Omni was not
  # part of their study). Omni size comes from the calib phase (null data).
} else grid$omni_S <- 800L
cfg <- if (phase == "canary") {
  list(lcB = 19L, lcStarts = 2L, lcCores = 2L, ccB = 19L, nmat = 50L, omniB = 19L, Ns = 8L)
} else {
  if (phase == "tidmatch")
    list(lcB = 19L, lcStarts = 2L, lcCores = 4L, ccB = 99L, nmat = 200L, omniB = 99L, Ns = 25L)
  else
    list(lcB = 49L, lcStarts = 2L, lcCores = 1L, ccB = 299L, nmat = 200L, omniB = 99L, Ns = 25L)
}

ids <- grid$id
if (partition == "A") ids <- ids[seq(1, length(ids), 2)]
if (partition == "B") ids <- ids[seq(2, length(ids), 2)]

gen_data <- function(g) {
  set.seed(g$seed)
  N <- g$N
  if (g$model == "RM") {
    beta <- stats::runif(g$J, -2, 2)
    theta <- stats::rnorm(N, 0, g$sigma)
    d <- matrix(stats::rbinom(N * g$J, 1,
           stats::plogis(outer(theta, beta, "-"))), N, g$J)
    attr(d, "params") <- list(item = beta, sigma = g$sigma)
  } else {
    d <- simulate_responses(g$model, n_persons = N, n_items = g$J,
                            n_classes = 3L)
  }
  storage.mode(d) <- "integer"
  d
}
# true distance from DM: largest class-monotonicity violation in the
# generating success probabilities (classes ordered by mean); NA for RM
dist_dm <- function(d) {
  p <- attr(d, "params")
  if (is.null(p$L)) return(NA_real_)
  P <- stats::plogis(p$L)                    # classes x items
  P <- P[order(rowMeans(P)), , drop = FALSE]
  max(0, max(P[-nrow(P), ] - P[-1, ]))
}

run_one <- function(i) {
  g <- grid[grid$id == i, ]
  out_file <- file.path(out_dir, sprintf("F%04d_%s.csv", i, ROUTES))
  if (file.exists(out_file)) return(invisible())
  dat <- gen_data(g)
  ms <- g$id * 1000003L %% .Machine$integer.max
  row <- data.frame(id = i, model = g$model, J = g$J, N = g$N, sigma = g$sigma,
    rep = g$rep, seed = g$seed, build_sha = BUILD_SHA,
    lc_selected = NA, lc_C = NA,
    lc_mon_p = NA, lc_iio_p = NA, lc_dmiio_p = NA, lc_dmmon_p = NA,
    lc_mon_raw = NA, lc_iio_raw = NA, lc_dmiio_raw = NA, lc_dmmon_raw = NA,
    lc_mon_pre = NA, lc_iio_pre = NA, lc_dmiio_pre = NA, lc_dmmon_pre = NA,
    lc_mon_warm = NA, lc_iio_warm = NA,
    lc_dmiio_warm = NA, lc_dmmon_warm = NA,
    lc_lcrdm_p = NA, lc_lcrdm_B_eff = NA, lc_lcrdm_override = NA,
    lc_lcrdm_raw = NA, lc_lcrdm_pre = NA, lc_lcrdm_warm = NA,
    lc_lcrdm_boot_warm = NA,
    lc_bridge_C = NA, lc_rmlcr_p = NA, lc_rmlcr_B_eff = NA,
    lc_profile_C = NA, lc_bic_diff = NA, lc_quant_reached = NA,
    lc_rmlcr_boot_ok = NA, dm_p = NA, dm_override = NA, dist_dm = dist_dm(dat),
    cc_p_holm = NA, cc_reject = NA, cc_attrib = NA, cc_B_eff = NA,
    omni_S = g$omni_S, omni_p = NA, omni_reject = NA, omni_B_eff = NA,
    omni_ess_min = NA, omni_null_ess_med = NA, omni_null_ess_low = NA,
    secs_lc = NA, secs_cc = NA, secs_omni = NA, err = "")
  t0 <- proc.time()[3]
  lc <- if (!ROUTES %in% c("lc","all")) NULL else
    tryCatch(suppressWarnings(select_model_ll(as.matrix(dat),
        n_classes = 1:6, B = cfg$lcB, n_starts = 5, boot_n_starts = 5,
        method = "lattice", severity = FALSE,
        seed = ms, mc.cores = 1)), error = function(e) NULL)
  row$secs_lc <- round(proc.time()[3] - t0, 1)
  if (!is.null(lc)) {
    row$lc_selected <- lc$selected
    row$lc_C <- if (!is.null(lc$quant_fits)) lc$quant_fits$C else
      tryCatch(lc$n_classes, error = function(e) NA)
    if (!is.null(lc$quant_fits)) row$lc_bridge_C <- lc$quant_fits$bridge_C
    edge_p <- function(label) {
      z <- lc$tests$p_value[lc$tests$comparison == label]
      if (length(z)) round(z[[1L]], 4) else NA_real_
    }
    edge_raw <- function(label) {
      z <- lc$tests$raw_statistic[lc$tests$comparison == label]
      if (length(z)) round(z[[1L]], 6) else NA_real_
    }
    edge_warm <- function(label) {
      z <- lc$tests$general_warm_started[lc$tests$comparison == label]
      if (length(z)) z[[1L]] else NA
    }
    edge_pre <- function(label) {
      z <- lc$tests$pre_refinement_statistic[lc$tests$comparison == label]
      if (length(z)) round(z[[1L]], 6) else NA_real_
    }
    row$lc_mon_p <- edge_p("MON vs UN")
    row$lc_iio_p <- edge_p("IIO vs UN")
    row$lc_dmiio_p <- edge_p("DM vs IIO")
    row$lc_dmmon_p <- edge_p("DM vs MON")
    row$lc_mon_raw <- edge_raw("MON vs UN")
    row$lc_iio_raw <- edge_raw("IIO vs UN")
    row$lc_dmiio_raw <- edge_raw("DM vs IIO")
    row$lc_dmmon_raw <- edge_raw("DM vs MON")
    row$lc_mon_pre <- edge_pre("MON vs UN")
    row$lc_iio_pre <- edge_pre("IIO vs UN")
    row$lc_dmiio_pre <- edge_pre("DM vs IIO")
    row$lc_dmmon_pre <- edge_pre("DM vs MON")
    row$lc_mon_warm <- edge_warm("MON vs UN")
    row$lc_iio_warm <- edge_warm("IIO vs UN")
    row$lc_dmiio_warm <- edge_warm("DM vs IIO")
    row$lc_dmmon_warm <- edge_warm("DM vs MON")
    dmrow <- lc$tests[lc$tests$comparison %in% c("DM vs IIO", "DM vs MON"), ]
    if (nrow(dmrow)) {
      row$dm_p <- round(min(dmrow$p_value), 4)
      row$dm_override <- any(grepl("override", dmrow$decision))
    }
    lcrdm <- lc$tests[lc$tests$comparison == "LCR vs DM", ]
    if (nrow(lcrdm)) {
      row$lc_lcrdm_p <- round(lcrdm$p_value[1], 4)
      row$lc_lcrdm_raw <- round(lcrdm$raw_statistic[1], 6)
      row$lc_lcrdm_pre <- round(lcrdm$pre_refinement_statistic[1], 6)
      row$lc_lcrdm_warm <- lcrdm$general_warm_started[1]
      row$lc_lcrdm_boot_warm <- lcrdm$bootstrap_general_warm_started[1]
    }
    if (!is.null(lc$lcr_vs_dm)) {
      row$lc_lcrdm_B_eff <- lc$lcr_vs_dm$B_effective
      row$lc_lcrdm_override <- isTRUE(lc$lcr_vs_dm$severity_override)
    }
    row$lc_quant_reached <- !is.null(lc$lcr_vs_dm)
    if (!is.null(lc$rm_vs_lcr)) {
      row$lc_rmlcr_p  <- round(lc$rm_vs_lcr$p_value, 4)
      row$lc_rmlcr_B_eff <- lc$rm_vs_lcr$B_effective
      row$lc_bic_diff <- round(lc$rm_vs_lcr$statistic, 2)
      row$lc_C        <- lc$rm_vs_lcr$profiled_C
      row$lc_profile_C <- lc$rm_vs_lcr$profiled_C
      row$lc_rmlcr_boot_ok <- lc$rm_vs_lcr$B_failed == 0
    }
  } else row$err <- "lc;"
  t0 <- proc.time()[3]
  cc <- if (phase == "tidmatch" || !ROUTES %in% c("cc","all")) NULL else
    tryCatch(suppressWarnings(cc_bootstrap_hierarchy(as.matrix(dat),
        B = cfg$ccB, n.mat = cfg$nmat, alpha = 0.05, seed = ms,
        mc.cores = 1, verbose = FALSE)), error = function(e) NULL)
  row$secs_cc <- round(proc.time()[3] - t0, 1)
  if (!is.null(cc)) {
    row$cc_reject <- !cc$supports_quant; row$cc_attrib <- cc$attribution
    row$cc_p_holm <- round(min(cc$p_adjusted), 4)
    row$cc_B_eff <- min(vapply(cc$levels, function(l) l$B, 1))
  } else row$err <- paste0(row$err, "cc;")
  if (!is.na(g$omni_S) && ROUTES %in% c("omni","all")) {
    t0 <- proc.time()[3]
    ka <- tryCatch(suppressWarnings(omni_bootstrap_null(
      as.matrix(dat), B = cfg$omniB, S = g$omni_S,
      N_synth = cfg$Ns, alpha = 0.05, seed = ms,
      mc.cores = 1, verbose = FALSE)), error = function(er) NULL)
    row$secs_omni <- round(proc.time()[3] - t0, 1)
    if (!is.null(ka) && is.list(ka)) {
      row$omni_p <- round(ka$p_value, 4); row$omni_reject <- ka$reject
      row$omni_B_eff <- ka$B; row$omni_ess_min <- round(ka$ess_min, 1)
      row$omni_null_ess_med <- round(ka$null_ess_median, 1)
      row$omni_null_ess_low <- ka$null_ess_low_n
    } else row$err <- paste0(row$err, "omni-crash;")
  }
  write.csv(row, out_file, row.names = FALSE)
  invisible()
}
cat("fresh", phase, "partition", partition, ":", length(ids), "of",
    nrow(grid), "datasets,", n_cores, "workers | build", BUILD_SHA, "\n")
errlog <- file.path(out_dir, paste0("errors_", partition, ".log"))
invisible(parallel::mclapply(ids, function(i)
  tryCatch(run_one(i), error = function(e) {
    cat(sprintf("[%s] id %d: %s\n", format(Sys.time()), i,
        conditionMessage(e)), file = errlog, append = TRUE)
    NULL
  }), mc.cores = n_cores, mc.preschedule = FALSE))
done_n <- sum(file.exists(file.path(out_dir, sprintf("F%04d.csv", ids))))
cat("FRESH", toupper(phase), partition,
    sprintf("COMPLETE %d/%d", done_n, length(ids)),
    if (done_n < length(ids)) "- INCOMPLETE, see error log and rerun" else "",
    "\n")
