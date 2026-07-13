# Validation of QuantFit against the Torres Irribarra & Diakow simulation
# archive: 1080 datasets with known generating models. All three routes, one
# result file per dataset (crash-proof, resumable, snapshotable).
suppressMessages(library(QuantFit))

args <- commandArgs(trailingOnly = TRUE)
partition <- if (length(args) >= 1) toupper(args[1]) else "ALL"   # A, B, or ALL
n_cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 2L)
rep_cap <- if (length(args) >= 3) as.integer(args[3]) else 3L
base <- Sys.getenv("TID_BASE", getwd())
data_dir <- file.path(base, "tid_data")
res_dir <- file.path(base, "tid_results")

gm <- read.csv(file.path(data_dir, "generatingModels.csv"))
gc_ <- read.csv(file.path(data_dir, "generatingConditions.csv"))
meta <- merge(gm[, c("id", "genM")], gc_, by = "id")
models <- c("UN", "MON", "IIO", "DM", "LCR", "RM")   # genM 0..5
scale_of <- c(UN = "nominal", MON = "ordinal", IIO = "ordinal",
              DM = "ordinal", LCR = "quant", RM = "quant")

# Order: clean unidimensional first (the paper's core recovery target),
# stratified across models so snapshots stay balanced; contaminated after.
set.seed(42)
meta$clean <- meta$dCor == 0 & meta$nId2 == 0
meta <- meta[meta$rep <= rep_cap, ]
ord <- with(meta, order(!clean, nI, ave(runif(nrow(meta)), genM, FUN = rank), runif(nrow(meta))))
ids <- meta$id[ord]
# Two-machine split: A takes odd positions of the stratified order, B even -
# both partitions stay clean-first and model-balanced, with zero overlap.
if (partition == "A") ids <- ids[seq(1, length(ids), 2)]
if (partition == "B") ids <- ids[seq(2, length(ids), 2)]

run_one <- function(id) {
  out_file <- file.path(res_dir, paste0("TA", id, ".csv"))
  if (file.exists(out_file)) return(invisible())
  m <- meta[meta$id == id, ]
  row <- data.frame(id = id, genM = models[m$genM + 1], nI = m$nI,
                    slope = m$slope, dCor = m$dCor, nId2 = m$nId2,
                    rep = m$rep, clean = m$clean,
                    lc_selected = NA, lc_scale = NA, C_used = NA,
                    cc_p = NA, cc_reject = NA, cc_attrib = NA,
                    kara_p = NA, kara_reject = NA,
                    secs_lc = NA, secs_cc = NA, secs_kara = NA, err = "")
  dat <- tryCatch({
    e <- new.env(); load(file.path(data_dir, paste0("TA", id, ".Rdata")), envir = e)
    get(ls(e)[1], e)$obsData
  }, error = function(er) NULL)
  if (is.null(dat)) { row$err <- "load"; write.csv(row, out_file, row.names = FALSE); return(invisible()) }
  storage.mode(dat) <- "integer"

  # LC route: package defaults (auto class count, alpha_quant 0.05, power check)
  t0 <- proc.time()[3]
  lc <- tryCatch(suppressWarnings(
    select_model_ll(dat, n_classes = 1:6, B = 59, n_starts = 5,
                    boot_n_starts = 3, seed = 1, mc.cores = 1)),
    error = function(er) NULL)
  row$secs_lc <- round(proc.time()[3] - t0, 1)
  if (!is.null(lc)) {
    row$lc_selected <- lc$selected
    row$lc_scale <- scale_of[lc$selected]
    row$C_used <- lc$n_classes
  } else row$err <- paste0(row$err, "lc;")

  # CC route: sequential hierarchy, empirical latent null
  t0 <- proc.time()[3]
  cc <- tryCatch(suppressWarnings(
    cc_bootstrap_hierarchy(dat, B = 40, n.mat = 25, seed = 1,
                           mc.cores = 1, verbose = FALSE)),
    error = function(er) NULL)
  row$secs_cc <- round(proc.time()[3] - t0, 1)
  if (!is.null(cc)) {
    row$cc_reject <- !cc$supports_quant
    row$cc_attrib <- cc$attribution
    last <- cc$levels[[cc$stopped_at]]
    row$cc_p <- round(last$p_value, 4)
  } else row$err <- paste0(row$err, "cc;")

  # Kara route: empirical latent null
  t0 <- proc.time()[3]
  ka <- tryCatch(suppressWarnings(
    kara_bootstrap_null(dat, B = 12, S = 3000, N_synth = 25, seed = 1,
                        mc.cores = 1, verbose = FALSE)),
    error = function(er) NULL)
  row$secs_kara <- round(proc.time()[3] - t0, 1)
  if (!is.null(ka)) {
    row$kara_p <- round(ka$p_value, 4)
    row$kara_reject <- ka$reject
  } else row$err <- paste0(row$err, "kara;")

  write.csv(row, out_file, row.names = FALSE)
  invisible()
}

cat("Starting TI&D validation: partition", partition, "-", length(ids), "datasets,", n_cores, "workers\n")
invisible(parallel::mclapply(ids, function(i)
  tryCatch(run_one(i), error = function(e) NULL), mc.cores = n_cores))
cat("ALL DONE\n")
