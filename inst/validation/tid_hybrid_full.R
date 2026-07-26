# HYBRID on the real TI&D data at FULL N=5000, on exactly the datasets the prior
# lattice run (tid_results/, tid_runner.R) already covers -- so the two pair by
# TA id and give the first REAL-DATA head-to-head.
#
# Baseline being paired against (lattice, full N=5000, B=59, n_classes=1:6):
#   clean=TRUE  (216): exact 81%  | UN 36/36 MON 36/36 IIO 15/36 DM 33/36
#                                   LCR 31/36 RM 23/36
#   clean=FALSE (108): exact 27%  | IIO 0/18 DM 1/18 LCR 1/18 RM 0/18
#     (non-clean = slope variation / correlated dimensions / extra structure,
#      i.e. genuine misspecification - the lattice collapses there)
#
# CAVEAT to carry into any writeup: the lattice arm is PRIOR-SESSION code at
# different settings (B=59, boot_n_starts=5). So the overall comparison is
# indicative, not same-code. The IIO comparison specifically is safe: the
# lattice's ~42% IIO on this data agrees with a current same-code run on
# simulated data (45%), the older audit (~40%), and TI&D's human raters (41%).
#
# Env: TIDF_NI (default 6,12,24 - nI=48 is excluded by default because its
#      Lindsay bridge grain is ceiling(49/2)=25 classes, which makes the quant
#      edge extremely expensive at N=5000), TIDF_CLEAN (1=clean only, 0=all),
#      TIDF_B, TIDF_CORES.
# Resumable: one CSV per TA id.
suppressMessages(library(QuantFit))
DD <- "tid_data"
models <- c("UN","MON","IIO","DM","LCR","RM")
scale_of <- c(UN="nominal",MON="ordinal",IIO="ordinal",DM="ordinal",
              LCR="quant",RM="quant")
nIs   <- as.integer(strsplit(Sys.getenv("TIDF_NI","6,12,24"), ",")[[1]])
clean_only <- Sys.getenv("TIDF_CLEAN","1") == "1"
B     <- as.integer(Sys.getenv("TIDF_B","49"))     # the hybrid.s validated config
cores <- as.integer(Sys.getenv("TIDF_CORES","8"))

# the exact datasets the lattice baseline covers. tid_results/ lives in the
# working directory on the original machine; on other machines (e.g. the uni
# Mac) fall back to the copy committed inside the package repo.
res_dir <- Sys.getenv("TIDF_BASELINE", "")
if (!nzchar(res_dir))
  res_dir <- if (dir.exists("tid_results")) "tid_results" else
    "QuantFit/inst/validation/tid_lattice_baseline"
prior <- do.call(rbind, lapply(list.files(res_dir, full.names = TRUE), read.csv))
prior <- prior[!is.na(prior$lc_selected) & prior$lc_selected != "", ]
# TIDF_SET: "clean"  = the clean conditions (the head-to-head vs the lattice)
#           "lat_un" = every dataset the lattice called UN (the partial-order
#                      question: does a better ordinal layer recover them?)
#           "union"  = both (default; one pass, resumable, no stacked jobs)
# NOTE: all of the lattice's WRONG UN calls are in the non-clean conditions -
# on clean data it never mislabelled anything UN. So the partial-order question
# lives entirely in the misspecified conditions.
set_which <- Sys.getenv("TIDF_SET", "union")
base <- prior[prior$nI %in% nIs, ]
s_clean <- base[base$clean %in% c(TRUE, "TRUE"), ]
s_latun <- base[base$lc_selected == "UN", ]
sel <- switch(set_which,
  clean    = s_clean,
  lat_un   = s_latun,
  # "nonclean": ALL misspecified-condition datasets (slope variation /
  # correlated dimensions; every one is nI=24). This is the regime where the
  # prior lattice run COLLAPSED (27% exact; IIO 0/18, RM 0/18, nearly all
  # verdicts MON or UN). Most resolve in the fast ordinal layer even at nI=24;
  # only 2x2-DM-reaching datasets pay the (locally unreachable) 13-class
  # bridge and may not complete - report coverage explicitly.
  nonclean = base[!(base$clean %in% c(TRUE, "TRUE")), ],
  base[base$id %in% union(s_clean$id, s_latun$id), ])
if (identical(set_which, "clean") && !clean_only) sel <- base

# Process cheap-and-decisive datasets FIRST. UN/MON/IIO resolve in the 2x2 and
# never touch the quant edge (~10-20s); DM/LCR/RM pay 10-21 min for it. Ordering
# by expected cost means the IIO comparison -- the load-bearing claim -- completes
# early and survives interruptions, instead of queueing behind the quant cases.
sel <- sel[order(match(sel$genM, c("IIO","MON","UN","DM","LCR","RM"))), ]

out <- Sys.getenv("TIDF_OUT","tid_hybrid_full_out"); dir.create(out, showWarnings=FALSE)
load_d <- function(id) { e <- new.env()
  load(file.path(DD, sprintf("TA%d.Rdata", id)), envir = e)
  d <- get(ls(e)[1], e)$obsData; storage.mode(d) <- "integer"; d }

run <- function(k) {
  cs <- sel[k, ]; f <- file.path(out, sprintf("h_%d.csv", cs$id))
  if (file.exists(f)) return(invisible())
  d <- load_d(cs$id)                      # FULL N = 5000, no subsampling
  t0 <- proc.time()[3]
  r <- tryCatch(select_model_hybrid(d, n_classes = 3L, B = B,
           lr_boot_n_starts = 2L, mc.cores = 1L, seed = 1, verbose = FALSE),
         error = function(e) NULL)
  s <- if (is.null(r)) NA_character_ else r$selected
  write.csv(data.frame(id = cs$id, truth = cs$genM, nI = cs$nI,
    clean = cs$clean, N = nrow(d),
    truth_scale = scale_of[cs$genM], hybrid = s,
    hybrid_scale = if (is.na(s)) NA else scale_of[s],
    lattice = cs$lc_selected,
    secs = round(proc.time()[3] - t0, 1)), f, row.names = FALSE)
}
cat("TI&D hybrid FULL N:", nrow(sel), "datasets (nI", paste(nIs, collapse=","),
    "| clean_only", clean_only, "| B", B, ")\n")
invisible(parallel::mclapply(seq_len(nrow(sel)),
  function(k) tryCatch(run(k), error = function(e) NULL), mc.cores = cores,
  mc.preschedule = FALSE))
cat("TIDF DONE\n")
