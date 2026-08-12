# EXTENSION runner: the unrun IIO and DM datasets of the TI&D grid at
# nI <= 24, BOTH selectors at full N = 5000 - tightening the cells where the
# human-vs-machine gap lives (rater IIO 2-10%, DM 9-22% vs hybrid 70%/89%).
#
# Unlike tid_hybrid_full.R (which pairs against the prior lattice baseline),
# this runner RUNS the lattice itself: the extension sets have no baseline.
# Bundle-mode provenance: no git checkout on the fleet machine, so the rows
# record the INSTALLED build's stamped sha (quantfit_build_info) and the
# runner refuses an unstamped or method-mismatched install; the drift/dirty
# checks that require git are documented as not applicable off-checkout.
#
# Env: TIDX_DATA (TA archive / subset), TIDX_OUT (results dir; resumable),
#      TIDX_DONE (file of already-covered ids, one per line),
#      TIDX_B (49), TIDX_CORES (8), TIDX_MODELS ("IIO,DM")
suppressMessages(library(QuantFit))
DD    <- Sys.getenv("TIDX_DATA", "tid_data")
out   <- Sys.getenv("TIDX_OUT", "tid_extend_out")
B     <- as.integer(Sys.getenv("TIDX_B", "49"))
cores <- as.integer(Sys.getenv("TIDX_CORES", "8"))
want  <- strsplit(Sys.getenv("TIDX_MODELS", "IIO,DM"), ",")[[1]]
nImax <- as.integer(Sys.getenv("TIDX_NIMAX", "24"))
dir.create(out, showWarnings = FALSE, recursive = TRUE)

bi <- quantfit_build_info()
if (is.na(bi$sha) || identical(bi$sha, "unstamped"))
  stop("installed QuantFit build is not stamped; rebuild via stamp-then-install")
METHOD <- "max-T-studentized-v4"
if (!identical(bi$po_method, METHOD))
  stop("installed poset method ", bi$po_method, " != ", METHOD)
SHA <- as.character(bi$sha)
CONFIG <- sprintf("B=%d;N=full;hybC=3;latC=2:6;alpha=0.05", B)

lv <- c("UN","MON","IIO","DM","LCR","RM")
scale_of <- c(UN="nominal",MON="ordinal",IIO="ordinal",DM="ordinal",
              LCR="quant",RM="quant")
gm <- read.csv(file.path(DD, "generatingModels.csv"))
gc_ <- read.csv(file.path(DD, "generatingConditions.csv"))
names(gc_)[1] <- "id"
gm$genM <- lv[gm$genM + 1L]
sel <- merge(gm[, c("id","genM")], gc_[, c("id","nI","slope","dCor","nId2")],
             by = "id")
done_f <- Sys.getenv("TIDX_DONE", "")
done <- if (nzchar(done_f) && file.exists(done_f))
  as.integer(readLines(done_f)) else integer(0)
sel <- sel[sel$genM %in% want & sel$nI <= nImax & !(sel$id %in% done), ]
sel$clean <- sel$dCor == 0 & sel$nId2 == 0
sel <- sel[order(match(sel$genM, c("IIO","DM")), !sel$clean, sel$nI), ]
cat("TI&D extension:", nrow(sel), "datasets (", paste(want, collapse="+"),
    ", nI<=", nImax, ", excluding", length(done), "covered ) | build", SHA, "\n")

load_d <- function(id) { e <- new.env()
  load(file.path(DD, sprintf("TA%d.Rdata", id)), envir = e)
  d <- get(ls(e)[1], e)$obsData; storage.mode(d) <- "integer"; d }
g <- function(x, fld) if (is.null(x) || is.null(x[[fld]])) NA else x[[fld]]

run <- function(k) {
  cs <- sel[k, ]; f <- file.path(out, sprintf("x_%d.csv", cs$id))
  if (file.exists(f)) {
    prev <- tryCatch(read.csv(f, nrows = 1), error = function(e) NULL)
    if (!is.null(prev) && identical(as.character(prev$sha), SHA) &&
        identical(as.character(prev$config), CONFIG)) return(invisible())
    unlink(f)
  }
  d <- load_d(cs$id); t0 <- proc.time()[3]
  rh <- select_model_hybrid(d, n_classes = 3L, B = B, lr_boot_n_starts = 2L,
                            mc.cores = 1L, seed = 1, verbose = FALSE)
  if (is.null(rh) || is.na(rh$selected)) stop("hybrid returned no verdict")
  rl <- select_model_ll(d, n_classes = 2:6, B = B, n_starts = 5,
                        boot_n_starts = 2L, method = "lattice", seed = 1,
                        mc.cores = 1L, verbose = FALSE)
  if (is.null(rl) || is.na(rl$selected)) stop("lattice returned no verdict")
  write.csv(data.frame(method = METHOD, sha = SHA, config = CONFIG,
    id = cs$id, truth = cs$genM, nI = cs$nI, clean = cs$clean, N = nrow(d),
    hybrid = rh$selected, lattice = rl$selected,
    h_class_shape = g(rh$poset$class, "shape"),
    h_item_shape = g(rh$poset$item, "shape"),
    l_class_shape = g(rl$poset$class, "shape"),
    l_poset_refusal = if (is.null(rl$poset_refusal)) "" else rl$poset_refusal,
    h_quant_failed = isTRUE(rh$quant_edge_failed),
    l_quant_failed = isTRUE(rl$quant_edge_failed),
    secs = round(proc.time()[3] - t0, 1)), f, row.names = FALSE)
}
fails <- parallel::mclapply(seq_len(nrow(sel)), function(k)
  tryCatch({ suppressWarnings(run(k)); NULL }, error = function(e)
    data.frame(id = sel$id[k], truth = sel$genM[k],
               error = conditionMessage(e))),
  mc.cores = cores, mc.preschedule = FALSE)
fails <- do.call(rbind, Filter(Negate(is.null), fails))
if (!is.null(fails) && nrow(fails)) {
  write.csv(fails, file.path(out, "FAILURES.csv"), row.names = FALSE)
  cat("TIDX FAILED:", nrow(fails), "datasets errored\n")
  quit(save = "no", status = 1L)
}
unlink(file.path(out, "FAILURES.csv"))
cat("TIDX DONE\n")
