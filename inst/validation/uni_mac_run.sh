#!/bin/bash
# =============================================================================
# QuantFit: stable-machine batch (uni Mac)
#
# Runs everything that was structurally unreachable on the laptop (its
# environment kills background jobs on an irregular ~15-50 min cycle, and
# datasets needing >50 min of CPU never complete there). All phases are
# RESUMABLE - one CSV per dataset, skip-if-exists - and results land inside the
# repo (QuantFit/inst/validation/tid_realdata_results/), which already carries
# the completed datasets as committed state, so only missing work runs.
#
# LAUNCH (from the QuantModelFitting directory, the one containing tid_data/):
#
#   git pull --ff-only
#   nohup caffeinate -dims bash QuantFit/inst/validation/uni_mac_run.sh \
#       > uni_mac_run.log 2>&1 & disown
#
# caffeinate -dims prevents display/idle/disk/system sleep for the lifetime of
# the run; nohup + disown lets you close the terminal. IMPORTANT: keep the
# machine PLUGGED IN and the LID OPEN (caffeinate cannot prevent lid-close
# sleep unless in clamshell mode with an external display).
#
# Watch progress:   tail -f uni_mac_run.log
# Stop:             pkill -f tid_hybrid_full
# Bring results back: cd QuantFit && git add inst/validation/tid_realdata_results \
#                     && git commit -m "uni mac: stable-machine batch results" && git push
#
# Env knobs: CORES (default 8), RUN48=1 to enable the optional nI=48 phase.
# Rough budget at 8 cores: Phase 1 ~3-4 h, Phase 2 ~1-2 h, Phase 3 ~8-12 h.
# =============================================================================
set -u
CORES="${CORES:-8}"
ROOT="$(pwd)"
if [ ! -d tid_data ] || [ ! -d QuantFit ]; then
  echo "ERROR: run from the QuantModelFitting directory (needs tid_data/ and QuantFit/)"
  exit 1
fi
OUT="$ROOT/QuantFit/inst/validation/tid_realdata_results"
mkdir -p "$OUT"
stamp() { date "+%Y-%m-%d %H:%M:%S"; }
phase() { echo ""; echo "=== [$(stamp)] $1 ==="; }

phase "Install current QuantFit"
R CMD INSTALL QuantFit --no-docs || { echo "INSTALL FAILED"; exit 1; }

# ---------------------------------------------------------------------------
# Phase 1 - HIGHEST VALUE: the 18 stalled misspecified quant-edge datasets
# (non-clean, nI=24, full N=5000: 2 DM + 2 LCR + 14 RM). The hybrid's 2x2
# reached DM on these and the quant edge never completed on the laptop.
# This resolves the misspecification finding's open question: does "fails
# better" (structure preserved to DM) become "partially survives" (RM correctly
# recovered where the lattice collapsed to MON/UN)?
# ---------------------------------------------------------------------------
phase "Phase 1: non-clean nI=24 quant-edge datasets (the stalled 18)"
TIDF_OUT="$OUT" TIDF_SET=nonclean TIDF_NI=24 TIDF_B=49 TIDF_CORES="$CORES" \
  Rscript QuantFit/inst/validation/tid_hybrid_full.R

# ---------------------------------------------------------------------------
# Phase 2 - finish the CLEAN real-data head-to-head quant cells at nI 6,12
# (LCR was 5/18 covered, RM 0/18; ~31 datasets, 10-21 min each).
# ---------------------------------------------------------------------------
phase "Phase 2: clean nI=6,12 remaining LCR/RM cells"
TIDF_OUT="$OUT" TIDF_SET=union TIDF_NI=6,12 TIDF_B=49 TIDF_CORES="$CORES" \
  Rscript QuantFit/inst/validation/tid_hybrid_full.R

# ---------------------------------------------------------------------------
# Phase 3 - the clean nI=24 cells (never attempted on the laptop; quant-edge
# datasets here carry the 13-class Lindsay bridge, ~45-90 min each).
# ---------------------------------------------------------------------------
phase "Phase 3: clean nI=24 cells"
TIDF_OUT="$OUT" TIDF_SET=union TIDF_NI=24 TIDF_B=49 TIDF_CORES="$CORES" \
  Rscript QuantFit/inst/validation/tid_hybrid_full.R

# ---------------------------------------------------------------------------
# Phase 4 (OPTIONAL, RUN48=1) - nI=48: Lindsay bridge = 25 classes. Very
# expensive; only attempt on a machine that can run for days.
# ---------------------------------------------------------------------------
if [ "${RUN48:-0}" = "1" ]; then
  phase "Phase 4 (optional): nI=48"
  TIDF_OUT="$OUT" TIDF_SET=union TIDF_NI=48 TIDF_B=49 TIDF_CORES="$CORES" \
    Rscript QuantFit/inst/validation/tid_hybrid_full.R
fi

# ---------------------------------------------------------------------------
# Summary: print the answer tables into the log.
# ---------------------------------------------------------------------------
phase "Summary"
Rscript - <<'RS'
res <- "QuantFit/inst/validation/tid_realdata_results"
d <- do.call(rbind, lapply(list.files(res, "^h_", full.names = TRUE), read.csv))
lv <- c("UN","MON","IIO","DM","LCR","RM")
cln <- d$clean %in% c(TRUE, "TRUE")

cat("\n--- Phase 1 answer: the misspecified DM-reaching datasets ---\n")
nc <- d[!cln & d$truth %in% c("DM","LCR","RM"), ]
if (nrow(nc)) {
  print(table(truth = factor(nc$truth, lv), hybrid = factor(nc$hybrid, lv)))
  rmn <- nc[nc$truth == "RM", ]
  cat(sprintf("misspecified RM recovered as RM: %d/%d (lattice baseline: 0/18)\n",
              sum(rmn$hybrid == "RM"), nrow(rmn)))
}

cat("\n--- Clean head-to-head coverage & per-model recovery (hybrid vs lattice) ---\n")
cl <- d[cln, ]
for (ni in sort(unique(cl$nI))) {
  s <- cl[cl$nI == ni, ]
  cat(sprintf("nI=%d (n=%d): ", ni, nrow(s)))
  for (m in lv) { t <- s[s$truth == m, ]
    if (nrow(t)) cat(sprintf("%s %d/%d|%d ", m, sum(t$hybrid == m), nrow(t),
                             sum(t$lattice == m))) }
  cat("\n")
}
cat("\n(pair with tid_results/ for full lattice comparison; analysis scripts in\n inst/validation/. Commit tid_realdata_results/ and push to bring these home.)\n")
RS

phase "ALL PHASES COMPLETE"
