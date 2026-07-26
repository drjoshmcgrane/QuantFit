#!/bin/bash
# =============================================================================
# QuantFit stable-machine batch - iCLOUD BUNDLE VARIANT (uni Mac)
#
# Everything needed travels in one iCloud folder (the established channel:
# freshB/tidB used the same pattern). Results are written INTO the bundle, so
# they sync back to the home machine automatically as they are produced.
#
# Bundle layout (created on the home machine):
#   QF_uniMacBatch/
#     run.sh                     <- this script
#     QuantFit_<ver>.tar.gz      <- current package build
#     tid_hybrid_full.R          <- runner
#     tid_lattice_baseline/      <- 324 lattice-baseline CSVs
#     tid_realdata_results/      <- resumable state (h_*.csv) + NEW RESULTS
#
# LAUNCH ON THE UNI MAC - fully self-contained, run from anywhere (the bundle
# carries a needed-files-only subset of tid_data, ~83 MB):
#
#   B=~/Library/Mobile\ Documents/com~apple~CloudDocs/QF_uniMacBatch
#   nohup caffeinate -dims bash "$B/run.sh" > ~/uni_mac_run.log 2>&1 & disown
#
# caffeinate -dims keeps display/idle/disk/system awake for the whole run;
# nohup+disown survives closing the terminal. KEEP IT PLUGGED IN, LID OPEN
# (lid-close sleep overrides caffeinate without clamshell+external display).
#
# Watch:  tail -f uni_mac_run.log      Stop:  pkill -f tid_hybrid_full
# Results appear on the home machine via iCloud sync of the bundle.
# Env knobs: CORES (default 8), RUN48=1 for the optional nI=48 phase.
# =============================================================================
set -u
CORES="${CORES:-8}"
BUNDLE="$(cd "$(dirname "$0")" && pwd)"
stamp() { date "+%Y-%m-%d %H:%M:%S"; }
phase() { echo ""; echo "=== [$(stamp)] $1 ==="; }

# TA archive: use a full local tid_data/ if the cwd has one, else the bundled
# needed-files-only subset (covers phases 1-3; phase 4 / nI=48 requires the
# full archive).
if [ -d tid_data ]; then export TIDF_DATA="$PWD/tid_data"
else export TIDF_DATA="$BUNDLE/tid_data"; fi
echo "TA archive: $TIDF_DATA"
[ -d "$TIDF_DATA" ] || { echo "ERROR: no tid_data found"; exit 1; }

phase "Materialise any iCloud-evicted files in the bundle"
if find "$BUNDLE" -name "*.icloud" | grep -q .; then
  echo "dataless files found - requesting download"
  brctl download "$BUNDLE" 2>/dev/null || true
  for i in $(seq 1 60); do
    find "$BUNDLE" -name "*.icloud" | grep -q . || break; sleep 10
  done
  find "$BUNDLE" -name "*.icloud" | grep -q . && {
    echo "ERROR: bundle still not fully downloaded"; exit 1; }
fi

phase "Install missing CRAN dependencies"
Rscript -e 'need <- c("Rcpp","alabama","nloptr","numDeriv")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")'

phase "Install packaged QuantFit"
TARBALL="$(ls "$BUNDLE"/QuantFit_*.tar.gz | sort | tail -1)"
R CMD INSTALL "$TARBALL" || { echo "INSTALL FAILED"; exit 1; }

OUT="$BUNDLE/tid_realdata_results"; mkdir -p "$OUT"
export TIDF_BASELINE="$BUNDLE/tid_lattice_baseline"
RUNNER="$BUNDLE/tid_hybrid_full.R"

phase "Phase 1: the 18 stalled misspecified quant-edge datasets (non-clean nI=24)"
TIDF_OUT="$OUT" TIDF_SET=nonclean TIDF_NI=24 TIDF_B=49 TIDF_CORES="$CORES" \
  Rscript "$RUNNER"

phase "Phase 2: remaining clean nI=6,12 LCR/RM cells"
TIDF_OUT="$OUT" TIDF_SET=union TIDF_NI=6,12 TIDF_B=49 TIDF_CORES="$CORES" \
  Rscript "$RUNNER"

phase "Phase 3: clean nI=24 cells"
TIDF_OUT="$OUT" TIDF_SET=union TIDF_NI=24 TIDF_B=49 TIDF_CORES="$CORES" \
  Rscript "$RUNNER"

if [ "${RUN48:-0}" = "1" ]; then
  phase "Phase 4 (optional): nI=48"
  TIDF_OUT="$OUT" TIDF_SET=union TIDF_NI=48 TIDF_B=49 TIDF_CORES="$CORES" \
    Rscript "$RUNNER"
fi

phase "Summary"
OUT="$OUT" Rscript - <<'RS'
res <- Sys.getenv("OUT")
d <- do.call(rbind, lapply(list.files(res, "^h_", full.names = TRUE), read.csv))
lv <- c("UN","MON","IIO","DM","LCR","RM"); cln <- d$clean %in% c(TRUE, "TRUE")
cat("\n--- Phase 1 answer: misspecified DM-reaching datasets ---\n")
nc <- d[!cln & d$truth %in% c("DM","LCR","RM"), ]
if (nrow(nc)) {
  print(table(truth = factor(nc$truth, lv), hybrid = factor(nc$hybrid, lv)))
  rmn <- nc[nc$truth == "RM", ]
  cat(sprintf("misspecified RM recovered as RM: %d/%d (lattice baseline: 0/18)\n",
              sum(rmn$hybrid == "RM"), nrow(rmn)))
}
cat("\n--- clean coverage & per-model recovery (hybrid n / total | lattice n) ---\n")
cl <- d[cln, ]
for (ni in sort(unique(cl$nI))) { s <- cl[cl$nI == ni, ]
  cat(sprintf("nI=%d (n=%d): ", ni, nrow(s)))
  for (m in lv) { t <- s[s$truth == m, ]
    if (nrow(t)) cat(sprintf("%s %d/%d|%d ", m, sum(t$hybrid == m), nrow(t),
                             sum(t$lattice == m))) }
  cat("\n") }
cat("\nResults sync home automatically via iCloud (this folder).\n")
RS

phase "ALL PHASES COMPLETE"
