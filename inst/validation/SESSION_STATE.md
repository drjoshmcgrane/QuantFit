# SESSION STATE SNAPSHOT (2026-08-16) - read this first in a new session

## The frozen validation pipeline
- **Build f413264 = THE frozen pipeline** (stamped, installed on laptop).
  All constants sim-derived or analytic. USER DECISION: TID is presented as
  SIMPLE VALIDATION - "calibrated on simulation, applied once to TID,
  compared against sealed human ratings" - made TRUE by the one-shot run
  below, not asserted over the adaptive history.
- v0.3.4 axis fixes (both in f413264):
  1. IIO-axis null class count BIC-selected 2:6 floored at routing C
     (fixed-C null falsely rejected: 95% on DM/C4/J48 -> 10%;
     iio_axis_null_finding.md; 640-dataset study committed).
  2. Adaptive Monte-Carlo precision: borderline axis p-values re-decided at
     B_refine=999; band is ANALYTIC: (1.5/(B+1), alpha+3*sqrt(a(1-a)/B)) =
     (0.030, 0.143) at B=49 - no data-tuned constants
     (axis_monte_carlo_finding.md; at B=49 only 2/9 borderline datasets
     reproduced across seeds).
- Design rule established (axis_fixed_C_finding.md, 1360 datasets): fix C
  where the fit makes the STATISTIC (MON axis: size 0/240 at fixed C=3,
  true-C would be 15-100% false rej); adapt C where the fit makes the NULL
  (IIO axis). MON axis and poset UNCHANGED, correctly so.
- Quant gate: IDENTICAL machinery in both selectors; only the route to DM
  differs (2x2 vs LR-edge walk). Gate refusals matched 9:9.

## RUNNING on laptop: /tmp/tid_final.R -> qf_evidence/tid_final/
One-shot uniform hybrid run, all 900 nI<=24 TID datasets, frozen f413264,
6 cores, resumable (v_<id>.csv per dataset), started 2026-08-15 night,
~2-3 days. THE FINAL TABLES COME FROM THIS RUN ONLY (hybrid arm; lattice
rows in tid_extension_results/ + tid_realdata_results/ are single-config
and stand). On completion: merge with lattice + rater data
(tid_rater_merged.csv), produce final five-way tables (overall + per-nI +
misspecified), update consolidated + findings, restore R CMD check, final
stamp. Crashers TA218/376 will appear as NA rows (record as crashes;
TA552 completed earlier - known_crashers.md is CORRECTED: intermittent,
not deterministic; 30+ retry attempts 0 successes for 218/376).

## RUNNING on uni Mac + laptop shares: nI=48 quant block
results48/ in the QF_reval/extend bundle. Ordinal 90/90 DONE and committed
(tid_ni48_finding.md: hybrid 94%, IIO 83% vs humans/minBIC 0/30 = 0%).
Quant (DM/LCR/RM x30): ~5/90 done, ~20h/cell, USER SAID LET IT FINISH
(weeks). Old-build DM artifact rows were cleared; runs on 5d3f8fc/697bcae
mix via run48fixed.sh (uni Mac) + /tmp/laptop48q.sh (laptop, may need
relaunch after tid_final finishes). NOTE: these predate f413264; when the
block completes decide whether to re-run its handful of DM cells under
f413264 or annotate (axis fix matters exactly at J=48 DM).

## Evidence landmarks (all committed)
- Delta scans prove build-equivalence: axis_delta_scan2 (898 datasets,
  old-vs-f413264ish config): 7 flips only; regen7/ has their new verdicts.
- Pre-freeze headline tables (superseded by tid_final when it lands):
  clean 537: hybrid 92/lattice 83/RPD 75/DTI 68/minBIC 65; scale 97/96/83;
  false-quant 0.0%/0.8% vs 24%; demotions 7.9%/11.2%. Misspecified 360:
  everyone collapses vs labels; selectors 118/120 refusals = designed.
  nI-split + rater findings in tid_rater_finding.md.
- PLAN_gate_improvement.md / PLAN_triangulation.md: the two queued studies
  (gate: severity port + info-scaled min_effect + MC hygiene; sim-first).
  Axis studies partly superseded the axis-size item there.

## Discipline reminders
- Stamp-then-install before ANY evidence run (validation_shared.R gate).
- Evidence dirs OUTSIDE the checkout; never /tmp for >1-day artifacts
  (macOS purge ate a fleet once).
- Tune on simulation; TID is validation-only now (user directive).
- No interim reads of cost-ordered runs; complete blocks only.
