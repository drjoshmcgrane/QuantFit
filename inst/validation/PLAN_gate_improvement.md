# STUDY PLAN 1: Quantitative-gate improvement (drafted 2026-08-13)

## Motivation (measured facts)
Clean TA grid (537/538, five-way table in tid_rater_finding.md): hybrid
quant->qual demotions 14/178 (7.9%); lattice 20/178 (11.2%). Decomposition:
- 9 hybrid demotions = GATE refusals (LCR-vs-DM rejects) -> land on DM.
  Concentrated at nI=24 (11.9% demotion rate there vs 5-7% at nI 6/12) -
  severity-at-scale: min_effect = 1 raw LL unit is ~zero at N=5000xJ=24.
- 5 hybrid demotions = IIO-AXIS false rejections on true-quant data -> land
  on MON (never reach the gate). Separate mechanism, needs its own fix.
- Only 5 datasets demoted by BOTH selectors => irreducible data-ambiguity
  core ~3% (matches the historical power-check estimate 2/60).
- Within-quant confusions (scale-correct): hybrid 5 (4 RM->LCR), lattice 10
  (9 RM->LCR; grain-profiling discreteness bias).

## Levers (in priority order)
1. severity=TRUE for the gate (ll_general_null power check; ALREADY
   IMPLEMENTED in select_ll, opt-in). Historical K=30 measurement: quant
   demotions 7/60 -> 2/60 with no false-quant cost. The hybrid edge does NOT
   currently implement severity - needs porting into .hybrid_quant_edge.
2. Information-scaled min_effect: replace the raw-LL margin (1.0) with an
   equivalence margin scaled per observation, e.g. min_effect = c * N * J
   with c calibrated (or expressed as a probability-scale discrepancy like
   the axes' eps). This is an EQUIVALENCE-TEST reformulation: H0 becomes
   "constrained adequate up to negligible deviation". Direction of effect:
   retention easier => fixes over-rejection at scale. CAUTION: do NOT apply
   blindly to the lattice's ordinal edges - their failure mode is
   OVER-retention (degenerate nulls, IIO ceiling); margins would worsen it.
   Uniform philosophy, per-boundary margins.
3. Monte Carlo hygiene: gate B 49 -> 199; lr_boot_n_starts 2 -> 3-5.
   Under-optimized general-model refits in replicates SHRINK null LRs and
   bias toward demotion.
4. alpha last (measured trade at K=30: 0.05->0.01 gives recovery 53->58/60,
   false-quant 1->4/90).
5. Separate small study: IIO-axis size on true-RM/LCR data (the 5 ->MON
   demotions). Check the axis's parametric IIO null calibration on
   Rasch-generated data; if size inflated, the axis null may need the
   equivalence margin too (its stat is already probability-scale).

## Design discipline (bias guard - user-endorsed)
TUNE ON SIMULATION, CONFIRM ON TA ONCE. The TA grid is development-adjacent
(documented in the honesty section); using its demotions to pick the config
and then reporting the same grid = tuning-on-test. Protocol:
- Fresh simulation grid (simulate_responses, NEW seed range): truths
  LCR/RM/DM/IIO (DM/IIO for false-quant side), N in {1500, 5000},
  J in {6, 12, 24}, K>=20/cell, B=199.
- Cross levers: severity {off,on} x min_effect {1, c*N*J for c in 2-3
  values} x boot_n_starts {2,5} x alpha {0.05, 0.01}. Record per-cell
  quant sensitivity + false-quant rate -> ROC; pick operating point per the
  documented loss function (zero-false-quant priority).
- Run the chosen config ONCE on the TA quant cells; report as confirmation.
- Rows must carry method/sha/config canaries (stamp-then-install workflow;
  see validation_shared.R; outputs OUTSIDE the checkout; never /tmp).

## Deliverables
- New finding: gate_operating_point_finding.md with the ROC + chosen point.
- Possible code changes: severity in .hybrid_quant_edge; min_effect
  N-scaling (new default, old behaviour via argument); NEWS + docs + tests
  (forced-failure mocks exist as templates in test-select-ll.R).
- Update CONSOLIDATED_RESULTS defaults table with evidence.
