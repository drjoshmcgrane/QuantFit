# STUDY PLAN 2: Formal triangulation of the three routes (drafted 2026-08-13)

## The construction (agreed with user)
Current quant_fit deliberately refuses to combine the LC / CC / Kara p's
(dependent, different nulls). Formalization via the SHARED null "data are
Rasch":
1. Fit Rasch once; B parametric-bootstrap replicates (reuse each route's
   mask/banding machinery; ONE shared replicate set).
2. Per replicate compute the TRIPLE: LC edge statistic, CC violation rate,
   Kara global KL -> joint null distribution WITH dependence (inherited from
   shared replicates, not assumed).
3. Disjunctive test: studentized max-statistic / min-p over the triple
   (identical construction to the poset max-T; one calibrated (1-alpha)
   quantile absorbs multiplicity + dependence). Exact familywise size for
   "reject interval scalability if ANY route objects".
4. Conjunctive CERTIFICATION (all routes consent) is not auto-calibrated:
   measure its operating characteristics by simulation over the truth space
   {Rasch, LCR, each ordinal, contaminated/multidim} -> P(certificate |
   truth) table. Severity reading: certificate informative because each
   route has demonstrated power to object.
5. Divergence-as-diagnosis lookup: patterned disagreement localizes the
   violated assumption (canonical: LC=quant-discrete + CC consent + Kara
   reject => latent SHAPE, not additivity - the LCR cross-route study).
   Build the pattern table from existing cross-route studies + the new sims.

## Caveats to state
- The shared null is Rasch specifically: "reject" = not-Rasch-scalable.
  Discrete quantity (LCR) needs a second certification target (LCR-null
  shared bootstrap) or explicit framing.
- Kara route cost dominates: reuse the shared replicates (B=99 KaraChecks
  runs ~ current omni_B anyway); banding per replicate as in quant_fit.
- MAR-by-ability limitation applies to all three routes' nulls
  (missing_data_finding.md): scope claims accordingly.

## Implementation sketch
- quant_fit(combine = "joint"): shared replicate loop, triple recording,
  max-T layer (mirror .poset_refine's decide machinery), report joint p +
  per-route contributions + concordance pattern + its calibrated
  interpretation.
- Validation grid: size on Rasch truth (target alpha); certificate rates
  across truths; K>=50 for the size rows. Provenance canaries as always.
- Deliverable: triangulation_finding.md + consolidated section update +
  NEWS; reviewer will demand resolution proofs (B vs quantile depth - the
  max-T argument transfers verbatim).

## CONTEXT SNAPSHOT (2026-08-13, for post-compaction recovery)
- Main TA grid: 627/628 synced (1 clean LCR cell in flight, uni Mac);
  2 documented crashers TA376/TA552 (known_crashers.md, segfault at
  LCR-vs-DM bridge, engine/seed-independent; reproducer committed).
- Full-grid analysis waiter: /tmp/waitrun.sh -> /tmp/fullgrid_analysis.txt
  (fires at 628; misspecified table + demotion autopsy).
- nI=48 block: running on laptop (singleton scheduler /tmp/singleton48.sh,
  3 slots, ordinal-first IIO->MON->UN->DM->LCR->RM), results ->
  QF_reval/extend/results48/ (bundle, synced); stamped build 5d3f8fc;
  ~30-60 min/ordinal cell. Quant 90 cells (25-class bridge) = explicit
  DECISION POINT when ordinal 90 done: run/sample/exclude-by-cost.
- Rater analysis committed: tid_rater_finding.md (five-way tables incl.
  per-nI split), tid_rater_merged.csv (full 1080).
- Key headline (clean 537): hybrid 92% / lattice 83% / RPD 75% / DTI 68% /
  minBIC 65%; scale 97/96/83; false-quant 0.0%/0.8% vs 24% all
  non-calibrated; demotions 7.9%/11.2%.
- Bias framing agreed: sims = validation, TA = development-adjacent
  confirmation, raters = sealed out-of-sample. To add to consolidated
  honesty section.
- Uni Mac bundle path: ~/Library/Mobile Documents/com~apple~CloudDocs/
  QF_reval/extend; laptop full tid_data local.
