# Partial order: generators, calibrated tests, and their validation

Supersedes the v1 (count/threshold-era) version of this file. The annotation
scheme it validated was replaced - at the user's challenge that partial order
must be **tested the way the kinds of ordinality are tested** - by calibrated
hypothesis tests (commit 6ee0046); raw v1 results remain in
`po_validation_results/`. **WITHDRAWN AND REDESIGNED AGAIN (2026-07-29, external review round 4):**
the v2 "continuous asymmetry" statistic algebraically collapsed to the
absolute difference of side-average probabilities (|mean(d+)| - |mean(d-)| =
|mean(d)|) - it measured average-profile alignment, not dominance - and its
permutation null tested exchangeability, not the antichain hypothesis. A
constructed crossing-antichain with unequal averages rejected at p = 0.006.
ALL v2 size/power claims are therefore withdrawn (raw files remain in
po_validation2_results/). The v2 PO_ITEMS generator also silently delivered
crossing margins far below request (achieving ~0.003-0.03 against a 0.10
target), so its populations carried more comparable pairs than specified -
the v2 item-side "power" numbers validated detection of unrequested
structure.

## The redesigned test (v3, current)

Per DIRECTED pair, the violation mass V(x >= y) = mean over the opposing side
of pmax(0, P_y - P_x) is zero iff x weakly dominates y. ONE parametric
bootstrap from the fitted UN table (B refits, class labels aligned by
best-permutation profile matching) gives the sampling distribution of every
directed V. A pair is COMPARABLE iff, at Bonferroni-corrected one-sided
levels alpha/(2m), the upper bound of one direction's violation is <= eps
(the axis tolerance) while the lower bound of the reverse exceeds eps.
"partial" = >= 1 demonstrated pair; "antichain" = none demonstrated (a
failure to demonstrate, not a certified absence). No permutation null;
conservative under ANY antichain including crossing profiles with unequal
averages. Failures: b_eff reported, hard refusal below max(20, B/2).

PO_ITEMS is now constructive for graded posets (band per level; one class
carries the within-band grid in order, another fully reversed, so every
within-band pair provably crosses), audits its achieved margin, and ERRORS
when the requested margin exceeds the geometric ceiling (~slope x band_width
/ (items_per_level - 1)): the mirror truth is only demonstrable at short
tests or C > 3, and the validation design reflects that (C = 4 arms at
feasible margins; the C = 3, J = 24 case is included as a documented
below-tolerance boundary).

## v4 CANONICAL grid results (458 datasets, stamped build 1454231, per-row
## method/SHA/config canaries; raw in po_validation4_results/)

Inference: STUDENTIZED SIMULTANEOUS MAX-STATISTIC bounds (one calibrated
(1 - alpha/2) quantile of the studentized max deviation per family;
per-direction bootstrap scale; poset bootstrap floor 99, b_eff = 99
everywhere). Zero cell failures; zero non-transitive relations observed.

**Familywise size (class side, alpha = 0.05):**

| antichain class | n | false partials |
|---|---|---|
| XANTI nI=8 (deep crossings, unequal means) | 100 | **0** (CI 0-3.6%) |
| XANTI nI=24 | 100 | **0** (CI 0-3.6%) |
| **NEARANTI nI=8 (min crossing depth 0.0150 vs eps 0.01 - the boundary)** | 50 | **0** (CI 0-7.1%) |
| UN control | 24 | 1 (~4%) |

The construction is conservative at and beyond the tolerance boundary - no
false demonstration in 250 antichain datasets.

**Class-side power: 8/8 in every cell** - PO at margins {.05,.10,.15} x nI
{6,12,24}, PO_INV (accept-best, achieved ~0.115-0.125) at all nI, and the
POLYTOMOUS arm (PO 4-category, categorical bootstrap) 8/8 with 8/8 correct V
types. V-type recovery 8/8 at nI >= 12 (2-6/8 at nI = 6, point-estimate
readout).

**Item-side power: 8/8 at every accept-best margin** (C=3 J=8 ach .016;
C=4 J=8 .032; C=4 J=12 .018; C=4 J=24 .008) with demonstrated-pair counts
always BELOW the population count (e.g. 104.9 vs 182; 5.4 vs 17) -
under-demonstration, never over-demonstration; per-row pair identities in
the CSVs make the subset property verifiable.

## Deep-B stability (32 anchors, B = 99 vs 499; deepB_results/)

Every anchor validated (expected side present, b_eff checks, warnings
escalated). **Shape decisions: 32/32 identical.** 13/32 anchors show
pair-SET differences of 1-3 near-threshold pairs (predominantly item-side),
in the conservative direction at deeper B - the expected behaviour of a
top-order-statistic bound. Verdicts are depth-stable; individual marginal
pair demonstrations are not, which is why `poset_B = 499` is the documented
recommendation for final analyses.

## K=20 audit grid under max-T (audit_k20_maxT_results/, build 936a7ac*)

*936a7ac differs from 1454231 only in display/metadata code (quant_fit flag
propagation, print methods, refusal recording) - `git diff 1454231..936a7ac`
shows no statistical-path change.*

Six-model: hybrid 111/120 vs lattice 103/120 (McNemar 11:3, p = 0.057);
IIO 18/20 vs 9/20 - the headline comparison reproduced under single-run
provenance. **PO arm: hybrid 20/20, lattice 19/20** (18/20 correct V types
each, mean demonstrated pairs 1.9 vs truth 2). 280/280 rows, zero failures;
every lattice-IIO poset absence (9/20) is recorded as STRUCTURAL - the
lattice's BIC selects C = 2 on IIO-truth data, where no class poset exists
(the hybrid runs at fixed C = 3 and always reports).

Notable structural finding from the first (round-12) run, reproduced here:
the LATTICE's poset refinement at its BIC-selected class count refuses on a
large fraction of IIO-truth datasets (round-12: 9/20) - bootstrap UN refits
collapse at C = 5-6, the b_eff guard triggers, and the refusal is now
recorded per-row (poset_refused). The hybrid (fixed C = 3) never refuses.
Poset refinement at high BIC-selected C is fragile by nature of the refit
bootstrap; the fixed-C hybrid is the reliable carrier of the poset verdict.


## Historical: v3 (per-pair Bonferroni era, superseded)

po_validation3_results/)

Zero fit failures; b_eff = 49/49 in every dataset; ~10x cheaper than v2
(one shared bootstrap per dataset; median 3-13 s).

**Size.**

- **XANTI (the reviewer's counterexample class - crossing antichain with
  unequal class averages, the case that broke v2): 0/16 false partials**
  (nI = 8 and 24), all correctly routed to UN.
- UN control, class side: 1/24 (~4%, at alpha = 5% familywise per side).
- IIO control, class side: 2/17 flagged, both at >= 2 POPULATION comparable
  pairs (row-sorting creates real dominance); 0 flags among true-antichain
  draws.

**Class-side power (PO truths): 8/8 in EVERY cell** - all margins
{.05,.10,.15} x all nI {6,12,24}, and PO_INV at all nI. The v2 nI = 6 power
gap (~50-60%) is GONE: demonstrating the two real dominance pairs of a V is
much easier than the collapsed statistic's task. V-type recovery: 8/8 at
nI >= 12, 3-6/8 at nI = 6 (the type readout is point-estimate based).
Routing: PO -> UN 72/72; PO_INV -> IIO 23/24.

**Item-side power (PO_ITEMS at feasible margins): 8/8 in every arm** -
C=3 J=8 (achieved margin .016), C=4 J=8 (.032), C=4 J=12 (.018), and even
the below-tolerance boundary arm C=4 J=24 (.008): the demonstrated pairs are
the big cross-layer dominances, which don't depend on the within-layer
crossing margin. Critically, demonstrated-pair counts are always BELOW the
population count (e.g. 135.8 vs 182 at J=24; 8.9 vs 17 at C=4 J=8) -
Bonferroni under-demonstrates, never over-demonstrates: no incomparable pair
was falsely declared anywhere in the grid. Routing: PO_ITEMS -> MON 32/32.

**Population-faithfulness (item side of UN/MON controls).** Random TI&D
draws genuinely order many item pairs (~28% of pairs on UN truth, ~67% on
MON truth), and the test demonstrates a subset of exactly those - item-side
"partial" on such data is correct, not a false positive. Consequence for
interpretation: in the UN cell the item-side BINARY verdict will usually be
"partial" on any realistic data; the information is in the demonstrated pair
set, not the binary.

## Historical: v2 machinery description (superseded)



## The tests (both selectors, identical code)

Each side of the fitted UN table gets a three-level verdict - chain / partial /
antichain - bracketed by two calibrated tests: the existing axis test (MON for
the class side, IIO crossing for the item side) rules on "chain", and the NEW
**antichain test** rules on "any dominance at all". H0 = no side-consistent
dominance, simulated by permutation nulls of the fitted table (class side: each
item's class-probability vector permuted; item side: each class's item vector
permuted), simulate -> refit UN -> statistic -> p = (1 + #{null >= obs})/(B+1)
at the selector's alpha. The statistic is the **continuous dominance-asymmetry
mass** (sum over pairs of |mean one-way violation - mean other-way|):
label-invariant, eps-free, non-saturating. eps survives only in the
descriptive comparable-pair count and the C = 3 isomorphism type (single / V /
Lambda). Cell routing: UN cell -> both sides tested, IIO -> class side, MON ->
item side, DM -> none; verdicts are capped at "partial" (a chain would
contradict the routing axis) and never change `selected` or the scale claim.

## The generators (TI&D's own idiom throughout)

- **PO** (free items): per-item random linear extensions of a class poset over
  the U(-4,4) column draws. Comparable pairs dominate exactly on every item;
  incomparable pairs cross by >= `po_margin` both ways (rejection loop). The
  chain poset reproduces MON's column sort byte-identically; the antichain is
  UN. Routes to the UN cell.
- **PO invariant** (`item_order = "invariant"`): class lines over one shared
  sorted item spine, so IIO holds while class dominance follows the poset -
  routes to the IIO cell. Blind rejection sampling was computationally
  infeasible here (~3 min/dataset); the generator is now CONSTRUCTIVE: extreme
  slope separation for the crossing pair, intercepts placed in topological
  order (dominance exact by construction), crossing at the item-mass median,
  logits recentred on the incomparable classes; capped best-of-400 with an
  achieved-margin warning. 0.1 s/dataset; worst observed margin 0.112 vs the
  0.12 target.
- **PO_ITEMS** (`"layers2"`): theta chain (MON holds exactly) with per-class
  item values along random linear extensions of an ITEM poset - routes to the
  MON cell. Capped best-of-2000; at J = 8 the achieved crossing margin often
  falls short (warning fires), which is visible in the power table below.

## Validation grid results (K = 8 per cell)

**Routing** (six-model verdict): PO free -> UN **72/72**; PO invariant -> IIO
**23/24** (1 -> UN at nI = 12, both poset tests still flagged partial);
PO_ITEMS -> MON **23/24** (1 -> DM at nI = 8, the capped-search margin
shortfall). Controls: UN 24/24, MON 23/24, IIO 15/20 (the known IIO<->DM
boundary; unrelated to the poset layer).

**Class antichain test.**

- *Size*: on population antichains it holds at alpha - UN truth 1/24 rejections
  (4.2%); among the 4 IIO-truth datasets whose class side is a true population
  antichain, 1 rejection. Overall 2/28 (~7%) at alpha = 0.05.
- *Population-faithfulness*: 16/20 IIO-truth datasets contain GENUINE class
  dominance in the population (sorting rows makes order-statistic differences
  one-signed), and rejections track it: 10/12 rejections sit at >= 2 population
  comparable pairs.
- *Power* (PO truths, class side genuinely partial):

  | nI | PO m=.05 | PO m=.10 | PO m=.15 | PO_INV m=.12 |
  |----|----------|----------|----------|--------------|
  | 6  | 4/8      | 6/8      | 3/8      | 4/8          |
  | 12 | **8/8**  | **8/8**  | **8/8**  | **8/8**      |
  | 24 | **8/8**  | **8/8**  | **8/8**  | **8/8**      |

  Every detection also recovered the correct isomorphism type (V). The floor
  is **J-governed** (J observations per class pair): nI = 6 is marginal
  (~50-60%), nI >= 12 is saturated (mean p = 0.02 everywhere). Crossing margin
  barely matters once J >= 12.

**Item antichain test.**

- *Power* (PO_ITEMS): 1/7 at nI = 8, 5/8 at nI = 12, **8/8** at nI = 24. Two
  forces: C = 3 observations per item pair keeps per-pair signal weak, but the
  statistic accumulates over pairs, so J still buys power; the nI = 8 floor is
  also partly the generator's own achieved-margin shortfall.
- *Population-faithfulness*: on PO-free data (item side unconstrained but
  genuinely ~40-50% comparable in population) rejections concentrate where the
  population comparable share is high - 13/22 when > 50% vs 3/48 below. On MON
  truth (item side ~70% comparable in population) it rejects 18/23, correctly:
  those draws really do carry heavy item dominance.

## What this settles

1. **Partial order is now tested the same way ordinality is tested**: same
   permutation-calibration pattern, same alpha, same fitted-UN reference, one
   test per claim (axis test for "chain", antichain test for "any order"),
   with measured size and power - not an annotation with a heuristic
   threshold.
2. **Both selectors carry the tests** (`select_model_hybrid` and
   `select_model_ll`), so the audit grid scores PO symmetrically: correct =
   selected UN + class test rejects.
3. **Power floors are data-limited and symmetric with the axes' own floors**:
   the class test needs J >= ~12 items; the item test needs either more
   categories per pair or enough pairs. These are the same physics that give
   the MON/IIO axes their floors, now quantified for the poset layer.
4. The six-model vocabulary still cannot express any of this - every PO
   dataset files as UN/IIO/MON - which is the framework gap the refinement
   exists to close, demonstrated by construction.

## Caveats

- PO_ITEMS at J = 8 is under-margined by construction (capped search warns);
  its 1/7 power there is partly a generator limit, not purely a test limit.
- C = 3 throughout; richer posets at C >= 4 (diamond, N-shape) untested.
- Polytomous PO untested (the generators are dichotomous-idiom).
