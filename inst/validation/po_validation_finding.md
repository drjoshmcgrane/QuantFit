# Partial order: generators, calibrated tests, and their validation

Supersedes the v1 (count/threshold-era) version of this file. The annotation
scheme it validated was replaced - at the user's challenge that partial order
must be **tested the way the kinds of ordinality are tested** - by calibrated
hypothesis tests (commit 6ee0046); raw v1 results remain in
`po_validation_results/`. This file documents the redesigned machinery and its
validation grid (`po_validation.R` v2, 192 datasets, B = 49, N = 1500, K = 8).

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
