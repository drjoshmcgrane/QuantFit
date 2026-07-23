# The IIO recovery ceiling is largely procedural: direct property tests roughly double it

## Summary

Invariant Item Ordering (IIO) has been the hardest of the six latent-structure
models to recover for every method attempted: Torres Irribarra & Diakow's
expert human graphical classification recovered it at **41%**, minimum-BIC at
**0%**, and our likelihood-ratio "lattice" (adjacent-edge) selector at
**~42%**. This near-universal ~40% floor has been treated as an intrinsic
identifiability limit.

It is largely **not**. Testing the two defining structural PROPERTIES directly
against the data — rather than comparing two near-equivalent constrained-model
fits — lifts IIO recovery to **87%** overall and **100% / 93%** at J = 12 / 24,
with the ordinal/nominal layer reaching **94.4% exact** classification. The
residual is concentrated at J = 6 (a short-test regime TI&D never tested) and a
small irreducible core of IIO data that genuinely satisfies monotonicity by
chance.

## Why the LR-edge approach fails on IIO

The lattice decides DM vs IIO with a parametric-bootstrap likelihood-ratio
test. On doubly-monotone data DM and IIO fit **identically** (observed LR ~0.06
= optimisation noise). The bootstrap null — simulated under the fitted DM — is
even "cleaner" (perfectly double-monotone), so null LRs cluster at ~0. The tiny
positive observed LR then exceeds the whole degenerate null, giving p = 1/(B+1)
and a **spurious rejection of a negligible effect**. Symmetrically, on true-IIO
data the misspecified DM fit inflates the null so the observed LR cannot exceed
it (false retention). The edge test has near-zero power in both tails.

Measured DM-vs-IIO null 95th percentile: **0.04** on true-DM data vs **12-18**
on true-IIO data — the test is operating in a degenerate regime on exactly the
data it must classify.

## The manifest 2x2 (TI&D's own logic, automated and calibrated)

TI&D classified by hand using **constraint presence** (their Figure 7/8 encode
each fit with a 2-digit dummy for MON-present, IIO-present), not model
comparison. Automated:

- **IIO axis** — invariant item ordering. Order persons by rest-score; do item
  response functions cross (an overall-harder item become easier within a
  group)? Model-free crossing magnitude, calibrated by a parametric null in
  which item ordering is imposed (simulate under fitted DM).
- **MON axis** — class/person monotonicity. Order the fitted unconstrained
  latent classes by mean success; how much does any item's class-probability
  DECREASE across the ordered classes? Calibrated by **data resampling** (the
  statistic's own sampling distribution). NOTE: a parametric simulate-and-refit
  null does NOT work here — refitting UN to the null data reintroduces the very
  unconstrained-fit noise the statistic avoids, collapsing power (it drove IIO
  recovery to 1/6 in testing). Data resampling recovers 23/24.

  IIO holds & MON holds -> DM ;  IIO holds & MON violated -> IIO
  IIO violated & MON holds -> MON ;  both violated -> UN

## Results (180 datasets: 4 ordinal/nominal models x J{6,12,24} x 15 reps, N=1500)

Confusion (manifest 2x2):

           selected
    true    UN  MON  IIO  DM
    UN      44    .    1   .     98%
    MON      2   43    .   .     96%
    IIO      1    .   39   5     87%
    DM       .    .    1  44     98%

    exact 94.4%

IIO recovery by J:  J6 10/15 (67%)   J12 15/15 (100%)   J24 14/15 (93%)

Comparison on IIO:  manifest 2x2 ~87%  |  LR-edge lattice ~42%  |
TI&D humans 41%  |  min-BIC 0%.

## Cost

180 datasets in ~2 minutes (model-free IIO axis + one UN refit with cheap data
resampling for the MON axis; no LCR bootstrap, no lattice climb) versus hours
for the LR-edge lattice.

## Implementation

`select_model_manifest()` (R/select_manifest.R) - a SEPARATE selector,
independent of `select_model_ll()`. Ordinal/nominal layer by the manifest 2x2;
when DM is reached it enters the DM -> LCR -> RM quantitative sequence reusing
the same calibrated machinery as the lattice. The ordinal layer is validated
(above); the quantitative sequence reuses `fit_rm`/`rm_vs_lcr_test` and
inherits their known GC-segfault and cost characteristics (orthogonal to this
finding; fork-isolate or fix separately).

## Caveats

- Single N (1500), single generator per model; the effect is large and
  J-consistent but wants replication across N and item designs.
- Thresholds (IIO p>0.05; MON eps=0.03) are reasonable but not exhaustively
  tuned; the data-resampling calibration is the substantive part.
- ~40% of the IIO floor is genuinely irreducible (IIO data that satisfies
  monotonicity by chance) - but far more of it than assumed was procedural.

## Claim for the paper

The IIO recovery ceiling that has stood since Torres Irribarra & Diakow (and
which our own LR-edge selector reproduced) is largely an artifact of
model-comparison testing in a degenerate regime. Testing the invariant-item-
ordering and class-monotonicity properties directly - which is what TI&D did by
eye - roughly doubles IIO recovery (41% -> 87%) and is far cheaper.
