# The nominal/ordinal boundary is too coarse: partial orders have no rung

## The objection (user)

The six-model hierarchy jumps straight from **nominal** (UN) to **ordinal**
(MON/IIO/DM). MON asks whether a **total order on classes** exists. Its negation
therefore lumps together two very different situations:

- **antichain** - no class dominates any other on all items -> genuinely nominal
- **partial order** - e.g. c1 dominates c2 and c3 on every item, but c2 and c3
  cross each other -> *not* a total order, but emphatically **not nominal**

A partial order is strictly more measurement structure than nominal, and neither
TI&D's model set nor QuantFit has anywhere to put it. This is a limitation of the
framework, inherited, not a bug in the implementation.

## Instrument

The **class-dominance poset**, computed from the already-fitted UN class x item
table (no new fit, no new bootstrap). Class a dominates class b iff the per-item
mean downward violation is negligible:

    mean_j max(0, P[j,b] - P[j,a]) <= eps

which reuses the calibrated `mon_eps` rather than inventing a threshold. For C
classes there are C(C-1)/2 pairs:

| all pairs comparable | none comparable | in between |
|----------------------|-----------------|------------|
| **chain** (= MON holds) | **antichain** (nominal) | **partial order** |

### Validation on simulated data (8 reps each, J=12, N=1500, eps=0.01)

| true model | chain | partial | antichain |
|------------|-------|---------|-----------|
| DM         | 8     |         |           |
| MON        | 8     |         |           |
| **IIO**    | 1     | **6**   | 1         |
| UN         |       |         | **8**     |

Total-order models give chains, unstructured data gives antichains, and IIO -
which violates MON by construction - lands in the **partial** category. The
instrument behaves as theory predicts.

## Result on REAL TI&D data (full N = 5000)

Restricted to the datasets **both** selectors classify as nominal, so this cannot
be dismissed as one method's weakness. (The hybrid recovers only **1 of 27** of
the lattice's UN verdicts - both bury them equally, so this is a framework
limitation, not a lattice one.)

| poset shape           | n  | true generating model |
|-----------------------|----|-----------------------|
| antichain (correct)   | 19 | **all 19 true UN**    |
| **partial order**     | **7 (27%)** | **4 IIO, 1 MON**, 2 UN |

**Every one of the 5 misclassified datasets shows partial-order structure; none
is an antichain.** 19 of the 21 correctly-nominal datasets *are* antichains. So
on this data the poset separates "correctly nominal" from "wrongly nominal"
almost perfectly.

## Why it matters

A partial-order rung would not merely relabel output - it would **flag exactly
the cases both current selectors get wrong**, at a cost of 2 false positives
among 21 true-UN datasets. For a package whose purpose is defensible scale-type
claims, reporting "partial order" instead of "nominal" is a strictly better
answer when it is true, and it is true surprisingly often.

It is also cheap: the poset is free from the existing UN fit, and would calibrate
like the IIO axis (parametric DM null, bootstrap p-value) rather than needing
another hand-tuned threshold.

## Caveats

- n = 26 with only 5 misclassified cases, so "5/5" is **suggestive, not
  established**.
- All of these are **misspecified (non-clean) conditions** - slope variation,
  correlated dimensions - which is where partial structure is most expected. On
  clean conditions the lattice never wrongly said UN at all.
- This affects the **nominal/ordinal boundary only**. It does not touch the
  quantitativeness claim, which is the package's headline purpose.
- Adding the rung is an **extension to TI&D's model set**, i.e. new research,
  not a correction of an implementation error.

## Files

`poset_diag.R` (diagnostic), `tid_hybrid_full.R` (the real-data hybrid run),
`tid_realdata_results/`.
