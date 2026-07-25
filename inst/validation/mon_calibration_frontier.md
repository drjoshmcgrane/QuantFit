# Can the MON axis's `mon_eps` threshold be replaced by a calibrated null?

## Motivation

A self-audit flagged a design asymmetry in the manifest 2x2: the IIO axis is
calibrated properly (parametric bootstrap null → p-value → compare to alpha,
so its arbitrary statistic scale cancels), whereas the MON axis compares a raw
resample q05 to a **hand-tuned scalar `mon_eps`**. That magic number had already
caused one catastrophic regression this session (`0.04` collapsed IIO recovery
to ~20%; see `manifest_mon_calibration_finding.md`). The obvious fix: give MON a
p-value too, and delete the threshold.

**This was tested and it does not work.** Documented here so it is not retried.

## What was tested

A parametric null in which MON holds EXACTLY: fit UN, find the best class
ordering, project the fitted class×item probability table onto the monotone
cone, simulate from the projection, refit UN, recompute the statistic →
`p = (1 + #{null >= obs}) / (B+1)`, `holds = p > alpha`. Two projections:

- **pava** — weighted isotonic (pool-adjacent-violators) per item. Verified
  against `stats::isoreg`. POOLS violating classes to equal values.
- **sort** — sort each item's class probabilities along the ordered classes.
  Enforces monotonicity while preserving each item's spread exactly, so the null
  retains the observed class separation.

(A `fit_mon` parametric null was tested earlier and is worse still: constrained
EM ABSORBS an IIO crossing, inflating the null; IIO power 6-19%.)

## Result — a frontier, not a fix (10 reps/condition, N=1500)

Hold-rate; HOLDS rows want ~1.0, VIOLATED rows want ~0.0:

| condition            | pava | sort | **eps (default)** |
|----------------------|------|------|-------------------|
| DM holds             | 1.0  | 1.0  | 1.0 |
| MON holds            | 1.0  | 0.9  | 1.0 |
| LCR holds            | 1.0  | 1.0  | 1.0 |
| RM holds (J=6)       | 1.0  | 1.0  | 1.0 |
| **DM sep=0.6 holds** | 0.9  | 0.9  | **0.3** |
| **IIO violated J=8** | 1.0 (0% power) | 0.7 (30% power) | **0.2 (80% power)** |
| **IIO violated J=12**| 0.9 (10% power) | 0.6 (40% power) | **0.0 (100% power)** |
| UN violated          | 0.0  | 0.0  | 0.0 |

Reading: the three methods sit on ONE power/size frontier and **none dominates**.
Moving from `eps` → `sort` → `pava` buys hold-rate in the artificially
low-separation corner and pays for it, roughly one-for-one, in IIO power.

Confirmed functionally on a single IIO dataset: `eps` → **IIO** (correct);
`null` → **DM** (missed, MON p = 0.95).

## Decision: `eps` remains the default

1. **IIO detection is the entire purpose of the manifest route** (it is where it
   beats the LR lattice). A calibration that drops IIO power to 30-40% forfeits
   the reason the approach exists.
2. The null's only advantage is at class separation < ~1 logit — a regime
   `simulate_responses` does not produce and which is absent from the TI&D data.
   It is also a regime where weak-IIO and collapsed-DM are genuinely
   indistinguishable (an identifiability limit, not a calibration target).
3. So the threshold is **doing real work**, not papering over a missing
   calibration. The audit's "magic number" criticism is answered: the number is
   not removable without losing the capability.

`mon_method = "null"` is exposed (sort projection, the better-balanced variant)
for the case where class separation is known to be low and false IIO calls cost
more than missed ones. Default `mon_method = "eps"`.

## Residual caveat (honest)

`mon_eps = 0.01` is still a hand-tuned constant validated at N=1500 on
`simulate_responses` plus spot checks (J=6/8/12, RM/LCR, TI&D). It is the most
fragile parameter in the manifest route and the one most likely to need
re-derivation for a materially different item pool. The frontier above is the
evidence for its value, not a proof of optimality.
