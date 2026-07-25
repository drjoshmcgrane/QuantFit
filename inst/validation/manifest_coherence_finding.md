# Is the "manifest" route internally consistent? No — and the fix is to go
# fully latent, not fully manifest

## The objection

The manifest route mixes three different kinds of test on three different
structures:

| step | structure actually tested | framework |
|------|---------------------------|-----------|
| IIO axis | rest-score groups × items (raw data) | manifest statistic, latent null (fitted DM) |
| **MON axis** | **fitted UN latent classes × items** | **fully model-based** |
| double cancellation | sum-score bands × items | manifest axiom (additive conjoint measurement) |
| DIP axis | score distribution | manifest |

Three consequences:

1. **The stated reason for dropping single cancellation does not follow.** The
   code justified testing double cancellation ALONE because "ordinality is
   already established (DM reached via the IIO+MON axes)". Single cancellation
   (independence) on a two-factor array is essentially double monotonicity — but
   the 2×2 established DM on the **latent class × item** array, while double
   cancellation is then evaluated on the **manifest sum-score band × item**
   array. Ordinality in one array does not license assuming it in the other.
   (The separate *statistical* argument — avoiding Holm dilution — is unaffected.)
2. **"Manifest" is a misnomer**: the MON axis fits an unconstrained latent class
   model and inspects its class × item table.
3. Testing additivity on sum-score bands leans mildly on the sum score being a
   meaningful aggregate — related to (though weaker than) the reason
   theta-ordering was rejected as a person-ordering option.

## The obvious fix was tried and FAILED

Proposal: make the MON axis manifest too, using **Mokken monotone homogeneity**
(MH): for each item, P(X_j = 1 | rest score) must be non-decreasing. Implemented
as `.manifest_mh_stat` / `.manifest_mh_holds` — the exact counterpart of the IIO
axis (same rest-score construction, same person weighting, same parametric DM
null, same bootstrap p-value, and **no `mon_eps`**). Available as
`mon_method = "mh"`.

Result (10 reps/condition, N=1500, hold-rates; HOLDS rows want ~1, VIOL want ~0):

| condition            | **MH (manifest)** | **eps (latent, default)** |
|----------------------|-------------------|---------------------------|
| DM holds             | 1.0 | 1.0 |
| MON holds            | 1.0 | 1.0 |
| RM holds (J=6)       | 1.0 | 1.0 |
| **LCR holds**        | **0.8** (false violations) | 1.0 |
| DM sep=0.6 holds     | 1.0 | 0.1 |
| **IIO violated J=8** | **0.6 (40% power)** | **0.1 (90% power)** |
| **IIO violated J=12**| **0.9 (10% power)** | **0.1 (90% power)** |
| UN violated          | 0.0 | 0.0 |

The manifest MH axis is **strictly worse where it matters**: near-blind to IIO
data (10–40% power vs 90%), and it additionally mis-fires on LCR (0.8). Its only
advantage is the artificial low-separation corner.

## Why it fails — the substantive point

**Mokken's MH and TI&D's MON are not the same property.**

- TI&D's **MON** is a constraint on the *latent* class × item table: there exists
  a class ordering making every item's probability monotone across classes.
- Mokken's **MH** is the *observable consequence* of a monotone latent model, and
  it is strictly weaker. The rest score orders persons by their **actual
  responses**, so when latent classes cross, persons from those classes are
  blended across rest-score groups and the conditional regressions come out
  monotone anyway.

So latent class-monotonicity violations do not reliably surface as manifest
rest-score non-monotonicity. There is **no faithful manifest proxy for MON**, and
therefore the six-model hierarchy (UN/MON/IIO/DM/LCR/RM) cannot be tested
manifestly with fidelity. IIO is different: item-order invariance *does* have a
faithful manifest counterpart (the crossing statistic), which is why that axis
works model-free.

## Conclusion: the coherent design is FULLY LATENT

The inconsistency is real, but the resolution runs the opposite way to the
obvious fix:

- TI&D's hierarchy is **intrinsically latent**. The latent MON axis is not the
  anomaly.
- The anomaly is **importing manifest double cancellation** — an axiom about an
  observable conjoint array — as the DM→quant step of a latent-model hierarchy.
- The **hybrid** (`dm_quant = "lr"`: 2×2 ordinal layer + likelihood-ratio
  LCR-vs-DM edge) keeps everything in one framework and is therefore the
  *conceptually* correct design, not merely the empirically better one.

This gives a principled reading of the empirical result that the hybrid fixes
the manifest route's DM→quant leak (DM 7/10 → 10/10): the pure-manifest route was
not just underpowered at that seam, it was testing a different object.

`mon_method = "mh"` is retained (documented as inferior) so this is not retried,
and because it is the correct axis if one ever wants a genuinely manifest,
Mokken-style scale analysis rather than TI&D's latent hierarchy.

## Status of `mon_eps`

Not removable. Two attempts have now failed for different reasons: a calibrated
monotone-cone null (`mon_method = "null"`, see `mon_calibration_frontier.md`)
trades away IIO power along a size/power frontier; the manifest MH axis tests a
different property. `mon_eps = 0.01` with the resample-q05 rule remains the
default, and remains the most fragile parameter in the package.
