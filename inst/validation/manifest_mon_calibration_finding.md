# MON-axis separation-robust calibration (and a correct model×separation generator)

> **CORRECTION (post-hoc).** The `mon_eps = 0.04` recommended below was
> calibrated ONLY on the artificial `gen_correct` generator at low class
> separation (sep=0.6). On the actual deployment generator
> (`simulate_responses`, well-separated classes) a real IIO class-monotonicity
> violation gives a per-item q05 of only ~0.01-0.06, so **eps=0.04 sits above it
> and silently classifies IIO data as DM** — collapsing IIO recovery from the
> validated ~87-100% (at J=12) to ~20% in the head-to-head. IIO detection is the
> entire purpose of the manifest route, so this was a serious regression.
> **The default is now `mon_eps = 0.01`**, which restores IIO recovery while
> keeping DM/MON/UN size at 100% on `simulate_responses` (verified: eps=0.01
> gives sr_IIO power 100%, sr_DM/MON size 100%, sr_UN power 100%). The low-sep
> `gen_correct` "size" the 0.04 value protected (sep<~1 logit) is a regime the
> real generator never produces AND where weak-IIO and collapsed-DM are
> genuinely indistinguishable — an identifiability limit, not a tuning target.
> Lesson: calibrate the threshold on the DEPLOYMENT generator, not an artificial
> stress corner. The perm-min statistic and resample-q05 machinery below are
> retained and correct; only the threshold was mis-set.

## The bug that motivated this

The earlier `gen_sep` generator (in `manifest_iio_mon_robustness.R`) built IIO/MON
data by *permuting* a doubly-monotone logit table — permuting classes per item
for IIO, items per class for MON. That breaks **both** invariances at once, so
the IIO/MON rows were artifacts (catastrophic 32%/56% "size"). Only the
unpermuted DM row was trustworthy — and it revealed a **real** MON-axis size
failure: on clean DM data the MON axis fired *false violations* at a rate of
34% overall, rising to **97% at low class separation** (sep=0.6 logits).

## Root cause of the MON size failure

The MON statistic ordered the fitted-UN classes by their mean success and summed
each item's downward class-probability movement. At low separation the classes
are near-flat, so ordering-by-mean flips close classes on sampling noise →
spurious downward movement. Two compounding scale errors made it catastrophic:

1. The statistic was **summed over J items** but compared to a fixed `eps=0.03`,
   i.e. an effective per-item threshold of `0.03/J ≈ 0.0025` — absurdly tight.
2. Mean-ordering injected noise the *definition* of MON does not require.

## The fix

MON holds iff **some** class order makes every item monotone, so the statistic
is now the **minimum total downward movement over all class orderings**
(exact for C≤5; mean-ordering fallback above that), **normalised per item**.
Calibration stays data-resampling but keys on the mechanism that actually
separates the classes: a genuine crossing is *structural* and survives
resampling (its q05 lower bound stays high); ordering noise on near-flat classes
*shakes out* (q05 collapses toward 0). Property holds iff the resample q05 is at
or below a per-item tolerance `mon_eps` (default **0.04**, N-scaled `√(1500/N)`).

The q05 separation is clean at every separation, including the corner where the
point estimates overlap (per-item stat, J=12, N=1500, mean of 6 reps):

| sep | DM/MON q05 (holds) | IIO/UN q05 (violated) |
|-----|--------------------|-----------------------|
| 0.6 | 0.021, 0.022       | 0.063, 0.108          |
| 1.2 | 0.000, 0.000       | 0.061, 0.203          |
| 2.0 | 0.000, 0.000       | 0.085, 0.295          |

A parametric **monotone null** was rejected: it fixes size but `fit_mon`
absorbs the IIO crossing, collapsing IIO-detection power to 6–19%.

## Correct generator

Each model is now a class×item logit table satisfying **exactly** its
invariances (verified on the population probabilities):

- **DM**  `L = θ_c − β_j`                          (classes ordered, items ordered)
- **IIO** `L = −a_c·β_j + b_c`, `a_c>0` distinct   (item order fixed; class curves cross → MON violated)
- **MON** `L = θ_c·a_j + b_j`, `a_j>0` distinct    (class order fixed; item curves cross → IIO violated)
- **UN**  both crossing

`sep` scales the signal (class spacing / slope spread). Property check confirms
DM = mon&iio, IIO = iio-not-mon, MON = mon-not-iio, UN = neither.

## Post-fix validation (432 datasets, model × J∈{6,12,24} × N∈{750,1500,3000} × sep∈{0.6,1.2,2.0})

**MON axis** — hold-rate (want DM/MON high = size ok; IIO/UN low = power):

| model | sep0.6 | sep1.2 | sep2.0 |
|-------|--------|--------|--------|
| DM    | 0.86   | 1.00   | 1.00   |
| MON   | 0.75   | 1.00   | 1.00   |
| IIO   | 0.11   | 0.06   | 0.03   |
| UN    | 0.00   | 0.00   | 0.03   |

Size 86–100% (only MON@sep=0.6 dips to 75% — classes 0.6 logits apart *and*
crossing items, the identifiability corner). Power 83–100%. Replaces the old
rule's **0%** size at sep=0.6.

**IIO axis** — size perfect (DM/IIO hold 92–100%, no false violations). Power to
detect IIO-violation is separation-dependent: strong at sep=2.0 (MON→violated
89%, UN 97%) but weak at sep=0.6, because the MON generator's item-slope crossing
is genuinely tiny there. This is a power (not size) limit, concentrated at weak
signal — consistent with the standing IIO identifiability floor.

## Files

- `R/select_manifest.R`: `.class_orderings`, perm-min per-item `.manifest_mon_stat`,
  `.manifest_mon_holds` (unchanged mechanism, new statistic), `mon_eps` default 0.04.
- `inst/validation/manifest_iio_mon_robustness.R`: correct `gen_sep`, `eps=0.04`.
