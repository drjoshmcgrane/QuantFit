# Monte-Carlo resolution of the IIO axis: borderline decisions were not
# reproducible

## Discovery

Regenerating the 9 TI&D datasets whose IIO-axis decision flipped under the
v0.3.4 null fix produced verdicts that did not match the scan's prediction
(three RM datasets predicted to flip reproduced their old verdict). Cause:
the axis p-value is a bootstrap proportion with standard error
~sqrt(p(1-p)/B) = 0.031 at the shipped B = 49, while the decision threshold
is 0.05. Borderline datasets are therefore decided by the seed.

## Measurement (15 datasets x B in {49, 199} x 5 seeds)

Decision agreement across seeds on IDENTICAL data:

| datasets | B = 49 | B = 199 |
|----------|--------|---------|
| borderline (the 9 flips) | **2/9 unanimous** | 6/9 unanimous |
| clearly decided (3 controls) | 3/3 | 3/3 |

p-value spread (max - min over 5 seeds) on borderline datasets: median
0.080 (max 0.240) at B = 49; median 0.045 (max 0.070) at B = 199. Matches
the binomial prediction; reaching a spread comfortably inside the 0.05
boundary needs B ~ 999.

## Fix: adaptive precision

`.manifest_iio_holds()` gains `B_refine` (default 999) and `refine_band`
(default 0.02-0.15). The axis draws B replicates as before; if the p-value
lands in the indeterminate band, the SAME null is extended to `B_refine`
draws and the decision is taken on the pooled distribution. Decisive
datasets pay nothing; borderline ones become reproducible. Exposed as
`iio_B_refine` / `iio_refine_band` on `select_model_hybrid()`.

## Honest consequence for the v0.3.4 attribution

The nine "flips" attributed to the class-count fix are partly Monte-Carlo
noise: at nI <= 24 the fix's footprint is inside the resolution of the
shipped B. The fix's unambiguous value is the large-bias regime
(95% -> 10% false rejection at J = 48, 15% -> 0% at J = 24); its effect on
individual borderline verdicts at shorter tests cannot be separated from
seed noise without the refinement above. All tables quoting hybrid
verdicts should be produced with adaptive precision enabled.
