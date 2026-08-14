# The nI = 48 block: where perception hits zero and calibration peaks

TI&D's longest test length (48 items, N = 5000, clean conditions only -
the nI=48 grid has no multidimensional contamination). Ordinal block
complete: 90/90 (UN, MON, IIO x 30). Quant block (DM/LCR/RM x 30) runs
separately; each of those pays the 25-class Lindsay bridge.

## Ordinal block (90 datasets)

| truth | n | DTI | RPD | min-BIC | **hybrid** | lattice |
|-------|---|-----|-----|---------|------------|---------|
| UN | 30 | 100% | 100% | 100% | 100% | 100% |
| MON | 30 | 93% | 87% | 87% | 100% | 100% |
| **IIO** | 30 | **0%** | **0%** | **0%** | **83%** | 37% |
| overall | 90 | 64% | 62% | 62% | **94%** | 79% |
| scale-type | | 73% | 70% | 71% | **100%** | 98% |
| false-quant (of 90 non-quant) | | 24.4% | 25.6% | 24.4% | **0.0%** | 1.1% |

**Both raters and min-BIC score exactly 0/30 on IIO at 48 items** - the
worst human IIO performance at any length (nI 6/12/24 gave 7/7/0% and
17/27/10%), against the hybrid's 83%, its best or near-best. Forty-eight
items of evidence add nothing to graphical judgement about invariant item
ordering, while calibrated testing converts the extra information into
recovery. The false-quantitative rate for perception stays at ~25%: on
long tests, IIO data look MORE Rasch-like, not less.

Runtime note (used for the quant-block cost decision): IIO cells take a
median 162 min because they route to DM and pay the 25-class bridge;
UN/MON resolve in ~2 min in the 2x2 without touching it.

## Reading

This is the cleanest statement of the instrument's case. At the length
where an expert has the most information to look at, the eye extracts
nothing about the hierarchy's hardest distinction, and the penalty-based
criterion is structurally blind to it (equal effective parameter counts).
The property-based 2x2 with calibrated axes recovers 5 of every 6 IIO
datasets, and neither selector issues a single false quantitative
certificate.
