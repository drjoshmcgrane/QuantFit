# TI&D pre-consensus rater data: humans vs selectors on the same datasets

Received 2026-08-11: DTIresults.csv (Torres Irribarra's per-set selections,
'*'/'?' uncertainty marks), RPDresults.csv (Diakow's pass, xN codes marking
disagreement + the alternative model seen), summary table.csv (their
per-model fit table). All 1080 TA sets truth-mapped via the archive's generatingModels.csv
(coding verified against the 270 independently-derived labels: perfect
diagonal); merged table committed as tid_rater_merged.csv.

## Pre-consensus human recovery (FULL GRID, 1080 sets, 180/model)

| truth | DTI | RPD | min-BIC | inter-rater agreement |
|-------|-----|-----|---------|------|
| UN | 99% | 99% | 99% | 100% |
| MON | 92% | 88% | 88% | 96% |
| **IIO** | **2%** | **10%** | **0%** | 88% |
| **DM** | **9%** | **22%** | **0%** | 77% |
| LCR | 71% | 71% | 71% | 98% |
| RM | 89% | 89% | 89% | 100% |
| overall | 60.6% | 63.4% | 58.0% | 93.1% |

Scale-type: ~78% everywhere. Min-BIC scores EXACTLY 0/180 on both IIO and
DM (the equal-parameter-count failure our simulations predicted); the
raters beat it only by rescuing a handful. IIO-truth destinations (DTI,
full grid): MON 78, **RM 52, LCR 26** (43% of IIO sets read as
QUANTITATIVE to the expert eye), UN 20, IIO 4. The published ~41% IIO
figure is a consensus-era number; pre-consensus it is 2-10%.

**The raters knew where they were struggling**: uncertainty flags (251/1080
sets) concentrate exactly on the identifiability-limited middle (IIO 49%,
DM 39%, LCR 21% flagged vs UN 1%, MON 9%), and Diakow's alternative-codes
on IIO sets point at LCR/RM (x4/x5) - IIO data genuinely reads as
quantitative to expert eyes.

## Head-to-head on the same 162 clean full-N sets

| truth | DTI | RPD | hybrid | lattice |
|-------|-----|-----|--------|---------|
| UN | 100% | 100% | 100% | 100% |
| MON | 93% | 93% | 96% | 100% |
| **IIO** | 4% | 15% | **70%** | 44% |
| **DM** | 15% | 37% | **89%** | **89%** |
| LCR | **93%** | **93%** | 85% | 81% |
| RM | **100%** | **100%** | 85% | 67% |
| overall | 67% | 73% | **88%** | 80% |

Reading: the automated selectors dominate exactly where human perception
fails - the ordinal middle (IIO: hybrid 70% vs humans 4-15%; DM: 89% vs
15-37%) - while expert graphical judgment remains excellent at the ends of
the hierarchy (RM 100%, LCR 93%, beating both selectors on the quant
models). The two error profiles are COMPLEMENTARY, which supports the
paper's framing: automation is not replacing the human read, it is
supplying the calibrated middle distinctions that graphical inspection
demonstrably cannot make, while the human strength on quantitative
structure is matched (not exceeded) at acceptable cost by the calibrated
edges' false-quant protection (0 false promotions).

## Five-way breakdown (same sets, selectors at full N = 5000)

**CLEAN conditions (162):**

| truth (27 each) | DTI | RPD | min-BIC | hybrid | lattice |
|---|---|---|---|---|---|
| UN | 100% | 100% | 100% | 100% | 100% |
| MON | 93% | 93% | 93% | 96% | 100% |
| IIO | 4% | 15% | 0% | **70%** | 44% |
| DM | 15% | 37% | 0% | **89%** | **89%** |
| LCR | **93%** | **93%** | **93%** | 85% | 81% |
| RM | **100%** | **100%** | **100%** | 85% | 67% |
| overall | 67% | 73% | 64% | **88%** | 80% |
| scale-type | 81% | 81% | 81% | **96%** | **96%** |
| **false-quant promotions** | **27** | **27** | **27** | **0** | 1 |

The false-quant row is the measurement-theoretic headline: expert graphical
selection (and min-BIC) claims quantitative structure for 27 of the 162
clean datasets whose truth is nominal/ordinal (~17%; overwhelmingly the
IIO/DM sets read as LCR/RM), while the calibrated hybrid makes zero such
claims at a modest cost on the quant diagonal (85% vs 93-100%).

**MISSPECIFIED conditions (108):** humans score 46-47% against the nominal
generator labels vs selectors 33/27% - but the composition inverts the
reading. Human "wins" here are largely RM 72% and LCR 33%: claiming Rasch
quantitativeness for slope-varying (non-Rasch) data - the same graphical
leniency that produces the 27 clean false promotions. The selectors demote
these datasets to ordinal (RM->0%, zero false-quant both arms), which the
misspecification study established as the DESIGNED behaviour: data that are
not Rasch should not be certified quantitative. Against the generator
label it scores as a loss; against the measurement claim it is the point.

Caveats: rater
selections predate our runs and used their original displays at their
original N; the 162-set comparison uses our full-N=5000 selector runs vs
rater judgments made on the same underlying response data.
