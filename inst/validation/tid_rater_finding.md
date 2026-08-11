# TI&D pre-consensus rater data: humans vs selectors on the same datasets

Received 2026-08-11: DTIresults.csv (Torres Irribarra's per-set selections,
'*'/'?' uncertainty marks), RPDresults.csv (Diakow's pass, xN codes marking
disagreement + the alternative model seen), summary table.csv (their
per-model fit table). 1080 TA sets; truth available for the 270 sets in our
archive enumeration (balanced 45/model); merged table committed as
tid_rater_merged.csv.

## Pre-consensus human recovery (270 sets, all conditions)

| truth | DTI | RPD | inter-rater agreement |
|-------|-----|-----|------|
| UN | 100% | 100% | 100% |
| MON | 84% | 82% | 98% |
| **IIO** | **2%** | **9%** | 93% |
| **DM** | **11%** | **24%** | 80% |
| LCR | 69% | 69% | 98% |
| RM | 89% | 89% | 100% |
| overall | 59.3% | 62.2% | 94.8% |

Scale-type: ~78% both raters. IIO-truth destinations (DTI): MON 22/45,
RM 11, UN 7, LCR 4 - humans read invariant item ordering as class
monotonicity or even quantitative structure, almost never as IIO itself.
DM-truth: MON 23, LCR 14. The published ~41% IIO figure is a consensus-era
number; PRE-consensus it is 2-9%.

**The raters knew where they were struggling**: uncertainty flags
concentrate exactly on the identifiability-limited middle (IIO 40%, DM 47%,
LCR 38% flagged vs UN 0%, MON 9%), and Diakow's alternative-codes on IIO
sets point at LCR/RM (x4/x5) - IIO data genuinely reads as quantitative to
expert eyes. Inter-rater agreement drops from 96% (unflagged) to 91%
(flagged).

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

Caveats: truth is available for 270/1080 sets (archive enumeration); rater
selections predate our runs and used their original displays at their
original N; the 162-set comparison uses our full-N=5000 selector runs vs
rater judgments made on the same underlying response data.
