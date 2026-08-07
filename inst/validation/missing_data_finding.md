# Missing-data validation

**STATUS (2026-07-28): the 240-run MCAR grid below was produced BEFORE the
audit fixes** (mean-based mask rank-matching; IIO-only axis null) **and is
being regenerated under the fixed code, extended with genuine MAR
mechanisms.** Its numbers are historical. Two of its original claims were
also wrong as written and are corrected here (external review, items 4-5).

## Corrected reading of the pre-fix MCAR grid

Numerical headlines (verified from `missing_data_results/`): exact 105/120
(10% MAR) and 104/120 (25%) vs 111/120 complete; scale-type 110-111 vs 113;
zero false quantitative promotions; zero fit failures.

**Corrections to the original interpretation:**

1. "The entire exact cost sits at the LCR <-> RM grain boundary" was FALSE.
   The actual error transitions (both rates pooled): LCR->IIO 8, LCR->MON 3,
   LCR->DM 1, LCR->RM 1, RM->LCR 4, RM->DM 3, RM->MON 1, IIO->DM 4,
   MON->DM 2, MON->UN 2, DM->MON 1, UN->IIO 1. Quant-truth errors dominate
   (21 of 31) but most leave the quant family entirely (demotions to
   IIO/MON/DM), and the ordinal cells contribute 10 errors of their own.
2. "Demotions only" was FALSE: UN->IIO at 25% is a nominal->ordinal
   promotion (one case).
3. The claim that SURVIVES: **no false promotion to quantitative structure
   in any of the 240 masked runs** - masking never fabricated quantity.

## Post-fix validation design (missing_data_validation.R v2)

Same paired audit-grid datasets (6 truths x K = 20, J = 8, N = 1500, B = 49),
three missingness mechanisms:

- `mcar`: iid cellwise Bernoulli (10% and 25%);
- `mar_ability`: genuine MAR - odd items are never masked (anchor half), even
  items are masked with probability logistic in the standardised anchor-half
  score (slope -1.2: low scorers skip more; ~18% overall);
- `mar_item`: item-dependent rates (0.05 to 0.45 across items, ~25% overall),
  person-independent.

`mar_ability` is the mechanism the rank-matched mask transfer exists for; the
pre-fix transfer would have CONFOUNDED it (it manufactured the same
dependence under MCAR), so only post-fix results are interpretable.

## Post-fix results (mdv2, 480 runs, uni Mac, build 8230019*)

*Provenance: produced under build 8230019 (mean-based mask transfer and
IIO-only axis null included; verdict-path code identical to the stamped
release for dichotomous masked data). Rows predate the per-row SHA canary;
the run log and bundle timestamps document the build.*

| mechanism | miss | exact | scale-type | stable vs complete | false quant |
|-----------|------|-------|------------|--------------------|-------------|
| MCAR 10% | 10% | 104/120 | 106/120 | 107 | 1 |
| MCAR 25% | 25% | 102/120 | 108/120 | 101 | 0 |
| item-dependent (.05-.45) | 25% | 105/120 | 106/120 | 103 | 0 |
| **MAR-by-ability (anchor-half)** | **13%** | **80/120** | **94/120** | 79 | **0** |

(Complete-data reference: 112/120 exact, 114/120 scale.)

**The headline is the MAR-by-ability row.** At only ~13% missingness, genuine
score-dependent missingness costs ~27 points of exact accuracy - far more
than 25% MCAR - and the error anatomy identifies the mechanism precisely:
every major error cell is an IIO-axis FALSE REJECTION (DM->MON 13/20,
IIO->UN 9/20, RM->MON 6/20, LCR->MON 8/20; UN/MON nearly untouched).
Ability-dependent missingness induces apparent item-response crossings; the
mask-transfer null is supposed to reproduce that artifact, but mean-based
rank matching ATTENUATES the score-missingness dependence in replicates
(induced r -0.73 observed vs ~-0.32 transferred - the documented limitation
of matching on a noisy observed-mean proxy), so the null under-represents
the artifact and the axis over-rejects.

What holds: **zero false quantitative promotions under every mechanism**
(0/480, one MCAR-10% exception at 1), and scale-type degradation is
demotion-shaped. What does not: ordinal-axis calibration under genuine
MAR-by-ability. Consequence for practice: with score-dependent missingness,
treat MON/UN verdicts with suspicion of over-demotion; the safe claims are
the quant/non-quant boundary and demonstrated dominances.

Fix direction (future work, explicitly not attempted here): replace
rank-transfer with a parametric missingness model - fit P(miss | observed
anchor score) and REGENERATE masks in each replicate from that model, which
preserves the full dependence instead of an attenuated copy.


