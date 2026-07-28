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

## Post-fix results

_[RE-RUN PENDING - table inserted on completion.]_
