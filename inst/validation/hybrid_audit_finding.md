# Hybrid selector — large-set audit (generalisation check at J=8)

The hybrid (`select_model_manifest(dm_quant="lr", lr_boot_n_starts=2)`) run on
the STANDARD large audit grid used for `select_model_ll`: 6 models × 30 reps,
N=1500, **J=8 dichotomous**, n_classes=3, B=49. This is J=8 vs the head-to-head's
J=12 — a generalisation test at a shorter test length. Resumable
(`hybrid_audit.R`, one CSV per dataset); raw CSVs in `hybrid_audit_results/`.

## Result — 180 datasets (lr_n_classes = 2:6)

**Exact-model recovery 92.2% · Scale-type accuracy 95.0%**

> Updated from an earlier run that used `lr_n_classes = 3` (fixed C) in the LR
> delegation — a handicap that forced the quant-edge fits at C=3 and cost LCR/RM
> recovery. Passing the lattice's own C-range (2:6) recovered LCR 23→26 and
> RM 23→25 with ZERO collateral (UN/MON/IIO/DM unchanged). See the
> `lr_n_classes` note below.

truth → selected:

|      | UN | MON | IIO | DM | LCR | RM | recovery |
|------|----|-----|-----|----|-----|----|----------|
| UN   | 30 |     |     |    |     |    | **30/30 (100%)** |
| MON  | 1  | 29  |     |    |     |    | 29/30 (97%) |
| IIO  |    |     | 26  | 4  |     |    | **26/30 (87%)** |
| DM   |    |     |     | 30 |     |    | **30/30 (100%)** |
| LCR  |    |     | 1   | 3  | 26  |    | 26/30 (87%) |
| RM   |    |     |     | 4  | 1   | 25 | 25/30 (83%) |

Scale-type confusion (all 9 errors conservative: quant→ordinal or MON→nominal):

| truth\got | nominal | ordinal | quant |
|-----------|---------|---------|-------|
| nominal   | 30      | 0       | 0     |
| ordinal   | 1       | 89      | 0     |
| quant     | 0       | 8       | 52    |

### The `lr_n_classes` fix

The first audit forced the LR delegation to a single C=3 (the 2x2's class count).
Isolation test on the LCR/RM→DM misses (vary ONLY n_classes, C=3 vs 2:6, all
else equal): 3/4 recovered to the correct quant model. Root cause: at C=3 the
quant-edge fits are misspecified for finely-graded data, so the LCR-vs-DM edge
retains DM. Default is now `lr_n_classes = 2:6` (the standalone lattice's range),
which recovered 5 of the 13 quant leaks at J=8, cleanly (no UN/MON/IIO/DM change),
lifting exact 89.4→92.2% and scale-type 92.2→95.0%. The remaining 8 quant→DM are
genuine J=8 short-test identifiability (additive LCR/RM indistinguishable from
doubly-monotone DM at 8 items) — the conservative/safe direction; J≥10 recovers
them.

## Reading

- **IIO 87% generalises** — matches the 80% at J=12 and the original manifest
  finding; confirms the earlier 20% was purely the `mon_eps=0.04` regression, now
  fixed. For contrast the prior *lattice* audit on this exact grid recovered
  IIO ≈ 40%; the hybrid **more than doubles it**, and beats the lattice on
  overall exact too (92.2% vs the lattice's 87.2%).
- **UN/DM 100%, MON 97%** — the ordinal/nominal layer is essentially perfect.
- **Quant recovery LCR 87% / RM 83%** after the `lr_n_classes` fix (was 77%/77%
  when the delegation forced C=3). Scale-type is now 95.0%, essentially matching
  the lattice's 96.1% while beating it on exact and IIO. The residual 8 errors
  are all quant→DM — the CONSERVATIVE direction (under-claiming quantitativeness,
  the safe error) — and are genuine J=8 short-test identifiability: at 8 items an
  additive LCR/RM is, for some draws, statistically indistinguishable from a
  doubly-monotone DM (`boot_n_starts` 2→5 flips only ~1/4; the C-range fix
  recovered ~5/13; J≥10 recovers the rest).

## Bottom line

At J=8 the hybrid is the best selector overall: exact 92.2% and scale-type 95.0%,
beating the lattice on exact (87.2%) and IIO (87% vs ≈40%) and matching it on
scale-type (96.1%). It trades a little quant sensitivity at very short tests
(conservative, safe direction) for a decisive IIO gain — consistent with the
J=12 head-to-head where it led on both axes. The manifest route's IIO advantage
is robust across test length.
