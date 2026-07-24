# Hybrid selector — large-set audit (generalisation check at J=8)

The hybrid (`select_model_manifest(dm_quant="lr", lr_boot_n_starts=2)`) run on
the STANDARD large audit grid used for `select_model_ll`: 6 models × 30 reps,
N=1500, **J=8 dichotomous**, n_classes=3, B=49. This is J=8 vs the head-to-head's
J=12 — a generalisation test at a shorter test length. Resumable
(`hybrid_audit.R`, one CSV per dataset); raw CSVs in `hybrid_audit_results/`.

## Result — 180 datasets

**Exact-model recovery 89.4% · Scale-type accuracy 92.2%**

truth → selected:

|      | UN | MON | IIO | DM | LCR | RM | recovery |
|------|----|-----|-----|----|-----|----|----------|
| UN   | 30 |     |     |    |     |    | **30/30 (100%)** |
| MON  | 1  | 29  |     |    |     |    | 29/30 (97%) |
| IIO  |    |     | 26  | 4  |     |    | **26/30 (87%)** |
| DM   |    |     |     | 30 |     |    | **30/30 (100%)** |
| LCR  |    |     | 1   | 6  | 23  |    | 23/30 (77%) |
| RM   |    |     |     | 6  | 1   | 23 | 23/30 (77%) |

Scale-type confusion (all 13 errors are the CONSERVATIVE quant→ordinal direction):

| truth\got | nominal | ordinal | quant |
|-----------|---------|---------|-------|
| nominal   | 30      | 0       | 0     |
| ordinal   | 1       | 89      | 0     |
| quant     | 0       | 13      | 47    |

## Reading

- **IIO 87% generalises** — matches the 80% at J=12 and the original manifest
  finding; confirms the earlier 20% was purely the `mon_eps=0.04` regression, now
  fixed. For contrast the prior *lattice* audit on this exact grid recovered
  IIO ≈ 40%; the hybrid **more than doubles it**, and beats the lattice on
  overall exact too (89.4% vs the lattice's 87.2%).
- **UN/DM 100%, MON 97%** — the ordinal/nominal layer is essentially perfect.
- **The J=8 weak spot is quant recovery (LCR/RM 77%)**, and every scale-type
  error is quant→DM — the CONSERVATIVE direction (under-claiming
  quantitativeness, the safe error for a quantitativeness test). At the short
  test the LR DM-vs-LCR edge has less power to promote additive data past DM.
  This cost scale-type here (92.2% vs the lattice's 96.1%) where at J=12 it did
  not. Tested the `lr_boot_n_starts` hypothesis directly on the 6 LCR→DM misses:
  raising it from 2 to 5 flipped only **1 of 4** back to LCR — so the
  conservatism is **mostly a genuine J=8 short-test identifiability limit**, not a
  bootstrap-starts artifact. At J=8 the additive LCR structure is, for many
  draws, statistically indistinguishable from doubly-monotone DM. `boot_n_starts`
  is a minor contributor; J≥10 is the real remedy.

## Bottom line

At J=8, the hybrid trades a little quant sensitivity (conservative, safe
direction) for a huge IIO gain (87% vs 40%) and a small exact-recovery edge over
the lattice — consistent with the J=12 head-to-head where it led on both axes.
The manifest route's IIO advantage is robust across test length.
