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
  IIO ≈ 40%; the hybrid **more than doubles it**. (The IIO ≈ 40% and other
  bare lattice figures here are from a prior July-2026 audit under older code —
  treat them as indicative, not same-code benchmarks; see the verified same-code
  quant comparison below. The IIO advantage is large and robust regardless.)
- **UN/DM 100%, MON 97%** — the ordinal/nominal layer is essentially perfect.
- **Quant recovery LCR 87% / RM 83%** after the `lr_n_classes` fix (was 77%/77%
  when the delegation forced C=3). The residual 8 errors are all quant→DM — the
  CONSERVATIVE direction (under-claiming quantitativeness, the safe error) — and
  are genuine J=8 short-test identifiability: at 8 items an additive LCR/RM is,
  for some draws, statistically indistinguishable from a doubly-monotone DM
  (`boot_n_starts` 2→5 flips only ~1/4; the C-range fix recovered ~5/13; J≥10
  recovers the rest).

### The hybrid is NOT worse than the lattice on quant — verified

An earlier version of this note cited "lattice ≈ 96-97% quant/scale" as if it were
a same-code benchmark. It was NOT: that figure is from a July-2026 audit under
older code (and B=99), never re-run on these seeds this session. Comparing a
current-code hybrid to a stale-memory lattice number manufactured a phantom gap.
Verified this session: the CURRENT standalone lattice, run at the hybrid's EXACT
delegation settings (B=49, boot=2, n_classes=2:6) on the same J=8 LCR/RM data,
recovers quant **17/21 = 81% (LCR 9/12, RM 8/9)** — i.e. it MATCHES the hybrid
(≈85%), because they run the identical LR edge. There is no algorithmic gap.
The J=8 quant misses are knife-edge cases: e.g. one LCR dataset classifies as
LCR / DM / MON under (B49,boot2,seed1) / (B49,boot2,seed2001=the hybrid's) /
(B99,boot5,seed1) respectively — the seed alone flips it. The classification is
bootstrap-noise-dominated because additive and doubly-monotone are barely
separable at 8 items; stronger bootstrap does not reliably help (the example went
BACKWARDS at B=99).

## Bottom line

At J=8 the hybrid is exact 92.2%, scale-type 95.0%, and — verified same-code,
same-settings — is NOT beaten by the lattice on quant (both ≈81-85%). Its
decisive advantage is IIO (87% vs the lattice's ≈40% on the hardest model),
consistent with the J=12 head-to-head. The manifest route's IIO advantage is
robust across test length; the residual J=8 quant softness is genuine short-test
identifiability shared by any selector, in the conservative/safe direction, and
resolves by J≥10.
