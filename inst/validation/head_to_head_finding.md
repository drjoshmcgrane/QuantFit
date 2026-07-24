# Six-model head-to-head: manifest 2×2 selector vs LR-edge lattice

Both selectors run on the SAME datasets (shared seeds), full generating
hierarchy, production settings (manifest DC n.mat=500, cc_B=49; lattice B=99).
Grid: J=12, N∈{1500,3000}, 5 reps → 60 datasets per selector. Run as two
separate passes (`SELECTOR=lattice|manifest`) because the lattice pass is ~50×
slower per dataset; join by dataset id. Script: `head_to_head.R`.

## MANIFEST selector — complete (60/60)

**Exact model 46/60 (77%). Scale-type 51/60 (85%).**

truth → selected:

|      | UN | MON | IIO | DM | LCR | RM | recovery |
|------|----|-----|-----|----|-----|----|----------|
| UN   | 9  |     | 1   |    |     |    | 9/10 |
| MON  |    | 10  |     |    |     |    | **10/10** |
| IIO  |    | 1   | 2   | 2  | 1   | 4  | **2/10** |
| DM   |    |     |     | 7  | 3   |    | 7/10 |
| LCR  |    |     |     |    | 9   | 1  | 9/10 |
| RM   |    |     |     |    | 1   | 9  | 9/10 |

Scale-type confusion (all errors in the LIBERAL direction; nothing demoted):

| truth\got | nominal | ordinal | quant |
|-----------|---------|---------|-------|
| nominal   | 9       | 1       | 0     |
| ordinal   | 0       | 22      | 8     |
| quant     | 0       | 0       | **20** |

### Reading
- **Strong**: MON 10/10, LCR/RM/UN 9/10. Quant is never demoted (20/20) — the
  conservative direction is airtight.
- **IIO is the weak spot (2/10)**: at the `simulate_responses` default class
  separation the IIO axis is underpowered (matches the separation-gated IIO
  power map — strong only at wide separation), so IIO data leaks, mostly to
  quant (5/10 → LCR/RM). This is the standing IIO identifiability floor.
- **DM→quant (3/10)**: draw-dependent. The sorted-random-logit DM draws vary in
  non-additivity strength; the weakly non-additive ones slip past double
  cancellation (the DM→quant power limit documented in
  `double_cancellation_finding.md`). Recovery by N was 5/5 @1500 vs 2/5 @3000
  here — noisy at 5 reps/cell; the leak is governed by per-draw additivity, not
  cleanly by N.
- The 8 ordinal→quant leaks (5 IIO + 3 DM) are the entire scale-type error
  budget; both are the hardest theoretical boundaries.

## LATTICE selector — complete

truth → selected:

|      | UN | MON | IIO | DM | LCR | RM | recovery |
|------|----|-----|-----|----|-----|----|----------|
| UN   | 10 |     |     |    |     |    | **10/10** |
| MON  |    | 10  |     |    |     |    | **10/10** |
| IIO  | 1  |     | 3   | 5  |     | 1  | 3/10 |
| DM   |    |     |     | 10 |     |    | **10/10** |
| LCR  |    | 1   |     |    | 9   |    | 9/10 |
| RM   |    |     |     |    | 2   | 8  | 8/10 |

Lattice B was reduced from production 99 to **49** with `boot_n_starts=2`:
production B=99 was attempted and confirmed INFEASIBLE on this machine — a single
quant-edge lattice dataset (LCR bootstrap at 7 classes) exceeds the ~20-min
background-job kill ceiling and never completes, so the LCR/RM rows could never
finish. B=49 brings each dataset under the ceiling; p-granularity 1/50 is
adequate at alpha=0.05.

## Head-to-head (60/60 shared datasets)

|                | MANIFEST 2×2 | LATTICE (LR-edge) |
|----------------|--------------|-------------------|
| exact model    | 46/60 (77%)  | **50/60 (83%)**   |
| scale-type     | 51/60 (85%)  | **57/60 (95%)**   |

Per-model recovery: UN 9/10 vs 10/10 · MON 10/10 vs 10/10 · IIO 2/10 vs 3/10 ·
DM 7/10 vs **10/10** · LCR 9/10 vs 9/10 · RM 9/10 vs 8/10.

### The decisive difference is WHERE the errors land
- **Lattice keeps errors within scale-type.** Its IIO misses go IIO→DM (both
  ordinal) — only 1 IIO error leaks to quant — so scale-type holds at 95%, and
  DM is recovered perfectly (10/10).
- **Manifest leaks ordinal→quant.** IIO→LCR/RM (5/10) and DM→LCR (3/10) are its
  entire scale-type error budget. The DM→quant leak is the double-cancellation
  power limit (see `double_cancellation_finding.md`): the lattice's LR-based
  DM-vs-quant edge is more powerful than the manifest's double-cancellation gate
  at N≤3000 / J=12.
- **Both** recover MON perfectly (10/10) and handle UN/LCR/RM well; **IIO is the
  hard model for both** (exact 2–3/10) — the standing identifiability floor.

### Verdict
At J=12, N≤3000, the LR-edge lattice is the stronger selector, chiefly because
its misclassifications respect scale-type boundaries. The manifest 2×2 is
competitive on the nominal/ordinal split and MUCH cheaper to run (fast axes
clear UN/MON/IIO in ~1s; only DM/LCR/RM pay the DC step, vs the lattice paying a
full multi-model bootstrap on every dataset), but its double-cancellation
DM→quant gate is the weak seam — the very seam the lattice tests by likelihood.
A hybrid — manifest 2×2 for the ordinal layer, LR edge for DM→quant — is the
natural follow-up.

## Note on settings vs cost

At production fidelity the head-to-head is a multi-hour grind, and the machine's
~20-min background-job ceiling forced the lattice down to B=49. The manifest
selector is the faster of the two by ~50× per quant dataset — itself a practical
finding: the manifest route is far cheaper to run at production DC fidelity.
