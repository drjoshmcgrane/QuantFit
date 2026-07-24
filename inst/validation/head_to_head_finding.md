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

## LATTICE selector — in progress

Lattice (`select_model_ll`, B=99) costs ~21 min/dataset on the quant/DM data
(LCR bootstrap at ceiling((J+1)/2)=7 classes), so its 60-dataset pass is being
accumulated across relaunches (the machine periodically kills background jobs).
Table to be filled once complete. Prior full validation of this selector
(selection_audit.R, dichotomous, B=99) reported exact 87.2% / scale-type 96.1%
under simulate_responses with the single-gate quant path.

## Note on settings vs cost

At production `n.mat=500` / `B=99` the head-to-head is a multi-hour grind. The
manifest selector is the faster of the two (fast axes clear UN/MON/IIO in ~1s;
only DM/LCR/RM pay the ~12-min DC), which is itself a practical finding: the
manifest route is cheaper to run than the LR-edge lattice at production fidelity.
