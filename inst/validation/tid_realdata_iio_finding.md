# IIO recovery on the REAL TI&D data (full N = 5000): hybrid vs lattice

The central claim of the 2x2 route is that it recovers IIO where the LR-edge
lattice cannot. A self-audit found that claim rested on (a) a lattice figure
remembered from prior-session code and (b) a single underpowered comparison.
This is the same-data, paired test on the **real TI&D archive** rather than
`simulate_responses`.

## Design

`select_model_hybrid()` on the TA datasets at **full N = 5000** (no
subsampling), paired by TA id against the prior lattice run stored in
`tid_results/` (`tid_runner.R`). nI in {6,12}; 9 IIO datasets each = 18 cells.
Runner: `tid_hybrid_full.R`.

## Result — all 18 IIO cells

|                    | IIO recovered |
|--------------------|---------------|
| **hybrid (2x2)**   | **13/18 (72%)** |
| lattice (LR edge)  | 7/18 (39%)      |
| McNemar            | 7 vs 1 discordant, **p = 0.070 — NOT significant** |

By test length: hybrid 8/9 at nI=6, 5/9 at nI=12.

Confusion (hybrid rows, lattice columns):

|            | lat DM | lat IIO | lat MON |
|------------|--------|---------|---------|
| **hyb DM** | 3      | 1       | 1       |
| **hyb IIO**| **7**  | 6       | 0       |

## Reading

- The **effect is large**: 72% vs 39%, near double. The mechanism is the
  predicted one - **7 of the 8 discordant cases are the lattice sending
  IIO -> DM**, the degenerate-null failure (DM and IIO fit near-identically on
  doubly-monotone data, so the LR test has no power to separate them). The 2x2's
  direct constraint check does not depend on likelihood separation.
- But it is **not significant at alpha = 0.05** on this data alone.

Across both experiments the direction is consistent and the effect large, while
neither is individually decisive:

| experiment                 | hybrid       | lattice     | McNemar |
|----------------------------|--------------|-------------|---------|
| simulated, K=20 (J=8), COMPLETE grid | **18/20 (90%)** | 9/20 (45%) | 9:0, **p = 0.0039** |
| **real TI&D, full N=5000** | 13/18 (72%)  | 7/18 (39%)  | 7:1, p = 0.070 |

(The simulated row was first reported from an interim 64-dataset read as "12/12
vs 6/12, p = 0.031"; the completed 120-dataset grid gives the numbers above -
larger n, same direction, and now decisively significant.) So the IIO advantage
IS established at alpha = 0.05 on the complete simulated grid, and the real-data
replication shows the same ~2x effect without individually reaching
significance. 16 of 17 discordant pairs across the two favour the hybrid.

The lattice's IIO ceiling is now measured four independent ways and agrees:
**39% here, 45% (simulated, same-code), 42% (prior full TI&D run), 41% (TI&D's
own human raters)**.

## FINAL full-coverage table (uni Mac batch, clean nI 6/12, full N = 5000, 108/108)

The quant cells were completed on an uninterrupted machine (Phase 2 of the
stable-machine batch). Per-model, of 18 each (hybrid / lattice):
UN 18/18 · MON 17/18 · **IIO 13/7** · DM 17/16 · **LCR 16/15** · **RM 16/15**.
**Exact: hybrid 97/108 (90%) vs lattice 89/108 (82%)**; scale-type 98% vs 95%;
McNemar overall 12:4, p = 0.077. False quantitativeness promotions: hybrid 0,
lattice 1. The real-data quant half confirms the shared-edge argument: hybrid
and lattice are statistically identical on LCR/RM, and the hybrid's overall
lead is carried by IIO.

## Coverage of the earlier laptop-only run (superseded by the table above)

| truth | covered | hybrid | lattice |
|-------|---------|--------|---------|
| UN    | 18/18   | 18     | 18      |
| MON   | 18/18   | 17     | 18      |
| **IIO** | **18/18** | **13** | **7** |
| DM    | 18/18   | 17     | 16      |
| LCR   | 5/18    | 4      | 4       |
| **RM**| **0/18**| — | — |

Over the four COMPLETE cells (UN/MON/IIO/DM, n = 72): hybrid 65/72 (90%),
lattice 59/72 (82%), McNemar 8 vs 2 discordant, **p = 0.109 - not significant**.
No aggregate is quoted over all rows, because the missing cells are not a random
subset: they are systematically the quantitative models, so any overall figure
would be biased.

**RM was not reachable at full N = 5000 in this environment** (0/18 across many
attempts) and LCR only partially (5/18). The cost is structural: an RM dataset
must fit DM and LCR at the bridge grain, run a 49-replicate equivalence
bootstrap, profile LCR across the grain grid, and then run `rm_vs_lcr_test`,
whose every bootstrap replicate refits the LCR profile - plausibly 1-3 hours per
dataset at N = 5000. This is a stated coverage gap, not an omission; it needs a
machine that can run uninterrupted for hours.

## METHODOLOGICAL WARNING: interim reads of this run are biased

Datasets were deliberately processed cheapest-first (UN/MON/IIO resolve in the
2x2 in ~10-20s; DM/LCR/RM pay 10-21 min for the quantitative edge) so the
decisive cells would survive interruptions. **That ordering is not random and it
is correlated with the outcome**: "cheap" means "the 2x2 settled it in the
ordinal layer", which is exactly when the hybrid gets IIO right.

The consequence was measurable. At 15/18 cells the interim read was 13/15 vs
6/15, p = 0.0156 (significant). The final 18/18 result is 13/18 vs 7/18,
p = 0.070 (not significant). **Any interim read of a cost-ordered run is
optimistically biased for whichever method resolves cases cheaply.** Quote final
numbers only.
