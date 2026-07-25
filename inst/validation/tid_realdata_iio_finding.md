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
| simulated, K=20 (J=8)      | 12/12        | 6/12        | p = 0.031 |
| **real TI&D, full N=5000** | 13/18 (72%)  | 7/18 (39%)  | p = 0.070 |
| combined discordant pairs  | 13           | 1           | -       |

13 of 14 discordant pairs across the two favour the hybrid. Pooling them would
give a very small p, but that combines two different data sources and is a
judgement call, so it is flagged rather than asserted.

The lattice's IIO ceiling is now measured four independent ways and agrees:
**39% here, 45% (simulated, same-code), 42% (prior full TI&D run), 41% (TI&D's
own human raters)**.

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
