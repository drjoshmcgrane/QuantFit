# Under misspecification, both selectors collapse — the hybrid somewhat less,
# and in a structurally different way

## Design

All 108 **non-clean** TI&D datasets (slope variation / correlated dimensions —
genuine model misspecification; every one is nI=24), full N = 5000, hybrid vs
the prior lattice run (`tid_results/`). This is the most realistic regime — real
data is always somewhat misspecified — and the regime where the lattice was
known to collapse (27% exact overall).

**Coverage: 90/108.** The 18 missing datasets are precisely those whose 2x2
reached DM and stalled at the locally unreachable nI=24 quant edge (bridge = 13
classes at N = 5000): 2 DM + 2 LCR + **14 RM**. Stated, not hidden — and see
below, because *which* datasets are missing is itself the key signal.

## Result on the covered 90 (paired)

| truth | hybrid | lattice |
|-------|--------|---------|
| UN    | **18/18** | **18/18** |
| MON   | **14/18** | 9/18 |
| IIO   | 2/18   | 0/18 |
| DM    | 0/16   | 0/16 |
| LCR   | 0/16   | 1/16 |
| RM    | 0/4    | 0/4 |
| overall | 34/90 (38%) | 28/90 (31%) |

Both selectors largely collapse. The hybrid degrades somewhat less (MON +5,
IIO +2) but this is a slower collapse, not robustness.

## The structural difference: what happens to RM data

| the 18 misspecified RM datasets | verdict |
|---------------------------------|---------|
| lattice | MON 10, UN 7, DM 1 — the quantitative structure is **lost in the ordinal layer** |
| hybrid  | 4 → MON; **14 reached DM** and proceeded to the quant edge (stalled locally) |

The hybrid's 2x2 preserves the double-monotone structure of misspecified Rasch
data on 14/18 where the lattice loses it entirely. Reaching DM is strictly
closer to the truth (RM satisfies DM's constraints) — but whether the quant edge
would then correctly promote to RM is **unresolved here**; those 18 runs need an
uninterrupted machine. (Directly tested here: a single caffeinated attempt got a
~50-min window - the longest observed - and still zero of the 8 in-flight
datasets finished, so each needs > 50 min; and the workers do NOT survive the
kill as orphans, so windows cannot be chained.) So the fair statement is: the hybrid keeps the question
alive that the lattice forecloses.

## The one clean positive: failures are all in the safe direction

Under misspecification neither selector EVER over-claims: UN is 18/18 for both
(no false structure), and every error moves *down* the hierarchy (quant/ordinal
data called MON/UN), never up toward quantitativeness. For a package whose
purpose is defensible quantitativeness claims, misspecification produces
under-claiming, not false positives.

## Caveats

- The covered-90 comparison EXCLUDES the hybrid's 18 most promising datasets
  (the DM-reaching ones), so the 38%-vs-31% gap is, if anything, an
  underestimate of the hybrid's relative position; but the converse resolution
  (quant edge demotes them all) is also possible.
- Single misspecification family (TI&D's slope/dCor conditions), nI=24 only.
- The lattice arm is prior-session code at B=59; indicative, not same-code.

## Bottom line

"Is the hybrid robust where the lattice collapses?" — **No, but it fails
better**: slower degradation, structure preserved further down the pipeline on
exactly the datasets where it matters (RM), and all failures in the
conservative direction. Completing the 18 stalled quant-edge runs on a stable
machine is the single measurement that would settle whether "fails better"
becomes "partially survives".
