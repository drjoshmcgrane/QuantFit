# Six-model head-to-head: manifest 2×2 selector vs LR-edge lattice

> **API STATUS (updated).** The selector called "manifest" below (2x2 ordinal
> layer closed with a **double-cancellation** additivity test) has since been
> **REMOVED**. Its 2x2 ordinal layer survives in **`select_model_hybrid()`**
> (formerly `select_model_manifest`), which closes DM->quant
> with the **LR edge** instead. Double cancellation is untouched as its own
> route: `cc_bootstrap_null()` / `cc_bootstrap_hierarchy()`. Reason for the
> removal: TI&D's hierarchy is defined on the *latent* class x item table and
> class monotonicity has no faithful manifest proxy, so mixing an
> observable-conjoint axiom into a latent hierarchy was incoherent as well as
> empirically worse -- see `manifest_coherence_finding.md`.

> **MAJOR CORRECTION (supersedes the numbers first recorded here).** The first
> manifest/hybrid runs used `mon_eps = 0.04`, which was an IIO-killing regression
> (calibrated on an artificial low-separation corner; it classified IIO data as
> DM on the real generator — see `manifest_mon_calibration_finding.md`). With the
> fix (`mon_eps = 0.01`) IIO recovery went 2/10 → 8/10 and the whole verdict
> flipped. **All tables below are the corrected `mon_eps = 0.01` results.** The
> lattice pass is unaffected (it does not use the MON axis) and is unchanged.

Both selectors run on the SAME datasets (shared seeds), full generating
hierarchy (manifest DC n.mat=500, cc_B=49; lattice B=49 - see the settings
note below, which supersedes an earlier B=99 claim).
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
| IIO  | 1  |     | 8   |    |     | 1  | **8/10** |
| DM   |    |     |     | 7  | 3   |    | 7/10 |
| LCR  |    |     |     |    | 9   | 1  | 9/10 |
| RM   |    |     |     |    | 1   | 9  | 9/10 |

**Manifest exact 52/60 (87%), scale-type 54/60 (90%).**

### Reading
- **IIO 8/10** — the manifest route does exactly what it was built for and
  **beats the lattice (3/10) decisively** on the hardest model. (The broken
  eps=0.04 had hidden this at 2/10.) The 2 misses are 1 IIO→UN (IIO-axis miss)
  and 1 IIO→RM.
- MON 10/10; UN/LCR/RM 9/10; quant never demoted.
- **DM→quant (3/10) is now the manifest's one real weak spot**: the
  sorted-random-logit DM draws vary in non-additivity strength; the weak ones
  slip past double cancellation (the DM→quant power limit in
  `double_cancellation_finding.md`). This is the seam the HYBRID closes with the
  LR edge.

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

Lattice B was reduced from production 99 to **49** with `boot_n_starts=2`, which
also MATCHES the hybrid's LR delegation so the paired comparison is at equal
settings. (An earlier version of this note claimed B=99 was "confirmed
infeasible on this machine" because of a ~20-min job-kill ceiling. That belief
was later tested and DISPROVEN — a probe ran 14+ min and an 8-worker job 16+ min
uninterrupted; the kills coincided with several concurrent jobs oversubscribing
the cores. B=49 is a matched-settings choice, not a forced one.)

## Three-way head-to-head (60/60 shared datasets, mon_eps = 0.01)

`select_model_hybrid(lr_boot_n_starts = 2)` is the HYBRID:
manifest 2×2 for the nominal/ordinal layer, the lattice's likelihood-ratio
LCR-vs-DM edge for DM→quant.

|                | **HYBRID**   | MANIFEST     | LATTICE      |
|----------------|--------------|--------------|--------------|
| exact model    | **55/60 (92%)** | 52/60 (87%) | 50/60 (83%)  |
| **scale-type** | **58/60 (97%)** | 54/60 (90%) | 57/60 (95%)  |

Per-model recovery (hybrid / manifest / lattice), of 10 each:
UN 9/9/10 · MON 10/10/10 · **IIO 8/8/3** · **DM 10/7/10** · LCR 10/9/9 · RM 8/9/8.

### Verdict — the hybrid wins on BOTH axes; the manifest beats the lattice
This REPLACES the earlier (broken-eps) verdict that "the lattice is the stronger
selector." With IIO detection working:

- **IIO: manifest & hybrid 8/10 vs lattice 3/10.** The manifest route does what
  it was built for and beats the LR lattice decisively on the hardest model —
  the LR DM-vs-IIO edge is degenerate (near-zero power both tails), the direct
  property test is not. This is the whole reason the manifest approach exists.
- **Manifest alone beats the lattice on exact (87% vs 83%)** — driven entirely by
  IIO (8 vs 3). The lattice still edges scale-type (95% vs 90%) because the
  manifest's DM→quant double-cancellation leak (DM 7/10) costs 3 ordinal→quant
  errors.
- **The hybrid is the clear best (92% exact, 97% scale-type)**: it keeps the
  manifest's IIO strength (8/10) AND fixes the DM seam with the LR edge
  (DM 7→10/10), so it leads on exact by ~9 points and ties/edges the top
  scale-type. Only one ordinal→quant leak; quant 20/20.
- It is also CHEAPER than the full lattice: UN/MON/IIO data is cleared by the
  ~1s manifest 2×2 and never invokes the LR machinery; only DM-reaching datasets
  pay the edge cost.

The hybrid is the recommended selector: manifest-2×2 economy + IIO power on the
ordinal layer, LR-edge rigor on the one boundary (DM→quant) that decides
quantitativeness.

## Note on settings vs cost

At production fidelity the head-to-head is a multi-hour grind. The 2×2 ordinal
layer is far cheaper than the lattice (~1s vs minutes per dataset when the data
resolves to UN/MON/IIO), which is why the hybrid — which only invokes the LR
machinery on DM-reaching data — costs less than the full lattice.
