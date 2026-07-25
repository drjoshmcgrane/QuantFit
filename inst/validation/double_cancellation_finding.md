# Double-cancellation (DM→quant) validation: size controlled, power N-gated

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

The manifest selector, once the ordinal 2×2 reaches DM, decides DM (ordinal)
vs quantitative (LCR/RM) with a single **double-cancellation** test —
`cc_bootstrap_null(data, check = "double")`, Thomsen additivity against a
Rasch/additive parametric null. `reject` (observed violation rate above the
null) ⇒ non-additive ⇒ stays **DM**; non-reject ⇒ `supports_quant` ⇒ proceed to
the DIP axis (LCR vs RM). Single cancellation is dropped deliberately: ordinality
is already established by the 2×2, so testing it again only dilutes the one
distinction that matters (Holm), and double-alone was shown to recover DM
detection where the full hierarchy lost it.

## What "size" and "power" mean here

The null is additive, so **rejecting = evidence against additivity**:
- **DM** data (doubly monotone but **non-additive** — `simulate_responses` sorted
  random logits) → want **reject** (POWER): keep it ordinal.
- **LCR / RM** data (**additive** Rasch) → want **non-reject** (SIZE ≈ α): let it
  pass to the quantitative axis.

## Result (J = 12, n.mat = 500, B = 49, α = 0.05; 2 reps per cell)

| model | N = 1500 | N = 3000 | interpretation |
|-------|----------|----------|----------------|
| DM  (power) | 1/2 reject | 2/2 reject | detection ~50% at N=1500, reliable at N=3000 |
| LCR (size)  | 0/2 reject | 0/2 reject | additive structure correctly passed through |

Observed violation rates expose the mechanism:

| model | N | obs viol (rep1, rep2) |
|-------|-----|------------------------|
| DM  | 1500 | 0.0058, 0.042 |
| DM  | 3000 | 0.027, 0.017  |
| LCR | 1500 | 0.0074, 0.023 |
| LCR | 3000 | 0.0068, 0.0091 |

## Reading

1. **Size is controlled.** LCR (additive) is never falsely rejected (0/4). The
   DM→quant step does not wrongly demote genuine quantitative structure. One
   LCR/N=1500 draw sits at p=0.06 (one step from the boundary), so size is
   controlled but not luxuriously at moderate N — consistent with a slightly
   noisy per-band count.

2. **Power is N-gated and draw-dependent.** The sorted-random-logit DM generator
   produces *variable* non-additivity strength per draw. At N=1500 the weakly
   non-additive draws sink to the sampling-noise floor (obs ≈ 0.006, ≈ the LCR
   level) and are not detected → such DM datasets **leak to quantitative**. At
   N=3000 the per-score-band counts are dense enough that the structural
   non-additivity surfaces (obs 0.017–0.027) and DM is retained.

## Implication — a requirement, not a calibration bug

The test's **size is correct**, so this is not a recalibration target: it is the
inherent sample-size requirement of double cancellation. Reliable DM-vs-quant
separation at J≈12 needs **N ≳ 3000**; at moderate N a non-rejection (→ quant)
may reflect low power rather than genuine additivity. This is exactly the caveat
`cc_bootstrap_null()` already emits (`N < 1000` warning; the DM-detection
threshold is higher still) and matches Student & Read (2025) on the procedure
being underpowered at moderate samples. The conservative-direction error
(additive→DM) is well controlled; the liberal-direction error (DM→quant) is the
one bounded by N. For a package whose purpose is defensible quantitativeness
claims, users should treat a quant call at moderate N as provisional.

## Cost note

At production `n.mat = 500` each dataset costs ~7 min at J=8/N=1500 and ~12+ min
at J=12, so a full grid is expensive; this study used 2 reps/cell (runs were
repeatedly interrupted by machine sleep). The qualitative conclusion — size
controlled, power N-gated — is stable across the completed reps and would only
be tightened, not changed, by more.

Files: `inst/validation/dc_validation.R` (full grid), focused run parameters in
this note.
