# The "PO" model: TI&D-idiom data whose ground truth is a partial order

## The generator

TI&D's models 0-3 are one C x J matrix of U(-4,4) logits with each column
(item) fully sorted for MON and left unsorted for UN. Sorting a column imposes
the CHAIN; not sorting is the ANTICHAIN. `simulate_responses("PO")` is the
natural interpolation inside the same idiom: each column's values are sorted
descending and assigned to classes along a **random linear extension of a
specified partial order**, drawn independently per item.

- Comparable pairs precede in *every* linear extension, so they dominate
  **exactly, on every item**.
- Incomparable pairs get random per-item relative order and therefore cross; a
  rejection loop (TI&D's own device - cf. their LCR separation loop) guarantees
  each incomparable pair crosses in both directions by >= `po_margin` (default
  0.05) on the probability scale.
- **The nesting is exact**: the chain poset reproduces TI&D's MON column sort
  byte-identically, and the antichain reduces to their UN. Verified across 600
  draws: zero dominance violations, zero crossing failures.
- The item side is left unconstrained (as in UN/MON), isolating the class
  structure as the signal. Posets: `"V"` (default), `"Lambda"`, `"single"`, or
  a user dominance matrix.

## Validation at TI&D spec (N = 5000, C = 3; poset x nI {6,12,24} x K = 10)

### Hybrid arm (120 datasets)

| metric | result |
|--------|--------|
| PO data routed to the UN cell | **90/90** (no IIO-cell leak) |
| **poset-rung sensitivity** (flags "partial") | **89/90 (99%)** |
| specificity on matched UN controls | 28/30 (93%); **20/20 at nI >= 12** |
| dominance-count recovery (floored q05 = true pair count) | **85/89** (V/Lambda -> 2: 28/30 each; single -> 1: 29/29) |

All three errors sit at nI = 6: the one missed PO case is the weakest structure
("single", one dominance pair) at the shortest test, and both false partials
are nI = 6 UN controls.

### Lattice arm (nI = 12, 24 datasets)

The six-model lattice files **every** PO dataset as plain `UN` - identical to
its verdict on true UN data. The partial order is silently conflated with
"nominal", exactly the gap identified in `partial_order_finding.md`; nothing in
the six-model vocabulary can express or flag it.

## What this settles

1. **The rung works when the structure is real.** Its earlier poor sensitivity
   (2/5) on TI&D's misspecified non-clean data is now explained: that dominance
   structure sat mostly at noise level and correctly failed to survive
   resampling. On genuine, margin-backed partial orders the calibrated rung is
   near-perfect (99%) AND recovers the number of dominance relations almost
   exactly. It distinguishes real from incidental partial order - the desired
   behaviour on both sides.
2. **The seventh model completes the validation family**: ground-truth PO data
   now exists inside TI&D's own generative idiom, the hybrid routes it
   correctly, and the refinement labels it correctly while leaving genuine UN
   data alone at realistic test lengths.
3. **The framework gap is now demonstrable by construction, not just by
   diagnosis**: a selector limited to the six models cannot even in principle
   report the truth for these data; the hybrid's poset refinement can and does.

## Caveats

- `po_margin = 0.05` guarantees detectable crossing; PO structure with shallower
  incomparability would grade into the low-sensitivity regime seen on the
  non-clean TI&D data. Sensitivity is a function of crossing depth (and nI = 6
  is marginal even at 0.05).
- The rung still refines only the UN cell; a partial order combined with
  invariant item ordering would land in the IIO cell, where no poset is
  computed. That combination does not arise from this generator (item side
  unconstrained) but is a possible extension.
- C = 3 throughout; larger C admits richer posets (untested).
