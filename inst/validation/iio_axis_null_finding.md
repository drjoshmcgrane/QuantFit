# The IIO axis's null class count (v0.3.4 fix and its validation)

## The defect

The IIO axis is a model-free crossing statistic calibrated against a null
simulated from an IIO model fitted to the data. Until v0.3.4 that null was
fitted at the 2x2's fixed routing class count (default C = 3). When the
data's latent heterogeneity exceeds C, the fitted null model cannot
reproduce the observed crossing level, its simulated statistics sit
systematically below the observed one, and the axis falsely rejects
invariant item ordering.

Discovered on the TI&D nI = 48 block: the first four DM-truth datasets all
returned hybrid = MON while the lattice returned DM. A direct diagnostic
(committed under axis_null_fix_evidence/) isolated the cause immediately -
the OBSERVED statistic is identical at every C because it is model-free;
only the null moves:

| dataset | statistic | null C=3 | null C=4 | null C=5 |
|---------|-----------|----------|----------|----------|
| TA173 | 3.2632 | p = 0.020 | p = 0.260 | p = 0.240 |
| TA220 | 0.5302 | p = 0.020 | p = 0.220 | p = 0.160 |

## The fix

`iio_null_C_range` (default 2:6, floored at the routing C): the null's
class count is chosen by BIC. ROUTING still uses the fixed C, so the
stability rationale for a small fixed C in the property tests is
untouched; only the reference distribution adapts. Legacy behaviour is
`iio_null_C_range = n_classes`. All four affected TI&D datasets hold under
the adaptive null (C = 4 selected in each; p = 0.22, 0.28, 0.48, 0.96).

## Size/power validation (640 datasets, both nulls on identical data)

4 truths x J in {6,12,24,48} x true C in {3,4} x K = 20.

**Size (IIO/DM truths; rejection = false rejection): 8.8% -> 1.6% overall.**
Concentrated entirely in DM truth with C = 4 > routing C:

| J | 6 | 12 | 24 | 48 |
|---|---|----|----|----|
| old | 0% | 0% | 15% | **95%** |
| new | 0% | 0% | **0%** | **10%** |

**Power (MON/UN truths; rejection = correct): 99.4% -> 98.8% overall.**
UN unaffected at every length; MON unaffected except MON/C4/J=6, where it
falls 95% -> 85% (BIC cannot detect the 4th class at 6 items, so a minority
of datasets get a richer null than the evidence supports).

Selected null counts confirm the mechanism: at J = 48 with 4 true classes
BIC picks exactly 4; at J = 6 it stays near 3.

## Honest residuals

- **10% false rejection remains at J = 48 with C = 4** - vastly better than
  95%, but above the nominal 5%. The axis is approximately, not exactly,
  calibrated for long tests with rich latent structure. Candidate
  refinements (deferred, see PLAN_gate_improvement.md): widen the null
  range beyond 6, or add a negligibility margin to the axis as planned for
  the quantitative gate.
- **10-point power loss at MON/C4/J = 6.** Short tests with undetectable
  extra classes now occasionally retain IIO that should be rejected.
- The bias existed at EVERY test length and was decisive only where the
  statistic aggregates over many item pairs (1128 at J = 48 vs 28 at
  J = 8). This is why fixed C passed validation at J <= 24 (DM recovery
  97/100/93%) and collapsed at J = 48.

## Consequence

Every pre-0.3.4 hybrid verdict was computed with the biased null. The
lattice arm and all rater comparisons are unaffected (they do not use this
axis). Hybrid results are being regenerated; pre-fix numbers are retained
as superseded, per the practice established for the poset redesign.
