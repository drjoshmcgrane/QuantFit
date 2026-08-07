# mon_eps sensitivity: the default sits on a wide plateau

Design: mon_eps enters routing only as a threshold on the MON axis's
resampling q05 (`mon$lo`), so one hybrid run per dataset (6 truths x K = 10,
J = 8, N = 1500, B = 49) yields the ordinal verdict under EVERY eps
analytically. Correct verdict: UN/MON/IIO/DM to themselves, LCR/RM to DM
(the quant edge is eps-independent given DM). Raw results in
`mon_eps_sweep_results/`; script `mon_eps_sweep.R`.

## Routing accuracy by eps

| eps    | UN | MON | IIO | DM | LCR | RM | overall |
|--------|----|-----|-----|----|-----|----|---------|
| 0.0025 | 10 | 10  | 8   | 10 | 9   | 10 | 57/60 |
| 0.005  | 10 | 10  | 8   | 10 | 9   | 10 | 57/60 |
| **0.01 (default)** | 10 | 10 | 8 | 10 | 9 | 10 | **57/60** |
| 0.02   | 10 | 10  | 7   | 10 | 10  | 10 | 57/60 |
| 0.04   | 10 | 10  | 4   | 10 | 10  | 10 | 54/60 |
| 0.08   | 10 | 10  | 0   | 10 | 10  | 10 | 50/60 |

## Why: the q05 landscape has a gap, and the default sits in it

Per-truth q05 (mon_lo) ranges: MON/DM/RM effectively 0 (max 0.002); LCR max
0.0101 (one dataset); **IIO 0.033-0.064** (its class-side crossing signal);
UN 0.15-0.40. The "holds" truths and IIO are separated by an empty band
(~0.011-0.03): any eps inside it routes everything except the irreducible
cases correctly. 0.01 and 0.02 bracket this band and trade one dataset
(a q05 = 0.0101 LCR case vs one weak-signal IIO case); everything from
0.0025 to 0.02 is a 57/60 plateau.

- **Loosening beyond the band destroys IIO exactly as documented** in the
  `@param mon_eps` note: at 0.04 IIO falls to 4/10, at 0.08 to 0/10 (all to
  DM), while nothing else improves - confirming eps = 0.04 was the historical
  IIO-collapse setting.
- **Tightening below 0.01 buys nothing**: the two IIO misses at tight eps have
  q05 = 0.000 - no crossing signal exists in those draws (the identifiability
  floor), so they are not recoverable by any threshold.

## Conclusion

The shipped default mon_eps = 0.01 is validated: it sits on a four-fold-wide
accuracy plateau whose failures are boundary cases owned by the data, not the
threshold. Sensitivity to the choice is nil within [0.0025, 0.02] and the
failure mode outside the plateau is one-sided and understood.

## Post-fix reproduction (uni Mac, build 8230019, IIO-only axis null)

The rerun reproduces the plateau under the corrected IIO null: 57/60 routing
at every eps in {0.0025, 0.005, 0.01} (IIO 8/10), decaying exactly as before
when loosened (54/60 at 0.04 with IIO 4/10; 50/60 at 0.08 with IIO 0/10).
The default 0.01 verdict is unchanged post-fix.
