# Missing-data validation: masked hybrid degrades gracefully and conservatively

240-run grid (`missing_data_validation.R`): the EXACT audit-grid datasets
(same seeds, 6 truths x K = 20, J = 8, N = 1500, B = 49) with iid MAR masks at
10% and 25%, paired verdict-by-verdict against the complete-data hybrid arm in
`audit_k20_results/`. Raw results in `missing_data_results/`.

## Results (masked-correct / complete-correct, of 20)

| truth | 10% MAR | stability | 25% MAR | stability |
|-------|---------|-----------|---------|-----------|
| UN    | 20/20   | 20/20     | 19/20   | 19/20 |
| MON   | 19/19   | 20/20     | 17/19   | 18/20 |
| IIO   | 17/18   | 19/20     | 19/18   | 19/20 |
| DM    | 19/20   | 19/20     | 20/20   | 20/20 |
| LCR   | 14/17   | 15/20     | 13/17   | 13/20 |
| RM    | 16/17   | 18/20     | 16/17   | 13/20 |

- **Exact**: 105/120 (87.5%) at 10%, 104/120 (86.7%) at 25%, vs 111/120
  (92.5%) complete.
- **Scale-type: essentially unimpaired** - 110/120 and 111/120 vs 113/120
  complete. Losing a quarter of all responses costs about two scale-type
  verdicts.
- **Zero false quantitativeness promotions in all 240 masked runs; zero fit
  failures.**

## Reading

1. The ordinal layer is nearly mask-invariant: UN/MON/IIO/DM lose ~1 verdict
   per cell at either rate (IIO even gains one at 25% - draw noise around the
   known floor, not signal).
2. The entire exact-accuracy cost sits at the LCR <-> RM grain boundary:
   with a quarter of the observations gone, discrete-vs-continuous latent
   shape blurs first (LCR 13-14/20; RM stability 13/20 at 25% - verdicts
   shuffle WITHIN the quant family, so scale claims survive).
3. Degradation is in the safe direction: masked misses are demotions
   (quant -> DM/ordinal), never fabricated quantity - matching the
   misspecification study's behaviour and the design goal.
4. Runtime is unaffected (median ~20 s/dataset ordinal cells; the masked
   likelihood adds no measurable cost).
