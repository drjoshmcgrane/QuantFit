# QuantFit: consolidated validation results

The authoritative summary of every final number, with pointers to the finding
documents that carry the detail. Rules of this document: FINAL numbers only -
interim reads of cost-ordered runs have flipped conclusions twice (documented
in `tid_realdata_iio_finding.md`) and are banned; every claim names its source
file; superseded (pre-fix / count-era) figures appear only where marked
historical. Last consolidated: 2026-07-28.

**AUDIT NOTE (2026-07-28, external review):** three correctness fixes landed
after the grids below were produced - mean-based mask rank-matching (the
total-based transfer manufactured ability-dependent missingness under MCAR),
a practical-equivalence guard on the degenerate-null retention (observed LR
must also be negligible), and an IIO-only null for the hybrid's IIO axis
(previously simulated from DM). Sections marked [RE-RUN PENDING/DONE] record
which grids have been regenerated under the fixed code; unmarked simulation
results predate the fixes and are historical until re-verified.

---

## 1. The framework

Torres Irribarra & Diakow's six-model hierarchy on the latent class x item
table - UN (nominal), MON / IIO / DM (ordinal), LCR / RM (quantitative) -
plus this package's refinement: calibrated PARTIAL-ORDER tests on whichever
side the six-model verdict leaves unordered (section 6).

Two selectors share one quantitative edge:
- **lattice** (`select_model_ll`): calibrated LR tests over the complete
  adjacent-edge ordinal lattice at a BIC-selected class count, then the
  published DM -> LCR -> RM succession with the LCR-vs-DM comparison at the
  fixed Lindsay bridge grain. The severity (estimated-power) check is OPT-IN
  (`severity = TRUE`, default FALSE); the historical single LCR-vs-UN
  `alpha_quant` gate describes the 2026-07-12 audit era, not current code.
- **hybrid** (`select_model_hybrid`): the 2x2 property layer (IIO crossing
  statistic + MON perm-min q05 vs `mon_eps`) routes to UN/MON/IIO/DM; from DM
  the same LR machinery (`ll_equivalence_test` at the Lindsay bridge grain
  ceiling((score_max+1)/2), then `rm_vs_lcr_test`) decides quantity.

## 2. Lattice selector on its own terms (post-fix audits)

Source: `selection_audit.R`, memory of 2026-07-12 audits; all pre-fix numbers
obsolete.

- **Dichotomous** (K=30/model, B=99, N=1500, J=8, C=3): exact 87.2%,
  scale-type 96.1%; nominal 100%, ordinal 94%, quant 97%.
- **Polytomous** (4-category, K=20, J=6): exact 89.9%, scale-type 97.5% -
  the single-gate fix carries cleanly.
- **With the estimated-power check** (K=30, commit 9e7143d): MON 30/30,
  ordinal->nominal 0, scale-type 96.7%; the 2 residual quant->DM demotions are
  data-level ambiguity (fitted-UN references genuinely separate).
- **alpha_quant = 0.05 (SHIPPED)**: false-quant 1/90 (1.1%), quant recovery
  53/60; 0.01 trades to 4/90 false vs 58/60 recovery - a loss function, not an
  accuracy difference. USER DECISION: 0.05.
- **IIO <-> DM is the irreducible boundary**: IIO exact recovery ~40% under
  the lattice, cf. TI&D's own human raters at 41%.

## 3. Hybrid vs lattice, same code, same data (paired grids)

Source: `audit_selector.R` + `audit_k20_results/`, `hybrid_audit` findings.

- **Complete K=20 paired grid** (J=8, B=49): **IIO 18/20 (90%) hybrid vs 9/20
  (45%) lattice, McNemar 9:0, p = 0.0039** - the load-bearing result. Overall
  111/120 vs 103/120 (p = 0.057). DM boundary held: DM 10/10 vs 7/10 for the
  removed double-cancellation variant (why the LR edge stays).
- **Test length** (`test-length-recovery-finding`): J=20 rescues LCR (+30
  points) but NOT IIO - the IIO limit is class-count-, not item-, bound.
- **IIO ceiling measured four independent ways**: 39% (real TI&D), 45%
  (simulated same-code), 42% (prior full TI&D run), 41% (TI&D's human raters).
  Convergent identifiability floor (`iio-ceiling-maxstat-negative`).

## 4. Real TI&D development data (full N = 5000)

Source: `tid_realdata_iio_finding.md`, `tid_hybrid_full.R`,
`tid_realdata_results/` (all phases complete 2026-07-27).

**DEFINITIVE clean table, 162/162 (nI 6/12/24, 9 datasets per cell):**

| nI | UN | MON | IIO | DM | LCR | RM |
|----|----|-----|-----|----|-----|----|
| 6  | 9/9\|9 | 8/9\|9 | **8/9\|4** | 8/9\|7 | 9/9\|8 | 9/9\|9 |
| 12 | 9/9\|9 | 9/9\|9 | **5/9\|3** | 9/9\|9 | 7/9\|7 | 7/9\|6 |
| 24 | 9/9\|9 | 9/9\|9 | **6/9\|5** | 7/9\|8 | 7/9\|7 | **7/9\|3** |

- **Exact: hybrid 142/162 (88%) vs lattice 130/162 (80%), McNemar 20:8,
  p = 0.036 - significant.** Scale-type 96% both. False quantitativeness:
  hybrid 0, lattice 1.
- Leads: IIO 19|12 (as everywhere), RM gap opens at nI=24 (7/9 vs 3/9);
  LCR a tie throughout (shared-edge argument confirmed on real data).
- **Misspecified conditions** (108/108; slope variation / non-Rasch):
  hybrid 33% vs lattice 27% exact; RM->RM 0/18 BOTH (7/18 RM reach DM under
  the hybrid - structure degrades to ordinal, never to false quantity);
  **zero false quant claims in every misspecified run**
  (`misspecification_finding.md`).

## 5. Partial order (the seventh structure)

Source: `po_validation_finding.md` (calibrated era; v1 count-era superseded),
`po_validation.R` v2 + `po_validation2_results/`, commit 6ee0046.

- Each unordered side gets chain/partial/antichain via two calibrated tests:
  the axis + a NEW antichain test (permutation null of the fitted UN table,
  continuous dominance-asymmetry statistic, eps-free). In BOTH selectors;
  never changes `selected`.
- **Size** ~alpha on true antichains (2/28); controls scored against
  POPULATION posets - row-sorting creates real class dominance and the test
  correctly finds it.
- **Power (class test): J-governed.** ~50-60% at J=6; **8/8 in every cell at
  J>=12, every margin, correct V type always** (mean p 0.02). Item test:
  8/8 at nI=24; nI=8 floor partly the generator's margin shortfall.
- **Audit arm** (K=20, J=8, margin 0.10): all 40 datasets route UN in both
  selectors; class detection lattice 17/20, hybrid 13/20 - both on the J=8
  floor; discordance 5:1, below K=20 resolution.
- Generators: "PO" (per-item linear extensions; chain = MON byte-identically),
  "PO invariant" (constructive, routes to IIO cell), "PO_ITEMS" (routes to
  MON cell). All in TI&D's U(-4,4) idiom.

## 6. Defaults and their evidence

| default | value | evidence |
|---------|-------|----------|
| alpha (edges/axes) | 0.05 | conventional; power-check guards rejections |
| alpha_quant (lattice gate) | 0.05 | 4x fewer false-quant claims vs 0.01 (sec. 2) |
| mon_eps (hybrid MON axis) | 0.01 | **plateau study 2026-07-28** (`mon_eps_finding.md`): 57/60 routing flat across {0.0025..0.02}; empty q05 band 0.011-0.03 separates holds-truths (<=0.01) from IIO signal (0.033-0.064); loosening to 0.04/0.08 collapses IIO to 4/10 / 0/10; tightening buys nothing (missed IIO have q05 = 0) |
| latent nulls (CC/Kara) | empirical (Bock-Aitkin) | normal null false-rejects 6/12 on LCR data vs 1/12 empirical (`karachecks` memories); normal KEPT in rm_vs_lcr where shape IS the hypothesis |
| Lindsay bridge grain | ceiling((score_max+1)/2), floor 2, no cap | removes grain selection from the scale test (sec. 1) |

## 7. Conjoint / axiom routes (non-selector evidence)

- **ConjointChecks bootstrap null** (`cc-bootstrap-null-student-read`):
  Rasch-bootstrap percentile p replaces the last CC heuristic; size study
  bimodal additive 0/24 rejections (sufficiency-protected); CC hierarchy
  diagnostic (single->double->triple with attribution).
- **KaraChecks**: validated against Karabatsos 2018 ACMtest.m; UNRELIABLE on
  sum-score/additive data (`kara-sumscore-unreliable`) - triangulate with
  CC + LC routes, never alone.
- Cross-route consistency on LCR data: CC 1/12, Kara(empirical) 1/12,
  Kara(normal) 6/12 false - hence the empirical default.

## 8. Missing data (MAR, masked likelihood)

Design (`missing-data-design`, user's principle): masked likelihood in all
fits; observation-weighted cells for conjoint; every bootstrap / permutation
null re-imposes the observed mask rank-matched. MAR assumed throughout.

- Lattice: NA-capable since the design landed.
- Hybrid: NA-capable 2026-07-28 (pairwise-complete IIO statistic, masked
  poset nulls, `.hyb_item_probs` shape fix). Smoke: MON/IIO/DM/RM verdicts
  identical complete vs 10% MAR.
- **Deep validation [RE-RUN PENDING]**: the first 240-run MCAR grid predates
  the audit fixes and two of its claims were misreported (see the corrected
  `missing_data_finding.md`: errors were NOT confined to the LCR <-> RM
  boundary, and one UN -> IIO promotion occurred). Its surviving claim: zero
  false QUANT promotions in 240 masked runs. The v2 grid (fixed code; MCAR +
  genuine MAR-by-ability + item-dependent rates) replaces it here on
  completion.

## 9. TI&D graded-180 pre-consensus rater data

**RESERVED.** Awaiting the pre-consensus rater-level dataset from Torres
Irribarra. Planned analysis, so the section is ready to fill:

- Run both selectors on the rater-level (pre-consensus) gradings; compare
  against the published consensus verdicts and TI&D's own graphical/human
  selections (their raters agreed 41% on IIO - the benchmark the IIO ceiling
  keeps reproducing).
- Poset refinement on every non-DM verdict: do human pre-consensus gradings
  carry partial-order structure that the consensus process erases?
- Missing-data machinery applies directly if rater coverage is incomplete.
- _Results, table, and reading: TO BE ADDED when data arrive._

## 10. Known limits (honest-boundary list)

- IIO <-> DM: ~40-50% exact IIO is an identifiability floor shared by humans
  and both selectors; more items do not fix it (class-count limited).
- Two-in-sixty quant->DM demotions are data-level (chance-close draws);
  100% scale recovery is impossible in principle for nested families.
- Class antichain test needs J >= ~12; item antichain test needs long tests
  (or C > 3) for full power - the axes' own data floors, inherited.
- Kara route unreliable on additive data; CC sufficiency-protected but
  low-powered at small nI.
- Polytomous PO generators and C >= 4 posets untested.
