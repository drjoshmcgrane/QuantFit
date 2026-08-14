# Fixed vs adaptive class counts in the 2x2 axes: a principled asymmetry

Two studies, 1360 datasets, prompted by the v0.3.4 IIO-axis defect.

## The rule that emerged

**Fix the class count where the fit produces the STATISTIC; adapt it where
the fit produces the NULL.**

| axis | what the fit produces | underfitting (C < true) | overfitting (C > warranted) | correct policy |
|------|----------------------|-------------------------|------------------------------|----------------|
| MON | the statistic (min downward movement over class orderings of the fitted UN table) | SMOOTHS violations away -> conservative: size perfect, power sags | inflates noise into apparent violations -> catastrophic false rejection | **fixed small C** |
| IIO | the null (simulate from a fitted IIO model; the statistic is model-free) | null too TIGHT -> observed looks extreme -> false rejection | null loose -> mild power loss | **BIC-selected C** |

## Evidence: MON axis (720 datasets; routing C = 3 vs C = truth, same data)

Rejection rates, routing C=3 -> C=truth:

| truth | trueC | J=8 | J=24 | J=48 |
|-------|-------|-----|------|------|
| MON | 3 | 0->0% | 0->0% | 0->0% |
| MON | 4 | 0->**40%** | 0->0% | 0->0% |
| MON | 5 | 0->**100%** | 0->0% | 0->0% |
| DM | 4 | 0->30% | 0->0% | 0->0% |
| DM | 5 | 0->100% | 0->5% | 0->0% |
| IIO | 3 | 70->70% | 95->95% | 100->100% |
| IIO | 4 | 65->100% | 85->100% | 95->100% |
| IIO | 5 | 70->100% | 75->100% | 85->100% |
| UN | any | 100->100% | 100->100% | 100->100% |

- **Size at fixed C = 3: 0.0% (0 of 240).** No masking-induced size problem.
- **Power at fixed C = 3: 91.1%**, declining with latent richness
  (94.2% / 90.8% / 88.3% at true C = 3/4/5) - the masking cost, modest.
- **Using the true C instead: size 15.3%**, up to 100% at J = 8 with C = 5
  (five classes on eight items overfits; the noise reads as downward
  movement). Power would rise to 98.1%, but at an unacceptable size cost.

## Evidence: IIO axis (640 datasets; see iio_axis_null_finding.md)

- Fixed C = 3 null: **95% false rejection** on DM truth with true C = 4 at
  J = 48 (15% at J = 24); overall size 8.8%.
- BIC-selected null: size 1.6% overall, 10% in that worst cell; power
  99.4% -> 98.8% (one cell affected: MON/C4/J=8, 95% -> 85%).

## Consequence for the design

The hybrid keeps a fixed routing C (default 3) for the property tests -
now with measured justification rather than a stability argument - and
adapts only the IIO axis's reference distribution. No change is warranted
for the MON axis, the poset refinement (fitted-UN statistic, same logic as
MON), or the quantitative edge (its grain is set by the Lindsay bridge,
not by the routing C).
