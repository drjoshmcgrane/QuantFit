# Crash-prone datasets (TI&D archive) - CORRECTED 2026-08-13 (evening)

## Status

| dataset | truth | nI | outcome |
|---------|-------|----|---------|
| TA552 | LCR | 12 | **RESOLVED** - completed 2026-08-13 (43 min, build 5d3f8fc): hybrid LCR, lattice LCR, `quant_edge_failed = TRUE`. A lattice-only run also completed independently (LCR, `rm_stage_failed = TRUE`). |
| TA376 | IIO | 24 | crashes (hybrid runner and lattice-only, repeatedly, both machines) |
| TA218 | LCR | 24 | crashes (hybrid runner and lattice-only, repeatedly, both machines) |

Coverage at nI <= 24: **898/900 verdicts** (628 extension + 270 original),
2 outstanding.

## CORRECTION to the earlier characterisation

The first version of this file called the crashes "deterministic in the
data, engine-independent, seed-independent" and labelled all three cells
IIO-truth. Both claims were wrong:

- **Not deterministic.** TA552 crashed repeatedly (forked and singleton,
  use_cpp TRUE and FALSE, seeds 1 and 2) and then COMPLETED normally on a
  later attempt with no relevant code change (697bcae -> 5d3f8fc differs
  only by an env-var nI ceiling in a validation runner). The honest
  description is an INTERMITTENT memory-corruption bug whose manifestation
  depends on allocation state - which varies with machine, load, and which
  cells ran before in the same process. Repeated failure on one machine is
  not proof of determinism; it is proof of a reproducible *context*.
- **Truth labels were wrong**: TA552 and TA218 are LCR-truth (nI 12 and
  24), TA376 is IIO-truth (nI 24). The earlier "all IIO-truth" statement
  came from an unchecked assumption, not the archive key.

What survives from the original analysis: the crash site (the LCR-vs-DM
bridge stage of the quantitative edge), the signature (SIGSEGV, 'invalid
permissions', wild address = heap corruption), and the fact that a
segfaulting fork produces no catchable R condition - which is why phase-1
stragglers vanished without a failure manifest, and why the
process-per-dataset scheduler (exit code 139) was needed to see them.

`crash_reproducer.R` remains useful, but should be run repeatedly: it
reproduces the crash *often* on TA218/TA376, not certainly.

## Consequence for the analyses

TA552's verdicts are included in the committed evidence. TA218 and TA376
are excluded and reported as crashes, never as verdicts. Analysis
denominators state 898/900 accordingly.
