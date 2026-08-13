# Known-crasher datasets (build 697bcae era)

TA376, TA552 (IIO-truth) and TA218 (clean LCR-truth, nI = 24; confirmed
2026-08-13 after masquerading as a perpetually-running "last cell" - it was
dying silently in uni Mac forks) deterministically
SEGFAULT the selector at the quantitative edge's LCR-vs-DM bridge stage -
heap corruption ('invalid permissions', wild address), reproduced across
forked/non-forked execution, both EM engines, and multiple seeds. See
`crash_reproducer.R` for the minimal reproducer and crash signature.

Consequences for the full-grid extension analysis: coverage is 1077/1080
at nI <= 24 (897 of 900 attempted; 627 extension + 270 original verdicts);
these three cells are EXCLUDED and reported as crashes, not as verdicts.
(They also explain the two phase-1 stragglers that died silently on the
uni Mac: a fork that segfaults produces no error row - mclapply loses the
worker without a catchable condition, which is why no failure manifest
appeared. The runners' failure accounting cannot see signal-death in a
fork; the process-per-dataset scheduler used on the laptop can, and did.)

The crash signature matches the long-unresolved fork-isolated segfault in
the Kara runner (attribute nodes of young objects GC'd while payload
lives). These datasets are the first DETERMINISTIC reproducers - the entry
point for the eventual bug hunt. Prime suspect: compiled code shared by
both engines (constrained-fit optimizers) rather than the C++ EM.
