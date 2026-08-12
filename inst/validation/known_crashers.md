# Known-crasher datasets (build 697bcae era)

TA376 and TA552 of the TI&D archive (both IIO-truth) deterministically
SEGFAULT the selector at the quantitative edge's LCR-vs-DM bridge stage -
heap corruption ('invalid permissions', wild address), reproduced across
forked/non-forked execution, both EM engines, and multiple seeds. See
`crash_reproducer.R` for the minimal reproducer and crash signature.

Consequences for the full-grid extension analysis: coverage is 1078/1080;
these two cells are EXCLUDED and reported as crashes, not as verdicts.
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
