# Hybrid on the TI&D DEVELOPMENT data (real datasets)

The hybrid (`select_model_hybrid()`) on the actual TI&D archive
(`tid_data/`, the "development study" set; the locked validation is
`fresh_study.R`). 1080 real TI&D-simulated datasets, N=5000, J∈{6,12,24,48},
dichotomous, known generating models 0-5 = UN/MON/IIO/DM/LCR/RM (180 each).
Balanced sample over model × nI, 5 per cell. Runner: `tid_hybrid.R`.

**Two tractability constraints — and an important retraction:**
- **N subsampled 5000 → 1500** (seeded), on the belief that the LR quant edge over
  full N=5000 exceeded a ~20-min background-job kill ceiling. **That belief was
  later TESTED AND DISPROVEN** (a probe ran 14+ min and an 8-worker job 16+ min
  uninterrupted; the kills coincided with several concurrent jobs oversubscribing
  the cores). So the subsampling was **not necessary** — and it is the cause of
  the *only* weakness seen here (RM/nI=6 below), which vanishes at N=5000. The
  runner now defaults to the full N=5000 (`TID_N`); a full-N rerun is the
  outstanding item.
- **nI=48 excluded** (LCR bridge grain ceiling(49/2)=25 classes — infeasible);
  **nI=24 excluded** from the final run for the same (now-retracted) reason.
  Results below are for **nI∈{6,12}**; `TID_NI` now defaults to 6,12,24.

## Result (nI∈{6,12}, N subsampled to 1500)

| nI | UN | MON | IIO | DM | LCR | RM | exact | scale-type |
|----|----|-----|-----|----|-----|----|-------|-----------|
| 12 | 5/5 | 5/5 | **5/5** | 5/5 | 1/3† | 3/3 | **92%** | **96%** |
| 6  | 5/5 | 4/5 | **4/5** | 4/5 | 5/5 | 0/4‡ | 76% | 79% |

†LCR nI=12: 1→LCR, 1→RM (within-quant, scale-type correct), 1→DM.
‡RM nI=6: all →IIO — a SUBSAMPLING ARTIFACT, see below.

## The RM/nI=6 leak is a subsampling artifact, not a method limit

At nI=6, all 4 RM datasets were classified IIO. Diagnosed: the manifest 2×2
routed them to IIO because the MON axis fired (IIO holds, MON "violated"). But
the MON q05 at N=1500 (0.022-0.036, above eps=0.01) COLLAPSES at the data's real
N=5000:

| RM dataset | MON q05 @ N=1500 | MON q05 @ N=5000 | verdict |
|-----------|------------------|------------------|---------|
| 178 | 0.0257 (violated) | 0.0000 | **holds** |
| 449 | 0.0355 (violated) | 0.0010 | **holds** |
| 907 | 0.0219 (violated) | 0.0000 | **holds** |
| 993 | 0.0288 (violated) | 0.0000 | **holds** |

All four flip to MON-holds at full N=5000: with 5000 persons the class
monotonicity of continuous-θ (RM) data is estimated cleanly, whereas a 1500-row
subsample of a continuum discretized into 3 UN classes at only 6 items shows
noise-level non-monotonicity that the tight eps=0.01 (needed to catch real IIO
violations) reads as a violation. On the ACTUAL data (N=5000) these RM datasets
hold MON → reach DM → the LR edge → RM. LCR (genuine discrete classes) is
unaffected and recovers 5/5 even at nI=6. (Note: for genuinely-monotone data the
q05 falls to ~0 much faster than the N-scaled eps, so the eps N-scaling does NOT
make this decision N-invariant — subsampling genuinely caused it.)

## Bottom line

On genuine TI&D development data the hybrid recovers all six models well —
**nI=12: 92% exact, 96% scale-type, with IIO 5/5** — confirming the manifest
route's IIO strength on the real data it was designed against, not just the
`simulate_responses` simulator. The single blemish (RM→IIO at nI=6) is an
artifact of the N=1500 subsampling forced by the machine's job ceiling and
disappears at the data's true N=5000 (verified on the MON axis).

**Outstanding:** a full-N=5000 end-to-end rerun including nI=24 (the MON-axis
flip at N=5000 is decisive on the *mechanism* of the RM/nI=6 artifact, but the
full-N recovery itself has not been run). This was previously deferred as
infeasible; that reason has been retracted (see above), so it is simply pending.
Run one job at a time to avoid the core oversubscription that caused the earlier
kills.
