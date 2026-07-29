# select_model_hybrid and its missing-data path (2026-07-28 audit coverage)

test_that(".impose_mask preserves masks without manufacturing MAR under MCAR", {
  set.seed(11); N <- 4000; J <- 8
  th <- rnorm(N)
  P <- stats::plogis(outer(th, seq(-1.5, 1.5, length.out = J), "-"))
  obs <- matrix(rbinom(N * J, 1, P), N, J)
  obs[matrix(runif(N * J) < 0.25, N, J)] <- NA        # genuine MCAR
  th2 <- rnorm(N)
  P2 <- stats::plogis(outer(th2, seq(-1.5, 1.5, length.out = J), "-"))
  sim <- matrix(rbinom(N * J, 1, P2), N, J)
  m <- QuantFit:::.impose_mask(sim, obs)
  # per-person mask multiset transfers exactly
  expect_identical(sum(is.na(m)), sum(is.na(obs)))
  expect_identical(sort(rowSums(is.na(m))), sort(rowSums(is.na(obs))))
  # MCAR must not become ability-dependent missingness (pre-fix defect:
  # total-based rank matching induced r ~ -0.26 here)
  expect_lt(abs(cor(th2, rowSums(is.na(m)))), 0.06)
})

test_that(".impose_mask preserves the direction of score-dependent missingness", {
  # NOTE: masking here depends on latent theta, which is MNAR relative to the
  # observed responses - deliberately, as a DIRECTION-PRESERVATION property
  # check for the rank matching (the strongest dependence available). The
  # genuinely-MAR anchor-half mechanism is exercised in
  # inst/validation/missing_data_validation.R.
  set.seed(12); N <- 4000; J <- 8
  th <- rnorm(N)
  P <- stats::plogis(outer(th, seq(-1.5, 1.5, length.out = J), "-"))
  obs <- matrix(rbinom(N * J, 1, P), N, J)
  pm <- stats::plogis(-1.1 - 0.9 * th)                # low scorers skip more
  obs[matrix(runif(N * J) < pm, N, J)] <- NA
  th2 <- rnorm(N)
  P2 <- stats::plogis(outer(th2, seq(-1.5, 1.5, length.out = J), "-"))
  sim <- matrix(rbinom(N * J, 1, P2), N, J)
  m <- QuantFit:::.impose_mask(sim, obs)
  # attenuated (observed mean is a noisy ability proxy) but clearly negative
  expect_lt(cor(th2, rowSums(is.na(m))), -0.15)
})

test_that(".hyb_item_probs normalises both engine formats", {
  P <- matrix(runif(12), 4, 3)
  expect_identical(QuantFit:::.hyb_item_probs(list(item_probs = P)), P)
  # masked/polytomous engine: per-item list of class x category matrices
  li <- lapply(1:4, function(j) {
    p1 <- runif(3); cbind(1 - p1, p1)                 # dichotomous K = 2
  })
  M <- QuantFit:::.hyb_item_probs(list(item_probs = li))
  expect_equal(dim(M), c(4L, 3L))
  expect_equal(M[2, ], li[[2]][, 2])                  # expected score = P(X=1)
})

test_that("select_model_hybrid runs on complete and masked data", {
  skip_on_cran()
  d <- simulate_responses("MON", n_persons = 600, n_items = 6, n_classes = 3,
                          seed = 31)
  d <- if (is.list(d)) d$data else d; storage.mode(d) <- "integer"
  r <- select_model_hybrid(d, n_classes = 3L, B = 9, n_starts = 2L,
                           lr_boot_n_starts = 1L, seed = 1, verbose = FALSE)
  expect_s3_class(r, "qlselect_hybrid")
  expect_true(r$selected %in% c("UN", "MON", "IIO", "DM", "LCR", "RM"))
  expect_false(r$quant_edge_failed)
  set.seed(5); dm <- d
  dm[matrix(runif(length(d)) < 0.10, nrow(d))] <- NA
  rm_ <- suppressMessages(
    select_model_hybrid(dm, n_classes = 3L, B = 9, n_starts = 2L,
                        lr_boot_n_starts = 1L, seed = 1, verbose = FALSE))
  expect_true(rm_$selected %in% c("UN", "MON", "IIO", "DM", "LCR", "RM"))
})

test_that("quant edge completes on DM truth and reports its evidence", {
  skip_on_cran()
  # At these smoke settings (N = 600, B = 9) the DM <-> LCR boundary is
  # legitimately ambiguous, so no verdict is asserted - calibration claims
  # live in the validation grids. The unit contract: the edge runs, reports
  # its test object, and the failure flag stays FALSE.
  d <- simulate_responses("DM", n_persons = 600, n_items = 6, n_classes = 3,
                          seed = 77)
  d <- if (is.list(d)) d$data else d; storage.mode(d) <- "integer"
  r <- select_model_hybrid(d, n_classes = 3L, B = 9, n_starts = 2L,
                           lr_boot_n_starts = 1L, seed = 1, verbose = FALSE)
  expect_true(r$selected %in% c("UN", "MON", "IIO", "DM", "LCR", "RM"))
  expect_false(r$quant_edge_failed)
  expect_s3_class(r$quant$lcr_vs_dm, "qleqtest")
  # a quantitative verdict must carry a completed RM stage (or the flag)
  if (r$selected %in% c("LCR", "RM")) {
    expect_false(isTRUE(r$quant$rm_stage_failed))
    expect_false(is.null(r$quant$rm_vs_lcr))
  }
})
