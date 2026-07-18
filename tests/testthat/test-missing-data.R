# Missing-data support: masked likelihood for the model routes,
# observation-weighted count matrices for the conjoint routes, and
# mask-matched bootstrap nulls. All under MAR.

gen_rasch_na <- function(n, J, rate, seed) {
  set.seed(seed)
  theta <- rnorm(n)
  d <- matrix(rbinom(n * J, 1, plogis(outer(theta, seq(-1.5, 1.5, length.out = J), "-"))), n, J)
  d[matrix(runif(n * J) < rate, n, J)] <- NA
  d
}

test_that("masked validators accept NA and reject bad values", {
  d <- matrix(c(0L, 1L, NA, 1L, 0L, 1L), 3)   # both categories observed per item
  expect_silent(QuantFit:::.validate_poly(d, allow_na = TRUE))
  expect_error(QuantFit:::.validate_poly(d), "complete data")
  expect_error(QuantFit:::.validate_poly(matrix(c(0.5, NA, 1, 0, 1, 0), 3),
                                         allow_na = TRUE), "integer")
})

test_that("binary data with NA routes through the masked poly engine", {
  skip_on_cran()
  d <- gen_rasch_na(600, 6, 0.15, seed = 1)
  f <- fit_un(d, 2, n_starts = 2, seed = 1)
  expect_true(isTRUE(f$polytomous))            # poly-engine representation
  expect_equal(f$n_par, (2 - 1) + 2 * 6)       # m = 1 reduction holds
  expect_true(is.finite(f$loglik))
  # masked C++ == masked R
  f0 <- fit_un(d, 2, n_starts = 2, seed = 1, use_cpp = FALSE)
  expect_lt(abs(f$loglik - f0$loglik), 1e-6)
})

test_that("masked RM matches mirt's native NA handling", {
  skip_on_cran()
  skip_if_not_installed("mirt")
  d <- gen_rasch_na(800, 8, 0.15, seed = 2)
  ours <- fit_rm(d)
  mm <- suppressMessages(mirt::mirt(as.data.frame(d), 1, itemtype = "Rasch",
                                    verbose = FALSE))
  expect_lt(abs(ours$loglik - as.numeric(mirt::extract.mirt(mm, "logLik"))), 0.5)
  # EAP works; ML/WLE refuse informatively (score sufficiency needs complete data)
  sc <- rm_scores(ours, "EAP")
  expect_true(all(is.finite(sc$theta)))
  expect_error(rm_scores(ours, "ML"), "complete")
})

test_that("person_order ladder governs missing-data conditioning", {
  d <- gen_rasch_na(1000, 8, 0.15, seed = 3)
  # default = complete-case: incomplete respondents dropped (with a message),
  # so every score-group row has constant N (the assumption-free frame)
  expect_message(p <- PrepareChecks(d, ss.lower = 10), "incomplete")
  expect_true(all(p$n <= p$N))
  expect_true(all(apply(p$N, 1, function(r) length(unique(r)) == 1)))
  # facility / adjusted keep everyone: cells weighted by observations,
  # so N varies within rows
  for (po in c("facility", "adjusted")) {
    pk <- PrepareChecks(d, ss.lower = 10, person_order = po)
    expect_true(all(pk$n <= pk$N))
    expect_true(any(apply(pk$N, 1, function(r) length(unique(r)) > 1)))
    expect_gte(nrow(pk$N), 3L)
  }
  # complete-data behaviour unchanged: constant N per score group, no message
  dc <- gen_rasch_na(1000, 8, 0, seed = 3)
  expect_silent(pc <- PrepareChecks(dc, ss.lower = 10))
  expect_true(all(apply(pc$N, 1, function(r) length(unique(r)) == 1)))
})

test_that("prepare_polytomous handles missing responses as out-of-play", {
  set.seed(4)
  n <- 1200; J <- 6; M <- 3
  theta <- rnorm(n); dl <- seq(-1, 1, length.out = J); tau <- c(-0.8, 0, 0.8)
  resp <- matrix(0L, n, J)
  for (j in seq_len(J)) {
    num <- cbind(0, t(apply(outer(theta, dl[j] + tau, "-"), 1, cumsum)))
    P <- exp(num) / rowSums(exp(num))
    resp[, j] <- rowSums(runif(n) > t(apply(P, 1, cumsum)))
  }
  resp[matrix(runif(n * J) < 0.15, n, J)] <- NA
  pp <- prepare_polytomous(resp, ss.lower = 10, cell.lower = 3)
  expect_true(all(pp$n <= pp$N))
  expect_false(any(pp$N == 0))
})

test_that("impose_mask rank-matches and is a no-op on complete data", {
  set.seed(5)
  obs <- gen_rasch_na(300, 6, 0.2, seed = 6)
  sim <- matrix(rbinom(300 * 6, 1, 0.5), 300, 6)
  out <- QuantFit:::.impose_mask(sim, obs)
  expect_identical(sum(is.na(out)), sum(is.na(obs)))       # same mask size
  # rank matching: per-person missing counts transfer by score rank
  obs_cnt <- rowSums(is.na(obs))[order(rowSums(obs, na.rm = TRUE))]
  out_cnt <- rowSums(is.na(out))[order(rowSums(sim))]
  expect_identical(obs_cnt, out_cnt)
  # complete observed data -> unchanged
  expect_identical(QuantFit:::.impose_mask(sim, sim), sim)
})

test_that("cc_bootstrap_null runs and calibrates on masked additive data", {
  skip_on_cran()
  d <- gen_rasch_na(1200, 10, 0.15, seed = 7)
  r <- suppressWarnings(cc_bootstrap_null(d, B = 12, n.mat = 12, seed = 3,
                                          verbose = FALSE))
  expect_s3_class(r, "ccnull")
  expect_false(r$reject)                     # additive data, masked pipeline
})
