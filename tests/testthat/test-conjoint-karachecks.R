# Regression tests for KaraChecks (Karabatsos, 2018) against the MATLAB
# ACMtest.m fixtures in tests/testthat/fixtures/matlab/.

matlab_fixture <- function(file) {
  test_path("fixtures", "matlab", file)
}

# Perline, Wright & Wainer (1979) parole data: 9 score groups x 9 items,
# the data set behind the 81-cell MATLAB fixtures.
perline_data <- function() {
  N <- matrix(c(15, 47, 61, 84, 82, 86, 60, 47, 8), 9, 9, byrow = FALSE)
  per <- structure(c(0, 0.06, 0.07, 0.18, 0.13, 0.13, 0.17, 0.17,
    1, 0, 0.04, 0.15, 0.24, 0.33, 0.28, 0.47, 0.85, 1, 0, 0.04, 0.08,
    0.12, 0.3, 0.64, 0.85, 1, 1, 0, 0.19, 0.39, 0.4, 0.51, 0.58,
    0.82, 0.98, 1, 0, 0.06, 0.18, 0.52, 0.73, 0.95, 1, 1, 1, 0,
    0.23, 0.33, 0.51, 0.68, 0.91, 0.93, 1, 1, 0.27, 0.51, 0.61,
    0.64, 0.68, 0.77, 0.9, 1, 1, 0, 0.21, 0.52, 0.68, 0.84, 0.97,
    0.97, 1, 1, 0.73, 0.64, 0.67, 0.7, 0.78, 0.78, 0.9, 1, 1),
    .Dim = c(9L, 9L))
  list(N = N, n = round(per * N))
}

# 4 score groups x 3 items simulated Rasch fixture (cells stored column-major,
# test score varying fastest within item)
simrasch_data <- function() {
  N_vec <- scan(matlab_fixture("simrasch_n.csv"), quiet = TRUE)
  n_vec <- scan(matlab_fixture("simrasch_r.csv"), quiet = TRUE)
  res <- read.csv(matlab_fixture("simrasch_matlab_results.csv"), header = FALSE)
  list(N = matrix(N_vec, 4, 3), n = matrix(n_vec, 4, 3), results = res)
}

test_that("deterministic two-stage GLM + isotonic pipeline matches MATLAB on the Perline data", {
  d <- perline_data()
  nr <- nrow(d$N)
  nc <- ncol(d$N)
  testscore <- rep(seq(-(nr - 1) / 2, (nr - 1) / 2, length.out = nr), nc)
  item <- rep(1:nc, each = nr)
  N_vec <- as.vector(d$N)
  n_vec <- as.vector(d$n)
  X <- cbind(testscore, model.matrix(~ factor(item) - 1))
  dat <- n_vec / N_vec

  fit1 <- suppressWarnings(glm(cbind(n_vec, N_vec - n_vec) ~ X - 1, family = binomial))
  xhat <- unname(predict(fit1, type = "response"))
  fit2 <- suppressWarnings(glm(cbind(n_vec, N_vec - n_vec) ~ xhat, family = binomial))
  xhat2 <- unname(predict(fit2, type = "response"))
  ty <- lsqisotonic(xhat2, dat, N_vec)

  xhat_matlab <- scan(matlab_fixture("matlab_xhat.csv"), sep = ",", quiet = TRUE)
  xhat2_matlab <- scan(matlab_fixture("matlab_xhat2.csv"), sep = ",", quiet = TRUE)
  ty_matlab <- scan(matlab_fixture("matlab_ty_obs.csv"), quiet = TRUE)

  # fixtures stored to ~5 significant digits; observed agreement is ~5e-6
  expect_lt(max(abs(xhat - xhat_matlab)), 1e-4)
  expect_lt(max(abs(xhat2 - xhat2_matlab)), 1e-4)
  expect_lt(max(abs(ty - ty_matlab)), 1e-4)
})

test_that("theta_0 is the exact Beta-Binomial posterior mean and matches MATLAB", {
  d <- simrasch_data()
  # deterministic part of the MATLAB results file: column 4 is theta_0
  thetau <- (0.5 + as.vector(d$n)) / (0.5 + 0.5 + as.vector(d$N))
  expect_lt(max(abs(thetau - d$results$V4)), 1e-4)
})

test_that("KaraChecks on the simrasch fixture: deterministic parts exact, KL correlates with MATLAB", {
  skip_on_cran()
  d <- simrasch_data()

  set.seed(101)
  out <- KaraChecks(d$N, d$n, S = 500, N_synth = 20, mc.cores = 1, verbose = FALSE)

  # deterministic: unrestricted posterior mean
  expect_equal(out$theta_0, matrix((0.5 + d$n) / (1 + d$N), 4, 3), tolerance = 1e-12)

  # deterministic: observed isotonic summary statistic recomputed independently
  nr <- 4; nc <- 3
  testscore <- rep(seq(-(nr - 1) / 2, (nr - 1) / 2, length.out = nr), nc)
  item <- rep(1:nc, each = nr)
  N_vec <- as.vector(d$N)
  n_vec <- as.vector(d$n)
  X <- cbind(testscore, model.matrix(~ factor(item) - 1))
  fit1 <- suppressWarnings(glm(cbind(n_vec, N_vec - n_vec) ~ X - 1, family = binomial))
  fit2 <- suppressWarnings(glm(cbind(n_vec, N_vec - n_vec) ~ unname(predict(fit1, type = "response")),
                               family = binomial))
  ty <- lsqisotonic(unname(predict(fit2, type = "response")), n_vec / N_vec, N_vec)
  expect_equal(out$t_obs, matrix(ty, 4, 3), tolerance = 1e-12)

  # structure
  expect_true(is.matrix(out$KL) && all(dim(out$KL) == c(4, 3)))
  expect_true(is.matrix(out$ZMn) && all(dim(out$ZMn) == c(4, 3)))
  expect_identical(out$violations, out$KL > 0.01)
  expect_identical(out$n_violations, sum(out$violations))
  expect_equal(out$global_KL, sum(out$KL))

  # stochastic: per-cell KL should track MATLAB's (column 13 of results file)
  expect_gt(cor(as.vector(out$KL), d$results$V13), 0.8)
})

test_that("KaraChecks vector input reproduces matrix input given testscore/item", {
  d <- simrasch_data()
  testscore <- rep(seq(-1.5, 1.5, 1), 3)
  item <- rep(1:3, each = 4)

  set.seed(9)
  o_mat <- KaraChecks(d$N, d$n, S = 40, N_synth = 4, mc.cores = 1, verbose = FALSE)
  set.seed(9)
  o_vec <- KaraChecks(as.vector(d$N), as.vector(d$n), S = 40, N_synth = 4,
                      mc.cores = 1, verbose = FALSE,
                      testscore = testscore, item = item)

  expect_identical(o_vec$KL, o_mat$KL)
  expect_identical(o_vec$theta_bar, o_mat$theta_bar)
  expect_identical(o_vec$t_obs, o_mat$t_obs)
  expect_true(is.matrix(o_vec$KL)) # complete grid -> reshaped to 4x3
})

test_that("KaraChecks vector input without testscore/item errors clearly", {
  d <- simrasch_data()
  expect_error(
    KaraChecks(as.vector(d$N), as.vector(d$n), S = 10, verbose = FALSE),
    "testscore.*item|item.*testscore"
  )
  expect_error(
    KaraChecks(as.vector(d$N), as.vector(d$n), S = 10, verbose = FALSE,
               testscore = c(1, 2), item = rep(1:3, each = 4)),
    "same length"
  )
  expect_warning(
    KaraChecks(d$N, d$n, S = 10, N_synth = 2, mc.cores = 1, verbose = FALSE,
               testscore = rep(1, 12), item = rep(1, 12)),
    "ignored"
  )
})

test_that("KaraChecks returns vectors when testscore/item do not form a complete grid", {
  d <- simrasch_data()
  testscore <- rep(seq(-1.5, 1.5, 1), 3)
  item <- rep(1:3, each = 4)
  keep <- -1 # drop one cell -> incomplete grid
  set.seed(9)
  out <- KaraChecks(as.vector(d$N)[keep], as.vector(d$n)[keep], S = 20, N_synth = 3,
                    mc.cores = 1, verbose = FALSE,
                    testscore = testscore[keep], item = item[keep])
  expect_false(is.matrix(out$KL))
  expect_length(out$KL, 11)
  expect_false(is.matrix(out$ZMn))
})
