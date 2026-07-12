# Helpers ---------------------------------------------------------------

# Partial-credit (PCM) data: theta ~ N(0,1), evenly spaced step difficulties.
gen_pcm <- function(nP, J, M, seed) {
  set.seed(seed)
  theta <- rnorm(nP)
  d <- seq(-1, 1, length.out = J); tau <- seq(-0.8, 0.8, length.out = M)
  resp <- matrix(0L, nP, J)
  for (j in seq_len(J)) {
    num <- cbind(0, t(apply(outer(theta, d[j] + tau, "-"), 1, cumsum)))
    P <- exp(num) / rowSums(exp(num))
    resp[, j] <- rowSums(runif(nP) > t(apply(P, 1, cumsum)))
  }
  resp
}

# Two-class polytomous LCA with class 2 stochastically dominating class 1.
gen_mon_poly <- function(nP, J, M, seed) {
  set.seed(seed)
  truth <- lapply(seq_len(J), function(j) {
    lo <- runif(M + 1, .3, 1.7); lo <- lo / sum(lo)
    hi <- lo * seq(0.5, 2, length.out = M + 1); hi <- hi / sum(hi)
    rbind(lo, hi)
  })
  cls <- sample(2, nP, TRUE); resp <- matrix(0L, nP, J)
  for (j in seq_len(J))
    resp[, j] <- sapply(cls, function(c) sample(0:M, 1, prob = truth[[j]][c, ]))
  resp
}

# Tests -----------------------------------------------------------------

test_that("polytomous detection and validation work", {
  expect_true(QuantFit:::.is_polytomous(matrix(c(0, 1, 2, 1), 2)))
  expect_false(QuantFit:::.is_polytomous(matrix(c(0, 1, 1, 0), 2)))
  # non-consecutive categories are rejected
  expect_error(QuantFit:::.validate_poly(matrix(c(0L, 2L, 0L, 2L), 2)),
               "consecutive")
})

test_that("polytomous parameter counts reduce to the dichotomous ones", {
  # every m_j = 1 (binary) -> sumM = J -> matches count_parameters()
  cc1 <- rep(1L, 6)
  expect_equal(QuantFit:::count_parameters_poly("UN", cc1, 3),
               QuantFit:::count_parameters("UN", 6, 3))
  expect_equal(QuantFit:::count_parameters_poly("LCR", cc1, 3),
               QuantFit:::count_parameters("LCR", 6, 3))
  expect_equal(QuantFit:::count_parameters_poly("RM", cc1, 3),
               QuantFit:::count_parameters("RM", 6, 3))
  # polytomous counts use total category steps sumM
  cc <- c(3L, 3L, 2L)  # sumM = 8
  expect_equal(QuantFit:::count_parameters_poly("UN", cc, 2), (2 - 1) + 2 * 8)
  expect_equal(QuantFit:::count_parameters_poly("LCR", cc, 2), 2 * 2 + 8 - 2)
  expect_equal(QuantFit:::count_parameters_poly("RM", cc, 2), 8 + 1)
})

test_that("all six models fit polytomous data with correct structure", {
  skip_on_cran()
  resp <- gen_pcm(1200, 5, 3, seed = 11)
  sumM <- sum(apply(resp, 2, max))
  fu <- fit_un(resp, 2, n_starts = 3, seed = 1)
  expect_true(isTRUE(fu$polytomous))
  expect_equal(fu$n_par, (2 - 1) + 2 * sumM)
  expect_length(fu$item_probs, 5)              # list of category-prob matrices
  expect_equal(nrow(fu$item_probs[[1]]), 2)    # rows index classes
  # rows are proper probability vectors
  expect_true(all(abs(rowSums(fu$item_probs[[1]]) - 1) < 1e-6))

  fl <- fit_lcr(resp, 2, n_starts = 3, seed = 1)
  expect_equal(fl$n_par, 2 * 2 + sumM - 2)
  fr <- fit_rm(resp)
  expect_equal(fr$n_par, sumM + 1)
  expect_true(isTRUE(fr$polytomous))
})

test_that("ordinal nesting holds on constraint-consistent data", {
  skip_on_cran()
  resp <- gen_mon_poly(2000, 5, 3, seed = 3)
  fu <- fit_un(resp, 2, n_starts = 4, seed = 1)
  fm <- fit_mon(resp, 2, n_starts = 4, seed = 1)
  fd <- fit_dm(resp, 2, n_starts = 4, seed = 1)
  # constrained LL cannot exceed the unconstrained LL
  expect_lte(fm$loglik, fu$loglik + 1e-3)
  expect_lte(fd$loglik, fm$loglik + 1e-3)
  # on monotone data the MON constraint barely binds
  expect_lt(fu$loglik - fm$loglik, 5)
  # fitted MON respects stochastic ordering across classes
  ok <- all(vapply(fm$item_probs, function(P) {
    G <- t(apply(P, 1, function(r) rev(cumsum(rev(r)))))
    all(apply(G[, -1, drop = FALSE], 2, diff) >= -1e-4)
  }, logical(1)))
  expect_true(ok)
})

test_that("quantitative models win BIC on PCM data", {
  skip_on_cran()
  resp <- gen_pcm(2000, 6, 3, seed = 7)
  fu <- fit_un(resp, 3, n_starts = 3, seed = 1)
  fl <- fit_lcr(resp, 3, n_starts = 3, seed = 1)
  fr <- fit_rm(resp)
  expect_lt(BIC(fr), BIC(fu))    # Rasch beats unconstrained on Rasch data
  expect_lt(BIC(fr), BIC(fl))    # continuous beats discrete-Rasch
})

test_that("simulate_from_qlfit produces valid polytomous data", {
  skip_on_cran()
  resp <- gen_pcm(800, 4, 3, seed = 5)
  fu <- fit_un(resp, 2, n_starts = 2, seed = 1)
  sim <- QuantFit:::simulate_from_qlfit(fu, 500)
  expect_equal(dim(sim), c(500, 4))
  expect_true(all(sim >= 0 & sim <= 3))
  expect_true(QuantFit:::.is_polytomous(sim))
})
