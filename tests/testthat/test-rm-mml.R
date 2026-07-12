# The Rasch / partial-credit model is estimated by the package's own marginal
# maximum-likelihood engine (Gauss-Hermite quadrature EM). These tests confirm
# it reproduces mirt (when installed) and that the C++ and R engines agree.

gen_rasch <- function(nP, J, seed) {
  set.seed(seed)
  matrix(rbinom(nP * J, 1, plogis(outer(rnorm(nP), seq(-1.5, 1.5, length.out = J), "-"))),
         nP, J)
}
gen_pcm <- function(nP, J, M, seed) {
  set.seed(seed)
  theta <- rnorm(nP) * 1.1; d <- seq(-1, 1, length.out = J)
  tau <- seq(-0.8, 0.8, length.out = M); resp <- matrix(0L, nP, J)
  for (j in seq_len(J)) {
    num <- cbind(0, t(apply(outer(theta, d[j] + tau, "-"), 1, cumsum)))
    P <- exp(num) / rowSums(exp(num))
    resp[, j] <- rowSums(runif(nP) > t(apply(P, 1, cumsum)))
  }
  resp
}

test_that("our Rasch/PCM MML matches mirt (LL, sigma)", {
  skip_on_cran()
  skip_if_not_installed("mirt")
  for (dat in list(binary = gen_rasch(2000, 6, 11),
                   poly   = gen_pcm(2000, 6, 3, 12))) {
    ours <- QuantFit:::em_rasch_mml(dat, n_quad = 61, seed = 1)
    mf <- suppressMessages(mirt::mirt(as.data.frame(dat), 1,
                                      itemtype = "Rasch", verbose = FALSE))
    expect_lt(abs(ours$loglik - as.numeric(mirt::extract.mirt(mf, "logLik"))), 0.5)
    expect_lt(abs(ours$sigma - sqrt(mirt::coef(mf, simplify = TRUE)$cov[1])), 0.02)
  }
})

test_that("fit_rm works without mirt for binary and polytomous data", {
  skip_on_cran()
  fb <- fit_rm(gen_rasch(1200, 6, 3))
  expect_equal(fb$model_type, "RM")
  expect_equal(fb$n_par, 6 + 1)                 # J + 1
  expect_false(isTRUE(fb$polytomous))
  fp <- fit_rm(gen_pcm(1200, 6, 3, 4))
  expect_true(isTRUE(fp$polytomous))
  expect_equal(fp$n_par, 6 * 3 + 1)             # sumM + 1
  # scores + fit stats compute without error
  sc <- rm_scores(fp, "EAP")
  expect_equal(nrow(sc), 1200)
  expect_true(all(is.finite(sc$theta)))
  expect_s3_class(rm_itemfit(fp), "data.frame")
})

test_that("rm_scores EAP and ML are sensible and monotone in score", {
  skip_on_cran()
  resp <- gen_rasch(1500, 8, 7)
  fit <- fit_rm(resp)
  eap <- rm_scores(fit, "EAP")$theta
  ml  <- rm_scores(fit, "ML")$theta
  score <- rowSums(resp)
  # both increase with the total score
  expect_gt(cor(eap, score), 0.95)
  expect_gt(cor(ml[is.finite(ml)], score[is.finite(ml)]), 0.95)
})

test_that("polytomous C++ and R engines agree to machine precision", {
  skip_on_cran()
  resp <- gen_pcm(600, 5, 3, 21); storage.mode(resp) <- "integer"
  C <- 3
  ip <- lapply(rep(3L, 5), function(m) { P <- matrix(runif(C * (m + 1), .5, 1.5), C, m + 1); P / rowSums(P) })
  cp <- c(.5, .3, .2)
  e1 <- QuantFit:::poly_estep(resp, ip, cp, use_cpp = TRUE)
  e0 <- QuantFit:::poly_estep(resp, ip, cp, use_cpp = FALSE)
  expect_lt(abs(e1$loglik - e0$loglik), 1e-8)
  expect_lt(max(abs(e1$posteriors - e0$posteriors)), 1e-12)
  f1 <- fit_un(resp, 2, n_starts = 2, seed = 1, use_cpp = TRUE)
  f0 <- fit_un(resp, 2, n_starts = 2, seed = 1, use_cpp = FALSE)
  expect_lt(abs(f1$loglik - f0$loglik), 1e-6)
})
