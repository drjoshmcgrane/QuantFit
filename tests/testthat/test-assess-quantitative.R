test_that("assess_quantitative returns a well-formed triangulated verdict", {
  skip_on_cran()
  set.seed(1)
  n <- 500; J <- 8
  dat <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n),
                       seq(-1.8, 1.8, length.out = J), "-"))), n, J)

  v <- assess_quantitative(dat, n_bands = 5, kara_S = 2000, kara_N_synth = 20,
                           cc_n_mat = 10, B = 19, mc.cores = 2, seed = 1,
                           verbose = FALSE)

  expect_s3_class(v, "quantverdict")
  expect_true(is.character(v$verdict) && nchar(v$verdict) > 0)
  expect_true(v$support >= 0 && v$support <= 3)
  expect_true(all(c("lc", "cc", "kara") %in% names(v)))
  expect_identical(v$n_bands, 5)

  # each route reports availability and a quantitative-support flag
  for (route in c("lc", "cc", "kara")) {
    expect_true("supports_quant" %in% names(v[[route]]))
  }
  # triple cancellation ran (n_bands >= 4)
  expect_true(isFALSE(v$cc$available) || "triple_violation" %in% names(v$cc))
  expect_output(print(v), "triangulated judgement")
  expect_output(print(v), "ability bands")
})

test_that("assess_quantitative respects triple = FALSE", {
  skip_on_cran()
  set.seed(2)
  n <- 400; J <- 8
  dat <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n),
                       seq(-1.5, 1.5, length.out = J), "-"))), n, J)
  v <- assess_quantitative(dat, n_bands = 5, triple = FALSE, kara_S = 1500,
                           kara_N_synth = 15, cc_n_mat = 8, B = 19,
                           mc.cores = 2, seed = 2, verbose = FALSE)
  expect_s3_class(v, "quantverdict")
  if (isTRUE(v$cc$available)) expect_null(v$cc$triple_violation)
})
