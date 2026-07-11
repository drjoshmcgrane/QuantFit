test_that("assess_quantitative returns a well-formed triangulated verdict", {
  skip_on_cran()
  set.seed(1)
  n <- 500; J <- 8
  dat <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n),
                       seq(-1.8, 1.8, length.out = J), "-"))), n, J)

  v <- assess_quantitative(dat, n_bands = 5, cc_B = 12, kara_S = 2000,
                           kara_N_synth = 20, cc_n_mat = 8, B = 19,
                           mc.cores = 2, seed = 1, verbose = FALSE)

  expect_s3_class(v, "quantverdict")
  expect_true(is.character(v$verdict) && nchar(v$verdict) > 0)
  expect_true(v$support >= 0 && v$support <= 3)
  expect_true(all(c("lc", "cc", "kara") %in% names(v)))
  expect_identical(v$n_bands, 5)

  # each route reports availability and a quantitative-support flag
  for (route in c("lc", "cc", "kara")) {
    expect_true("supports_quant" %in% names(v[[route]]))
  }
  # CC route is bootstrap-null-calibrated (reports a percentile p-value)
  expect_true(isFALSE(v$cc$available) || "double_percentile" %in% names(v$cc))
  # triple ran (default TRUE)
  expect_true(isFALSE(v$cc$available) || "triple_percentile" %in% names(v$cc))
  expect_output(print(v), "triangulated judgement")
  expect_output(print(v), "Student & Read")
})

test_that("assess_quantitative respects triple = FALSE", {
  skip_on_cran()
  set.seed(2)
  n <- 400; J <- 8
  dat <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n),
                       seq(-1.5, 1.5, length.out = J), "-"))), n, J)
  v <- assess_quantitative(dat, n_bands = 5, triple = FALSE, cc_B = 10,
                           kara_S = 1500, kara_N_synth = 15, cc_n_mat = 8,
                           B = 19, mc.cores = 2, seed = 2, verbose = FALSE)
  expect_s3_class(v, "quantverdict")
  if (isTRUE(v$cc$available)) expect_null(v$cc$triple_percentile)
})

test_that("cc_bootstrap_null returns a well-formed ccnull object", {
  skip_on_cran()
  set.seed(3)
  n <- 600; J <- 12
  dat <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n),
                       seq(-2, 2, length.out = J), "-"))), n, J)
  res <- cc_bootstrap_null(dat, B = 15, n.mat = 12, mc.cores = 2, seed = 1,
                           verbose = FALSE)
  expect_s3_class(res, "ccnull")
  expect_true(res$observed >= 0 && res$observed <= 1)
  expect_true(res$percentile >= 0 && res$percentile <= 1)
  expect_equal(res$p_value, 1 - res$percentile, tolerance = 1e-12)
  expect_length(res$null, res$B)
  expect_true(all(diff(res$null) >= 0))       # sorted
  expect_type(res$reject, "logical")
  expect_output(print(res), "Student & Read")

  # parallel matches serial (independent per-replicate seeds)
  skip_on_os("windows")
  a <- cc_bootstrap_null(dat, B = 12, n.mat = 10, mc.cores = 1, seed = 7,
                         verbose = FALSE)
  b <- cc_bootstrap_null(dat, B = 12, n.mat = 10, mc.cores = 2, seed = 7,
                         verbose = FALSE)
  expect_identical(a$null, b$null)
  expect_identical(a$observed, b$observed)
})
