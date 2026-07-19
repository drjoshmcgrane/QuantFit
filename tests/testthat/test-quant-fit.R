test_that("quant_fit returns a well-formed triangulated verdict", {
  skip_on_cran()
  set.seed(1)
  n <- 500; J <- 8
  dat <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n),
                       seq(-1.8, 1.8, length.out = J), "-"))), n, J)

  # N < 1000 trips the CC under-power warning; not the object under test here
  v <- suppressWarnings(quant_fit(dat, n_bands = 5, cc_B = 10, omni_B = 8,
                           omni_S = 1500, omni_N_synth = 15, cc_n_mat = 8,
                           B = 19, mc.cores = 2, seed = 1, verbose = FALSE))

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

test_that("quant_fit respects triple = FALSE", {
  skip_on_cran()
  set.seed(2)
  n <- 400; J <- 8
  dat <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n),
                       seq(-1.5, 1.5, length.out = J), "-"))), n, J)
  v <- suppressWarnings(quant_fit(dat, n_bands = 5, triple = FALSE, cc_B = 8,
                           omni_B = 6, omni_S = 1200, omni_N_synth = 12,
                           cc_n_mat = 8, B = 19, mc.cores = 2, seed = 2,
                           verbose = FALSE))
  expect_s3_class(v, "quantverdict")
  if (isTRUE(v$cc$available)) expect_null(v$cc$triple_percentile)
})

test_that("kara_bootstrap_null returns a well-formed ccnull object", {
  skip_on_cran()
  set.seed(4)
  n <- 800; J <- 12
  dat <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n),
                       seq(-2, 2, length.out = J), "-"))), n, J)
  res <- kara_bootstrap_null(dat, n_bands = 5, B = 8, S = 1500, N_synth = 15,
                             mc.cores = 2, seed = 1, verbose = FALSE)
  expect_s3_class(res, "ccnull")
  expect_identical(res$check, "omni-KL")
  expect_true(res$observed >= 0)
  expect_true(res$percentile >= 0 && res$percentile <= 1)
  expect_true(all(c("kl_median", "kl_max") %in% names(res)))
  expect_output(print(res), "Omnibus cancellation-hierarchy")
})

test_that("assess_quantitative() is a deprecated alias for quant_fit()", {
  skip_on_cran()
  set.seed(5)
  n <- 300; J <- 6
  dat <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n),
                       seq(-1.5, 1.5, length.out = J), "-"))), n, J)
  # Run the (expensive) call once: record whether the deprecation warning fired
  # and muffle every warning (incl. the N < 1000 under-power note) so we can also
  # capture the return value, which expect_warning() would not hand back.
  saw_deprecation <- FALSE
  v <- withCallingHandlers(
    assess_quantitative(dat, n_bands = 5, cc_B = 6, omni_B = 4,
                        omni_S = 1000, omni_N_synth = 10, cc_n_mat = 6,
                        B = 15, mc.cores = 2, seed = 1, verbose = FALSE),
    warning = function(w) {
      if (grepl("deprecat", conditionMessage(w))) saw_deprecation <<- TRUE
      invokeRestart("muffleWarning")
    })
  expect_true(saw_deprecation)
  expect_s3_class(v, "quantverdict")
})

test_that("rm_vs_lcr_test prefers RM on Rasch data", {
  skip_on_cran()
  set.seed(6)
  n <- 800; J <- 12; C <- 3
  dat <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n),
                       seq(-2, 2, length.out = J), "-"))), n, J)
  rm_fit  <- fit_rm(dat, verbose = FALSE)
  lcr_fit <- fit_lcr(dat, n_classes = C, n_starts = 2, verbose = FALSE)
  res <- QuantFit:::rm_vs_lcr_test(dat, rm_fit, lcr_fit, n_classes = C,
                                   B = 20, mc.cores = 2, seed = 1)
  expect_true(all(c("statistic", "p_value", "available", "select_lcr") %in%
                    names(res)))
  expect_true(res$available)
  expect_true(res$p_value > 0 && res$p_value <= 1)
  # Rasch-generated data -> continuous RM preferred, discreteness not supported
  expect_false(res$select_lcr)
})

test_that("cc_bootstrap_null returns a well-formed ccnull object", {
  skip_on_cran()
  set.seed(3)
  n <- 600; J <- 12
  dat <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n),
                       seq(-2, 2, length.out = J), "-"))), n, J)
  # N < 1000 trips the under-power warning; not the object under test here
  res <- suppressWarnings(cc_bootstrap_null(dat, B = 15, n.mat = 12,
                                            mc.cores = 2, seed = 1,
                                            verbose = FALSE))
  expect_s3_class(res, "ccnull")
  expect_true(res$observed >= 0 && res$observed <= 1)
  expect_true(res$percentile >= 0 && res$percentile <= 1)
  # p_value uses the (1 + #{null >= obs}) / (B + 1) continuity correction, so it
  # is strictly positive and matches the exact formula (no longer 1 - percentile)
  expect_true(res$p_value > 0 && res$p_value <= 1)
  expect_equal(res$p_value,
               (1 + sum(res$null >= res$observed)) / (length(res$null) + 1),
               tolerance = 1e-12)
  # res$B counts successful null draws; B + n_failed == requested B (15)
  expect_equal(res$B + res$n_failed, 15)
  expect_length(res$null, res$B)
  expect_true(all(diff(res$null) >= 0))       # sorted
  expect_type(res$reject, "logical")
  expect_output(print(res), "Student & Read")

  # parallel matches serial (independent per-replicate seeds)
  skip_on_os("windows")
  a <- suppressWarnings(cc_bootstrap_null(dat, B = 12, n.mat = 10, mc.cores = 1,
                                          seed = 7, verbose = FALSE))
  b <- suppressWarnings(cc_bootstrap_null(dat, B = 12, n.mat = 10, mc.cores = 2,
                                          seed = 7, verbose = FALSE))
  expect_identical(a$null, b$null)
  expect_identical(a$observed, b$observed)
})

test_that("cc_bootstrap_hierarchy attributes failure to the first rejecting level", {
  skip_on_cran()
  set.seed(1)
  # Rasch data: no level should reject
  n <- 1000; J <- 10
  r <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n), seq(-1.5, 1.5, length.out = J), "-"))), n, J)
  h <- suppressWarnings(cc_bootstrap_hierarchy(r, B = 12, n.mat = 12, seed = 3,
                                               verbose = FALSE))
  expect_s3_class(h, "cchier")
  expect_identical(h$attribution, "none")
  expect_true(h$supports_quant)
  expect_named(h$levels, c("single", "double", "triple"))

  # unstructured LCA data: fails, and at the single level (ordering breaks)
  u <- simulate_responses("UN", n_persons = 1000, n_items = 10, n_classes = 3,
                          seed = 5)
  h2 <- suppressWarnings(cc_bootstrap_hierarchy(u, B = 12, n.mat = 12, seed = 3,
                                                verbose = FALSE))
  expect_false(h2$supports_quant)
  expect_identical(h2$attribution, "single")
  # sequential early stop: deeper levels not run after a rejection
  expect_identical(names(h2$levels), "single")
  expect_output(print(h2), "single-cancellation")
})
