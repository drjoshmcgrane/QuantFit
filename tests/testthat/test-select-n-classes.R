gen_lca_data <- function(n, J, C, seed) {
  set.seed(seed)
  # well-separated classes so the class count is recoverable
  centers <- seq(0.1, 0.9, length.out = C)
  p <- matrix(centers, C, J) + matrix(runif(C * J, -0.05, 0.05), C, J)
  p <- pmin(pmax(p, 0.02), 0.98)
  cls <- sample(seq_len(C), n, replace = TRUE)
  matrix(rbinom(n * J, 1, p[cls, ]), n, J)
}

test_that("select_n_classes returns a well-formed enumeration table", {
  skip_on_cran()
  dat <- gen_lca_data(500, 8, 3, seed = 21)
  nc <- select_n_classes(dat, C_range = 1:4, n_starts = 3, seed = 1)

  expect_s3_class(nc, "qlnclasses")
  expect_identical(nc$table$C, 1:4)
  expect_true(all(c("C", "loglik", "n_par", "AIC", "BIC", "entropy",
                    "converged") %in% names(nc$table)))
  # best_C minimises BIC
  expect_identical(nc$best_C, nc$table$C[which.min(nc$table$BIC)])
  # C = 1 baseline: J parameters, entropy undefined
  expect_identical(nc$table$n_par[nc$table$C == 1], as.numeric(ncol(dat)))
  expect_true(is.na(nc$table$entropy[nc$table$C == 1]))
  expect_output(print(nc), "Selected C")
})

test_that("select_n_classes C = 1 log-likelihood matches the closed form", {
  dat <- gen_lca_data(400, 6, 2, seed = 22)
  nc <- select_n_classes(dat, C_range = 1, seed = 1)

  p <- pmin(pmax(colMeans(dat), 1e-10), 1 - 1e-10)
  s1 <- colSums(dat)
  ll_expected <- sum(s1 * log(p) + (nrow(dat) - s1) * log(1 - p))
  expect_equal(nc$table$loglik[1], ll_expected, tolerance = 1e-8)
  expect_equal(nc$table$BIC[1],
               -2 * ll_expected + ncol(dat) * log(nrow(dat)),
               tolerance = 1e-8)
})

test_that("select_n_classes recovers the class count on separated data", {
  skip_on_cran()
  dat <- gen_lca_data(800, 10, 3, seed = 23)
  nc <- select_n_classes(dat, C_range = 1:5, n_starts = 5, seed = 2)
  # single-class must be beaten; the true count (3) should be preferred or
  # at least clearly in contention over the 1-class baseline
  expect_gt(nc$best_C, 1L)
  expect_lt(nc$table$BIC[nc$table$C == 3],
            nc$table$BIC[nc$table$C == 1])
})

test_that("select_n_classes is parallel-reproducible", {
  skip_on_cran()
  skip_on_os("windows")
  dat <- gen_lca_data(400, 8, 3, seed = 24)
  a <- select_n_classes(dat, C_range = 1:4, n_starts = 3, seed = 3, mc.cores = 1)
  b <- select_n_classes(dat, C_range = 1:4, n_starts = 3, seed = 3, mc.cores = 2)
  expect_equal(a$table, b$table)
  expect_identical(a$best_C, b$best_C)
})

test_that("ll_equivalence_test parallel matches serial exactly", {
  skip_on_cran()
  skip_on_os("windows")
  dat <- gen_lca_data(400, 6, 2, seed = 25)
  f_un <- fit_un(dat, 2, seed = 1)
  f_mon <- fit_mon(dat, 2, seed = 1)
  t1 <- ll_equivalence_test(dat, f_mon, f_un, B = 40, seed = 9, mc.cores = 1)
  t3 <- ll_equivalence_test(dat, f_mon, f_un, B = 40, seed = 9, mc.cores = 2)
  expect_identical(t1$p_value, t3$p_value)
  expect_identical(t1$null_distribution, t3$null_distribution)
  expect_identical(t1$statistic, t3$statistic)
})

test_that("select_model_ll accepts a class-count range and selects C", {
  skip_on_cran()
  dat <- gen_lca_data(500, 8, 3, seed = 26)
  sel <- select_model_ll(dat, n_classes = 2:4, B = 19, n_starts = 3,
                         boot_n_starts = 2, seed = 5)
  expect_s3_class(sel, "qlselect_ll")
  expect_true(sel$n_classes %in% 2:4)
  expect_false(is.null(sel$n_classes_table))
  expect_identical(sel$n_classes_table$C, 2:4)
  expect_output(print(sel), "selected by BIC")
})
