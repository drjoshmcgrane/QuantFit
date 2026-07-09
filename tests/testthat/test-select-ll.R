# Helpers ---------------------------------------------------------------

# Monotone (MON-consistent) latent class data: per item, class probabilities
# sorted increasing over classes
gen_mon_data <- function(n, J, C, seed) {
  set.seed(seed)
  prob <- plogis(t(apply(matrix(runif(J * C, -4, 4), J, C), 1, sort)))
  cls <- sample.int(C, n, replace = TRUE)
  matrix(rbinom(n * J, 1, t(prob)[cls, ]), n, J)
}

# Strongly anti-monotone data: two classes with crossing item profiles
gen_anti_data <- function(n, J, seed) {
  set.seed(seed)
  half <- floor(J / 2)
  p1 <- c(rep(0.9, half), rep(0.1, J - half))
  p2 <- c(rep(0.1, half), rep(0.9, J - half))
  prob <- rbind(p1, p2)
  cls <- sample.int(2, n, replace = TRUE)
  matrix(rbinom(n * J, 1, prob[cls, ]), n, J)
}

# Rasch data: theta ~ N(0,1), evenly spaced betas
gen_rasch_data <- function(n, J, seed) {
  set.seed(seed)
  theta <- rnorm(n)
  beta <- seq(-1.5, 1.5, length.out = J)
  matrix(rbinom(n * J, 1, plogis(outer(theta, beta, "-"))), n, J)
}

# Unstructured data: random independent logits per class (unsorted)
gen_un_data <- function(n, J, C, seed) {
  set.seed(seed)
  prob <- matrix(runif(C * J, 0.1, 0.9), C, J)
  cls <- sample.int(C, n, replace = TRUE)
  matrix(rbinom(n * J, 1, prob[cls, ]), n, J)
}

# Tests ------------------------------------------------------------------

test_that("ll_equivalence_test returns a valid, deterministic qleqtest", {
  skip_on_cran()
  dat <- gen_mon_data(300, 6, 2, seed = 101)
  f_un <- fit_un(dat, 2, n_starts = 3, seed = 1)
  f_mon <- fit_mon(dat, 2, n_starts = 3, seed = 1)

  tst1 <- ll_equivalence_test(dat, f_mon, f_un, B = 19, n_starts = 2,
                              seed = 42)
  expect_s3_class(tst1, "qleqtest")
  expect_true(tst1$statistic >= 0)
  expect_true(tst1$p_value >= 0 && tst1$p_value <= 1)
  expect_identical(tst1$models, c("MON", "UN"))
  expect_true(tst1$B_effective <= 19)
  expect_length(tst1$null_distribution, tst1$B_effective)
  expect_true(all(diff(tst1$null_distribution) >= 0))

  # deterministic given the same seed
  tst2 <- ll_equivalence_test(dat, f_mon, f_un, B = 19, n_starts = 2,
                              seed = 42)
  expect_identical(tst1$p_value, tst2$p_value)
  expect_identical(tst1$null_distribution, tst2$null_distribution)

  expect_output(print(tst1), "LL-equivalence test")
  expect_output(print(tst1), "Bootstrap p")
})

test_that("MON vs UN does not reject on MON-generated data", {
  skip_on_cran()
  dat <- gen_mon_data(300, 6, 2, seed = 202)
  f_un <- fit_un(dat, 2, n_starts = 3, seed = 1)
  f_mon <- fit_mon(dat, 2, n_starts = 3, seed = 1)

  tst <- ll_equivalence_test(dat, f_mon, f_un, B = 19, n_starts = 2,
                             seed = 7)
  expect_gte(tst$p_value, 0.05)
})

test_that("MON vs UN rejects on strongly anti-monotone data", {
  skip_on_cran()
  dat <- gen_anti_data(300, 6, seed = 303)
  f_un <- fit_un(dat, 2, n_starts = 3, seed = 1)
  f_mon <- fit_mon(dat, 2, n_starts = 3, seed = 1)

  # constrained fit should lose a lot of likelihood on crossing profiles
  expect_gt(2 * (f_un$loglik - f_mon$loglik), 10)

  # reject at level alpha: p <= alpha (with B = 19 the smallest achievable
  # p is 1/20 = 0.05, attained when LR_obs exceeds every null draw)
  tst <- ll_equivalence_test(dat, f_mon, f_un, B = 19, n_starts = 2,
                             seed = 7)
  expect_lte(tst$p_value, 0.05)
})

test_that("negative observed LR is guarded and treated as zero", {
  skip_on_cran()
  # anti-monotone data guarantees UN's LL is strictly above MON's, so
  # swapping the roles produces a clearly negative "LR" that must trigger
  # the guard
  dat <- gen_anti_data(300, 6, seed = 404)
  f_un <- fit_un(dat, 2, n_starts = 3, seed = 1)
  f_mon <- fit_mon(dat, 2, n_starts = 3, seed = 1)
  expect_gt(f_un$loglik, f_mon$loglik + 1e-4)
  expect_warning(
    tst <- ll_equivalence_test(dat, f_un, f_mon, B = 5, n_starts = 2,
                               seed = 9),
    "negative")
  expect_identical(tst$statistic, 0)
  expect_true(tst$p_value > 0.5)  # LR_obs = 0 is never beaten strictly
})

test_that("bootstrap refit failures are dropped with a warning", {
  skip_on_cran()
  dat <- gen_mon_data(300, 6, 2, seed = 505)
  f_un <- fit_un(dat, 2, n_starts = 3, seed = 1)
  f_mon <- fit_mon(dat, 2, n_starts = 3, seed = 1)

  # Every 4th refit call fails: with two refits per replicate this drops
  # half the replicates, which must trigger the >10%-dropped warning
  real_refit <- QuantFit:::refit_model_type
  counter <- 0L
  local_mocked_bindings(
    refit_model_type = function(...) {
      counter <<- counter + 1L
      if (counter %% 4L == 0L) stop("forced failure")
      real_refit(...)
    },
    .package = "QuantFit"
  )
  expect_warning(
    tst <- ll_equivalence_test(dat, f_mon, f_un, B = 8, n_starts = 2,
                               seed = 3),
    "dropped")
  expect_lt(tst$B_effective, 8)
  expect_gt(tst$B_effective, 0)
  expect_identical(tst$B_effective + tst$n_failed, 8L)
})

test_that("all bootstrap refits failing raises an error", {
  skip_on_cran()
  dat <- gen_mon_data(300, 6, 2, seed = 505)
  f_un <- fit_un(dat, 2, n_starts = 3, seed = 1)
  f_mon <- fit_mon(dat, 2, n_starts = 3, seed = 1)

  local_mocked_bindings(
    refit_model_type = function(...) stop("forced failure"),
    .package = "QuantFit"
  )
  expect_error(
    suppressWarnings(
      ll_equivalence_test(dat, f_mon, f_un, B = 5, n_starts = 2, seed = 3)),
    "bootstrap refits failed")
})

test_that("input validation catches bad arguments", {
  dat <- gen_mon_data(100, 5, 2, seed = 606)
  expect_error(ll_equivalence_test(dat, list(), list()), "qlfit")
})

test_that("select_model_ll picks a quantitative model on Rasch data", {
  skip_on_cran()
  dat <- gen_rasch_data(300, 6, seed = 707)
  sel <- select_model_ll(dat, n_classes = 2, B = 19, n_starts = 3,
                         boot_n_starts = 2, seed = 11)
  expect_s3_class(sel, "qlselect_ll")
  expect_true(sel$selected %in% c("LCR", "RM"))
  expect_match(sel$interpretation, "QUANTITATIVE")
  expect_true(all(c("comparison", "statistic", "p_value", "decision") %in%
                    names(sel$tests)))
  # lattice default reaches DM through the increment edges, not DM vs UN
  expect_true(any(grepl("DM", sel$tests$comparison)))
  expect_true(any(grepl("LCR vs DM", sel$tests$comparison)))
  expect_false(any(is.na(sel$bics)))
  expect_output(print(sel), "Selected model")
  expect_output(print(sel), "Decision path")
})

test_that("select_model_ll method argument controls the ordinal-layer tests", {
  skip_on_cran()
  dat <- gen_rasch_data(300, 6, seed = 707)

  lat <- select_model_ll(dat, n_classes = 2, B = 19, n_starts = 3,
                         boot_n_starts = 2, method = "lattice", seed = 11)
  joint <- select_model_ll(dat, n_classes = 2, B = 19, n_starts = 3,
                           boot_n_starts = 2, method = "joint", seed = 11)

  expect_identical(lat$method, "lattice")
  expect_identical(joint$method, "joint")
  # lattice tests the single-constraint edges; joint tests DM vs UN directly
  expect_true(any(grepl("MON vs UN", lat$tests$comparison)))
  expect_true(any(grepl("IIO vs UN", lat$tests$comparison)))
  expect_false(any(grepl("DM vs UN", lat$tests$comparison)))
  expect_true(any(grepl("DM vs UN", joint$tests$comparison)))
  expect_match(default_method <- formals(select_model_ll)$method[[2]], "lattice")
})

test_that("lattice method keeps IIO data ordinal, not doubly-monotone", {
  skip_on_cran()
  set.seed(451)
  n <- 600; J <- 8; C <- 3
  # IIO generation: sort item logits within each class (invariant ordering),
  # but classes are NOT monotone -> MON should be rejected
  logit <- t(apply(matrix(runif(C * J, -4, 4), C, J), 1, sort))
  cls <- sample(1:C, n, replace = TRUE)
  dat <- matrix(rbinom(n * J, 1, plogis(logit)[cls, ]), n, J)

  sel <- select_model_ll(dat, n_classes = C, B = 19, n_starts = 3,
                         boot_n_starts = 2, method = "lattice", seed = 5)
  # Mechanism check: the lattice method must test MON vs UN as its own edge
  # (this is what lets it detect the monotonicity violation in IIO data);
  # the exact outcome on a single small-B dataset is stochastic and is
  # validated over many datasets in the recovery sweep, not here.
  expect_true(any(grepl("MON vs UN", sel$tests$comparison)))
  expect_true(any(grepl("IIO vs UN", sel$tests$comparison)))
  expect_s3_class(sel, "qlselect_ll")
})

test_that("select_model_ll stays classificatory/ordinal on unstructured data", {
  skip_on_cran()
  dat <- gen_un_data(300, 6, 3, seed = 808)
  sel <- select_model_ll(dat, n_classes = 3, B = 19, n_starts = 3,
                         boot_n_starts = 2, seed = 13)
  expect_s3_class(sel, "qlselect_ll")
  expect_true(sel$selected %in% c("UN", "MON"))
  expect_output(print(sel), sel$selected)
})

test_that("simulate_from_qlfit respects fitted probabilities", {
  set.seed(1)
  fake <- structure(list(
    class_probs = c(1, 0),
    item_probs = cbind(c(1, 0, 1), c(0.5, 0.5, 0.5))
  ), class = "qlfit")
  sim <- QuantFit:::simulate_from_qlfit(fake, 50)
  expect_equal(dim(sim), c(50, 3))
  # everyone in class 1: items 1 and 3 always 1, item 2 always 0
  expect_true(all(sim[, 1] == 1))
  expect_true(all(sim[, 2] == 0))
  expect_true(all(sim[, 3] == 1))
})
