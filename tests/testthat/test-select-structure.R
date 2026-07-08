test_that("check_mon_constraint detects monotone and non-monotone matrices", {
  mon <- rbind(c(0.2, 0.3), c(0.5, 0.6), c(0.8, 0.9))
  expect_true(check_mon_constraint(mon, tolerance = 0))
  # order-invariant: shuffled class rows are re-ordered by row mean
  expect_true(check_mon_constraint(mon[c(3, 1, 2), ], tolerance = 0))
  anti <- rbind(c(0.8, 0.2), c(0.5, 0.5), c(0.2, 0.8))
  expect_false(check_mon_constraint(anti, tolerance = 0.2))
  expect_true(check_mon_constraint(matrix(0.5, 1, 4)))
})

test_that("check_iio_constraint detects invariant and crossing item orders", {
  iio <- rbind(c(0.2, 0.4, 0.6), c(0.3, 0.5, 0.7))
  expect_true(check_iio_constraint(iio, tolerance = 0))
  crossing <- rbind(c(0.2, 0.8), c(0.8, 0.2))
  expect_false(check_iio_constraint(crossing, tolerance = 0.2))
  # sub-0.05 differences are ignored
  near_tie <- rbind(c(0.50, 0.52), c(0.52, 0.50))
  expect_true(check_iio_constraint(near_tie, tolerance = 0))
})

test_that("check_rasch_constraint accepts additive and rejects non-additive logits", {
  theta <- c(-1, 0, 1); beta <- c(-1, -0.3, 0.4, 1)
  rasch_prob <- plogis(outer(theta, beta, "-"))
  expect_true(check_rasch_constraint(rasch_prob, tolerance = 0.25))
  set.seed(7)
  noise_prob <- matrix(runif(12, 0.1, 0.9), 3, 4)
  expect_false(check_rasch_constraint(noise_prob, tolerance = 0.1))
})

test_that("select_model_constraint returns a qlselect object and selects sensibly", {
  skip_on_cran()
  set.seed(11)
  n <- 400; J <- 8
  theta <- rnorm(n); beta <- seq(-1.5, 1.5, length.out = J)
  rasch_dat <- matrix(rbinom(n * J, 1, plogis(outer(theta, beta, "-"))), n, J)

  sel <- select_model_constraint(rasch_dat, n_classes = 3, n_starts = 3, seed = 1)
  expect_s3_class(sel, "qlselect")
  expect_true(sel$selected %in% c("UN", "MON", "IIO", "DM", "LCR", "RM"))
  # Rasch data should be diagnosed as at least ordinal, plausibly quantitative
  expect_true(sel$constraints$mon)
  expect_output(print(sel), "Selected model")

  # unstructured data: random independent logits per class
  set.seed(12)
  cls <- sample(1:3, n, replace = TRUE)
  p_un <- matrix(runif(3 * J, 0.1, 0.9), 3, J)
  un_dat <- matrix(rbinom(n * J, 1, p_un[cls, ]), n, J)
  sel_un <- select_model_constraint(un_dat, n_classes = 3, n_starts = 3, seed = 1)
  expect_s3_class(sel_un, "qlselect")
  expect_named(sel_un$bics, c("UN", "LCR", "RM"))
})
