# Tests for compute_se() standard errors (hessian and bootstrap)

# --- helpers ----------------------------------------------------------------

# Generate 2-class LCA data with class-specific item probabilities
gen_lca2 <- function(n, p1, p2, pi1 = 0.5, seed = 1) {
  set.seed(seed)
  J <- length(p1)
  cls <- 1L + rbinom(n, 1, 1 - pi1)          # 1 or 2
  probs <- rbind(p1, p2)[cls, , drop = FALSE] # n x J
  matrix(rbinom(n * J, 1, probs), n, J)
}

# Generate Rasch data
gen_rasch <- function(n, delta, seed = 1) {
  set.seed(seed)
  J <- length(delta)
  theta <- rnorm(n)
  p <- plogis(outer(theta, delta, "-"))
  matrix(rbinom(n * J, 1, p), n, J)
}

# --- input validation -------------------------------------------------------

test_that("compute_se rejects non-qlfit objects", {
  expect_error(compute_se(list(a = 1)), "qlfit")
})

# --- analytic check: hand-rolled Hessian on a tiny UN case ------------------

test_that("hessian SEs match a hand-rolled observed-information computation", {
  skip_if_not_installed("numDeriv")

  J <- 4
  dat <- gen_lca2(300,
                  p1 = c(0.2, 0.25, 0.3, 0.2),
                  p2 = c(0.8, 0.75, 0.7, 0.8),
                  seed = 101)

  fit <- fit_un(dat, n_classes = 2, n_starts = 5, seed = 101)
  se_h <- compute_se(fit, method = "hessian", data = dat)

  expect_s3_class(se_h, "qlse")
  expect_equal(se_h$method, "hessian")

  # Hand-rolled observed-data log-likelihood in the same eta parameterization:
  # eta = (log(pi1/pi2), qlogis(item_probs) column-major)
  ll_hand <- function(eta) {
    cp1 <- plogis(eta[1])
    ip <- matrix(plogis(eta[-1]), J, 2)
    s <- 0
    for (i in seq_len(nrow(dat))) {
      x <- dat[i, ]
      l1 <- prod(ip[, 1]^x * (1 - ip[, 1])^(1 - x))
      l2 <- prod(ip[, 2]^x * (1 - ip[, 2])^(1 - x))
      s <- s + log(cp1 * l1 + (1 - cp1) * l2)
    }
    s
  }

  cp <- fit$class_probs
  ip <- fit$item_probs
  eta_hat <- c(log(cp[1] / cp[2]), qlogis(as.vector(ip)))

  H <- numDeriv::hessian(ll_hand, eta_hat)
  V <- solve(-H)

  # Delta method back to the probability scale
  se_class_hand <- rep(sqrt(V[1, 1]) * cp[1] * (1 - cp[1]), 2)
  se_items_hand <- matrix(sqrt(diag(V)[-1]), J, 2) * ip * (1 - ip)

  expect_equal(unname(se_h$se$class_probs), se_class_hand, tolerance = 1e-3)
  expect_equal(unname(se_h$se$item_probs), unname(se_items_hand),
               tolerance = 1e-3)
})

# --- hessian vs bootstrap agreement -----------------------------------------

test_that("hessian and bootstrap SEs agree within 30% on well-separated data", {
  skip_on_cran()

  dat <- gen_lca2(400,
                  p1 = c(0.15, 0.2, 0.25, 0.2, 0.15),
                  p2 = c(0.85, 0.8, 0.75, 0.8, 0.85),
                  seed = 202)

  fit <- fit_un(dat, n_classes = 2, n_starts = 5, seed = 202)

  se_h <- compute_se(fit, method = "hessian", data = dat)
  se_b <- suppressWarnings(
    compute_se(fit, method = "bootstrap", B = 80, n_starts = 2,
               seed = 202, data = dat)
  )

  expect_equal(se_b$method, "bootstrap")
  expect_gte(se_b$B_effective, 2)

  ratio_items <- se_b$se$item_probs / se_h$se$item_probs
  ratio_class <- se_b$se$class_probs / se_h$se$class_probs

  expect_true(all(is.finite(ratio_items)))
  expect_true(all(abs(ratio_items - 1) < 0.3))
  expect_true(all(abs(ratio_class - 1) < 0.3))
})

# --- calibration against the empirical sampling distribution ----------------

test_that("hessian SEs are calibrated against the empirical SD over sims", {
  skip_on_cran()

  n_sims <- 30
  J <- 4
  p1 <- c(0.2, 0.25, 0.3, 0.2)
  p2 <- c(0.8, 0.75, 0.7, 0.8)
  ip_true <- cbind(p1, p2)
  perms2 <- QuantFit:::all_permutations(2)

  est <- array(NA_real_, c(J, 2, n_sims))
  ses <- array(NA_real_, c(J, 2, n_sims))

  for (s in seq_len(n_sims)) {
    dat <- gen_lca2(300, p1, p2, seed = 300 + s)
    fit <- fit_un(dat, n_classes = 2, n_starts = 3, seed = 300 + s)
    se_s <- suppressWarnings(compute_se(fit, method = "hessian", data = dat))

    pm <- QuantFit:::best_permutation(fit$item_probs, ip_true, perms2)
    est[, , s] <- fit$item_probs[, pm]
    ses[, , s] <- se_s$se$item_probs[, pm]
  }

  emp_sd <- apply(est, c(1, 2), stats::sd)
  mean_se <- apply(ses, c(1, 2), mean, na.rm = TRUE)
  ratio <- mean_se / emp_sd

  # Median calibration ratio should be near 1; individual cells looser
  # (30 sims -> the empirical SD itself has ~13% Monte Carlo error)
  expect_true(median(ratio) > 0.65 && median(ratio) < 1.5)
  expect_true(all(ratio > 0.4 & ratio < 2.5))
})

# --- RM: SEs must match mirt's own ------------------------------------------

test_that("RM standard errors match mirt's own SEs", {
  skip_if_not_installed("mirt")

  dat <- gen_rasch(500, delta = seq(-1, 1, length.out = 5), seed = 404)

  fit <- fit_rm(dat, verbose = FALSE)
  se_rm <- compute_se(fit, method = "hessian")

  expect_named(se_rm$se, "delta")
  expect_equal(names(se_rm$se$delta), names(fit$delta))

  # Reference: mirt fit with SE = TRUE and the same settings as compute_se
  mref <- mirt::mirt(as.data.frame(dat), 1, itemtype = "Rasch",
                     SE = TRUE, quadpts = 61, verbose = FALSE)
  co <- suppressWarnings(suppressMessages(mirt::coef(mref, printSE = TRUE)))
  item_names <- setdiff(names(co), c("GroupPars", "lr.betas"))
  se_ref <- vapply(item_names, function(nm) co[[nm]]["SE", "d"], numeric(1))

  # delta = -d, so SE(delta) = SE(d)
  expect_equal(unname(se_rm$se$delta), unname(se_ref), tolerance = 1e-4)
})

# --- LCR: structure and delta_1 delta-method consistency --------------------

test_that("LCR SEs have the right structure and delta_1 is delta-method
           consistent with the vcov of the free deltas", {
  set.seed(31)
  n <- 500
  theta_t <- c(-1.5, 1.5)
  delta_t <- c(-0.9, -0.3, 0.3, 0.9)
  J <- length(delta_t)
  cls <- sample(1:2, n, replace = TRUE)
  p <- plogis(outer(theta_t, delta_t, "-"))
  dat <- matrix(rbinom(n * J, 1, p[cls, ]), n, J)

  fit <- fit_lcr(dat, n_classes = 2, n_starts = 5, seed = 31)
  se_lcr <- suppressWarnings(compute_se(fit, method = "hessian", data = dat))

  expect_equal(se_lcr$model, "LCR")
  expect_length(se_lcr$se$theta, 2)
  expect_length(se_lcr$se$delta, J)
  expect_length(se_lcr$se$class_probs, 2)
  expect_true(all(is.finite(se_lcr$se$theta) & se_lcr$se$theta > 0))
  expect_true(all(is.finite(se_lcr$se$delta) & se_lcr$se$delta > 0))

  # The deltas are sum-to-zero identified
  expect_equal(sum(fit$delta), 0, tolerance = 1e-6)

  # delta_1 = -sum(delta_2..delta_J)  =>  Var(delta_1) = 1' V_dfree 1
  idx <- se_lcr$eta_index$delta_free
  V_dfree <- se_lcr$vcov_eta[idx, idx, drop = FALSE]
  expect_equal(unname(se_lcr$se$delta[1]), sqrt(sum(V_dfree)),
               tolerance = 1e-10)
  expect_equal(unname(se_lcr$se$delta[-1]), unname(sqrt(diag(V_dfree))),
               tolerance = 1e-10)
})

# --- bootstrap label-switching alignment ------------------------------------

test_that("best_permutation resolves label switching", {
  set.seed(9)
  ref <- matrix(runif(10), 5, 2)

  perms2 <- QuantFit:::all_permutations(2)
  expect_equal(QuantFit:::best_permutation(ref[, 2:1], ref, perms2), c(2L, 1L))
  expect_equal(QuantFit:::best_permutation(ref, ref, perms2), c(1L, 2L))

  ref3 <- matrix(runif(15), 5, 3)
  perms3 <- QuantFit:::all_permutations(3)
  expect_equal(nrow(perms3), 6)
  pm <- QuantFit:::best_permutation(ref3[, c(3, 1, 2)], ref3, perms3)
  expect_equal(ref3[, c(3, 1, 2)][, pm], ref3)
})

test_that("bootstrap SEs are not inflated by label switching", {
  skip_on_cran()

  dat <- gen_lca2(300,
                  p1 = c(0.15, 0.2, 0.25, 0.2),
                  p2 = c(0.85, 0.8, 0.75, 0.8),
                  seed = 505)

  fit <- fit_un(dat, n_classes = 2, n_starts = 5, seed = 505)
  se_b <- suppressWarnings(
    compute_se(fit, method = "bootstrap", B = 25, n_starts = 2,
               seed = 505, data = dat)
  )

  # Classes differ by ~0.6 in every item; unresolved label switching would
  # produce SDs around 0.3, while sampling noise alone is ~0.03
  expect_true(all(se_b$se$item_probs < 0.15))
  expect_true(all(se_b$se$class_probs < 0.15))
})

# --- MON: active constraints produce NA SEs with a warning ------------------

test_that("active MON constraints yield NA SEs and a warning", {
  set.seed(606)
  dat <- matrix(rbinom(250 * 4, 1, 0.5), 250, 4)

  fit <- suppressWarnings(fit_mon(dat, n_classes = 2, n_starts = 3, seed = 606))

  # Force an active constraint: adjacent class probabilities exactly equal
  fit$item_probs[2, 2] <- fit$item_probs[2, 1]

  w <- testthat::capture_warnings(
    se_m <- compute_se(fit, method = "hessian", data = dat)
  )

  expect_true(any(grepl("active", w)))
  expect_gte(se_m$active_constraints$n_active, 1)
  expect_true(all(is.na(se_m$se$item_probs[2, 1:2])))
  # Class probabilities are unaffected by the item-level constraint handling
  expect_true(all(is.finite(se_m$se$class_probs)))
})

# --- print method -----------------------------------------------------------

test_that("print.qlse produces informative output", {
  dat <- gen_lca2(150,
                  p1 = c(0.2, 0.25, 0.3),
                  p2 = c(0.8, 0.75, 0.7),
                  seed = 707)

  fit <- fit_un(dat, n_classes = 2, n_starts = 3, seed = 707)
  se_h <- suppressWarnings(compute_se(fit, method = "hessian", data = dat))

  expect_output(print(se_h), "Standard errors for UN model")
  expect_output(print(se_h), "observed information")
  expect_output(print(se_h), "Class probabilities")
  expect_output(print(se_h), "Item response probabilities")
  printed <- capture.output(vis <- withVisible(print(se_h)))
  expect_false(vis$visible)
  expect_identical(vis$value, se_h)

  # Bootstrap branch of the print method (formatting only)
  se_fake <- se_h
  se_fake$method <- "bootstrap"
  se_fake$B <- 50
  se_fake$B_effective <- 48
  se_fake$n_failed <- 2
  expect_output(print(se_fake), "nonparametric bootstrap")
  expect_output(print(se_fake), "48 of 50", fixed = TRUE)
  expect_output(print(se_fake), "2 bootstrap refit(s) failed", fixed = TRUE)
})
