# Regression tests for the conjoint-measurement Gibbs samplers and their
# C++ helpers (omni.check, lsqisotonic, ksdensity, ConjointChecks driver).

matlab_fixture <- function(file) {
  test_path("fixtures", "matlab", file)
}

# Rasch-like counts with no 0/1 proportions
make_rasch_counts <- function(nr, per_cell = 60L) {
  ab <- seq(-1.5, 1.5, length.out = nr)
  dif <- seq(-1, 1, length.out = nr)
  p <- 1 / (1 + exp(-outer(ab, dif, "-")))
  N <- matrix(as.integer(per_cell), nr, nr)
  n <- round(p * N)
  n[n == 0] <- 1
  n[n == N] <- N[n == N] - 1
  list(N = N, n = n)
}

# Perline, Wright & Wainer (1979) parole data (Karabatsos 2018, Table 1);
# same data as in the ConjointChecks() example and the MATLAB fixtures.
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

test_that("pure-R omni.check equals the C++ sampler exactly (single and double, 3x3)", {
  d <- make_rasch_counts(3)
  for (chk in c("single", "double")) {
    set.seed(42)
    cpp <- omni.check(d$N, d$n, n.iter = 200, burn = 50, CR = c(.025, .975),
                      check = chk, use_cpp = TRUE)
    set.seed(42)
    rrr <- omni.check(d$N, d$n, n.iter = 200, burn = 50, CR = c(.025, .975),
                      check = chk, use_cpp = FALSE)
    expect_equal(rrr, cpp, tolerance = 0, info = paste("check =", chk))
  }
})

test_that("pure-R omni.check equals the C++ sampler exactly (triple, 4x4)", {
  d <- make_rasch_counts(4)
  set.seed(42)
  cpp <- omni.check(d$N, d$n, n.iter = 200, burn = 50, CR = c(.025, .975),
                    check = "triple", use_cpp = TRUE)
  set.seed(42)
  rrr <- omni.check(d$N, d$n, n.iter = 200, burn = 50, CR = c(.025, .975),
                    check = "triple", use_cpp = FALSE)
  expect_equal(rrr, cpp, tolerance = 0)
})

test_that("omni.check enforces submatrix sizes and chain-length preconditions", {
  d3 <- make_rasch_counts(3)
  d4 <- make_rasch_counts(4)
  expect_error(
    omni.check(d4$N, d4$n, n.iter = 200, burn = 50, CR = c(.025, .975), check = "double"),
    "3x3"
  )
  expect_error(
    omni.check(d3$N, d3$n, n.iter = 200, burn = 50, CR = c(.025, .975), check = "triple"),
    "4x4"
  )
  expect_error(
    omni.check(d3$N, d3$n, n.iter = 100, burn = 100, CR = c(.025, .975), check = "single"),
    "n.iter"
  )
})

test_that("lsqisotonic equals stats::isoreg for unit weights", {
  set.seed(11)
  x <- seq_len(40)
  y <- cumsum(rnorm(40)) / 5 + 0.1 * x
  got <- lsqisotonic(x, y, rep(1, 40))
  expect_equal(got, isoreg(x, y)$yf, tolerance = 1e-12)
})

test_that("lsqisotonic matches a hand-computed weighted case", {
  # PAVA on (1,0,2) with weights (1,3,1): first two pool to (1*1+3*0)/4 = 0.25
  got <- lsqisotonic(c(1, 2, 3), c(1, 0, 2), c(1, 3, 1))
  expect_equal(got, c(0.25, 0.25, 2), tolerance = 1e-12)

  # returns fits in the original (unsorted) order of x
  got2 <- lsqisotonic(c(3, 1, 2), c(2, 1, 0), c(1, 1, 3))
  expect_equal(got2, c(2, 0.25, 0.25), tolerance = 1e-12)
})

test_that("lsqisotonic reproduces the MATLAB isotonic stage on the Perline data", {
  d <- perline_data()
  N_vec <- as.vector(d$N)
  dat <- as.vector(d$n) / N_vec
  # matlab_xhat2.csv holds MATLAB's second-stage GLM fitted values (the x of
  # the isotonic regression); matlab_ty_obs.csv holds MATLAB's isotonic fit.
  xhat2_matlab <- scan(matlab_fixture("matlab_xhat2.csv"), sep = ",", quiet = TRUE)
  ty_matlab <- scan(matlab_fixture("matlab_ty_obs.csv"), quiet = TRUE)
  expect_length(xhat2_matlab, 81)
  ty <- lsqisotonic(xhat2_matlab, dat, N_vec)
  # fixture values are stored to ~5 significant digits
  expect_lt(max(abs(ty - ty_matlab)), 1e-4)
})

test_that("ksdensity matches an R reference implementation (MAD/0.6745 sigma, Scott bandwidth)", {
  ks_ref <- function(data, point) {
    n <- length(data)
    sigma <- median(abs(data - median(data))) / 0.6745
    if (sigma == 0 || is.na(sigma)) {
      r <- max(data) - min(data)
      sigma <- if (r > 0) r else 0
    }
    h <- if (sigma == 0) 1 else sigma * (4 / (3 * n))^0.2
    mean(dnorm((point - data) / h)) / h
  }
  set.seed(21)
  for (n in c(10, 100, 101)) {
    dat <- rbeta(n, 2, 5)
    for (pt in c(0.01, 0.2, 0.5, 0.99)) {
      expect_equal(ksdensity(dat, pt), ks_ref(dat, pt), tolerance = 1e-12)
    }
  }
  # degenerate (constant) data uses the h = 1 fallback in both implementations
  expect_equal(ksdensity(rep(0.3, 50), 0.3), ks_ref(rep(0.3, 50), 0.3), tolerance = 1e-12)
})

test_that("ksdensity bandwidth formula reproduces matlab_bandwidths.csv from the TYSTAR fixture", {
  # matlab_TYSTAR.csv: 100 synthetic isotonic statistics x 81 Perline cells
  # (iteration 1); matlab_bandwidths.csv: MATLAB's per-cell ksdensity bandwidth.
  TYSTAR <- as.matrix(read.csv(matlab_fixture("matlab_TYSTAR.csv"), header = FALSE))
  bw_matlab <- scan(matlab_fixture("matlab_bandwidths.csv"), quiet = TRUE)
  expect_equal(dim(TYSTAR), c(100L, 81L))
  expect_length(bw_matlab, 81)
  mad_col <- apply(TYSTAR, 2, function(x) median(abs(x - median(x))))
  expect_true(all(mad_col > 0))
  h <- (mad_col / 0.6745) * (4 / (3 * nrow(TYSTAR)))^0.2
  # fixture stored to ~5 significant digits
  expect_lt(max(abs(h - bw_matlab) / bw_matlab), 1e-3)
})

test_that("ConjointChecks errors informatively when every sampled submatrix has a 0/1 cell", {
  d <- make_rasch_counts(3, per_cell = 50L)
  d$n[1, 1] <- 0
  set.seed(5)
  expect_error(
    ConjointChecks(d$N, d$n, n.mat = 2, check = "double", adjust_extremes = FALSE),
    "No checkable 3x3 submatrices"
  )
})

test_that("ConjointChecks with adjust_extremes=TRUE runs and does not auto-flag the 0 cell", {
  d <- make_rasch_counts(3, per_cell = 50L)
  d$n[1, 1] <- 0
  set.seed(5)
  out <- ConjointChecks(d$N, d$n, n.mat = 2, check = "double", adjust_extremes = TRUE)
  expect_s4_class(out, "checks")
  # the 0 cell was checked (denominator counted) ...
  expect_true(out@check.counts[1, 1] >= 1)
  # ... and the adjusted observed proportion (0.001) is not automatically a violation
  expect_equal(out@tab[1, 1], 0)
  expect_false(is.na(out@means$weighted))
})

test_that("DoubleCancel/TripleCancel terminate with an informative error when no valid submatrix exists", {
  # regression: these sampling loops previously had no max_tries and hung forever
  d3 <- make_rasch_counts(3, per_cell = 50L)
  d3$n[2, 2] <- 0
  set.seed(1)
  expect_error(DoubleCancel(d3$N, d3$n, n.3mat = 1), "No checkable 3x3 submatrices")

  d4 <- make_rasch_counts(4, per_cell = 50L)
  d4$n[2, 2] <- 0
  set.seed(1)
  expect_error(TripleCancel(d4$N, d4$n, n.4mat = 1), "No checkable 4x4 submatrices")
})

test_that("weighted violation mean excludes never-checked (NA) cells from the denominator", {
  d <- make_rasch_counts(6)
  set.seed(6)
  out <- ConjointChecks(d$N, d$n, n.mat = 1, check = "double")
  checked <- !is.na(out@tab)
  expect_true(any(!checked)) # only one 3x3 out of a 6x6 was checked
  m2_expected <- sum(out@tab * d$N, na.rm = TRUE) / sum(d$N[checked])
  expect_equal(out@means$weighted, m2_expected)
})
