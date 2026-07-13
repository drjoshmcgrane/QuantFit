# Helpers ---------------------------------------------------------------

# Simulate partial-credit (PCM) data: theta ~ N(0,1), evenly spaced step
# difficulties. Returns an integer person x item matrix scored 0..M.
gen_pcm <- function(nP, J, M, seed) {
  set.seed(seed)
  theta <- rnorm(nP)
  delta <- matrix(seq(-1.2, 1.2, length.out = J), J, M) +
    matrix(seq(-1, 1, length.out = M), J, M, byrow = TRUE)
  resp <- matrix(0L, nP, J)
  for (j in seq_len(J)) {
    num <- cbind(0, t(apply(outer(theta, delta[j, ], "-"), 1, cumsum)))
    P <- exp(num) / rowSums(exp(num))
    resp[, j] <- rowSums(runif(nP) > t(apply(P, 1, cumsum)))
  }
  resp
}

# Tests -----------------------------------------------------------------

test_that("recode_adjacent produces correct adjacent-category coding", {
  # one item, categories 0..3 -> 3 sub-items
  x <- matrix(c(0L, 1L, 2L, 3L), ncol = 1)
  rec <- recode_adjacent(x)
  expect_equal(ncol(rec), 3L)
  expect_equal(attr(rec, "item"), c(1L, 1L, 1L))
  expect_equal(attr(rec, "step"), c(1L, 2L, 3L))
  # step 1 (1 vs 0): person scoring 0 -> 0, 1 -> 1, 2/3 -> NA
  expect_equal(rec[, 1], c(0L, 1L, NA, NA))
  # step 2 (2 vs 1): 1 -> 0, 2 -> 1, else NA
  expect_equal(rec[, 2], c(NA, 0L, 1L, NA))
  # step 3 (3 vs 2): 2 -> 0, 3 -> 1, else NA
  expect_equal(rec[, 3], c(NA, NA, 0L, 1L))
})

test_that("recode_adjacent rejects malformed input and passes NA through", {
  # missing responses are allowed: they recode to NA (out-of-play) on every
  # sub-item, so they simply do not count toward any cell's N
  rec <- recode_adjacent(matrix(c(0L, NA, 1L, 2L), ncol = 1))
  expect_true(all(is.na(rec[2, ])))
  expect_error(recode_adjacent(matrix(c(0, 1.5, 2), ncol = 1)), "integer")
  expect_error(recode_adjacent(matrix(c(0L, 0L, 0L), ncol = 1)),
               "single category")
})

test_that("prepare_polytomous returns a dense, valid N/n matrix", {
  resp <- gen_pcm(1500, 8, 3, seed = 11)
  pp <- prepare_polytomous(resp, ss.lower = 10, cell.lower = 5)
  expect_named(pp, c("N", "n"))
  expect_false(any(pp$N == 0))               # no empty cells
  expect_true(all(pp$n <= pp$N))             # counts are valid
  expect_true(all(pp$N >= 5))                # cell.lower honoured
  expect_gte(nrow(pp$N), 3)
  expect_gte(ncol(pp$N), 3)
  expect_equal(dim(pp$N), dim(pp$n))
  # attributes record the pruning and the surviving sub-item map
  expect_true(is.list(attr(pp, "dropped")))
  expect_equal(nrow(attr(pp, "sub_items")), ncol(pp$N))
})

test_that("prepare_polytomous output feeds ConjointChecks", {
  skip_on_cran()
  resp <- gen_pcm(1500, 8, 3, seed = 12)
  pp <- prepare_polytomous(resp, ss.lower = 10, cell.lower = 5)
  cc <- ConjointChecks(N = pp$N, n = pp$n, check = "double", use_cpp = TRUE)
  expect_s4_class(cc, "checks")
})

test_that("columns are difficulty-ordered when requested", {
  resp <- gen_pcm(1500, 8, 3, seed = 13)
  pp <- prepare_polytomous(resp, ss.lower = 10, cell.lower = 5,
                           order.columns = TRUE)
  prop <- colSums(pp$n) / colSums(pp$N)
  expect_false(is.unsorted(prop))
})

test_that("prepare_polytomous errors when the dense core is too small", {
  # a tiny, sparse data set cannot yield a 3x3 dense block
  resp <- gen_pcm(40, 3, 3, seed = 14)
  expect_error(
    prepare_polytomous(resp, ss.lower = 10, cell.lower = 5),
    "score groups|3x3|Dense-core"
  )
})
