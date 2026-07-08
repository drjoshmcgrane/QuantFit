# Regression tests for PrepareChecks, in particular the collapse.columns
# bug where tied column sums caused one group to be collapsed twice and
# another to be dropped (cs[i] vs cs.index[i]).

test_that("PrepareChecks produces consistent N and n matrices", {
  set.seed(1)
  resp <- matrix(rbinom(300 * 5, 1, 0.5), 300, 5)
  out <- PrepareChecks(resp, ss.lower = 2)
  expect_equal(dim(out$N), dim(out$n))
  expect_true(all(out$n <= out$N))
  expect_true(all(out$n >= 0))
  # every retained row of N is constant (same number of respondents per score group)
  expect_true(all(apply(out$N, 1, function(x) length(unique(x)) == 1)))
})

test_that("collapse.columns preserves totals when column sums tie", {
  # Build data in which several items have exactly tied column sums:
  # items 1 and 2 are identical, items 3 and 4 are identical.
  set.seed(42)
  n_resp <- 400
  base1 <- rbinom(n_resp, 1, 0.3)
  base2 <- rbinom(n_resp, 1, 0.6)
  other <- rbinom(n_resp, 1, 0.5)
  resp <- cbind(base1, base1, base2, base2, other)

  uncollapsed <- PrepareChecks(resp, ss.lower = 2, collapse.columns = FALSE)
  cs <- colSums(uncollapsed$n)
  expect_true(any(duplicated(cs))) # the regression scenario: tied column sums

  collapsed <- PrepareChecks(resp, ss.lower = 2, collapse.columns = TRUE)

  # totals must be preserved: before the fix, one tied group was collapsed
  # twice and another dropped entirely
  expect_equal(sum(collapsed$n), sum(uncollapsed$n))
  expect_equal(sum(collapsed$N), sum(uncollapsed$N))
  expect_equal(rowSums(collapsed$n), rowSums(uncollapsed$n))
  expect_equal(rowSums(collapsed$N), rowSums(uncollapsed$N))

  # one output column per distinct column sum
  expect_equal(ncol(collapsed$n), length(unique(cs)))

  # each collapsed column is the sum of the tied original columns
  for (lev in unique(cs)) {
    tmp <- uncollapsed$n[, cs == lev, drop = FALSE]
    expect_true(any(apply(collapsed$n, 2, function(col) isTRUE(all.equal(col, rowSums(tmp), check.attributes = FALSE)))),
                info = paste("collapsed group with column sum", lev))
  }
  expect_true(all(collapsed$n <= collapsed$N))
})

test_that("collapse.columns is a no-op partition when all column sums are distinct", {
  set.seed(7)
  resp <- matrix(rbinom(500 * 4, 1, c(0.2, 0.4, 0.6, 0.8)), 500, 4, byrow = TRUE)
  uncollapsed <- PrepareChecks(resp, ss.lower = 2, collapse.columns = FALSE)
  cs <- colSums(uncollapsed$n)
  skip_if(any(duplicated(cs)), "column sums happened to tie")
  collapsed <- PrepareChecks(resp, ss.lower = 2, collapse.columns = TRUE)
  expect_equal(ncol(collapsed$n), ncol(uncollapsed$n))
  expect_equal(sum(collapsed$n), sum(uncollapsed$n))
})
