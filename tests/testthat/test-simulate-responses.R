test_that("simulate_responses returns valid dichotomous data for every model", {
  for (mod in c("UN", "MON", "IIO", "DM", "LCR", "RM")) {
    d <- simulate_responses(mod, n_persons = 200, n_items = 6, n_classes = 3,
                            n_cat = 2, seed = 1)
    expect_equal(dim(d), c(200, 6))
    expect_true(all(d %in% c(0L, 1L)))
    expect_identical(attr(d, "model"), mod)
    expect_false(QuantFit:::.is_polytomous(d))
  }
})

test_that("simulate_responses returns valid polytomous data", {
  d <- simulate_responses("RM", n_persons = 300, n_items = 5, n_cat = 4, seed = 2)
  expect_equal(dim(d), c(300, 5))
  expect_true(all(d >= 0 & d <= 3))
  expect_true(QuantFit:::.is_polytomous(d))
  # consecutive categories 0..3 used (given enough n)
  expect_true(all(0:3 %in% unique(as.vector(d))))
})

test_that("MON data is stochastically ordered across classes; UN is not", {
  # Recover the generating class means and check the ordering property directly
  set.seed(10)
  dm <- simulate_responses("MON", n_persons = 3000, n_items = 8, n_classes = 3,
                          n_cat = 2, seed = 3)
  cls <- attr(dm, "params")$class
  # class-conditional P(X=1) per item, classes ordered as generated
  pm <- t(sapply(1:3, function(c) colMeans(dm[cls == c, , drop = FALSE])))
  # every item's success probability is non-decreasing across the ordered classes
  expect_true(mean(apply(pm, 2, function(col) all(diff(col) >= -0.05))) > 0.9)

  du <- simulate_responses("UN", n_persons = 3000, n_items = 8, n_classes = 3,
                          n_cat = 2, seed = 3)
  clsu <- attr(du, "params")$class
  pu <- t(sapply(1:3, function(c) colMeans(du[clsu == c, , drop = FALSE])))
  # unconstrained data should NOT be uniformly monotone across items
  expect_lt(mean(apply(pu, 2, function(col) all(diff(col) >= -0.05))), 0.9)
})

test_that("IIO data has an invariant item ordering across classes", {
  di <- simulate_responses("IIO", n_persons = 3000, n_items = 8, n_classes = 3,
                          n_cat = 2, seed = 5)
  cls <- attr(di, "params")$class
  pm <- t(sapply(1:3, function(c) colMeans(di[cls == c, , drop = FALSE])))  # C x J
  # the item ordering by success prob should agree across classes (columns
  # already sorted so each class is increasing across items)
  ranks <- t(apply(pm, 1, rank))
  # rank of each item is nearly identical across classes
  expect_true(mean(apply(ranks, 2, function(r) diff(range(r)) <= 1)) > 0.75)
})

test_that("input validation", {
  expect_error(simulate_responses("RM", n_cat = 1), "n_cat")
  expect_error(simulate_responses("XX"), "arg")
})
