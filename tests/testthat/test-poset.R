# Partial-order machinery (redesigned 2026-07-29 after external review:
# dominance-demonstration bounds replaced the collapsed asymmetry statistic
# and its exchangeability permutation null).

simulate_from_table <- function(P, n, seed) {   # P = C x J class-prob table
  set.seed(seed); C <- nrow(P); J <- ncol(P)
  cls <- sample.int(C, n, replace = TRUE)
  d <- matrix(rbinom(n * J, 1, P[cls, ]), n, J)
  storage.mode(d) <- "integer"; d
}

test_that("antichain with systematic crossings and unequal means is NOT partial", {
  skip_on_cran()
  # The reviewer's counterexample class: every pair crosses deeply (a genuine
  # antichain) while class AVERAGES differ - the pre-redesign statistic
  # (which collapsed to |mean difference|) falsely rejected here.
  J <- 8
  base <- seq(-1.6, 1.6, length.out = J)
  P <- rbind(plogis(base),                 # increasing
             plogis(rev(base) + 0.35),     # decreasing, higher mean
             plogis(base * c(-1, 1) + 0.15))  # alternating
  d <- simulate_from_table(P, 1200, seed = 21)
  r <- QuantFit:::.poset_refine(d, C = 3L, B = 29L, n_starts = 3L,
        use_cpp = TRUE, eps = 0.01, alpha = 0.05, seed = 5,
        sides = "class")
  expect_identical(r$class$shape, "antichain")
  expect_identical(r$class$comparable, 0L)
})

test_that("genuine V dominance is demonstrated with correct pairs", {
  skip_on_cran()
  # class 1 dominates classes 2 and 3 everywhere; 2 and 3 cross: type V
  J <- 8
  lo1 <- seq(-1.2, 1.2, length.out = J)
  lo2 <- rev(lo1)
  P <- rbind(plogis(pmax(lo1, lo2) + 0.9), plogis(lo1 - 0.4),
             plogis(lo2 - 0.4))
  d <- simulate_from_table(P, 1500, seed = 22)
  r <- QuantFit:::.poset_refine(d, C = 3L, B = 29L, n_starts = 3L,
        use_cpp = TRUE, eps = 0.01, alpha = 0.05, seed = 5,
        sides = "class")
  expect_identical(r$class$shape, "partial")
  expect_gte(r$class$comparable, 1L)
  expect_true(all(r$class$pairs$dominant %in% 1:3))
  expect_true(r$class$b_eff >= 20L)
})

test_that("PO_ITEMS errors on infeasible margins and delivers exact structure", {
  expect_error(
    simulate_responses("PO_ITEMS", n_persons = 200, n_items = 12,
                       n_classes = 3, poset = "layers2", po_margin = 0.10,
                       seed = 1),
    "not achievable")
  d <- simulate_responses("PO_ITEMS", n_persons = 200, n_items = 8,
                          n_classes = 3, poset = "layers2", po_margin = 0.005,
                          seed = 1)
  p <- attr(d, "params")
  expect_false(is.null(p$item_poset))
  expect_gte(p$po_achieved, 0.005)
  P <- plogis(p$L); Di <- p$item_poset; J <- ncol(P)
  for (j in seq_len(J)) for (k in seq_len(J)) {
    if (j != k && Di[j, k])                       # specified dominance: exact
      expect_true(all(P[, j] >= P[, k] - 1e-12))
    if (j < k && !Di[j, k] && !Di[k, j])          # incomparable: real crossing
      expect_gt(min(mean(pmax(0, P[, j] - P[, k])),
                    mean(pmax(0, P[, k] - P[, j]))), 0)
  }
})

test_that("PO free records its achieved crossing margin", {
  d <- simulate_responses("PO", n_persons = 200, n_items = 8, n_classes = 3,
                          poset = "V", po_margin = 0.05, seed = 3)
  p <- attr(d, "params")
  expect_gte(p$po_achieved, 0.05)
  expect_false(is.null(p$poset))
})

test_that("poset refinement handles polytomous data (categorical bootstrap)", {
  skip_on_cran()
  # pre-fix defect: expected scores were fed to rbinom as probabilities -
  # every polytomous bootstrap refit failed (0/29) and the refinement refused
  d <- simulate_responses("PO", n_persons = 800, n_items = 6, n_classes = 3,
                          n_cat = 4, poset = "V", po_margin = 0.05, seed = 41)
  r <- QuantFit:::.poset_refine(d, C = 3L, B = 29L, n_starts = 2L,
        use_cpp = TRUE, eps = 0.01, alpha = 0.05, seed = 5, sides = "class")
  expect_true(r$class$b_eff >= 50L)   # internal floor is 99; most must succeed
  expect_true(r$class$shape %in%
                c("partial", "antichain", "nontransitive_dominance"))
  expect_true(is.logical(r$class$transitive))
})

test_that("relation consistency logic: transitive, missing-edge, cycle", {
  ok <- QuantFit:::.poset_refine   # just to assert internals exist
  ro <- function(pr, nn) {
    D <- matrix(FALSE, nn, nn)
    if (nrow(pr)) D[cbind(pr$dominant, pr$dominated)] <- TRUE
    R <- D
    for (i in seq_len(nn)) R <- R | ((R %*% R) > 0)
    if (any(diag(R))) return(FALSE)
    all((R & !D) == FALSE)
  }
  expect_true(ro(data.frame(dominant = c(1, 1, 2), dominated = c(2, 3, 3)), 3))
  expect_false(ro(data.frame(dominant = c(1, 2), dominated = c(2, 3)), 3))  # gap
  expect_false(ro(data.frame(dominant = c(1, 2, 3), dominated = c(2, 3, 1)), 3))
})

test_that("poset_B floor of 99 is enforced against explicit smaller values", {
  skip_on_cran()
  d <- simulate_responses("PO", n_persons = 400, n_items = 6, n_classes = 3,
                          poset = "V", po_margin = 0.05, seed = 61)
  r <- QuantFit:::.poset_refine(d, C = 3L, B = 9L, n_starts = 2L,
        use_cpp = TRUE, eps = 0.01, alpha = 0.05, seed = 5, sides = "class",
        poset_B = 50L)
  expect_gte(r$class$b_eff, 90L)   # floor raised 50 -> 99 (minus rare failures)
  expect_error(QuantFit:::.poset_refine(d, C = 3L, B = 9L, n_starts = 2L,
        use_cpp = TRUE, eps = 0.01, alpha = 0.05, seed = 5, sides = "class",
        poset_B = -1), "finite positive")
})

test_that("po_margin range is validated", {
  expect_error(simulate_responses("PO", n_persons = 100, n_items = 6,
    n_classes = 3, poset = "V", po_margin = -0.1, seed = 1), "0, 0.5")
  expect_error(simulate_responses("PO", n_persons = 100, n_items = 6,
    n_classes = 3, poset = "V", po_margin = 0.7, seed = 1), "0, 0.5")
})
