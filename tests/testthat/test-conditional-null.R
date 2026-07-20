# Conditional-CML null generator: correctness and routing (external review).
test_that("CML recovers difficulties; conditional draws preserve footprints", {
  set.seed(5); th <- rnorm(1200); b <- seq(-1.5, 1.5, length.out = 10)
  d <- matrix(rbinom(1200 * 10, 1, plogis(outer(th, b, "-"))), 1200, 10)
  storage.mode(d) <- "integer"
  dm <- d; dm[sample(length(dm), 0.1 * length(dm))] <- NA
  dl <- QuantFit:::.cml_fit_general(dm)
  expect_gt(cor(vapply(dl, sum, 1), b), 0.99)
  nd <- QuantFit:::.conditional_null_general(dm, dl)
  expect_identical(is.na(nd), is.na(dm))
  expect_identical(rowSums(nd, na.rm = TRUE), rowSums(dm, na.rm = TRUE))
})

test_that("extreme scorers carry no information (margins consistency)", {
  set.seed(5); th <- rnorm(800); b <- seq(-1, 1, length.out = 8)
  d <- matrix(rbinom(800 * 8, 1, plogis(outer(th, b, "-"))), 800, 8)
  storage.mode(d) <- "integer"
  d[sample(length(d), 0.1 * length(d))] <- NA
  dl0 <- QuantFit:::.cml_fit_general(d)
  d2 <- rbind(d, matrix(0L, 30, 8))         # deterministic respondents
  dl1 <- QuantFit:::.cml_fit_general(d2)
  expect_lt(max(abs(unlist(dl0) - unlist(dl1))), 0.02)
})

test_that("null_method routing: hierarchy forwards; poly uses conditional; no MML dependency", {
  skip_on_cran()
  set.seed(6); n <- 500; J <- 6
  d <- matrix(rbinom(n * J, 1, plogis(outer(rnorm(n), seq(-1, 1, length.out = J), "-"))), n, J)
  storage.mode(d) <- "integer"
  env <- asNamespace("QuantFit")
  # hierarchy forwards empirical_mml (CML must NOT be called)
  called <- FALSE
  orig <- get(".cml_fit_general", env)
  unlockBinding(".cml_fit_general", env)
  assign(".cml_fit_general", function(...) { called <<- TRUE; orig(...) }, env)
  h <- suppressWarnings(cc_bootstrap_hierarchy(d, levels = "single", B = 9,
        n.mat = 10, null_method = "empirical_mml", seed = 1, verbose = FALSE))
  assign(".cml_fit_general", orig, env); lockBinding(".cml_fit_general", env)
  expect_false(called)
  # conditional path independent of MML (fit_rm mocked to fail)
  origf <- get("fit_rm", env)
  unlockBinding("fit_rm", env)
  assign("fit_rm", function(...) stop("mml down"), env)
  r <- tryCatch(suppressWarnings(cc_bootstrap_null(d, B = 9, n.mat = 10,
        seed = 1, verbose = FALSE)), error = function(e) NULL)
  assign("fit_rm", origf, env); lockBinding("fit_rm", env)
  expect_false(is.null(r))
  expect_identical(r$reject, r$p_value <= r$alpha)
  # polytomous routes through the conditional generator
  resp <- matrix(sample(0:3, n * J, replace = TRUE), n, J); storage.mode(resp) <- "integer"
  called2 <- FALSE
  orign <- get(".conditional_null_general", env)
  unlockBinding(".conditional_null_general", env)
  assign(".conditional_null_general",
         function(data, dl) { called2 <<- TRUE; orign(data, dl) }, env)
  rp <- tryCatch(suppressWarnings(cc_bootstrap_null(resp, B = 5, n.mat = 8,
        seed = 1, verbose = FALSE)), error = function(e) NULL)
  assign(".conditional_null_general", orign, env)
  lockBinding(".conditional_null_general", env)
  expect_true(called2)
})

test_that("CML objective is consistent on non-centred random difficulties", {
  # regression: the clipped-denominator seam drove the objective negative and
  # failed CML on ordinary data (UN/MON/IIO/DM/RM all erred pre-fix)
  set.seed(41); th <- rnorm(1500); b <- runif(10, -2, 2) + 0.7
  d <- matrix(rbinom(1500 * 10, 1, plogis(outer(th, b, "-"))), 1500, 10)
  storage.mode(d) <- "integer"
  dl <- QuantFit:::.cml_fit_general(d)          # must not error
  est <- vapply(dl, sum, numeric(1))
  expect_gt(cor(est, b), 0.99)
  expect_lt(max(abs(est)), 20)                  # non-boundary
  expect_lt(abs(sum(unlist(dl))), 1e-6)         # exact sum-to-zero
})

test_that("conditional CC completes on the previously failing UN design", {
  skip_on_cran()
  u <- simulate_responses("UN", n_persons = 1200, n_items = 10,
                          n_classes = 3, seed = 77)
  u <- if (is.list(u)) u$data else u; storage.mode(u) <- "integer"
  r <- suppressWarnings(cc_bootstrap_null(u, B = 19, n.mat = 30,
                                          seed = 1, verbose = FALSE))
  expect_s3_class(r, "ccnull")
  expect_gt(r$B, 15)                            # replicates actually ran
})
