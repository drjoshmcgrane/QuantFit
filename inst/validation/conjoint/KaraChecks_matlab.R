#' @export
KaraChecks_matlab <- function(N, n, matlab_thetas_file, matlab_rstar_file,
                               N_synth = 100, a = 0.5, b = 0.5,
                               mc.cores = 1, verbose = TRUE) {

  if (verbose) cat("Reading MATLAB random values...\n")
  thetas_matlab <- as.matrix(read.csv(matlab_thetas_file, header = FALSE))
  rstar_matlab <- as.matrix(read.csv(matlab_rstar_file, header = FALSE))

  S <- nrow(thetas_matlab)
  IJ <- ncol(thetas_matlab)

  if (verbose) {
    cat("  thetas:", nrow(thetas_matlab), "x", ncol(thetas_matlab), "\n")
    cat("  rstar:", nrow(rstar_matlab), "x", ncol(rstar_matlab), "\n")
    cat("  S =", S, ", N_synth =", N_synth, ", IJ =", IJ, "\n")
  }

  if (is.matrix(N)) {
    nr <- nrow(N)
    nc <- ncol(N)
    ability_values <- seq(-(nr - 1) / 2, (nr - 1) / 2, length.out = nr)
    testscore <- rep(ability_values, nc)
    item <- rep(1:nc, each = nr)
    N_vec <- as.vector(N)
    n_vec <- as.vector(n)
  } else {
    N_vec <- N
    n_vec <- n
    nr <- length(unique(testscore))
    nc <- length(unique(item))
  }

  item_dummies <- model.matrix(~ factor(item) - 1)
  X <- cbind(testscore, item_dummies)

  if (verbose) cat("Computing observed PAVA (ty)...\n")

  dat <- n_vec / N_vec
  fit1 <- suppressWarnings(glm(cbind(n_vec, N_vec - n_vec) ~ X - 1,
                                family = binomial(link = "logit")))
  xhat <- predict(fit1, type = "response")
  fit2 <- suppressWarnings(glm(cbind(n_vec, N_vec - n_vec) ~ xhat,
                                family = binomial(link = "logit")))
  xhat2 <- predict(fit2, type = "response")
  ty <- lsqisotonic(xhat2, dat, N_vec)

  if (verbose) cat("Running importance sampling with MATLAB randoms...\n")

  ws <- matrix(NA, S, IJ)
  TYSTAR <- matrix(NA, N_synth, IJ)

  for (s in 1:S) {
    theta_s <- thetas_matlab[s, ]

    start_idx <- (s - 1) * N_synth + 1

    for (m in 1:N_synth) {
      rstar <- rstar_matlab[start_idx + m - 1, ]

      fit_synth <- tryCatch({
        suppressWarnings(glm(cbind(rstar, N_vec - rstar) ~ X - 1,
                              family = binomial(link = "logit")))
      }, error = function(e) NULL)

      if (!is.null(fit_synth)) {
        xhat_synth <- predict(fit_synth, type = "response")
        TYSTAR[m, ] <- lsqisotonic(xhat_synth, rstar / N_vec, N_vec)
      } else {
        TYSTAR[m, ] <- rstar / N_vec
      }
    }

    for (j in 1:IJ) {
      ws[s, j] <- ksdensity(TYSTAR[, j], ty[j]) + .Machine$double.eps
    }

    if (verbose && s %% 5000 == 0) cat("  Iteration", s, "\n")
  }

  sumw <- colSums(ws)
  W <- ws / matrix(sumw, S, IJ, byrow = TRUE)
  thetaMn <- colSums(thetas_matlab * W)

  thetau <- (a + n_vec) / (a + b + N_vec)

  thetaMn[thetaMn < 1e-10] <- 1e-10
  thetaMn[thetaMn > 1 - 1e-10] <- 1 - 1e-10
  thetau[thetau < 1e-10] <- 1e-10
  thetau[thetau > 1 - 1e-10] <- 1 - 1e-10

  KL <- thetau * log(thetau / thetaMn) + (1 - thetau) * log((1 - thetau) / (1 - thetaMn))
  violations <- KL > 0.01

  if (exists("nr") && exists("nc")) {
    KL_mat <- matrix(KL, nr, nc)
    violations_mat <- matrix(violations, nr, nc)
    thetaMn_mat <- matrix(thetaMn, nr, nc)
    thetau_mat <- matrix(thetau, nr, nc)
  } else {
    KL_mat <- KL
    violations_mat <- violations
    thetaMn_mat <- thetaMn
    thetau_mat <- thetau
  }

  if (verbose) {
    cat("\nKaraChecks (MATLAB replication) complete.\n")
    cat("Global KL:", round(sum(KL), 4), "\n")
    cat("Violations:", sum(violations), "of", IJ, "cells\n")
  }

  list(
    KL = KL_mat,
    global_KL = sum(KL),
    theta_bar = thetaMn_mat,
    theta_0 = thetau_mat,
    violations = violations_mat,
    n_violations = sum(violations)
  )
}
