# Robustness / calibration of the manifest Silverman-critical-bandwidth DIP axis
# (LCR discrete vs RM continuous) with a fully MANIFEST RM-null bootstrap:
#   item difficulties beta_j = -logit(mean X_j)  (item margins)
#   latent SD sigma          = moment-match to the observed score variance
#                              (numeric root-find; no IRT model fit)
#   null: simulate continuous-theta Rasch at (beta, sigma), recompute h_crit.
# Observed h_crit above the null (1-alpha) quantile => discrete => LCR.
#
# SIZE  : RM data across J x sigma       (false-discrete should be ~ alpha)
# POWER : LCR data across J x n_classes  (want high; identifiability-limited)
suppressMessages(library(QuantFit))

n_modes <- function(x, h) {
  g <- seq(min(x) - 1, max(x) + 1, length.out = 400)
  d <- vapply(g, function(t) mean(stats::dnorm((t - x) / h)) / h, numeric(1))
  sum(diff(sign(diff(d))) < 0)
}
h_crit <- function(x) {                          # Silverman critical bandwidth / sd
  s <- stats::sd(x); if (s == 0) return(0)
  lo <- s / 60; hi <- s * 2.5
  for (i in 1:40) { mid <- (lo + hi) / 2; if (n_modes(x, mid) <= 1) hi <- mid else lo <- mid }
  hi / s
}
# manifest moment estimate of sigma: find sigma so simulated Var(score) matches
sim_scorevar <- function(beta, sigma, N = 4000L) {
  th <- stats::rnorm(N, 0, sigma)
  stats::var(rowSums(matrix(stats::rbinom(N * length(beta), 1,
    stats::plogis(outer(th, beta, "-"))), N)))
}
est_sigma <- function(data, beta) {
  target <- stats::var(rowSums(data))
  lo <- 0.05; hi <- 4
  for (i in 1:25) { mid <- (lo + hi) / 2
    if (sim_scorevar(beta, mid) < target) lo <- mid else hi <- mid }
  (lo + hi) / 2
}
dip_silverman <- function(data, B = 99L, alpha = 0.05, seed = 1) {
  x <- rowSums(data); obs <- h_crit(x)
  p <- pmin(pmax(colMeans(data), 1e-3), 1 - 1e-3); beta <- -stats::qlogis(p)
  sigma <- est_sigma(data, beta)
  n <- nrow(data); J <- ncol(data)
  set.seed(seed)
  null <- vapply(seq_len(B), function(b) {
    th <- stats::rnorm(n, 0, sigma)
    h_crit(rowSums(matrix(stats::rbinom(n * J, 1, stats::plogis(outer(th, beta, "-"))), n)))
  }, numeric(1))
  list(obs = obs, null95 = stats::quantile(null, 0.95, names = FALSE),
       p = (1 + sum(null >= obs)) / (B + 1), sigma = sigma)
}

gen <- function(m, J, N, seed, nclass = 3L, sigma = 1) { set.seed(seed)
  if (m == "RM") { b <- runif(J,-2,2); th <- rnorm(N, 0, sigma)
    d <- matrix(rbinom(N*J,1,plogis(outer(th,b,"-"))),N,J) }
  else { d <- simulate_responses(m, n_persons=N, n_items=J, n_classes=nclass, seed=seed)
    d <- if(is.list(d)) d$data else d }
  storage.mode(d) <- "integer"; d }

# grid: RM (size) across J x sigma ; LCR (power) across J x nclass
gRM  <- expand.grid(model="RM",  J=c(6L,12L,24L), sigma=c(0.5,1,2), nclass=3L, rep=1:6, stringsAsFactors=FALSE)
gLCR <- expand.grid(model="LCR", J=c(6L,12L,24L), sigma=1,          nclass=c(2L,3L,4L), rep=1:6, stringsAsFactors=FALSE)
cases <- rbind(gRM, gLCR); cases$N <- 1500L; cases$id <- seq_len(nrow(cases))
set.seed(55501); cases$seed <- sample.int(.Machine$integer.max, nrow(cases))
out <- Sys.getenv("DIPROB_OUT", "dip_silverman_out"); dir.create(out, showWarnings=FALSE)
run <- function(k) { cs <- cases[k,]; f <- file.path(out, sprintf("s%03d.csv", k))
  if (file.exists(f)) return(invisible())
  d <- gen(cs$model, cs$J, cs$N, cs$seed, cs$nclass, cs$sigma)
  r <- tryCatch(dip_silverman(d, B=99L, seed=cs$id), error=function(e) NULL)
  if (is.null(r)) return(invisible())
  write.csv(data.frame(model=cs$model, J=cs$J, sigma=cs$sigma, nclass=cs$nclass, rep=cs$rep,
    obs=round(r$obs,4), null95=round(r$null95,4), p=round(r$p,4),
    call=if(r$p<=0.05) "LCR" else "RM", sigma_est=round(r$sigma,3)), f, row.names=FALSE) }
cat("dip Silverman robustness:", nrow(cases), "datasets\n")
invisible(parallel::mclapply(seq_len(nrow(cases)),
  function(k) tryCatch(run(k), error=function(e) NULL), mc.cores=8))
cat("DIP SILVERMAN ROBUST DONE\n")
