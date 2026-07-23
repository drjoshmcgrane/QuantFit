# Corrected-build robustness for the IIO and MON axes: SIZE and POWER across
# N x J x class-separation. MON uses the N-scaled eps (eps ~ 1/sqrt(N)).
#   IIO holds : IIO, DM      (violated: MON, UN)   -> size on holds, power on violated
#   MON holds : MON, DM      (violated: IIO, UN)
# Separation lever: for the ordered-class models we vary how far apart the
# classes sit (sep in {0.6, 1.2, 2.0} logits) via a custom generator, so the
# cutoffs are stressed at weak and strong signal.
suppressMessages(library(QuantFit))
env <- asNamespace("QuantFit")
iio_h <- get(".manifest_iio_holds", env); mon_h <- get(".manifest_mon_holds", env)

# Generator that spans model x class-separation via the class x item logit
# table L[c,j], each model satisfying exactly its defining invariances (verified
# on the population probs: DM mon&iio, IIO iio-not-mon, MON mon-not-iio, UN
# neither). `sep` scales the signal (class spacing / slope spread).
#   DM  : L = theta_c - beta_j                (classes ordered, items ordered)
#   IIO : L = -a_c*beta_j + b_c , a_c>0 vary  (item order fixed; class curves cross)
#   MON : L = theta_c*a_j + b_j , a_j>0 vary  (class order fixed; item curves cross)
#   UN  : both crossing
gen_sep <- function(model, J, N, sep, seed) {
  set.seed(seed); C <- 3L
  beta  <- sort(runif(J, -1.5, 1.5))
  theta <- (seq_len(C) - (C + 1) / 2) * sep
  if (model == "DM") {
    L <- outer(theta, beta, "-")
  } else if (model == "IIO") {
    a <- 0.5 + (seq_len(C) - 1) * (sep / 2); b <- rnorm(C, 0, 0.3)
    L <- -outer(a, beta) + b
  } else if (model == "MON") {
    a <- 0.5 + (seq_len(J) - 1) / J * sep
    L <- outer(theta, a) + matrix(beta, C, J, byrow = TRUE)
  } else {                                          # UN
    a <- 0.5 + (seq_len(C) - 1) * (sep / 2)
    L <- -outer(a, beta) + matrix(rnorm(C * J, 0, sep), C, J)
  }
  P <- plogis(L); cl <- sample.int(C, N, replace = TRUE)
  d <- matrix(rbinom(N * J, 1, P[cl, ]), N, J); storage.mode(d) <- "integer"; d
}

B <- 49L; C <- 3L; ns <- 5L
cases <- expand.grid(model=c("UN","MON","IIO","DM"), J=c(6L,12L,24L),
                     N=c(750L,1500L,3000L), sep=c(0.6,1.2,2.0), rep=1:4,
                     stringsAsFactors=FALSE)
cases$id <- seq_len(nrow(cases)); set.seed(70701)
cases$seed <- sample.int(.Machine$integer.max, nrow(cases))
out <- Sys.getenv("IIOMON_OUT", "iiomon_robust_out"); dir.create(out, showWarnings=FALSE)
run <- function(k) { cs <- cases[k,]; f <- file.path(out, sprintf("r%04d.csv", k))
  if (file.exists(f)) return(invisible())
  d <- gen_sep(cs$model, cs$J, cs$N, cs$sep, cs$seed); ms <- cs$id
  eps_N <- 0.04 * sqrt(1500 / cs$N)   # per-item tolerance (perm-min q05)
  iio <- tryCatch(iio_h(d, C, B, ns, TRUE, ms),            error=function(e) list(holds=NA))
  mon <- tryCatch(mon_h(d, C, max(B%/%2L,20L), ns, TRUE, eps_N, ms+7L), error=function(e) list(holds=NA))
  write.csv(data.frame(model=cs$model, J=cs$J, N=cs$N, sep=cs$sep, rep=cs$rep,
    iio_holds=iio$holds, mon_holds=mon$holds), f, row.names=FALSE) }
cat("IIO/MON robustness:", nrow(cases), "datasets\n")
invisible(parallel::mclapply(seq_len(nrow(cases)),
  function(k) tryCatch(run(k), error=function(e) NULL), mc.cores=8))
cat("IIOMON ROBUST DONE\n")
