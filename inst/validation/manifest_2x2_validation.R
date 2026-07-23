# Full-scale validation of the manifest 2x2 ordinal classifier.
# 4 models (UN/MON/IIO/DM) x J in {6,12,24} x 15 reps = 180 datasets.
# IIO axis: manifest crossing, parametric-DM calibrated (worked well).
# MON axis: fitted-UN class-monotonicity, DATA-RESAMPLING calibrated.
# One CSV per dataset (resumable, snapshotable).
suppressMessages(library(QuantFit))
rmt <- QuantFit:::refit_model_type

iio_stat <- function(data, min_group = 20L) {
  data <- as.matrix(data); J <- ncol(data); N <- nrow(data)
  p <- colMeans(data); ord <- order(p, decreasing = TRUE)
  total <- rowSums(data); mag <- 0
  for (ai in seq_len(J-1L)) for (bi in (ai+1L):J) {
    a <- ord[ai]; b <- ord[bi]; rest <- total - data[,a] - data[,b]
    for (r in sort(unique(rest))) { g <- rest==r; if (sum(g) < min_group) next
      d <- mean(data[g,b]) - mean(data[g,a]); if (d>0) mag <- mag + d*sum(g) } }
  mag / N
}
iio_holds <- function(data, C=3L, B=40L, seed=1) {
  obs <- iio_stat(data); dm <- rmt("DM", as.matrix(data), C, 5L, TRUE)
  set.seed(seed); n <- nrow(data)
  null <- vapply(seq_len(B), function(b) iio_stat(QuantFit:::simulate_from_qlfit(dm,n)), numeric(1))
  (1 + sum(null >= obs))/(B+1) > 0.05
}
mon_stat <- function(data, C=3L, n_starts=8L) {
  un <- rmt("UN", as.matrix(data), C, n_starts, TRUE)
  P <- un$item_probs; P <- P[, order(colMeans(P)), drop=FALSE]
  sum(pmax(0, -t(apply(P,1,diff))))
}
mon_holds <- function(data, C=3L, B=25L, eps=0.03, seed=1) {
  set.seed(seed); n <- nrow(data)
  boot <- vapply(seq_len(B), function(b)
    mon_stat(data[sample.int(n,n,replace=TRUE),,drop=FALSE], C), numeric(1))
  stats::quantile(boot, 0.05, names=FALSE) <= eps      # holds if lower bound negligible
}
classify <- function(data, C=3L, seed=1) {
  ih <- iio_holds(data, C, 40L, seed); mh <- mon_holds(data, C, 25L, 0.03, seed+7L)
  if (ih && mh) "DM" else if (ih && !mh) "IIO" else if (!ih && mh) "MON" else "UN"
}
gen <- function(model,J,seed){ set.seed(seed)
  d <- simulate_responses(model, n_persons=1500L, n_items=J, n_classes=3L)
  d <- if(is.list(d))d$data else d; storage.mode(d)<-"integer"; d}

out <- "manifest_full_out"; dir.create(out, showWarnings=FALSE)
cases <- expand.grid(model=c("UN","MON","IIO","DM"), J=c(6L,12L,24L), rep=1:15,
                     stringsAsFactors=FALSE)
cases$id <- seq_len(nrow(cases)); set.seed(20260723)
cases$seed <- sample.int(.Machine$integer.max, nrow(cases))
run <- function(k){ cs <- cases[k,]; f <- file.path(out, sprintf("f%03d.csv", k))
  if(file.exists(f)) return(invisible())
  d <- gen(cs$model, cs$J, cs$seed)
  sel <- tryCatch(classify(d, seed=cs$id), error=function(e) NA)
  write.csv(data.frame(model=cs$model, J=cs$J, rep=cs$rep, selected=sel), f, row.names=FALSE) }
cat("manifest 2x2 full-scale:", nrow(cases), "datasets\n")
invisible(parallel::mclapply(seq_len(nrow(cases)),
  function(k) tryCatch(run(k), error=function(e) NULL), mc.cores=8))
cat("MANIFEST FULL DONE\n")
