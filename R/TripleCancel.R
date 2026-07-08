#' Gibbs sampler for triple cancellation in one 4x4 submatrix
#'
#' Internal function. Checks ONLY the triple cancellation bits in one 4x4
#' submatrix, using the 14 coherent triple cancellation tests of Kyngdon and
#' Richards (2006).
#'
#' @param N Matrix containing the total number of responses.
#' @param n Matrix containing the number of correct responses.
#' @param n.iter Total number of samples.
#' @param burn Number of initial samples that should be discarded.
#' @param CR Width of the credible region taken from the posterior, e.g. a
#'   95% credible region is `c(.025,.975)`.
#'
#' @return A list with elements `low`, `high`, and `mean`: matrices giving the
#'   lower and upper limits of the credible region and the posterior mean for
#'   each cell.
#'
#' @references
#' Perline, R., Wright, B. D., & Wainer, H. (1979). The Rasch model as
#' additive conjoint measurement. \emph{Applied Psychological Measurement},
#' 3(2), 237-255.
#'
#' @author Ben Domingue \email{ben.domingue@@gmail.com}
#'
#' @keywords internal
omni.check_triple <- function(N, n, n.iter, burn=1000, CR) {
  dat <- n/N
  chain <- list()
  inits <- dat
  inits -> foo1 -> foo2
  for (i in 1:nrow(foo1)) foo1[i,] <- i
  for (i in 1:ncol(foo2)) foo2[,i] <- i
  foo1 + foo2 -> index
  sort(as.numeric(dat)) -> hold
  counter <- 1
  for (i in unique(as.numeric(index))) {
    grep(paste("^",i,"$",sep=""), index) -> index2
    stop.here <- length(index2)
    counter:(counter+stop.here-1) -> replace
    inits[index2] <- hold[replace]
    counter <- max(replace)+1
  }
  old <- inits
  like <- function(theta, N, n) {
    n*log(theta) + (N-n)*log(1-theta)
  }
  old.ll <- inits
  for (i in 1:nrow(old.ll)) for (j in 1:ncol(old.ll)) like(inits[i,j], N[i,j], n[i,j]) -> old.ll[i,j]
  # 14 coherent triple cancellation tests from Kyngdon & Richards (2006)
  # Each entry: list(ant1, ant2, ant3, conseq) where each is c(r1,c1,r2,c2)
  tc_tests <- list(
    # Form 1
    list(ant1=c(2,1,1,2), ant2=c(3,2,2,3), ant3=c(4,3,3,4), conseq=c(4,1,1,4)),
    list(ant1=c(2,1,1,2), ant2=c(3,2,2,3), ant3=c(4,2,3,3), conseq=c(4,1,1,3)),
    list(ant1=c(2,1,1,2), ant2=c(3,1,2,2), ant3=c(4,3,3,4), conseq=c(4,1,1,4)),
    list(ant1=c(2,1,1,2), ant2=c(3,1,2,2), ant3=c(4,2,3,3), conseq=c(4,1,1,3)),
    # Form 2
    list(ant1=c(2,1,1,3), ant2=c(3,3,2,4), ant3=c(4,2,3,3), conseq=c(4,1,1,4)),
    list(ant1=c(2,1,1,3), ant2=c(3,3,2,4), ant3=c(4,1,3,2), conseq=c(4,2,1,4)),
    list(ant1=c(2,1,1,3), ant2=c(3,2,2,3), ant3=c(4,3,3,4), conseq=c(4,1,1,4)),
    list(ant1=c(2,1,1,3), ant2=c(3,2,2,3), ant3=c(4,1,3,2), conseq=c(4,2,1,4)),
    list(ant1=c(2,1,1,3), ant2=c(3,1,2,2), ant3=c(4,3,3,4), conseq=c(4,1,1,4)),
    # Form 3
    list(ant1=c(2,2,1,3), ant2=c(3,3,2,4), ant3=c(4,1,3,2), conseq=c(4,1,1,4)),
    list(ant1=c(2,2,1,3), ant2=c(3,1,2,2), ant3=c(4,3,3,4), conseq=c(4,1,1,4)),
    list(ant1=c(2,2,1,3), ant2=c(3,1,2,2), ant3=c(4,2,3,3), conseq=c(4,1,1,4)),
    list(ant1=c(2,1,1,2), ant2=c(3,3,2,4), ant3=c(4,2,3,3), conseq=c(4,1,1,4)),
    list(ant1=c(2,1,1,2), ant2=c(3,2,2,4), ant3=c(4,3,3,4), conseq=c(4,1,1,4))
  )
  for (I in 2:n.iter) {
    for (i in 1:nrow(dat)) for (j in 1:ncol(dat)) {
      lh1 <- 0
      lh2 <- 0
      rh1 <- 1
      rh2 <- 1
      lh3 <- 0
      rh3 <- 1
      for (test in tc_tests) {
        ant1 <- test$ant1
        ant2 <- test$ant2
        ant3 <- test$ant3
        conseq <- test$conseq
        t1 <- as.logical(old[ant1[1],ant1[2]] < old[ant1[3],ant1[4]])
        t2 <- as.logical(old[ant2[1],ant2[2]] < old[ant2[3],ant2[4]])
        t3 <- as.logical(old[ant3[1],ant3[2]] < old[ant3[3],ant3[4]])
        if (t1 & t2 & t3) {
          if (i==conseq[3] & j==conseq[4]) lh3 <- max(lh3, old[conseq[1],conseq[2]])
          if (i==conseq[1] & j==conseq[2]) rh3 <- min(rh3, old[conseq[3],conseq[4]])
        }
        if (!t1 & !t2 & !t3) {
          if (i==conseq[1] & j==conseq[2]) lh3 <- max(lh3, old[conseq[3],conseq[4]])
          if (i==conseq[3] & j==conseq[4]) rh3 <- min(rh3, old[conseq[1],conseq[2]])
        }
      }
      lh <- max(lh1, lh2, lh3)
      if (rh3 > lh) rh <- min(rh1, rh2, rh3) else rh <- min(rh1, rh2)
      if (rh < lh) rh <- 1
      draw <- runif(1, lh, rh)
      ar <- 2
      new.ll <- like(draw, N[i,j], n[i,j])
      if (!(old[i,j] %in% 0:1)) ar <- exp(new.ll - old.ll[i,j])
      if (ar > runif(1)) {
        old[i,j] <- draw
        old.ll[i,j] <- new.ll
      }
    }
    if (I > burn & I%%4==0) old -> chain[[as.character(I)]]
  }
  hi <- lo <- M <- chain[[1]]
  for (i in 1:4) for (j in 1:4) {
    post <- sapply(chain, function(x) x[i,j])
    lo[i,j] <- quantile(post, CR[1])
    hi[i,j] <- quantile(post, CR[2])
    M[i,j] <- mean(post)
  }
  list(low=lo, high=hi, mean=M)
}

#' Check triple cancellation in a sample of 4-matrices
#'
#' Checks ONLY the triple cancellation axiom of additive conjoint measurement
#' in a sample of 4x4 submatrices of `n`/`N`. The check uses the 14 coherent
#' triple cancellation tests within a 4x4 matrix.
#'
#' Note that this standalone function deliberately tests the triple
#' cancellation axiom IN ISOLATION: the sampler drops the single-cancellation
#' (row/column monotonicity) bounds. This differs from
#' [ConjointChecks()] with `check="triple"`, which imposes the single,
#' double, and triple cancellation axioms jointly. Results from the two
#' approaches are therefore not comparable.
#'
#' @param N Matrix containing the total number of responses.
#' @param n Matrix containing the number of correct responses.
#' @param n.4mat Number of 4-matrices to sample or the string `"adjacent"`
#'   if all adjacently formed 4-matrices are to be checked.
#' @param CR Width of the credible region taken from the posterior. Defaults
#'   to a 95% credible region (`c(.025,.975)`).
#' @param mc.cores The number of cores to parallelize over.
#'
#' @return An object of class [`checks`][checks-class] summarizing the
#'   detected violations.
#'
#' @references
#' Perline, R., Wright, B. D., & Wainer, H. (1979). The Rasch model as
#' additive conjoint measurement. \emph{Applied Psychological Measurement},
#' 3(2), 237-255.
#'
#' @author Ben Domingue \email{ben.domingue@@gmail.com}
#'
#' @seealso [ConjointChecks()], [SingleCancel()], [DoubleCancel()],
#'   [HiConjointChecks()]
#'
#' @examples
#' \dontrun{
#' # simulated Rasch example
#' n.items <- 20
#' n.respondents <- 2000
#' diff <- rnorm(n.items)
#' abil <- rnorm(n.respondents)
#' kern <- outer(abil, diff, "-")
#' pv <- exp(kern)/(1+exp(kern))
#' resp <- ifelse(pv > runif(n.items*n.respondents), 1, 0)
#' tmp <- PrepareChecks(resp)
#' out <- TripleCancel(tmp$N, tmp$n, n.4mat = 100)
#' summary(out)
#' }
#'
#' @export
TripleCancel <- function(N, n, n.4mat=1, CR=c(.025,.975), mc.cores=1) {
  proc.fun <- function(dummy, arg.list) {
    N <- arg.list[[1]]
    n <- arg.list[[2]]
    lof <- arg.list[[3]]
    CR <- arg.list[[4]]
    omni.check_triple <- lof[[1]]
    test <- 1
    nrows <- nrow(N)
    ncols <- ncol(N)
    max_tries <- 1000
    tries <- 0
    while (test > 0 & tries < max_tries) {
      rows <- sort(sample(1:nrows, 4, replace=FALSE))
      cols <- sort(sample(1:ncols, 4, replace=FALSE))
      nt <- N[rows,cols]
      nc <- n[rows,cols]
      dat <- nc/nt
      test <- sum(dat==1|dat==0)
      tries <- tries + 1
    }
    if (test > 0) return(NULL)
    out <- omni.check_triple(nt, nc, n.iter=3000, CR=CR)
    list(rows, cols, out)
  }
  proc.fun_adjacent <- function(dummy, arg.list) {
    r1 <- dummy[1]
    rows <- r1:(r1+3)
    c1 <- dummy[2]
    cols <- c1:(c1+3)
    N <- arg.list[[1]]
    n <- arg.list[[2]]
    lof <- arg.list[[3]]
    CR <- arg.list[[4]]
    omni.check_triple <- lof[[1]]
    nt <- N[rows,cols]
    nc <- n[rows,cols]
    dat <- nc/nt
    test <- sum(dat==1|dat==0)
    if (test > 0) {
      NULL
    } else {
      out <- omni.check_triple(nt, nc, n.iter=3000, CR=CR)
      list(rows, cols, out)
    }
  }
  dat <- n/N
  test <- ifelse(abs(dat-.5) <= .5, TRUE, FALSE)
  if (!all(test)) stop("There is a problem with n/N, values not between 0 and 1 (inclusive)")
  lof <- list(omni.check_triple)
  arg.list <- list(N, n, lof, CR)
  dummy <- list()
  if (n.4mat=="adjacent") {
    nr <- nrow(N)
    nc <- ncol(N)
    if (nr < 4 | nc < 4) stop("Matrix must have at least 4 rows and 4 columns for triple cancellation")
    for (i in 1:(nr-3)) for (j in 1:(nc-3)) c(i,j) -> dummy[[paste(i,j)]]
    out <- parallel::mclapply(dummy, proc.fun_adjacent, arg.list=arg.list, mc.cores=mc.cores)
  } else {
    nr <- nrow(N)
    nc <- ncol(N)
    if (nr < 4 | nc < 4) stop("Matrix must have at least 4 rows and 4 columns for triple cancellation")
    for (i in 1:n.4mat) dummy[[i]] <- i
    out <- parallel::mclapply(dummy, proc.fun, arg.list=arg.list, mc.cores=mc.cores)
  }
  destroy <- vapply(out, is.null, logical(1))
  out <- out[!destroy]
  if (length(out) == 0) {
    stop("No checkable 4x4 submatrices were found: every candidate contained ",
         "observed proportions of exactly 0 or 1. ",
         "Consider collapsing sparse cells (see PrepareChecks).")
  }
  compare <- function(dat, lim) {
    ifelse(dat < lim[[2]] & dat > lim[[1]], 0, 1)
  }
  dat <- n/N
  mat.num <- matrix(0, nrow(N), ncol(N))
  mat.den <- matrix(-1, nrow(N), ncol(N))
  for (i in 1:length(out)) {
    x <- out[[i]]
    ro <- x[[1]]
    co <- x[[2]]
    comp <- compare(dat[ro,co], x[[3]])
    for (ii in 1:4) for (jj in 1:4) {
      mat.den[ro[ii],co[jj]] <- mat.den[ro[ii],co[jj]] + 1
      mat.num[ro[ii],co[jj]] <- mat.num[ro[ii],co[jj]] + comp[ii,jj]
    }
  }
  mat.den <- ifelse(mat.den < 0, NA, mat.den)
  mat.den <- mat.den + 1
  tab <- mat.num/mat.den
  m1 <- mean(tab, na.rm=TRUE)
  weight <- N
  m2 <- sum(tab*weight, na.rm=TRUE)/sum(weight[!is.na(tab)])
  new("checks", N=N, n=n, Checks=out, tab=tab, means=list(unweighted=m1, weighted=m2), check.counts=mat.den)
}
