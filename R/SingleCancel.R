#' Gibbs sampler for single cancellation in one 3x3 submatrix
#'
#' Internal function. Checks ONLY single cancellation, in a single direction
#' (rows or columns), in one 3x3 submatrix.
#'
#' @param N Matrix containing the total number of responses.
#' @param n Matrix containing the number of correct responses.
#' @param n.iter Total number of samples.
#' @param burn Number of initial samples that should be discarded.
#' @param CR Width of the credible region taken from the posterior, e.g. a
#'   95% credible region is `c(.025,.975)`.
#' @param single Direction of the single cancellation constraint: `"R"` for
#'   rows or `"C"` for columns.
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
omni.check_single <- function(N, n, n.iter, burn=1000, CR, single) {
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
  for (I in 2:n.iter) {
    for (i in 1:nrow(dat)) for (j in 1:ncol(dat)) {
      lh1 <- if (j==1) 0 else old[i,j-1]
      lh2 <- if (i==1) 0 else old[i-1,j]
      rh1 <- if (j==ncol(dat)) 1 else old[i,j+1]
      rh2 <- if (i==nrow(dat)) 1 else old[i+1,j]
      if (single=="R") {
        lh <- lh1
        rh <- rh1
      }
      if (single=="C") {
        lh <- lh2
        rh <- rh2
      }
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
  for (i in 1:3) for (j in 1:3) {
    post <- sapply(chain, function(x) x[i,j])
    lo[i,j] <- quantile(post, CR[1])
    hi[i,j] <- quantile(post, CR[2])
    M[i,j] <- mean(post)
  }
  list(low=lo, high=hi, mean=M)
}

#' Check single cancellation in all adjacent 3x3 submatrices
#'
#' Checks ONLY the single cancellation axiom of additive conjoint measurement
#' in all adjacently formed 3x3 submatrices of `n`/`N`, in a single direction
#' (across rows or across columns).
#'
#' @param N Matrix containing the total number of responses.
#' @param n Matrix containing the number of correct responses.
#' @param CR Width of the credible region taken from the posterior. Defaults
#'   to a 95% credible region (`c(.025,.975)`).
#' @param single Direction of the single cancellation constraint: `"R"` for
#'   rows or `"C"` for columns.
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
#' @seealso [ConjointChecks()], [DoubleCancel()], [TripleCancel()]
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
#' out <- SingleCancel(tmp$N, tmp$n, single = "R")
#' summary(out)
#' }
#'
#' @export
SingleCancel <- function(N, n, CR=c(.025,.975), single, mc.cores=1) {
  proc.fun_adjacent <- function(dummy, arg.list) {
    r1 <- dummy[1]
    rows <- r1:(r1+2)
    c1 <- dummy[2]
    cols <- c1:(c1+2)
    N <- arg.list[[1]]
    n <- arg.list[[2]]
    lof <- arg.list[[3]]
    CR <- arg.list[[4]]
    single <- arg.list[[5]]
    omni.check_single <- lof[[1]]
    nt <- N[rows,cols]
    nc <- n[rows,cols]
    dat <- nc/nt
    test <- sum(dat==1|dat==0)
    if (test > 0) {
      NULL
    } else {
      out <- omni.check_single(nt, nc, n.iter=3000, CR=CR, single=single)
      list(rows, cols, out)
    }
  }
  dat <- n/N
  test <- ifelse(abs(dat-.5) <= .5, TRUE, FALSE)
  if (!all(test)) stop("There is a problem with n/N, values not between 0 and 1 (inclusive)")
  lof <- list(omni.check_single)
  arg.list <- list(N, n, lof, CR, single)
  dummy <- list()
  nr <- nrow(N)
  nc <- ncol(N)
  for (i in 1:(nr-2)) for (j in 1:(nc-2)) c(i,j) -> dummy[[paste(i,j)]]
  out <- parallel::mclapply(dummy, proc.fun_adjacent, arg.list=arg.list, mc.cores=mc.cores)
  destroy <- vapply(out, is.null, logical(1))
  out <- out[!destroy]
  if (length(out) == 0) {
    stop("No checkable 3x3 submatrices were found: every candidate contained ",
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
    for (i in 1:3) for (j in 1:3) {
      mat.den[ro[i],co[j]] <- mat.den[ro[i],co[j]] + 1
      mat.num[ro[i],co[j]] <- mat.num[ro[i],co[j]] + comp[i,j]
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
