#' Gibbs sampler for cancellation checks in a single submatrix
#'
#' Internal function. Runs the Metropolis-within-Gibbs sampler that checks the
#' single, double, or triple cancellation axiom in a single 3x3 (or 4x4 for
#' triple cancellation) submatrix of proportions.
#'
#' @param N Matrix containing the total number of responses.
#' @param n Matrix containing the number of correct responses.
#' @param n.iter Total number of samples.
#' @param burn Number of initial samples that should be discarded.
#' @param CR Width of the credible region taken from the posterior, e.g. a
#'   95% credible region is `c(.025,.975)`.
#' @param check Which cancellation axiom to impose: `"single"`, `"double"`,
#'   or `"triple"`.
#' @param use_cpp Use the compiled C++ sampler (`TRUE`, the default) or the
#'   pure R implementation.
#' @param adjust_extremes Nudge cells with observed proportions of exactly 0
#'   or 1 slightly inward so they can be checked. This imposes a 0.001/0.999
#'   pseudo-count adjustment: cells with `n == 0` are replaced by `0.001*N`
#'   and cells with `n == N` by `0.999*N`.
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
omni.check <- function(N, n, n.iter, burn=1000, CR, check, use_cpp=TRUE, adjust_extremes=FALSE) {
  if (check == "double" && !(nrow(N) == 3 && ncol(N) == 3)) {
    stop("check='double' requires a 3x3 submatrix; got ", nrow(N), "x", ncol(N))
  }
  if (check == "triple" && !(nrow(N) == 4 && ncol(N) == 4)) {
    stop("check='triple' requires a 4x4 submatrix; got ", nrow(N), "x", ncol(N))
  }
  if (n.iter <= burn + 4) {
    stop("n.iter (", n.iter, ") is too small relative to burn (", burn,
         "); no post-burn-in samples would be retained. Require n.iter > burn + 4.")
  }
  if (adjust_extremes) {
    n[n == 0] <- 0.001 * N[n == 0]
    n[n == N] <- 0.999 * N[n == N]
  }
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

  if (check=="single" && use_cpp) {
    chain <- CCIterateSingle(n.iter, old, old.ll, burn, N, n)
    mat_size <- nrow(dat)
  } else if (check=="double" && use_cpp) {
    chain <- CCIterateDouble(n.iter, old, old.ll, burn, N, n)
    mat_size <- nrow(dat)
  } else if (check=="triple" && use_cpp) {
    chain <- CCIterateTriple(n.iter, old, old.ll, burn, N, n)
    mat_size <- nrow(dat)
  } else if (check=="single" && !use_cpp) {
    for (I in 2:n.iter) {
      for (i in 1:nrow(dat)) for (j in 1:ncol(dat)) {
        lh1 <- if (j==1) 0 else old[i,j-1]
        lh2 <- if (i==1) 0 else old[i-1,j]
        rh1 <- if (j==ncol(dat)) 1 else old[i,j+1]
        rh2 <- if (i==nrow(dat)) 1 else old[i+1,j]
        lh <- max(lh1, lh2)
        rh <- min(rh1, rh2)
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
    mat_size <- nrow(dat)
  } else if (check=="double" && !use_cpp) {
    for (I in 2:n.iter) {
      for (i in 1:nrow(dat)) for (j in 1:ncol(dat)) {
        lh1 <- if (j==1) 0 else old[i,j-1]
        lh2 <- if (i==1) 0 else old[i-1,j]
        rh1 <- if (j==ncol(dat)) 1 else old[i,j+1]
        rh2 <- if (i==nrow(dat)) 1 else old[i+1,j]
        lh3 <- 0
        rh3 <- 1
        test.1 <- as.logical(old[2,1] < old[1,2])
        test.2 <- as.logical(old[3,2] < old[2,3])
        if (test.1 & test.2) {
          if (i==1 & j==3) lh3 <- old[3,1]
          if (i==3 & j==1) rh3 <- old[1,3]
        }
        if (!test.1 & !test.2) {
          if (i==3 & j==1) lh3 <- old[1,3]
          if (i==1 & j==3) rh3 <- old[3,1]
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
    mat_size <- nrow(dat)
  } else if (check=="triple" && !use_cpp) {
    list(c(1,2,3),c(1,2,4),c(1,3,4),c(2,3,4)) -> row_combos
    list(c(1,2,3),c(1,2,4),c(1,3,4),c(2,3,4)) -> col_combos
    dc_submats <- list()
    for (r in row_combos) for (co in col_combos) dc_submats[[length(dc_submats)+1]] <- list(rows=r, cols=co)
    tc_tests <- list(
      list(ant1=c(2,1,1,2), ant2=c(3,2,2,3), ant3=c(4,3,3,4), conseq=c(4,1,1,4)),
      list(ant1=c(2,1,1,2), ant2=c(3,2,2,3), ant3=c(4,2,3,3), conseq=c(4,1,1,3)),
      list(ant1=c(2,1,1,2), ant2=c(3,1,2,2), ant3=c(4,3,3,4), conseq=c(4,1,1,4)),
      list(ant1=c(2,1,1,2), ant2=c(3,1,2,2), ant3=c(4,2,3,3), conseq=c(4,1,1,3)),
      list(ant1=c(2,1,1,3), ant2=c(3,3,2,4), ant3=c(4,2,3,3), conseq=c(4,1,1,4)),
      list(ant1=c(2,1,1,3), ant2=c(3,3,2,4), ant3=c(4,1,3,2), conseq=c(4,2,1,4)),
      list(ant1=c(2,1,1,3), ant2=c(3,2,2,3), ant3=c(4,3,3,4), conseq=c(4,1,1,4)),
      list(ant1=c(2,1,1,3), ant2=c(3,2,2,3), ant3=c(4,1,3,2), conseq=c(4,2,1,4)),
      list(ant1=c(2,1,1,3), ant2=c(3,1,2,2), ant3=c(4,3,3,4), conseq=c(4,1,1,4)),
      list(ant1=c(2,2,1,3), ant2=c(3,3,2,4), ant3=c(4,1,3,2), conseq=c(4,1,1,4)),
      list(ant1=c(2,2,1,3), ant2=c(3,1,2,2), ant3=c(4,3,3,4), conseq=c(4,1,1,4)),
      list(ant1=c(2,2,1,3), ant2=c(3,1,2,2), ant3=c(4,2,3,3), conseq=c(4,1,1,4)),
      list(ant1=c(2,1,1,2), ant2=c(3,3,2,4), ant3=c(4,2,3,3), conseq=c(4,1,1,4)),
      list(ant1=c(2,1,1,2), ant2=c(3,2,2,4), ant3=c(4,3,3,4), conseq=c(4,1,1,4))
    )
    for (I in 2:n.iter) {
      for (i in 1:nrow(dat)) for (j in 1:ncol(dat)) {
        lh1 <- if (j==1) 0 else old[i,j-1]
        lh2 <- if (i==1) 0 else old[i-1,j]
        rh1 <- if (j==ncol(dat)) 1 else old[i,j+1]
        rh2 <- if (i==nrow(dat)) 1 else old[i+1,j]
        lh3 <- 0
        rh3 <- 1
        for (submat in dc_submats) {
          rs <- submat$rows
          cs <- submat$cols
          test.1 <- as.logical(old[rs[2],cs[1]] < old[rs[1],cs[2]])
          test.2 <- as.logical(old[rs[3],cs[2]] < old[rs[2],cs[3]])
          if (test.1 & test.2) {
            if (i==rs[1] & j==cs[3]) lh3 <- max(lh3, old[rs[3],cs[1]])
            if (i==rs[3] & j==cs[1]) rh3 <- min(rh3, old[rs[1],cs[3]])
          }
          if (!test.1 & !test.2) {
            if (i==rs[3] & j==cs[1]) lh3 <- max(lh3, old[rs[1],cs[3]])
            if (i==rs[1] & j==cs[3]) rh3 <- min(rh3, old[rs[3],cs[1]])
          }
        }
        lh4 <- 0
        rh4 <- 1
        for (test in tc_tests) {
          ant1 <- test$ant1
          ant2 <- test$ant2
          ant3 <- test$ant3
          conseq <- test$conseq
          t1 <- as.logical(old[ant1[1],ant1[2]] < old[ant1[3],ant1[4]])
          t2 <- as.logical(old[ant2[1],ant2[2]] < old[ant2[3],ant2[4]])
          t3 <- as.logical(old[ant3[1],ant3[2]] < old[ant3[3],ant3[4]])
          if (t1 & t2 & t3) {
            if (i==conseq[3] & j==conseq[4]) lh4 <- max(lh4, old[conseq[1],conseq[2]])
            if (i==conseq[1] & j==conseq[2]) rh4 <- min(rh4, old[conseq[3],conseq[4]])
          }
          if (!t1 & !t2 & !t3) {
            if (i==conseq[1] & j==conseq[2]) lh4 <- max(lh4, old[conseq[3],conseq[4]])
            if (i==conseq[3] & j==conseq[4]) rh4 <- min(rh4, old[conseq[1],conseq[2]])
          }
        }
        lh <- max(lh1, lh2, lh3, lh4)
        rh_dc <- min(rh3, rh4)
        if (rh_dc > lh) rh <- min(rh1, rh2, rh_dc) else rh <- min(rh1, rh2)
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
    mat_size <- nrow(dat)
  }
  if (length(chain) == 0) {
    stop("No post-burn-in samples were retained (n.iter=", n.iter, ", burn=", burn,
         "); increase n.iter relative to burn.")
  }
  hi <- lo <- M <- chain[[1]]
  for (i in 1:mat_size) for (j in 1:mat_size) {
    post <- sapply(chain, function(x) x[i,j])
    lo[i,j] <- quantile(post, CR[1])
    hi[i,j] <- quantile(post, CR[2])
    M[i,j] <- mean(post)
  }
  list(low=lo, high=hi, mean=M)
}

#' Check cancellation axioms in a sample of submatrices
#'
#' Given two matrices, `n` and `N` (which contain the number of correct
#' responses and the number of total responses for each cell), a Bayesian
#' check of the requested cancellation axiom (single, double, or triple
#' cancellation of additive conjoint measurement) is performed in `n.mat`
#' sampled submatrices. Double cancellation is checked in 3x3 submatrices;
#' triple cancellation in 4x4 submatrices. To check large numbers of
#' submatrices (to see why, see Domingue, 2014), parallel options help.
#'
#' @param N Matrix containing the total number of responses.
#' @param n Matrix containing the number of correct responses.
#' @param n.mat Number of submatrices to sample or the string `"adjacent"`
#'   if all adjacently formed submatrices are to be checked.
#' @param CR Width of the credible region taken from the posterior. Defaults
#'   to a 95% credible region (`c(.025,.975)`).
#' @param check Which cancellation axiom to test: `"single"`, `"double"`
#'   (the default), or `"triple"`.
#' @param mc.cores The number of cores to parallelize over. When
#'   `mc.cores > 1`, set `RNGkind("L'Ecuyer-CMRG")` (and a seed) before
#'   calling for reproducible results across parallel workers.
#' @param use_cpp Use the compiled C++ sampler (`TRUE`, the default) or the
#'   pure R implementation.
#' @param adjust_extremes Nudge cells with observed proportions of exactly 0
#'   or 1 slightly inward so that submatrices containing them can still be
#'   checked. If `FALSE` (the default), such submatrices are skipped. When
#'   `TRUE`, a 0.001/0.999 pseudo-count adjustment is imposed: cells with
#'   `n == 0` are treated as having proportion 0.001 and cells with `n == N`
#'   as having proportion 0.999, both inside the sampler and when the
#'   observed proportions are compared against the posterior credible
#'   regions.
#'
#' @return An object of class [`checks`][checks-class] summarizing the
#'   detected violations.
#'
#' @references
#' Domingue, B. (2014). Evaluating the equal-interval hypothesis with test
#' score scales. \emph{Psychometrika}, 79(1), 1-19.
#' \doi{10.1007/s11336-013-9342-4}
#'
#' Perline, R., Wright, B. D., & Wainer, H. (1979). The Rasch model as
#' additive conjoint measurement. \emph{Applied Psychological Measurement},
#' 3(2), 237-255.
#'
#' @author Ben Domingue \email{ben.domingue@@gmail.com}
#'
#' @seealso [PrepareChecks()], [HiConjointChecks()], [summary.checks()],
#'   [plot.checks()]
#'
#' @examples
#' \dontrun{
#' ######################################################
#' # parole data
#' # page 244 (table 2) of Perline, Wright, and Wainer
#' # about 9% were bad in perline
#' matrix(c(15,47,61,84,82,86,60,47,8),9,9,byrow=FALSE)->N
#' per <-structure(c(0, 0.06, 0.07, 0.18, 0.13, 0.13, 0.17, 0.17,
#'  1, 0, 0.04, 0.15, 0.24, 0.33, 0.28, 0.47, 0.85, 1, 0, 0.04, 0.08,
#'  0.12, 0.3, 0.64, 0.85, 1, 1, 0, 0.19, 0.39, 0.4, 0.51, 0.58,
#'  0.82, 0.98, 1, 0, 0.06, 0.18, 0.52, 0.73, 0.95, 1, 1, 1, 0,
#'  0.23, 0.33, 0.51, 0.68, 0.91, 0.93, 1, 1, 0.27, 0.51, 0.61,
#'  0.64, 0.68, 0.77, 0.9, 1, 1, 0, 0.21, 0.52, 0.68, 0.84, 0.97,
#'  0.97, 1, 1, 0.73, 0.64, 0.67, 0.7, 0.78, 0.78, 0.9, 1, 1),
#'  .Dim = c(9L, 9L) )
#' round(per*N)->n
#' ConjointChecks(N,n,n.mat=1)->out
#'
#' ######################################################
#' # Data from Rasch (1960)
#' # page 250 (table 5) of Perline, Wright, and Wainer
#' # about 4% showed violations
#' matrix(c(49,112,32,76,82,102,119,133,123,94,61,17,10),13,7,byrow=FALSE)->N
#' per <-structure(c(0, 0, 0, 0, 0.02, 0.01, 0.02, 0.03, 0.06, 0.09,
#'  0.23, 0.35, 0.7, 0.01, 0, 0.04, 0.05, 0.09, 0.09, 0.16, 0.28, 0.39,
#'  0.66, 0.8, 0.91, 0.85, 0, 0.02, 0.07, 0.07, 0.24, 0.28, 0.45, 0.59,
#'  0.76, 0.87, 0.9, 1, 0.85, 0.01, 0.04, 0.12, 0.21, 0.42, 0.62, 0.73,
#'  0.83, 0.9, 0.93, 0.98, 1, 1, 0.06, 0.11, 0.4, 0.7, 0.7, 0.79, 0.84,
#'  0.88, 0.94, 0.95, 0.98, 1, 1, 0.48, 0.84, 0.84, 0.86, 0.86, 0.9,
#'  0.95, 0.96, 0.98, 0.99, 0.99, 1, 1, 0.92, 0.98, 0.98, 0.99, 0.98,
#'  0.99, 0.99, 1, 1, 1, 1, 1, 1), .Dim = c(13L, 7L))
#' round(per*N)->n
#' ConjointChecks(N,n,n.mat=1)->out
#'
#' ###########
#' # simulated Rasch example
#' n.mat<-1000
#' n.items<-20
#' n.respondents<-2000
#' # simulate data
#' rnorm(n.items)->diff
#' rnorm(n.respondents)->abil
#' matrix(abil,n.respondents,n.items,byrow=FALSE)->m1
#' matrix(diff,n.respondents,n.items,byrow=TRUE)->m2
#' m1-m2 -> kern
#' exp(kern)/(1+exp(kern))->pv
#' runif(n.items*n.respondents)->test
#' ifelse(pv>test,1,0)->resp
#' ## now check
#' PrepareChecks(resp)->tmp
#' ConjointChecks(tmp$N,tmp$n,n.mat=n.mat,mc.cores=1)->rasch1000
#' }
#'
#' @export
ConjointChecks <- function(N, n, n.mat=1, CR=c(.025,.975), check="double", mc.cores=1, use_cpp=TRUE, adjust_extremes=FALSE) {
  if (!(check %in% c("single","double","triple"))) stop("check must be 'single', 'double', or 'triple'")
  if (check=="triple") mat_size <- 4 else mat_size <- 3
  proc.fun <- function(dummy, arg.list) {
    N <- arg.list[[1]]
    n <- arg.list[[2]]
    lof <- arg.list[[3]]
    CR <- arg.list[[4]]
    check <- arg.list[[5]]
    mat_size <- arg.list[[6]]
    use_cpp <- arg.list[[7]]
    adjust_extremes <- arg.list[[8]]
    omni.check <- lof[[1]]
    test <- 1
    nrows <- nrow(N)
    ncols <- ncol(N)
    max_tries <- 1000
    tries <- 0
    while (test > 0 & tries < max_tries) {
      rows <- sort(sample(1:nrows, mat_size, replace=FALSE))
      cols <- sort(sample(1:ncols, mat_size, replace=FALSE))
      nt <- N[rows,cols]
      nc <- n[rows,cols]
      dat <- nc/nt
      if (adjust_extremes) test <- 0 else test <- sum(dat==1|dat==0)
      tries <- tries+1
    }
    if (tries >= max_tries & !adjust_extremes) return(NULL)
    out <- omni.check(nt, nc, n.iter=3000, CR=CR, check=check, use_cpp=use_cpp, adjust_extremes=adjust_extremes)
    list(rows, cols, out)
  }
  proc.fun_adjacent <- function(dummy, arg.list) {
    r1 <- dummy[1]
    rows <- r1:(r1+arg.list[[6]]-1)
    c1 <- dummy[2]
    cols <- c1:(c1+arg.list[[6]]-1)
    N <- arg.list[[1]]
    n <- arg.list[[2]]
    lof <- arg.list[[3]]
    CR <- arg.list[[4]]
    check <- arg.list[[5]]
    use_cpp <- arg.list[[7]]
    adjust_extremes <- arg.list[[8]]
    omni.check <- lof[[1]]
    nt <- N[rows,cols]
    nc <- n[rows,cols]
    dat <- nc/nt
    test <- sum(dat==1|dat==0)
    if (test > 0 & !adjust_extremes) {
      NULL
    } else {
      out <- omni.check(nt, nc, n.iter=3000, CR=CR, check=check, use_cpp=use_cpp, adjust_extremes=adjust_extremes)
      list(rows, cols, out)
    }
  }
  proc.fun_exact <- function(dummy, arg.list) {
    N <- arg.list[[1]]; n <- arg.list[[2]]; lof <- arg.list[[3]]
    CR <- arg.list[[4]]; check <- arg.list[[5]]
    use_cpp <- arg.list[[7]]; adjust_extremes <- arg.list[[8]]
    rows <- dummy$rows; cols <- dummy$cols
    nt <- N[rows, cols]; ncnt <- n[rows, cols]
    dat <- ncnt/nt
    if (!adjust_extremes && sum(dat == 1 | dat == 0) > 0) return(NULL)
    out <- lof[[1]](nt, ncnt, n.iter = 3000, CR = CR, check = check,
                    use_cpp = use_cpp, adjust_extremes = adjust_extremes)
    list(rows, cols, out)
  }
  dat <- n/N
  test <- ifelse(abs(dat-.5) <= .5, TRUE, FALSE)
  if (!all(test)) stop("There is a problem with n/N, values not between 0 and 1 (inclusive)")
  nr <- nrow(N)
  nc <- ncol(N)
  if (nr < mat_size | nc < mat_size) stop(paste("Matrix must have at least", mat_size, "rows and", mat_size, "columns for", check, "cancellation"))
  lof <- list(omni.check)
  arg.list <- list(N, n, lof, CR, check, mat_size, use_cpp, adjust_extremes)
  dummy <- list()
  total_combos <- choose(nr, mat_size) * choose(nc, mat_size)
  exhaustive <- identical(n.mat, "all") ||
    (is.numeric(n.mat) && n.mat >= total_combos)
  if (identical(n.mat, "all") && total_combos > 2e5) {
    warning("n.mat='all' requested but the table has ", total_combos,
            " submatrices (> 2e5); sampling 5000 instead")
    exhaustive <- FALSE; n.mat <- 5000
  }
  if (n.mat=="adjacent") {
    for (i in 1:(nr-mat_size+1)) for (j in 1:(nc-mat_size+1)) c(i,j) -> dummy[[paste(i,j)]]
    out <- parallel::mclapply(dummy, proc.fun_adjacent, arg.list=arg.list, mc.cores=mc.cores)
  } else if (exhaustive) {
    # EXHAUSTIVE mode: every mat_size x mat_size submatrix once - the exact
    # population violation rate of the table, no Monte Carlo noise. Reached
    # explicitly via n.mat='all', or automatically whenever the requested
    # sample size meets or exceeds the number of distinct submatrices
    # (sampling more than the population is wasted duplication).
    rcomb <- utils::combn(nr, mat_size); ccomb <- utils::combn(nc, mat_size)
    k <- 0L
    for (i in seq_len(ncol(rcomb))) for (j in seq_len(ncol(ccomb))) {
      k <- k + 1L
      dummy[[k]] <- list(rows = rcomb[, i], cols = ccomb[, j])
    }
    out <- parallel::mclapply(dummy, proc.fun_exact, arg.list=arg.list, mc.cores=mc.cores)
  } else {
    for (i in 1:n.mat) dummy[[i]] <- i
    out <- parallel::mclapply(dummy, proc.fun, arg.list=arg.list, mc.cores=mc.cores)
  }
  destroy <- vapply(out, is.null, logical(1))
  out <- out[!destroy]
  if (length(out) == 0) {
    stop("No checkable ", mat_size, "x", mat_size, " submatrices were found: every candidate ",
         "contained observed proportions of exactly 0 or 1. ",
         "Consider adjust_extremes=TRUE or collapsing sparse cells (see PrepareChecks).")
  }
  compare <- function(dat, lim) {
    ifelse(dat < lim[[2]] & dat > lim[[1]], 0, 1)
  }
  dat <- n/N
  if (adjust_extremes) {
    # apply the same 0.001/0.999 pseudo-count adjustment used inside omni.check
    # so that 0/1 cells are not automatically flagged as violations
    dat[dat == 0] <- 0.001
    dat[dat == 1] <- 0.999
  }
  mat.num <- matrix(0, nrow(N), ncol(N))
  mat.den <- matrix(-1, nrow(N), ncol(N))
  for (k in 1:length(out)) {
    x <- out[[k]]
    ro <- x[[1]]
    co <- x[[2]]
    comp <- compare(dat[ro,co], x[[3]])
    for (i in 1:mat_size) for (j in 1:mat_size) {
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
  new("checks", N=N, n=n, Checks=out, tab=tab,
      means=list(unweighted=m1, weighted=m2,
                 coverage=c(checked=length(out), total=total_combos,
                            exhaustive=as.numeric(exhaustive))),
      check.counts=mat.den)
}
