#' Hierarchical cancellation checks in 4x4 submatrices
#'
#' Performs a hierarchical sequence of cancellation checks in sampled (or all
#' adjacent) 4x4 submatrices of `n`/`N`. Within each 4x4 submatrix the
#' function first tests single cancellation on the full 4x4 matrix; if that
#' passes, it tests double cancellation in each of the 16 embedded 3x3
#' submatrices; if all of those pass, it tests triple cancellation on the
#' full 4x4 matrix. Each submatrix is thus classified according to the
#' highest cancellation condition it satisfies.
#'
#' @param N Matrix containing the total number of responses.
#' @param n Matrix containing the number of correct responses.
#' @param n.mat Number of 4x4 submatrices to sample or the string
#'   `"adjacent"` if all adjacently formed 4x4 submatrices are to be checked.
#' @param CR Width of the credible region taken from the posterior. Defaults
#'   to a 95% credible region (`c(.025,.975)`).
#' @param mc.cores The number of cores to parallelize over.
#' @param use_cpp Use the compiled C++ sampler (`TRUE`, the default) or the
#'   pure R implementation.
#' @param adjust_extremes Nudge cells with observed proportions of exactly 0
#'   or 1 slightly inward so that submatrices containing them can still be
#'   checked. If `FALSE` (the default), such submatrices are skipped.
#'
#' @return A list with components:
#' \describe{
#'   \item{`N`, `n`}{The input matrices.}
#'   \item{`summary`}{A list giving the number of submatrices tested and the
#'     counts with status `PASSED_ALL`, `FAILED_SINGLE`, `FAILED_DOUBLE`,
#'     `FAILED_TRIPLE`, `SKIPPED`, and `ERROR` (a child test that failed
#'     with an error is reported with status `ERROR` rather than aborting
#'     the whole run).}
#'   \item{`results`}{A list with one element per tested submatrix containing
#'     its rows, columns, status, and the detailed results of each stage of
#'     the hierarchical test.}
#' }
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
#' @seealso [ConjointChecks()], [TripleCancel()], [PrepareChecks()]
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
#' out <- HiConjointChecks(tmp$N, tmp$n, n.mat = 10)
#' out$summary
#' }
#'
#' @export
HiConjointChecks <- function(N, n, n.mat=1, CR=c(.025,.975), mc.cores=1, use_cpp=TRUE, adjust_extremes=FALSE) {
  dat <- n/N
  test <- ifelse(abs(dat-.5) <= .5, TRUE, FALSE)
  if (!all(test)) stop("There is a problem with n/N, values not between 0 and 1 (inclusive)")
  nr <- nrow(N)
  nc <- ncol(N)
  if (nr < 4 | nc < 4) stop("Matrix must have at least 4 rows and 4 columns for hierarchical check")

  list(c(1,2,3),c(1,2,4),c(1,3,4),c(2,3,4)) -> row_combos_3x3
  list(c(1,2,3),c(1,2,4),c(1,3,4),c(2,3,4)) -> col_combos_3x3
  submat_3x3 <- list()
  for (r in row_combos_3x3) {
    for (co in col_combos_3x3) {
      paste(paste(r,collapse=","), paste(co,collapse=","), sep="_") -> key
      list(rows=r, cols=co) -> submat_3x3[[key]]
    }
  }

  check_pass <- function(obs, low, high) {
    all(obs >= low & obs <= high)
  }

  run_hierarchical_test <- function(N_sub, n_sub, rows_4x4, cols_4x4, CR, use_cpp, submat_3x3, adjust_extremes) {
    dat_sub_raw <- n_sub/N_sub
    if (any(dat_sub_raw == 0 | dat_sub_raw == 1) & !adjust_extremes) {
      return(list(
        rows=rows_4x4,
        cols=cols_4x4,
        status="SKIPPED",
        message="Contains 0 or 1 proportions - cannot test",
        single_result=NULL,
        double_results=NULL,
        double_failures=NULL,
        triple_result=NULL
      ))
    }

    if (adjust_extremes) {
      dat_sub <- dat_sub_raw
      dat_sub[dat_sub == 0] <- 0.001
      dat_sub[dat_sub == 1] <- 0.999
    } else {
      dat_sub <- dat_sub_raw
    }

    single_result <- NULL
    double_results <- list()
    double_failures <- list()
    triple_result <- NULL

    omni.check(N_sub, n_sub, n.iter=3000, CR=CR, check="single", use_cpp=use_cpp, adjust_extremes=adjust_extremes) -> out
    check_pass(dat_sub, out$low, out$high) -> single_pass
    list(result=out, pass=single_pass) -> single_result

    if (!single_pass) {
      return(list(
        rows=rows_4x4,
        cols=cols_4x4,
        status="FAILED_SINGLE",
        message="Failed at single cancellation on 4x4",
        single_result=single_result,
        single_failure=list(observed=dat_sub, low=out$low, high=out$high),
        double_results=NULL,
        double_failures=NULL,
        triple_result=NULL
      ))
    }

    TRUE -> double_pass_all
    for (key in names(submat_3x3)) {
      submat_3x3[[key]] -> sm
      sm$rows -> rows_3x3
      sm$cols -> cols_3x3
      N_sub[rows_3x3,cols_3x3] -> N_3x3
      n_sub[rows_3x3,cols_3x3] -> n_3x3
      n_3x3/N_3x3 -> dat_3x3
      omni.check(N_3x3, n_3x3, n.iter=3000, CR=CR, check="double", use_cpp=use_cpp, adjust_extremes=adjust_extremes) -> out
      check_pass(dat_3x3, out$low, out$high) -> pass
      list(rows=rows_3x3, cols=cols_3x3, result=out, pass=pass) -> double_results[[key]]
      if (!pass) {
        FALSE -> double_pass_all
        list(rows=rows_3x3, cols=cols_3x3, observed=dat_3x3, low=out$low, high=out$high) -> double_failures[[key]]
      }
    }

    if (!double_pass_all) {
      return(list(
        rows=rows_4x4,
        cols=cols_4x4,
        status="FAILED_DOUBLE",
        message=paste0("Failed double cancellation on ", length(double_failures), " of 16 embedded 3x3 matrices"),
        single_result=single_result,
        single_failure=NULL,
        double_results=double_results,
        double_failures=double_failures,
        triple_result=NULL
      ))
    }

    omni.check(N_sub, n_sub, n.iter=3000, CR=CR, check="triple", use_cpp=use_cpp, adjust_extremes=adjust_extremes) -> out
    check_pass(dat_sub, out$low, out$high) -> triple_pass
    list(result=out, pass=triple_pass) -> triple_result

    if (!triple_pass) {
      return(list(
        rows=rows_4x4,
        cols=cols_4x4,
        status="FAILED_TRIPLE",
        message="Passed single and double, failed at triple cancellation",
        single_result=single_result,
        single_failure=NULL,
        double_results=double_results,
        double_failures=NULL,
        triple_result=triple_result,
        triple_failure=list(observed=dat_sub, low=out$low, high=out$high)
      ))
    }

    list(
      rows=rows_4x4,
      cols=cols_4x4,
      status="PASSED_ALL",
      message="Passed all cancellation tests (single, double, triple)",
      single_result=single_result,
      single_failure=NULL,
      double_results=double_results,
      double_failures=NULL,
      triple_result=triple_result,
      triple_failure=NULL
    )
  }

  proc.fun <- function(dummy, arg.list) {
    N <- arg.list[[1]]
    n <- arg.list[[2]]
    CR <- arg.list[[3]]
    use_cpp <- arg.list[[4]]
    submat_3x3 <- arg.list[[5]]
    run_hierarchical_test <- arg.list[[6]]
    omni.check <- arg.list[[7]]
    adjust_extremes <- arg.list[[8]]
    nrows <- nrow(N)
    ncols <- ncol(N)
    test <- 1
    max_tries <- 100
    tries <- 0
    while (test > 0 & tries < max_tries) {
      rows <- sort(sample(1:nrows, 4, replace=FALSE))
      cols <- sort(sample(1:ncols, 4, replace=FALSE))
      nt <- N[rows,cols]
      nc <- n[rows,cols]
      dat <- nc/nt
      if (adjust_extremes) test <- 0 else test <- sum(dat==1|dat==0)
      tries <- tries+1
    }
    if (test > 0 & !adjust_extremes) {
      return(list(rows=rows, cols=cols, status="SKIPPED", message="Could not find 4x4 without 0/1 proportions"))
    }
    tryCatch(
      run_hierarchical_test(nt, nc, rows, cols, CR, use_cpp, submat_3x3, adjust_extremes),
      error = function(e) list(rows=rows, cols=cols, status="ERROR",
                               message=paste("Hierarchical test failed:", conditionMessage(e)))
    )
  }

  proc.fun_adjacent <- function(dummy, arg.list) {
    r1 <- dummy[1]
    rows <- r1:(r1+3)
    c1 <- dummy[2]
    cols <- c1:(c1+3)
    N <- arg.list[[1]]
    n <- arg.list[[2]]
    CR <- arg.list[[3]]
    use_cpp <- arg.list[[4]]
    submat_3x3 <- arg.list[[5]]
    run_hierarchical_test <- arg.list[[6]]
    omni.check <- arg.list[[7]]
    adjust_extremes <- arg.list[[8]]
    nt <- N[rows,cols]
    nc <- n[rows,cols]
    tryCatch(
      run_hierarchical_test(nt, nc, rows, cols, CR, use_cpp, submat_3x3, adjust_extremes),
      error = function(e) list(rows=rows, cols=cols, status="ERROR",
                               message=paste("Hierarchical test failed:", conditionMessage(e)))
    )
  }

  arg.list <- list(N, n, CR, use_cpp, submat_3x3, run_hierarchical_test, omni.check, adjust_extremes)
  dummy <- list()
  if (n.mat=="adjacent") {
    for (i in 1:(nr-3)) for (j in 1:(nc-3)) c(i,j) -> dummy[[paste(i,j)]]
    results <- parallel::mclapply(dummy, proc.fun_adjacent, arg.list=arg.list, mc.cores=mc.cores)
  } else {
    for (i in 1:n.mat) dummy[[i]] <- i
    results <- parallel::mclapply(dummy, proc.fun, arg.list=arg.list, mc.cores=mc.cores)
  }

  statuses <- vapply(results, function(x) {
    if (is.list(x) && !is.null(x$status) && is.character(x$status)) x$status else "ERROR"
  }, character(1))
  summary <- list(
    n_tested=length(results),
    n_passed_all=sum(statuses=="PASSED_ALL"),
    n_failed_single=sum(statuses=="FAILED_SINGLE"),
    n_failed_double=sum(statuses=="FAILED_DOUBLE"),
    n_failed_triple=sum(statuses=="FAILED_TRIPLE"),
    n_skipped=sum(statuses=="SKIPPED"),
    n_error=sum(statuses=="ERROR")
  )

  list(
    N=N,
    n=n,
    summary=summary,
    results=results
  )
}
