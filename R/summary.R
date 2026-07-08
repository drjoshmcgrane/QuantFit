#' Summarize checks produced by ConjointChecks
#'
#' Takes output from [ConjointChecks()] and produces summary measures of the
#' reported violations.
#'
#' @param object Object returned by [ConjointChecks()] of class
#'   [`checks`][checks-class].
#' @param ... further arguments passed to or from other methods.
#'
#' @return A list with the weighted and unweighted mean violation rates
#'   (`Means`) and the mean violation rate per item (`items`).
#'
#' @author Ben Domingue \email{ben.domingue@@gmail.com}
#'
#' @seealso [ConjointChecks()], [plot.checks()]
#'
#' @examples
#' \dontrun{
#' summary(rasch1000)
#' }
#'
#' @method summary checks
#' @export
summary.checks<-function(object, ...) {
  list(Means=object@means,items=colMeans(object@tab,na.rm=TRUE))
}
