#' Class "list.null"
#'
#' The formal S4 class union for `list.null`. This class contains either a
#' list or `NULL` and is used internally to hold the results of the
#' individual cancellation checks.
#'
#' @details Objects of this class are used internally.
#'
#' @author Ben Domingue \email{ben.domingue@@gmail.com}
#'
#' @seealso [ConjointChecks()], [summary.checks()], [plot.checks()]
#'
#' @keywords classes
#' @name list.null-class
#' @aliases list.null
#' @exportClass list.null
setClassUnion("list.null", c("list", "NULL"))

#' Class "checks"
#'
#' The formal S4 class for checks. This class contains transformed versions
#' of the raw response data as well as summaries of the checks.
#'
#' @details Objects of class `checks` contain all information returned by
#' [ConjointChecks()] (and by the related wrappers [SingleCancel()],
#' [DoubleCancel()], and [TripleCancel()]).
#'
#' @slot N matrix containing the number of respondents at each item/ability
#'   intersection.
#' @slot n matrix containing the number of correct responses at each
#'   item/ability intersection.
#' @slot Checks list containing information about each checked submatrix.
#' @slot tab matrix containing information about the detected violations at
#'   each item/ability intersection.
#' @slot means vector containing weighted and unweighted means for the
#'   detected violations (where weights are the number of individuals at each
#'   ability level).
#' @slot check.counts matrix giving the number of times an item/ability cell
#'   was sampled.
#'
#' @author Ben Domingue \email{ben.domingue@@gmail.com}
#'
#' @seealso [ConjointChecks()], [summary.checks()], [plot.checks()]
#'
#' @keywords classes
#' @name checks-class
#' @aliases checks
#' @exportClass checks
setClass("checks", representation(N = "matrix", n = "matrix", Checks = "list.null", tab = "matrix", means = "list", check.counts = "matrix"))
