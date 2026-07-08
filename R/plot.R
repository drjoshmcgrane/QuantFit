#' Plot checks produced by ConjointChecks
#'
#' Takes output from [ConjointChecks()] and produces a
#' [matplot][graphics::matplot] showing the percentage of reported violations
#' at each cell.
#'
#' @param x Object returned by [ConjointChecks()] of class
#'   [`checks`][checks-class].
#' @param items Vector of item numbers to include in a single plot. Defaults
#'   to all, but this is less helpful for diagnostic purposes.
#' @param item.labels Should item numbers be included? Defaults to `TRUE`. If
#'   length of `items` is unity (perhaps if the small multiple format of
#'   Tufte, 2001 is going to be used), then the item number gets printed
#'   below the x-axis. If the length of `items` is more than unity, the item
#'   number gets plotted in the figure above the largest proportion of
#'   violations for each item.
#' @param ... further arguments passed to or from other methods.
#'
#' @return No return value, called for side effects.
#'
#' @references
#' Tufte, E. R. (2001). The visual display of quantitative information
#' (2nd ed.). Cheshire, CT: Graphics Press.
#'
#' @author Ben Domingue \email{ben.domingue@@gmail.com}
#'
#' @seealso [ConjointChecks()], [summary.checks()]
#'
#' @examples
#' \dontrun{
#' opar <- par()
#' par(mfrow=c(3,2))
#' plot(rasch1000)
#' plot(rasch1000,items=c(5,10,15))
#' for (i in c(3,9,13,18)) plot(rasch1000,items=i)
#' par(opar)
#' }
#'
#' @importFrom graphics matplot par
#' @method plot checks
#' @export
plot.checks<-function(x, items=NULL, item.labels=TRUE, ...) {
  #graphical args
  modify.args<-function(...,items) {
    list(...)->args
    if (is.null(args$ylim)) args$ylim<-c(0,1)
    if (is.null(args$xlab)) args$xlab<-""
    if (is.null(args$ylab)) args$ylab<-"Proportion Violations"
    if (is.null(args$type)) args$type<-"l"
    if (is.null(args$lty)) args$lty<-rep(1,length(items))
    if (is.null(args$col)) args$col<-"black"
    if (is.null(args$xaxt)) args$xaxt<-"n"
    args
  }
  #
  if (is.null(items)) items <- 1:ncol(x@tab)
  #
  modify.args(...,items=items)->args
  x@tab[,items]->dat
  c(list(y=dat),args)->list.of.args
  do.call("matplot",list.of.args)
  mtext(side=1,line=1,paste("Increasing Ability"))
  if (item.labels) {
    if (length(items)==1) mtext(side=1,line=2,paste("Item",items)) else {
      apply(x@tab,2,which.max)->maxes.x
      apply(x@tab,2,max,na.rm=TRUE)->maxes.y
      for (i in 1:length(items)) text(maxes.x[items],maxes.y[items],items,pos=3)
    }
  }
}
