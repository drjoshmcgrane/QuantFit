#' Prepare raw response data for ConjointChecks
#'
#' Takes raw dichotomously coded response data and produces the `n` and `N`
#' matrices required by [ConjointChecks()] and related functions. Items
#' (columns) are ordered by total score and respondents are grouped by sum
#' score; sum-score groups with fewer than `ss.lower` respondents are
#' dropped.
#'
#' @param resp Raw dichotomously coded response data. Columns represent items
#'   and rows represent individuals.
#' @param ss.lower Only sum scores that have at least this many distinct
#'   individuals with that sum score will be used.
#' @param collapse.columns Sum over columns (items) that have the same total
#'   score.
#'
#' @return Returns a list with elements `N` and `n`, containing respectively
#'   the number of total responses and the number of correct responses for
#'   each cell.
#'
#' @author Ben Domingue \email{ben.domingue@@gmail.com}
#'
#' @seealso [ConjointChecks()], [HiConjointChecks()], [KaraChecks()]
#'
#' @examples
#' # simulated Rasch example
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
#' # now check
#' PrepareChecks(resp)->obj
#'
#' @export
PrepareChecks<-function(resp,ss.lower=10,collapse.columns=FALSE) {
  if (any(is.na(resp))) stop("Checks will only work with complete data. Suggestion: remove respondents with missing responses.")
  if (ss.lower==1) {
    message("ss.lower must be greater than 1, setting to 2.")
    ss.lower<-2
  }
  ncol(resp)->n.items
  #reorder columns
  colSums(resp)->cs
  resp[,order(cs)]->resp
  #group by sum scores
  rowSums(resp)->rs
  n<-N<-list()
  table(rs)->tab
  as.numeric(names(tab))[tab>=ss.lower]->lev
  for (s in lev) {
    resp[rs==s,]->tmp
    rep(nrow(tmp),n.items)->N[[as.character(s)]]
    colSums(tmp)->n[[as.character(s)]]
  }
  do.call("rbind",N)->N
  do.call("rbind",n)->n
  if (collapse.columns) {
    colSums(n)->cs
    sort(unique(cs))->cs.index
    n2<-N2<-list()
    for (i in 1:length(cs.index)) {
      cs.index[i]->lev
      n[,cs==lev]->tmp
      if (is(tmp,"matrix")) rowSums(tmp)->tmp
      tmp->n2[[i]]
      N[,cs==lev]->tmp
      if (is(tmp,"matrix")) rowSums(tmp)->tmp
      tmp->N2[[i]]
    }      
    do.call("cbind",n2)->n
    do.call("cbind",N2)->N
  }
  list(N=N,n=n)
}
