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
  if (ss.lower==1) {
    message("ss.lower must be greater than 1, setting to 2.")
    ss.lower<-2
  }
  ncol(resp)->n.items
  # Missing responses are allowed: every cell's count N(s, j) is the number of
  # respondents in score group s who ANSWERED item j (n = number correct among
  # them), so each cell proportion is weighted by its own number of
  # observations - the same weighting the adjacent-category polytomous matrix
  # uses for out-of-play cells. The score is the number correct among answered
  # items. Valid under MAR; the bootstrap null must impose the same
  # missingness pattern so pipeline effects cancel.
  has_na <- any(is.na(resp))
  #reorder columns
  colSums(resp, na.rm=TRUE)->cs
  resp[,order(cs)]->resp
  #group by sum scores (number correct among answered items)
  rowSums(resp, na.rm=TRUE)->rs
  n<-N<-list()
  table(rs)->tab
  as.numeric(names(tab))[tab>=ss.lower]->lev
  for (s in lev) {
    resp[rs==s,,drop=FALSE]->tmp
    if (has_na) {
      colSums(!is.na(tmp))->N[[as.character(s)]]
    } else {
      rep(nrow(tmp),n.items)->N[[as.character(s)]]
    }
    colSums(tmp, na.rm=TRUE)->n[[as.character(s)]]
  }
  do.call("rbind",N)->N
  do.call("rbind",n)->n
  if (has_na && any(N==0)) {
    keep <- rowSums(N==0)==0
    if (sum(keep) < 3L) stop("Too many empty cells after missing-data ",
                             "weighting; not enough complete score groups.")
    N <- N[keep,,drop=FALSE]; n <- n[keep,,drop=FALSE]
  }
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
