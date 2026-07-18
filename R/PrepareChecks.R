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
#' @param person_order How respondents are ordered/grouped when responses are
#'   missing. `"complete"` (default) drops incomplete respondents - the
#'   assumption-free frame. `"facility"` keeps everyone, grouping by
#'   proportion correct among answered items (difficulty-blind).
#'   `"adjusted"` groups by the observed count standardized against the
#'   facility-implied expectation for the answered items (difficulty-aware,
#'   imports a metric commensuration). Complete data are identical under all
#'   three; groups never split tied respondents.
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
PrepareChecks<-function(resp,ss.lower=10,collapse.columns=FALSE,
                        person_order=c("complete","facility","adjusted")) {
  person_order <- match.arg(person_order)
  if (ss.lower==1) {
    message("ss.lower must be greater than 1, setting to 2.")
    ss.lower<-2
  }
  ncol(resp)->n.items
  # Missing responses: the DEFAULT conditioning frame is complete cases - the
  # only assumption-free frame, since any person ordering across different
  # answered item sets imports extra-ordinal structure. person_order
  # "facility" (proportion correct among answered; difficulty-blind) and
  # "adjusted" (facility-anchored standardized count; difficulty-aware but
  # metric-flavoured) keep all respondents as documented approximations.
  # Cells stay observation-weighted throughout: N(s, j) counts respondents in
  # group s who ANSWERED item j. Valid under MAR; the bootstrap null must
  # impose the same missingness pattern so pipeline effects cancel.
  if (anyNA(resp) && person_order=="complete") {
    cc <- stats::complete.cases(resp)
    message("PrepareChecks: dropping ", sum(!cc), " of ", nrow(resp),
            " incomplete respondents (person_order=\"complete\"); use ",
            "person_order=\"facility\" or \"adjusted\" to keep them.")
    resp <- resp[cc,,drop=FALSE]
    if (nrow(resp) < 3L*ss.lower) stop("Too few complete cases; consider ",
                                       "person_order=\"facility\"/\"adjusted\".")
  }
  has_na <- any(is.na(resp))
  #reorder columns by weighted facility (identical order on complete data);
  #stamp original item ids first so the reordered columns stay identifiable
  if (is.null(colnames(resp))) colnames(resp) <- paste0("I", seq_len(n.items))
  fac <- colSums(resp,na.rm=TRUE)/pmax(colSums(!is.na(resp)),1L)
  resp[,order(fac)]->resp
  #person ordering
  if (!has_na) {
    rowSums(resp)->rs
  } else if (person_order=="facility") {
    rowMeans(resp,na.rm=TRUE)->rs
  } else {                                # adjusted
    p <- colMeans(resp,na.rm=TRUE)
    obs <- !is.na(resp)
    e <- as.vector(obs %*% p)
    v <- as.vector(obs %*% (p*(1-p)))
    rs <- (rowSums(resp,na.rm=TRUE)-e)/sqrt(pmax(v,1e-12))
  }
  n<-N<-list()
  if (!has_na) {
    # sum-score groups (discrete): keep value-groups with >= ss.lower persons
    table(rs)->tab
    as.numeric(names(tab))[tab>=ss.lower]->lev
    grp <- rs; glev <- lev; glab <- as.character(lev)
  } else {
    # continuous ordering: tie-preserving adjacent value-bands, targeting the
    # granularity a sum-score grouping would give (~J+1 bands) subject to the
    # ss.lower minimum - NOT bands of exactly ss.lower, which would fragment
    # n respondents into hundreds of noisy slivers
    vals<-sort(unique(rs)); cnt<-as.integer(table(factor(rs,levels=vals)))
    tgt <- max(ss.lower, ceiling(sum(cnt)/(n.items+1)))
    gid<-integer(length(vals)); g<-1L; acc<-0L
    for (k in seq_along(vals)) {
      gid[k]<-g; acc<-acc+cnt[k]
      if (acc>=tgt && k<length(vals)) { g<-g+1L; acc<-0L }
    }
    if (acc<tgt && acc>0L && g>1L) gid[gid==g]<-g-1L  # fold thin tail
    grp <- gid[match(rs,vals)]
    glev <- sort(unique(grp))
    glab <- vapply(glev, function(s) as.character(round(mean(rs[grp==s]),4)),
                   character(1))
  }
  for (k in seq_along(glev)) {
    s <- glev[k]
    resp[grp==s,,drop=FALSE]->tmp
    if (has_na) {
      colSums(!is.na(tmp))->N[[glab[k]]]
    } else {
      rep(nrow(tmp),n.items)->N[[glab[k]]]
    }
    colSums(tmp, na.rm=TRUE)->n[[glab[k]]]
  }
  do.call("rbind",N)->N
  do.call("rbind",n)->n
  colnames(N)<-colnames(n)<-colnames(resp)
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
