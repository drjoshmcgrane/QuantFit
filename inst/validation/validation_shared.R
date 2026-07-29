# Shared infrastructure for ALL evidence runners (round 10): provenance
# verification, canary comparison, and the antichain generators - one
# definition, sourced everywhere, so runners cannot drift apart.

qf_provenance_check <- function(method_expected) {
  bi <- QuantFit::quantfit_build_info()
  if (is.na(bi$sha) || identical(bi$sha, "unstamped"))
    stop("installed QuantFit build is not stamped (GitSHA); use the ",
         "stamp-then-install workflow")
  if (!file.exists("DESCRIPTION"))
    stop("no DESCRIPTION here - evidence runners must run from the package ",
         "root (fail-closed)")
  dcf <- trimws(unname(read.dcf("DESCRIPTION", fields = "GitSHA")[1, 1]))
  if (is.na(dcf) || !nzchar(dcf) || identical(dcf, "unstamped"))
    stop("checkout DESCRIPTION carries no GitSHA stamp")
  sha <- trimws(as.character(bi$sha))
  if (!identical(sha, as.character(dcf)))
    stop("installed build (", sha, ") != checkout stamp (", dcf, "); reinstall")
  if (!identical(bi$po_method, method_expected))
    stop("installed poset method (", bi$po_method, ") != expected (",
         method_expected, ")")
  head_sha <- suppressWarnings(system("git rev-parse --short HEAD",
                intern = TRUE, ignore.stderr = TRUE))
  if (length(head_sha) != 1L || !nzchar(head_sha))
    stop("cannot resolve HEAD (not a git checkout?) - fail-closed")
  st <- suppressWarnings(system(paste("git diff --quiet", paste0(sha, "..HEAD"),
          "-- R src inst tests NAMESPACE"), ignore.stderr = TRUE))
  if (!identical(st, 0L))
    stop("tracked package/validation files changed since stamped build ",
         sha, " (HEAD ", head_sha, "); re-stamp and reinstall")
  # DIRTY-STATE check: staged, unstaged, and untracked files in relevant
  # paths all invalidate provenance (the runners themselves are read from
  # the working tree)
  dirty <- suppressWarnings(system(
    "git status --porcelain -- R src inst tests NAMESPACE DESCRIPTION",
    intern = TRUE, ignore.stderr = TRUE))
  if (length(dirty))
    stop("working tree is dirty in provenance-relevant paths:\n  ",
         paste(dirty, collapse = "\n  "), "\ncommit or remove before ",
         "generating evidence")
  list(sha = sha, head = head_sha, method = bi$po_method)
}

qf_canary_match <- function(prev, sha, method, config) {
  !is.null(prev) &&
    identical(as.character(prev$method), method) &&
    identical(as.character(prev$sha), sha) &&
    identical(as.character(prev$config), config)
}

# Antichain generators (XANTI = deep crossings; NEARANTI = profile scale
# TUNED by bisection so the minimum pairwise crossing depth is ~0.015, just
# above the eps = 0.01 tolerance). Deterministic given (tr, nI, rep, N).
qf_anti_seed <- function(tr, nI, rep, C = 3L)
  51000L + 977L * rep + 13L * nI + 7L * C +
    match(tr, c("PO","PO_INV","PO_ITEMS","UN","MON","IIO",
                "XANTI","NEARANTI","PO_POLY"))

qf_gen_anti <- function(tr, nI, NR, rep, C = 3L) {
  set.seed(qf_anti_seed(tr, nI, rep, C))
  if (tr == "XANTI") { sc <- 1.6; sh <- c(0.35, 0.15) } else {
    depth_at <- function(scl) {
      b <- seq(-scl, scl, length.out = nI)
      Pt <- rbind(plogis(b), plogis(rev(b) + 0.24 * scl),
                  plogis(b * rep_len(c(-1, 1), nI) + 0.12 * scl))
      m <- Inf
      for (a in 1:2) for (bb in (a + 1):3)
        m <- min(m, min(mean(pmax(0, Pt[a, ] - Pt[bb, ])),
                        mean(pmax(0, Pt[bb, ] - Pt[a, ]))))
      m
    }
    sc <- stats::uniroot(function(x) depth_at(x) - 0.015, c(0.02, 1.5))$root
    sh <- c(0.24, 0.12) * sc
  }
  base <- seq(-sc, sc, length.out = nI)
  P <- rbind(plogis(base), plogis(rev(base) + sh[1]),
             plogis(base * rep_len(c(-1, 1), nI) + sh[2]))
  cls <- sample.int(3L, NR, replace = TRUE)
  d <- matrix(rbinom(NR * nI, 1, P[cls, ]), NR, nI)
  storage.mode(d) <- "integer"
  attr(d, "params") <- list(L = qlogis(P))
  d
}
