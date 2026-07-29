# Shared provenance verification for all evidence runners (round 9).
# FAIL-CLOSED: every check that cannot be performed stops the run.
qf_provenance_check <- function(method_expected) {
  bi <- QuantFit::quantfit_build_info()
  if (is.na(bi$sha) || identical(bi$sha, "unstamped"))
    stop("installed QuantFit build is not stamped (GitSHA); use the ",
         "stamp-then-install workflow")
  if (!file.exists("DESCRIPTION"))
    stop("no DESCRIPTION here - evidence runners must run from the package ",
         "root (fail-closed; the checkout comparison cannot be skipped)")
  dcf <- trimws(unname(read.dcf("DESCRIPTION", fields = "GitSHA")[1, 1]))
  if (is.na(dcf) || !nzchar(dcf) || identical(dcf, "unstamped"))
    stop("checkout DESCRIPTION carries no GitSHA stamp")
  sha <- trimws(as.character(bi$sha))
  if (!identical(sha, as.character(dcf)))
    stop("installed build (", sha, ") != checkout stamp (", dcf,
         "); reinstall")
  if (!identical(bi$po_method, method_expected))
    stop("installed poset method (", bi$po_method, ") != expected (",
         method_expected, ")")
  head_sha <- suppressWarnings(system("git rev-parse --short HEAD",
                intern = TRUE, ignore.stderr = TRUE))
  if (length(head_sha) != 1L || !nzchar(head_sha))
    stop("cannot resolve HEAD (not a git checkout?) - fail-closed")
  # no tracked package/validation file may have changed since the stamped
  # code commit (the stamp commit itself only touches DESCRIPTION)
  st <- suppressWarnings(system(paste("git diff --quiet", paste0(sha, "..HEAD"),
          "-- R src inst tests NAMESPACE"), ignore.stderr = TRUE))
  if (!identical(st, 0L))
    stop("tracked package/validation files changed since stamped build ",
         sha, " (HEAD ", head_sha, "); re-stamp and reinstall before ",
         "generating evidence")
  list(sha = sha, head = head_sha, method = bi$po_method)
}

# resume comparison must survive read.csv type coercion (a numeric-only SHA
# imports as integer)
qf_canary_match <- function(prev, sha, method, config) {
  !is.null(prev) &&
    identical(as.character(prev$method), method) &&
    identical(as.character(prev$sha), sha) &&
    identical(as.character(prev$config), config)
}
