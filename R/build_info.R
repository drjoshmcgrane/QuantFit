# Partial-order method version: bump whenever the poset inference changes.
.qf_po_method <- "max-T-studentized-v4"

#' Installed-build metadata
#'
#' Returns the identity of the INSTALLED QuantFit build: package version, the
#' git SHA stamped into DESCRIPTION at build time (`GitSHA:` field; NA when
#' the build was not stamped), and the partial-order method version constant.
#' Validation runners compare this against the expected checkout state before
#' generating or reusing results, so a stale installed package can never
#' produce rows labelled with a newer commit.
#'
#' @return A list with `version`, `sha`, and `po_method`.
#' @export
quantfit_build_info <- function() {
  d <- utils::packageDescription("QuantFit")
  sha <- d[["GitSHA"]]
  list(version = as.character(d[["Version"]]),
       sha = if (is.null(sha)) NA_character_ else sha,
       po_method = .qf_po_method)
}
