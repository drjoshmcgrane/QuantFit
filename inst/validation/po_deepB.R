# Deep-B decision-stability anchors (review round 6): rerun a fixed anchor
# subset of the v4 grid with poset_B = 499 and compare every demonstrated-pair
# decision against the B = 99 run. The simultaneous quantile is a top-order
# statistic; this measures how much the decisions move as it deepens.
suppressMessages(library(QuantFit))
`%||%` <- function(a, b) if (is.null(a)) b else a
out <- Sys.getenv("PODEEP_OUT", "po_deepB_out"); dir.create(out, showWarnings = FALSE)
METHOD <- "max-T-studentized-v4"
source("inst/validation/validation_shared.R")
prov <- qf_provenance_check(METHOD)
SHA <- prov$sha
HEAD_SHA <- prov$head
anchors <- expand.grid(truth = c("PO", "XANTI", "NEARANTI", "PO_ITEMS"),
                       rep = 1:8, stringsAsFactors = FALSE)
gen <- function(tr, rep) {
  # SAME seeds as the v4 grid for the corresponding cells
  NR <- 1500L
  if (tr == "PO") return(simulate_responses("PO", n_persons = NR, n_items = 12,
    n_classes = 3, poset = "V", po_margin = 0.10,
    seed = 51000L + 977L*rep + 13L*12L + 10L + 21L + 1L))
  if (tr == "PO_ITEMS") return(simulate_responses("PO_ITEMS", n_persons = NR,
    n_items = 8, n_classes = 4, poset = "layers2", po_margin = NULL,
    seed = 51000L + 977L*rep + 13L*8L + 28L + 3L))
  # XANTI / NEARANTI come from the SHARED generator (identical seeds and
  # construction to the v4 grid, including the tuned 0.015 boundary)
  qf_gen_anti(tr, 8L, NR, rep, 3L)
}
res <- parallel::mclapply(seq_len(nrow(anchors)), function(k) tryCatch({
  cs <- anchors[k, ]
  d <- gen(cs$truth, cs$rep)
  C <- if (cs$truth == "PO_ITEMS") 4L else 3L
  pstr <- function(x) if (is.null(x) || !nrow(x$pairs)) "" else
    paste(sprintf("%d>%d", x$pairs$dominant, x$pairs$dominated), collapse = ";")
  one <- function(pb) {
    r <- withCallingHandlers(
      select_model_hybrid(d, n_classes = C, B = 49, poset_B = pb,
                          mc.cores = 1L, seed = 1, verbose = FALSE),
      warning = function(w) {
        if (grepl("poset refinement failed|quantitative edge|RM-vs-LCR",
                  conditionMessage(w)))
          stop("selector warning treated as failure: ", conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    if (is.null(r) || is.na(r$selected)) stop("selector returned no verdict")
    # the anchor's EXPECTED side must be present with adequate depth -
    # unavailable comparisons must never count as stable
    side <- if (cs$truth == "PO_ITEMS") r$poset$item else r$poset$class
    if (is.null(side))
      stop("expected poset side missing (truth ", cs$truth, ", selected ",
           r$selected, ")")
    if (is.null(side$b_eff) || side$b_eff < max(99L, floor(0.9 * pb)))
      stop("poset b_eff ", side$b_eff, " below requirement for poset_B ", pb)
    if (is.na(side$shape)) stop("poset shape unavailable")
    c(cls_shape = r$poset$class$shape %||% NA,
      cls_pairs = pstr(r$poset$class),
      it_shape  = r$poset$item$shape %||% NA,
      it_pairs  = pstr(r$poset$item))
  }
  a <- one(99L); b <- one(499L)
  data.frame(method = "max-T-studentized-v4", sha = SHA, head = HEAD_SHA,
             config = "B=49;posetB=99v499;N=1500", truth = cs$truth,
             rep = cs$rep,
             shape99 = a["cls_shape"], shape499 = b["cls_shape"],
             pairs99 = a["cls_pairs"], pairs499 = b["cls_pairs"],
             ishape99 = a["it_shape"], ishape499 = b["it_shape"],
             ipairs99 = a["it_pairs"], ipairs499 = b["it_pairs"],
             stable = identical(a, b))
}, error = function(e) sprintf("anchor %d (%s rep %d): %s", k,
       anchors$truth[k], anchors$rep[k], conditionMessage(e))),
  mc.cores = 8, mc.preschedule = FALSE)
ok <- vapply(res, function(x) is.data.frame(x), logical(1))
fails <- res[!ok]
if (length(fails)) {
  write.csv(data.frame(error = vapply(fails, as.character, character(1))),
            file.path(out, "FAILURES.csv"), row.names = FALSE)
  cat("DEEPB FAILED:", length(fails), "anchors errored\n")
  quit(save = "no", status = 1L)
}
d <- do.call(rbind, res)
write.csv(d, file.path(out, "deepB_anchors.csv"), row.names = FALSE)
cat(sprintf("deep-B anchors: %d/%d decisions identical at B 99 vs 499\n",
            sum(d$stable), nrow(d)))
print(d[!d$stable, c("truth","rep","shape99","shape499","pairs99","pairs499")])
cat("DEEPB DONE\n")
