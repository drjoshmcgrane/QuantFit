# Deep-B decision-stability anchors (review round 6): rerun a fixed anchor
# subset of the v4 grid with poset_B = 499 and compare every demonstrated-pair
# decision against the B = 99 run. The simultaneous quantile is a top-order
# statistic; this measures how much the decisions move as it deepens.
suppressMessages(library(QuantFit))
`%||%` <- function(a, b) if (is.null(a)) b else a
out <- Sys.getenv("PODEEP_OUT", "po_deepB_out"); dir.create(out, showWarnings = FALSE)
METHOD <- "max-T-studentized-v4"
# BUILD VERIFICATION (round 8): the row SHA is the INSTALLED package's
# stamped build SHA, and generation refuses to run when the installed build
# does not match this checkout's stamp or the expected method version.
bi <- quantfit_build_info()
dcf <- tryCatch(trimws(unname(read.dcf("DESCRIPTION",
                fields = "GitSHA")[1, 1])), error = function(e) NA_character_)
if (is.na(bi$sha) || identical(bi$sha, "unstamped"))
  stop("installed QuantFit build is not stamped (GitSHA); rebuild via the ",
       "stamp-then-install workflow before generating evidence")
if (!is.na(dcf) && !identical(dcf, "unstamped") &&
    !identical(trimws(unname(bi$sha)), dcf))
  stop("installed QuantFit build (", bi$sha, ") does not match this ",
       "checkout's stamp (", dcf, "); reinstall before generating evidence")
if (!identical(bi$po_method, METHOD))
  stop("installed QuantFit poset method (", bi$po_method, ") != expected (",
       METHOD, ")")
SHA <- bi$sha
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
  sd <- 51000L + 977L*rep + 13L*8L + 21L + match(tr, c("PO","PO_INV","PO_ITEMS",
        "UN","MON","IIO","XANTI","NEARANTI","PO_POLY"))
  set.seed(sd)
  sc <- if (tr == "XANTI") 1.6 else 0.5
  sh <- if (tr == "XANTI") c(0.35, 0.15) else c(0.12, 0.06)
  base <- seq(-sc, sc, length.out = 8)
  P <- rbind(plogis(base), plogis(rev(base) + sh[1]),
             plogis(base * rep_len(c(-1, 1), 8) + sh[2]))
  cls <- sample.int(3L, NR, replace = TRUE)
  d <- matrix(rbinom(NR * 8, 1, P[cls, ]), NR, 8)
  storage.mode(d) <- "integer"; d
}
res <- parallel::mclapply(seq_len(nrow(anchors)), function(k) {
  cs <- anchors[k, ]
  d <- gen(cs$truth, cs$rep)
  C <- if (cs$truth == "PO_ITEMS") 4L else 3L
  pstr <- function(x) if (is.null(x) || !nrow(x$pairs)) "" else
    paste(sprintf("%d>%d", x$pairs$dominant, x$pairs$dominated), collapse = ";")
  one <- function(pb) {
    r <- suppressWarnings(select_model_hybrid(d, n_classes = C, B = 49,
           poset_B = pb, mc.cores = 1L, seed = 1, verbose = FALSE))
    c(cls_shape = r$poset$class$shape %||% NA,
      cls_pairs = pstr(r$poset$class),
      it_shape  = r$poset$item$shape %||% NA,
      it_pairs  = pstr(r$poset$item))
  }
  a <- one(99L); b <- one(499L)
  data.frame(sha = SHA, truth = cs$truth, rep = cs$rep,
             shape99 = a["cls_shape"], shape499 = b["cls_shape"],
             pairs99 = a["cls_pairs"], pairs499 = b["cls_pairs"],
             ishape99 = a["it_shape"], ishape499 = b["it_shape"],
             ipairs99 = a["it_pairs"], ipairs499 = b["it_pairs"],
             stable = identical(a, b))
}, mc.cores = 8, mc.preschedule = FALSE)
d <- do.call(rbind, res)
write.csv(d, file.path(out, "deepB_anchors.csv"), row.names = FALSE)
cat(sprintf("deep-B anchors: %d/%d decisions identical at B 99 vs 499\n",
            sum(d$stable), nrow(d)))
print(d[!d$stable, c("truth","rep","shape99","shape499","pairs99","pairs499")])
cat("DEEPB DONE\n")
