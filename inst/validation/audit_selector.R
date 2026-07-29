# Same-code, same-data audit of ANY of the three selectors on the standard
# large audit grid (6 models x K reps, N=1500, J=8 dichotomous, n_classes=3).
# Seeds match selection_audit.R / hybrid_audit.R exactly, so all selectors see
# IDENTICAL datasets and results are paired.
#
#   SELECTOR=lattice   select_model_ll (LR-edge lattice)
#   SELECTOR=hybrid    select_model_hybrid (2x2 ordinal layer + LR quant edge)
#
# Purpose: replace the stale-code lattice figures with properly powered,
# same-code paired comparisons (the IIO claim is the load-bearing one).
# Resumable: one CSV per (selector, model, rep).
suppressMessages(library(QuantFit))
SEL     <- Sys.getenv("SELECTOR", "lattice")
K       <- as.integer(Sys.getenv("AUDIT_K", "30"))
B       <- as.integer(Sys.getenv("AUDIT_B", "49"))
cores   <- as.integer(Sys.getenv("AUDIT_CORES", "8"))
n_items <- 8L
# "PO" is the seventh truth: partial class order in TI&D's own idiom
# (simulate_responses("PO"), free items, margin 0.10 - a middle value; the
# margin sweep in po_validation.R maps the detection floor). SCORING for PO
# truth: correct iff selected == "UN" AND the class poset test rejects
# (class_shape == "partial") - the six-model `selected` cannot express PO, the
# calibrated refinement carries it, and BOTH selectors now run that test.
models  <- c("UN","MON","IIO","DM","LCR","RM","PO")
scale_of <- c(UN="nominal",MON="ordinal",IIO="ordinal",DM="ordinal",
              LCR="quant",RM="quant",PO="partial")
out <- Sys.getenv("AUD_OUT","audit_out"); dir.create(out, showWarnings=FALSE)
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
CONFIG <- sprintf("B=%d;K=%d;J=%d;N=1500;alpha=0.05;posetfloor=99",
                  B, K, n_items)
grid <- expand.grid(model=models, rep=seq_len(K), stringsAsFactors=FALSE)

run <- function(k) {
  cs <- grid[k,]; f <- file.path(out, sprintf("%s_%s_%03d.csv", SEL, cs$model, cs$rep))
  if (file.exists(f)) {
    prev <- tryCatch(read.csv(f, nrows = 1), error = function(e) NULL)
    if (!is.null(prev) && identical(prev$method, METHOD) &&
        identical(prev$sha, SHA) && identical(prev$config, CONFIG))
      return(invisible())
    unlink(f)                      # stale method/build: regenerate
  }
  d <- if (cs$model == "PO")
    simulate_responses("PO", n_persons=1500, n_items=n_items, n_classes=3,
                       poset="V", po_margin=0.10, seed=7000*cs$rep + 7L)
  else simulate_responses(cs$model, n_persons=1500, n_items=n_items,
                          n_classes=3, seed=7000*cs$rep + match(cs$model, models))
  d <- if (is.list(d)) d$data else d; storage.mode(d) <- "integer"
  r <- tryCatch(switch(SEL,
      lattice  = select_model_ll(d, n_classes=2:6, B=B, n_starts=5,
                   boot_n_starts=2L, method="lattice", seed=1, mc.cores=1L,
                   verbose=FALSE),
      hybrid   = select_model_hybrid(d, n_classes=3L, B=B,
                   lr_boot_n_starts=2L, mc.cores=1L, seed=1, verbose=FALSE)),
    error=function(e) NULL)
  sel <- if (is.null(r)) NA_character_ else r$selected
  g <- function(x, fld) if (is.null(x) || is.null(x[[fld]])) NA else x[[fld]]
  write.csv(data.frame(method=METHOD, sha=SHA, config=CONFIG,
    selector=SEL, truth=cs$model, truth_scale=scale_of[cs$model],
    rep=cs$rep, selected=sel,
    selected_scale=if (is.na(sel)) NA else scale_of[sel],
    class_shape=g(r$poset$class,"shape"), class_comp=g(r$poset$class,"comparable"),
    class_type=g(r$poset$class,"type"),
    item_shape=g(r$poset$item,"shape"), item_comp=g(r$poset$item,"comparable")),
    f, row.names=FALSE)
}
cat("Audit [", SEL, "]:", nrow(grid), "datasets (J=8, N=1500, K=", K, ", B=", B, ")\n")
invisible(parallel::mclapply(seq_len(nrow(grid)),
  function(k) tryCatch(run(k), error=function(e) NULL), mc.cores=cores))
cat("AUDIT DONE\n")
