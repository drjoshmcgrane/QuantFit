# Deep missing-data validation of select_model_hybrid (masked likelihood, MAR).
#
# Reuses the EXACT audit-grid datasets (same generator seeds as
# audit_selector.R, K = 20, all six truths, J = 8, N = 1500), applies MAR masks
# at 10% and 25%, and runs the hybrid on the masked data. Because the
# complete-data hybrid arm already exists in audit_k20_results/, every masked
# verdict pairs with its own complete-data verdict - the masking effect is
# isolated dataset by dataset, not confounded with draw noise.
#
# MAR mechanism: iid Bernoulli(rate) over cells (covariate-free MAR = MCAR,
# the clean baseline; the machinery assumes MAR, and rank-matched mask
# re-imposition in nulls is exercised identically either way).
suppressMessages(library(QuantFit))
K     <- as.integer(Sys.getenv("MDV_K", "20"))
B     <- as.integer(Sys.getenv("MDV_B", "49"))
cores <- as.integer(Sys.getenv("MDV_CORES", "8"))
out   <- Sys.getenv("MDV_OUT", "missing_data_out"); dir.create(out, showWarnings = FALSE)
models <- c("UN", "MON", "IIO", "DM", "LCR", "RM")
rates  <- c(0.10, 0.25)
grid <- expand.grid(model = models, rep = seq_len(K), rate = rates,
                    stringsAsFactors = FALSE)
# cost order: cheap ordinal cells first, quant-edge truths last
grid <- grid[order(match(grid$model, c("IIO", "MON", "UN", "DM", "LCR", "RM"))), ]

run <- function(k) {
  cs <- grid[k, ]
  f <- file.path(out, sprintf("m_%s_%02d_%02.0f.csv", cs$model, cs$rep,
                              100 * cs$rate))
  if (file.exists(f)) return(invisible())
  d <- simulate_responses(cs$model, n_persons = 1500, n_items = 8,
                          n_classes = 3, seed = 7000 * cs$rep +
                            match(cs$model, models))   # audit-grid seeds
  d <- if (is.list(d)) d$data else d; storage.mode(d) <- "integer"
  set.seed(880000L + k)
  d[matrix(stats::runif(length(d)) < cs$rate, nrow(d))] <- NA
  t0 <- proc.time()[3]
  r <- tryCatch(suppressMessages(
         select_model_hybrid(d, n_classes = 3L, B = B, lr_boot_n_starts = 2L,
                             mc.cores = 1L, seed = 1, verbose = FALSE)),
       error = function(e) NULL)
  g <- function(x, fld) if (is.null(x)) NA else x[[fld]]
  write.csv(data.frame(truth = cs$model, rep = cs$rep, rate = cs$rate,
    selected = if (is.null(r)) NA_character_ else r$selected,
    class_shape = g(r$poset$class, "shape"), item_shape = g(r$poset$item, "shape"),
    secs = round(proc.time()[3] - t0, 1)), f, row.names = FALSE)
}
cat("missing-data validation:", nrow(grid), "runs (B =", B, ")\n")
invisible(parallel::mclapply(seq_len(nrow(grid)),
  function(k) tryCatch(run(k), error = function(e) NULL),
  mc.cores = cores, mc.preschedule = FALSE))
cat("MDV DONE\n")
