# Missing-data validation of select_model_hybrid, v2 (post-audit-fix code:
# mean-based mask rank-matching, IIO-only axis null, practical-equivalence
# retention guard).
#
# Reuses the EXACT audit-grid datasets (same generator seeds as
# audit_selector.R, K = 20, six truths, J = 8, N = 1500) so every masked
# verdict pairs with its complete-data counterpart. THREE mechanisms:
#
#   mcar        iid cellwise Bernoulli at 10% and 25%
#   mar_ability GENUINE MAR: odd items never masked (anchor half); even items
#               masked with prob plogis(-1.35 - 1.2 * z(anchor score)) - low
#               scorers skip more; ~18% overall. Depends only on observed
#               responses, so MAR holds by construction. This is the mechanism
#               the rank-matched mask transfer exists to serve.
#   mar_item    item-dependent rates seq(0.05, 0.45) across items (~25%
#               overall), person-independent.
suppressMessages(library(QuantFit))
K     <- as.integer(Sys.getenv("MDV_K", "20"))
B     <- as.integer(Sys.getenv("MDV_B", "49"))
cores <- as.integer(Sys.getenv("MDV_CORES", "8"))
out   <- Sys.getenv("MDV_OUT", "missing_data_out2"); dir.create(out, showWarnings = FALSE)
models <- c("UN", "MON", "IIO", "DM", "LCR", "RM")
grid <- rbind(
  expand.grid(model = models, rep = seq_len(K), mech = "mcar",
              rate = c(0.10, 0.25), stringsAsFactors = FALSE),
  expand.grid(model = models, rep = seq_len(K), mech = "mar_ability",
              rate = NA, stringsAsFactors = FALSE),
  expand.grid(model = models, rep = seq_len(K), mech = "mar_item",
              rate = NA, stringsAsFactors = FALSE))
grid <- grid[order(match(grid$model, c("IIO", "MON", "UN", "DM", "LCR", "RM"))), ]

apply_mask <- function(d, mech, rate, seed) {
  set.seed(seed); N <- nrow(d); J <- ncol(d)
  if (mech == "mcar") {
    d[matrix(stats::runif(N * J) < rate, N, J)] <- NA
  } else if (mech == "mar_ability") {
    anchor <- seq(1L, J, by = 2L); target <- seq(2L, J, by = 2L)
    z <- scale(rowSums(d[, anchor, drop = FALSE]))[, 1]
    pm <- stats::plogis(-1.35 - 1.2 * z)
    for (j in target) d[stats::runif(N) < pm, j] <- NA
  } else {
    rates <- seq(0.05, 0.45, length.out = J)
    for (j in seq_len(J)) d[stats::runif(N) < rates[j], j] <- NA
  }
  d
}

run <- function(k) {
  cs <- grid[k, ]
  f <- file.path(out, sprintf("m_%s_%02d_%s%s.csv", cs$model, cs$rep, cs$mech,
                 ifelse(is.na(cs$rate), "", sprintf("_%02.0f", 100 * cs$rate))))
  if (file.exists(f)) return(invisible())
  d <- simulate_responses(cs$model, n_persons = 1500, n_items = 8,
                          n_classes = 3, seed = 7000 * cs$rep +
                            match(cs$model, models))
  d <- if (is.list(d)) d$data else d; storage.mode(d) <- "integer"
  d <- apply_mask(d, cs$mech, cs$rate, 880000L + k)
  t0 <- proc.time()[3]
  r <- tryCatch(suppressMessages(suppressWarnings(
         select_model_hybrid(d, n_classes = 3L, B = B, lr_boot_n_starts = 2L,
                             mc.cores = 1L, seed = 1, verbose = FALSE))),
       error = function(e) NULL)
  g <- function(x, fld) if (is.null(x)) NA else x[[fld]]
  write.csv(data.frame(truth = cs$model, rep = cs$rep, mech = cs$mech,
    rate = cs$rate, miss_frac = round(mean(is.na(d)), 3),
    selected = if (is.null(r)) NA_character_ else r$selected,
    edge_failed = if (is.null(r)) NA else isTRUE(r$quant_edge_failed),
    class_shape = g(r$poset$class, "shape"), item_shape = g(r$poset$item, "shape"),
    secs = round(proc.time()[3] - t0, 1)), f, row.names = FALSE)
}
cat("missing-data validation v2:", nrow(grid), "runs (B =", B, ")\n")
invisible(parallel::mclapply(seq_len(nrow(grid)),
  function(k) tryCatch(run(k), error = function(e) NULL),
  mc.cores = cores, mc.preschedule = FALSE))
cat("MDV2 DONE\n")
