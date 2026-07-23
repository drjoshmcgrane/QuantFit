# Targeted long-test anchor for the automated TI&D quantitative lattice.
# Usage: Rscript bridge_anchor.R MODEL J [B] [INNER_CORES] [N]
suppressMessages(library(QuantFit))
a <- commandArgs(trailingOnly = TRUE)
model <- toupper(a[1])
J <- as.integer(a[2])
B <- if (length(a) >= 3L) as.integer(a[3]) else 19L
inner_cores <- if (length(a) >= 4L) as.integer(a[4]) else 4L
N <- if (length(a) >= 5L) as.integer(a[5]) else 1500L
stopifnot(model %in% c("LCR", "RM"), is.finite(J), J >= 6L,
          is.finite(B), B >= 19L, is.finite(inner_cores), inner_cores >= 1L,
          is.finite(N), N >= 100L)

out <- "bridge_anchor_tid032"
dir.create(out, showWarnings = FALSE)
f <- file.path(out, sprintf("A_%s_J%02d_B%03d_N%04d.csv", model, J, B, N))
if (file.exists(f)) {
  cat("ANCHOR EXISTS", f, "\n")
  quit(save = "no", status = 0L)
}

seed_data <- 880000L + J * 101L + match(model, c("LCR", "RM")) * 10007L
set.seed(seed_data)
if (model == "RM") {
  beta <- stats::runif(J, -2, 2)
  theta <- stats::rnorm(N)
  d <- matrix(stats::rbinom(N * J, 1,
    stats::plogis(outer(theta, beta, "-"))), N, J)
} else {
  d <- simulate_responses("LCR", n_persons = N, n_items = J,
                          n_classes = 3L, seed = seed_data)
  if (is.list(d)) d <- d$data
}
storage.mode(d) <- "integer"

t0 <- proc.time()[3]
sel <- select_model_ll(d, n_classes = 1:5, B = B, n_starts = 5,
  boot_n_starts = 5, method = "lattice", severity = FALSE,
  seed = 990000L + J * 103L + match(model, c("LCR", "RM")) * 10009L,
  mc.cores = inner_cores)
edge_p <- function(label) {
  z <- sel$tests$p_value[sel$tests$comparison == label]
  if (length(z)) z[[1L]] else NA_real_
}
edge_raw <- function(label) {
  z <- sel$tests$raw_statistic[sel$tests$comparison == label]
  if (length(z)) z[[1L]] else NA_real_
}
edge_pre <- function(label) {
  z <- sel$tests$pre_refinement_statistic[sel$tests$comparison == label]
  if (length(z)) z[[1L]] else NA_real_
}
edge_warm <- function(label) {
  z <- sel$tests$general_warm_started[sel$tests$comparison == label]
  if (length(z)) z[[1L]] else NA
}
edge_boot_warm <- function(label) {
  z <- sel$tests$bootstrap_general_warm_started[
    sel$tests$comparison == label]
  if (length(z)) z[[1L]] else NA_integer_
}
row <- data.frame(
  model = model, J = J, N = N, B = B, selected = sel$selected,
  mon_p = edge_p("MON vs UN"), iio_p = edge_p("IIO vs UN"),
  dmiio_p = edge_p("DM vs IIO"), dmmon_p = edge_p("DM vs MON"),
  lcrdm_p = edge_p("LCR vs DM"),
  mon_raw = edge_raw("MON vs UN"), iio_raw = edge_raw("IIO vs UN"),
  dmiio_raw = edge_raw("DM vs IIO"), dmmon_raw = edge_raw("DM vs MON"),
  lcrdm_raw = edge_raw("LCR vs DM"),
  mon_pre = edge_pre("MON vs UN"), iio_pre = edge_pre("IIO vs UN"),
  dmiio_pre = edge_pre("DM vs IIO"), dmmon_pre = edge_pre("DM vs MON"),
  lcrdm_pre = edge_pre("LCR vs DM"),
  mon_warm = edge_warm("MON vs UN"), iio_warm = edge_warm("IIO vs UN"),
  dmiio_warm = edge_warm("DM vs IIO"),
  dmmon_warm = edge_warm("DM vs MON"),
  lcrdm_warm = edge_warm("LCR vs DM"),
  mon_boot_warm = edge_boot_warm("MON vs UN"),
  iio_boot_warm = edge_boot_warm("IIO vs UN"),
  dmiio_boot_warm = edge_boot_warm("DM vs IIO"),
  dmmon_boot_warm = edge_boot_warm("DM vs MON"),
  lcrdm_boot_warm = edge_boot_warm("LCR vs DM"),
  lcrdm_B_eff = if (is.null(sel$lcr_vs_dm)) NA else sel$lcr_vs_dm$B_effective,
  bridge_C = if (is.null(sel$quant_fits)) NA else sel$quant_fits$bridge_C,
  rmlcr_p = if (is.null(sel$rm_vs_lcr)) NA else sel$rm_vs_lcr$p_value,
  rmlcr_B_eff = if (is.null(sel$rm_vs_lcr)) NA else sel$rm_vs_lcr$B_effective,
  prof_C = if (is.null(sel$rm_vs_lcr)) NA else sel$rm_vs_lcr$profiled_C,
  secs = proc.time()[3] - t0)
write.csv(row, f, row.names = FALSE)
print(row, row.names = FALSE)
