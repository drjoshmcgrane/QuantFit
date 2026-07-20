# Conditional-CML null generator (dichotomous, complete data).
#
# Estimates item difficulties by CONDITIONAL maximum likelihood (total score
# is sufficient, so the ability distribution drops out entirely) and generates
# null response patterns CONDITIONAL on each respondent's observed total
# score. For statistics defined over score groups (ConjointChecks,
# KaraChecks banding) this is the cleanest Rasch null: score-group sizes are
# preserved exactly in every replicate, nothing about the person distribution
# is estimated, and the item estimates inherit CML's distribution-free
# character (external methodological review; cf. Student & Read 2025, who use
# CML item estimation with plug-in person locations).

# Elementary symmetric functions of eps = exp(-beta), with running rescaling
# for numerical stability. Returns gamma[r+1] = esf of order r, and the
# rescale log-factor.
.esf <- function(eps) {
  J <- length(eps)
  g <- c(1, rep(0, J)); logC <- 0
  for (j in seq_len(J)) {
    g[2:(j + 1)] <- g[2:(j + 1)] + eps[j] * g[1:j]
    m <- max(g)
    if (m > 1e250) { g <- g / m; logC <- logC + log(m) }
  }
  list(gamma = g, logC = logC)
}

# CML estimation of item difficulties (Newton via BFGS on the conditional
# log-likelihood). Identification: mean(beta) = 0. Complete data only.
.cml_fit <- function(data, max_iter = 200L) {
  J <- ncol(data); r <- rowSums(data)
  keep <- r > 0 & r < J                       # extreme scores carry no info
  s_j <- colSums(data[keep, , drop = FALSE])  # item margins
  n_r <- tabulate(r[keep] + 1L, nbins = J + 1L)  # persons per score 0..J
  negll <- function(bfree) {
    if (any(!is.finite(bfree)) || max(abs(bfree)) > 30) return(1e12)
    b <- c(bfree, -sum(bfree))
    es <- .esf(exp(-pmin(pmax(b, -30), 30)))
    if (any(!is.finite(es$gamma))) return(1e12)
    lg <- log(es$gamma) + es$logC
    sum(s_j * b) + sum(n_r * lg)
  }
  o <- stats::optim(rep(0, J - 1L), negll, method = "BFGS",
                    control = list(maxit = max_iter))
  b <- c(o$par, -sum(o$par))
  names(b) <- colnames(data)
  b
}

# Suffix ESF table: G[[j]][r+1] = esf of order r over items j..J.
.esf_suffix <- function(eps) {
  J <- length(eps)
  G <- vector("list", J + 1L)
  G[[J + 1L]] <- c(1, rep(0, J))
  for (j in J:1) {
    g <- G[[j + 1L]]
    out <- g
    out[2:(J + 1L)] <- out[2:(J + 1L)] + eps[j] * g[1:J]
    G[[j]] <- out
  }
  G
}

# Draw one response pattern of total score r under Rasch(beta), sequentially:
# P(x_j = 1 | need s more from items j..J) = eps_j G[j+1][s] / G[j][s+1].
.r_conditional_pattern <- function(r, eps, G) {
  J <- length(eps)
  x <- integer(J); s <- r
  for (j in seq_len(J)) {
    if (s == 0L) break
    if (J - j + 1L == s) { x[j:J] <- 1L; break }
    p1 <- eps[j] * G[[j + 1L]][s] / G[[j]][s + 1L]
    if (stats::runif(1) < p1) { x[j] <- 1L; s <- s - 1L }
  }
  x
}

# Generate a full null dataset conditional on the observed score vector.
.conditional_null_dataset <- function(scores, beta) {
  eps <- exp(-beta)
  G <- .esf_suffix(eps)
  t(vapply(scores, .r_conditional_pattern, integer(length(beta)),
           eps = eps, G = G))
}

# ---- Generalized (polytomous + missing-data) conditional-CML machinery ----
# PCM conditional pattern probability given (answered set O, score r):
#   P(x | r, O) = prod_j eta_j(x_j) / gamma_r(O),  eta_j(c) = exp(-cumsum(d_j)[c])
# Dichotomous is the m_j = 1 special case. Persons are grouped by missingness
# pattern; each pattern gets its own suffix-ESF table.

.eta_list <- function(delta_list)
  lapply(delta_list, function(d) c(1, exp(-cumsum(pmin(pmax(d, -30), 30)))))

.gesf_suffix <- function(etas) {           # G[[j]][s+1] over items j..J
  J <- length(etas); maxr <- sum(vapply(etas, length, 1L)) - J
  G <- vector("list", J + 1L); G[[J + 1L]] <- c(1, rep(0, maxr))
  for (j in J:1) {
    g <- G[[j + 1L]]; out <- numeric(maxr + 1L); e <- etas[[j]]
    for (c in seq_along(e)) {
      k <- c - 1L
      out[(1 + k):(maxr + 1L)] <- out[(1 + k):(maxr + 1L)] +
        e[c] * g[1:(maxr + 1L - k)]
    }
    G[[j]] <- out
  }
  G
}

.cml_fit_general <- function(data, max_iter = 300L) {
  J <- ncol(data); m <- apply(data, 2, max, na.rm = TRUE)
  obs <- !is.na(data); r <- rowSums(data, na.rm = TRUE)
  pat <- apply(obs, 1, function(z) paste(as.integer(z), collapse = ""))
  upat <- unique(pat); pidx <- match(pat, upat)
  nstep <- sum(m); off <- c(0L, cumsum(m))
  # INFORMATIVE respondents only: pattern-extreme scorers (r = 0 or r = max
  # attainable on their answered set) have conditionally deterministic
  # patterns and must be excluded from margins AND denominators alike -
  # including them in margins only biases the difficulties under missingness
  # (external review: ~30-logit distortions from deterministic respondents).
  maxr_pat <- vapply(upat, function(u) {
    keep <- which(strsplit(u, "")[[1]] == "1"); sum(m[keep])
  }, 1)
  informative <- r > 0 & r < maxr_pat[pidx]
  di <- data[informative, , drop = FALSE]
  Sjk <- unlist(lapply(seq_len(J), function(j)
    vapply(seq_len(m[j]), function(k) sum(di[, j] >= k, na.rm = TRUE), 1)))
  negll <- function(dfree) {
    if (any(!is.finite(dfree)) || max(abs(dfree)) > 30) return(1e12)
    d <- c(dfree, -sum(dfree))
    dl <- lapply(seq_len(J), function(j) d[(off[j] + 1L):off[j + 1L]])
    tot <- sum(Sjk * d)
    for (u in seq_along(upat)) {
      keep <- which(strsplit(upat[u], "")[[1]] == "1")
      if (length(keep) < 2) next
      es <- .esf_gen_scaled(.eta_list(dl[keep]))
      rs <- r[pidx == u & informative]
      maxr <- sum(m[keep])
      if (!length(rs)) next
      lg <- log(es$g) + es$logC
      if (any(!is.finite(lg[rs + 1L]))) return(1e12)
      tot <- tot + sum(lg[rs + 1L])
    }
    tot
  }
  o <- stats::optim(rep(0, nstep - 1L), negll, method = "BFGS",
                    control = list(maxit = max_iter))
  if (o$convergence != 0L || !is.finite(o$value) || o$value >= 1e12 ||
      max(abs(o$par)) > 25)
    stop("CML estimation failed (convergence=", o$convergence,
         ", value=", signif(o$value, 4), ", max|par|=",
         signif(max(abs(o$par)), 3), "); the conditional null cannot be ",
         "constructed - check item co-observation connectivity")
  d <- c(o$par, -sum(o$par))
  lapply(seq_len(J), function(j) d[(off[j] + 1L):off[j + 1L]])
}

.esf_gen_scaled <- function(etas) {        # full-product ESF with rescaling
  maxr <- sum(vapply(etas, length, 1L)) - length(etas)
  g <- c(1, rep(0, maxr)); logC <- 0
  for (e in etas) {
    out <- numeric(maxr + 1L)
    for (c in seq_along(e)) {
      k <- c - 1L
      out[(1 + k):(maxr + 1L)] <- out[(1 + k):(maxr + 1L)] +
        e[c] * g[1:(maxr + 1L - k)]
    }
    g <- out; mx <- max(g)
    if (!is.finite(mx)) return(list(g = g, logC = logC))
    if (mx > 1e250) { g <- g / mx; logC <- logC + log(mx) }
  }
  list(g = g, logC = logC)
}

.conditional_null_general <- function(data, delta_list) {
  J <- ncol(data); obs <- !is.na(data)
  r <- rowSums(data, na.rm = TRUE)
  pat <- apply(obs, 1, function(z) paste(as.integer(z), collapse = ""))
  out <- matrix(NA_integer_, nrow(data), J)
  for (u in unique(pat)) {
    rows <- which(pat == u)
    keep <- which(strsplit(u, "")[[1]] == "1")
    if (!length(keep)) next
    etas <- .eta_list(delta_list[keep])
    G <- .gesf_suffix(etas)
    for (i in rows) {
      s <- r[i]; x <- integer(length(keep))
      for (jj in seq_along(keep)) {
        e <- etas[[jj]]
        if (s == 0L) break
        p <- vapply(seq_along(e), function(c) {
          k <- c - 1L
          if (k > s) 0 else e[c] * G[[jj + 1L]][s - k + 1L]
        }, 1)
        x[jj] <- sample.int(length(e), 1L, prob = p) - 1L
        s <- s - x[jj]
      }
      out[i, keep] <- x
    }
  }
  out
}
