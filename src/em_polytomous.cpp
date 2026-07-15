// Polytomous latent-class EM engine (C++ / RcppArmadillo).
//
// Companion to em_lca.cpp for the dichotomous engine. Each item-by-class
// combination carries a full category-probability vector; the E-step is a
// multinomial mixture in the log domain with log-sum-exp stabilisation. The
// same E-step drives the latent-class models and the Gauss-Hermite quadrature
// used for the partial-credit / Rasch marginal likelihood.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Row-wise log-sum-exp of an (n x C) matrix.
static inline arma::vec row_log_sum_exp(const arma::mat& M) {
  arma::vec mx = arma::max(M, 1);
  arma::mat e = arma::exp(M.each_col() - mx);
  return mx + arma::log(arma::sum(e, 1));
}

//' Polytomous multinomial E-step (internal)
//'
//' @param data Integer response matrix (n x J), entries 0..m_j.
//' @param item_probs List of length J; element j is a C x (m_j + 1) matrix of
//'   category probabilities (rows index classes, sum to 1).
//' @param class_probs Class probability vector (length C).
//' @return List with \code{posteriors} (n x C) and \code{loglik}.
//' @noRd
// [[Rcpp::export]]
List cpp_poly_estep(const arma::imat& data, const List& item_probs,
                    const arma::vec& class_probs) {
  const arma::uword n = data.n_rows, J = data.n_cols, C = class_probs.n_elem;
  arma::mat ll(n, C);
  arma::rowvec lcp = arma::log(class_probs.t());
  ll.each_row() = lcp;

  for (arma::uword j = 0; j < J; ++j) {
    // Access the item's probability matrix WITHOUT touching its attributes.
    // The per-iteration matrices from cpp_pcm_probs can lose their attribute
    // chain to a GC-lifetime anomaly (dim pairlist collected while the data
    // payload survives; see the segfault investigation in the TID validation).
    // as<arma::mat> walks the dim attribute and segfaults on such a node, so
    // we validate the header (safe reads) and index the column-major payload
    // manually: element (c, x) of the C x (m+1) matrix lives at [x*C + c].
    SEXP pj = item_probs[j];
    if (TYPEOF(pj) != REALSXP)
      stop("cpp_poly_estep: item_probs[%d] is not numeric (TYPEOF=%d)",
           (int)j + 1, TYPEOF(pj));
    const R_xlen_t len = Rf_xlength(pj);
    if (len < (R_xlen_t)C || (len % (R_xlen_t)C) != 0)
      stop("cpp_poly_estep: item_probs[%d] length %d not a multiple of n_classes %d",
           (int)j + 1, (int)len, (int)C);
    const arma::uword ncat = (arma::uword)(len / (R_xlen_t)C);
    const double* P = REAL(pj);
    // per-column logs, computed once
    arma::mat lP(C, ncat);
    for (arma::uword x = 0; x < ncat; ++x)
      for (arma::uword c = 0; c < C; ++c) {
        double p = P[x * C + c];
        if (!(p > 1e-12)) p = 1e-12;   // clamp; also maps NaN to floor
        if (p > 1.0) p = 1.0;
        lP(c, x) = std::log(p);
      }
    for (arma::uword i = 0; i < n; ++i) {
      const int x = data(i, j);
      // NA (INT_MIN) / negative codes, or a category beyond what this item's
      // probability matrix covers (can occur when a refit is scored on data
      // whose category range exceeds the fitted model's): treat as missing.
      if (x < 0 || static_cast<arma::uword>(x) >= ncat) continue;
      for (arma::uword c = 0; c < C; ++c) ll(i, c) += lP(c, x);
    }
  }
  arma::vec lrs = row_log_sum_exp(ll);
  arma::mat post = arma::exp(ll.each_col() - lrs);
  return List::create(_["posteriors"] = post,
                      _["loglik"] = arma::accu(lrs));
}

//' Expected category counts per item (internal)
//'
//' @param data Integer response matrix (n x J).
//' @param posteriors Posterior class-membership matrix (n x C).
//' @param cat_counts Integer vector of m_j (categories run 0..m_j).
//' @return List of length J; element j is a C x (m_j + 1) matrix of expected
//'   counts.
//' @noRd
// [[Rcpp::export]]
List cpp_poly_expected_counts(const arma::imat& data,
                              const NumericVector& posteriors,
                              const IntegerVector& cat_counts) {
  // posteriors is the n x C posterior matrix, taken as a plain NumericVector
  // so the conversion never reads its dim attribute (same GC-lifetime
  // anomaly guard as cpp_poly_estep; the matrix is recreated every EM
  // iteration and its attribute chain can be lost). Column-major: (i, c) is
  // posteriors[c * n + i].
  const arma::uword n = data.n_rows, J = data.n_cols;
  if (n == 0 || (posteriors.size() % (R_xlen_t)n) != 0)
    stop("cpp_poly_expected_counts: posteriors length %d not a multiple of n %d",
         (int)posteriors.size(), (int)n);
  const arma::uword C = (arma::uword)(posteriors.size() / (R_xlen_t)n);
  const double* post = posteriors.begin();
  List out(J);
  for (arma::uword j = 0; j < J; ++j) {
    const int m = cat_counts[j];
    arma::mat ec(C, m + 1, arma::fill::zeros);
    for (arma::uword i = 0; i < n; ++i) {
      const int x = data(i, j);
      // missing, or a category beyond the m_j this item was set up for: skip
      // (mirrors the bounds guard in cpp_poly_estep to avoid out-of-range col).
      if (x < 0 || x > m) continue;
      for (arma::uword c = 0; c < C; ++c) ec(c, x) += post[c * n + i];
    }
    out[j] = ec;
  }
  return out;
}

//' Partial-credit category probabilities (internal)
//'
//' For class/quadrature locations \code{theta} (length C) and item step
//' parameters \code{delta} (length m), returns the C x (m+1) matrix with
//' P(X = x | c) proportional to exp(sum_{k<=x} (theta_c - delta_k)).
//'
//' @param theta Numeric vector of person/class locations (length C).
//' @param delta Numeric vector of item step parameters (length m).
//' @return C x (m+1) matrix of category probabilities.
//' @noRd
// [[Rcpp::export]]
arma::mat cpp_pcm_probs(const arma::vec& theta, const arma::vec& delta) {
  const arma::uword C = theta.n_elem, m = delta.n_elem;
  arma::mat num(C, m + 1, arma::fill::zeros);          // column 0 = 0
  for (arma::uword c = 0; c < C; ++c) {
    double cum = 0.0;
    for (arma::uword k = 0; k < m; ++k) {
      cum += theta[c] - delta[k];
      num(c, k + 1) = cum;
    }
  }
  num.each_col() -= arma::max(num, 1);
  arma::mat e = arma::exp(num);
  return e.each_col() / arma::sum(e, 1);
}
