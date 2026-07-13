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
    arma::mat P = as<arma::mat>(item_probs[j]);          // C x (m+1)
    arma::mat lP = arma::log(arma::clamp(P, 1e-12, 1.0));
    for (arma::uword i = 0; i < n; ++i) {
      const int x = data(i, j);
      if (x < 0) continue;   // NA (INT_MIN) / negative codes: missing, skip
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
                              const arma::mat& posteriors,
                              const IntegerVector& cat_counts) {
  const arma::uword n = data.n_rows, J = data.n_cols, C = posteriors.n_cols;
  List out(J);
  for (arma::uword j = 0; j < J; ++j) {
    const int m = cat_counts[j];
    arma::mat ec(C, m + 1, arma::fill::zeros);
    for (arma::uword i = 0; i < n; ++i) {
      const int x = data(i, j);
      if (x < 0) continue;   // missing response contributes no expected count
      ec.col(x) += posteriors.row(i).t();
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
