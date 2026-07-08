// EM computational core for QuantFit latent class models (C++ port).
//
// This file is a line-faithful RcppArmadillo port of the reference R
// implementations in R/em_algorithm.R and R/constraints.R (which remain the
// authoritative, validated reference and are still used when use_cpp = FALSE).
//
// Numerical-equivalence notes:
//  * All reductions that R performs with sum()/colSums() are accumulated in
//    `long double` here, matching R's LDOUBLE accumulation.
//  * Matrix-vector products use the same BLAS that R's %*% uses (Armadillo is
//    configured by RcppArmadillo to call R's BLAS), with the same operation
//    order as the R code.
//  * Probability bounding uses the same eps = 1e-10 as bound_probs() and is
//    applied at exactly the same points as in the R code.
//  * The convergence check replicates check_convergence() (min_iter = 3,
//    relative change, and the "GEM decrease" rule that a decrease larger than
//    tol never counts as convergence). Warnings for GEM decreases and
//    degenerate classes are emitted by the R wrappers from the returned
//    ll_history / degenerate flag.
//  * No RNG is used: the EM is deterministic given the initial values, which
//    are always produced in R.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

using namespace Rcpp;

static const double PROB_EPS = 1e-10;

// R's sum()/colSums() accumulate in long double (LDOUBLE). On platforms
// where long double has no extra precision over double (e.g. Apple arm64 and
// other aarch64 ABIs with 64-bit long double), an in-order double
// accumulation -- such as a reference-BLAS dgemm -- is bit-identical to R's
// colSums, so the M-step numerators can be computed with a single dgemm.
// Where long double is wider (x86 80-bit, aarch64-linux quad), we keep the
// explicit long double loops to match R exactly.
static const bool kLongDoubleIsDouble = (LDBL_MANT_DIG == DBL_MANT_DIG);

// bound_probs() for a scalar: pmax(pmin(p, 1 - eps), eps)
static inline double bound_prob(double p) {
  return std::max(std::min(p, 1.0 - PROB_EPS), PROB_EPS);
}

// bound_probs() for a matrix
static arma::mat bound_probs_mat(arma::mat p) {
  for (arma::uword k = 0; k < p.n_elem; ++k) p(k) = bound_prob(p(k));
  return p;
}

// ---------------------------------------------------------------------------
// E-step (mirrors e_step() in R/em_algorithm.R)
// ---------------------------------------------------------------------------

// Fills `posteriors` and returns the observed-data log-likelihood.
// `one_minus_data` is (1 - data), precomputed once per EM run.
static double e_step_core(const arma::mat& data,
                          const arma::mat& one_minus_data,
                          const arma::mat& item_probs,
                          const arma::vec& class_probs,
                          arma::mat& posteriors) {
  const arma::uword n_obs = data.n_rows;
  const arma::uword n_items = data.n_cols;
  const arma::uword n_classes = class_probs.n_elem;

  arma::mat log_lik_mat(n_obs, n_classes);
  arma::vec log_p(n_items), log_1mp(n_items);

  for (arma::uword c = 0; c < n_classes; ++c) {
    for (arma::uword j = 0; j < n_items; ++j) {
      const double p = item_probs(j, c);
      log_p(j) = std::log(bound_prob(p));
      log_1mp(j) = std::log(bound_prob(1.0 - p));
    }
    // Same arithmetic order as R:
    // data %*% log_p + (1 - data) %*% log_1mp + log(class_probs[c])
    log_lik_mat.col(c) = data * log_p + one_minus_data * log_1mp
      + std::log(class_probs(c));
  }

  // Row-wise log-sum-exp with max shift (log_sum_exp() in R/utils.R);
  // the inner sum uses long double accumulation to match R's sum()
  arma::vec log_row_sums(n_obs);
  for (arma::uword i = 0; i < n_obs; ++i) {
    double max_x = log_lik_mat(i, 0);
    for (arma::uword c = 1; c < n_classes; ++c) {
      if (log_lik_mat(i, c) > max_x) max_x = log_lik_mat(i, c);
    }
    if (std::isinf(max_x)) {
      log_row_sums(i) = max_x;
      continue;
    }
    long double s = 0.0L;
    for (arma::uword c = 0; c < n_classes; ++c) {
      s += std::exp(log_lik_mat(i, c) - max_x);
    }
    log_row_sums(i) = max_x + std::log(static_cast<double>(s));
  }

  posteriors.set_size(n_obs, n_classes);
  for (arma::uword c = 0; c < n_classes; ++c) {
    for (arma::uword i = 0; i < n_obs; ++i) {
      posteriors(i, c) = std::exp(log_lik_mat(i, c) - log_row_sums(i));
    }
  }

  long double total = 0.0L;
  for (arma::uword i = 0; i < n_obs; ++i) total += log_row_sums(i);
  return static_cast<double>(total);
}

// ---------------------------------------------------------------------------
// Unconstrained M-step (mirrors m_step() in R/em_algorithm.R)
// ---------------------------------------------------------------------------

// class_counts receives colSums(posteriors) (used as PAVA weights by the
// constrained M-step); degenerate is any(class_counts < 1).
static void m_step_core(const arma::mat& data,
                        const arma::mat& posteriors,
                        arma::mat& item_probs,
                        arma::vec& class_probs,
                        arma::vec& class_counts,
                        bool& degenerate) {
  const arma::uword n_obs = data.n_rows;
  const arma::uword n_items = data.n_cols;
  const arma::uword n_classes = posteriors.n_cols;

  class_counts.set_size(n_classes);
  for (arma::uword c = 0; c < n_classes; ++c) {
    long double s = 0.0L;
    for (arma::uword i = 0; i < n_obs; ++i) s += posteriors(i, c);
    class_counts(c) = static_cast<double>(s);
  }

  class_probs = class_counts / static_cast<double>(n_obs);
  degenerate = arma::any(class_counts < 1.0);

  // Numerators of the weighted item means: colSums(data * weights) per class,
  // i.e. the (I x C) matrix data' * posteriors with in-order accumulation
  item_probs.set_size(n_items, n_classes);
  if (kLongDoubleIsDouble) {
    // Single dgemm: same products, same in-order accumulation as R's
    // colSums when long double == double (see note at top of file)
    const arma::mat num = data.t() * posteriors;
    for (arma::uword c = 0; c < n_classes; ++c) {
      const double w_sum = class_counts(c);  // equals R's sum(weights)
      for (arma::uword j = 0; j < n_items; ++j) {
        const double p = (w_sum > 0.0) ? num(j, c) / w_sum : 0.5;
        item_probs(j, c) = bound_prob(p);
      }
    }
  } else {
    // Portable path: explicit long double accumulation matching R's LDOUBLE
    for (arma::uword c = 0; c < n_classes; ++c) {
      const double w_sum = class_counts(c);  // equals R's sum(weights)
      for (arma::uword j = 0; j < n_items; ++j) {
        double p;
        if (w_sum > 0.0) {
          long double num = 0.0L;
          for (arma::uword i = 0; i < n_obs; ++i) {
            num += data(i, j) * posteriors(i, c);
          }
          p = static_cast<double>(num) / w_sum;
        } else {
          p = 0.5;
        }
        item_probs(j, c) = bound_prob(p);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Weighted PAVA (mirrors pava_increasing() / pava_decreasing() in
// R/constraints.R): classical block-merging with (weighted-mean, weight,
// count) blocks, identical merge arithmetic
// ---------------------------------------------------------------------------

static arma::vec pava_increasing_core(const arma::vec& x, const arma::vec& w) {
  const int n = static_cast<int>(x.n_elem);
  if (n <= 1) return x;

  std::vector<double> vals(n), wts(n);
  std::vector<int> cnts(n);
  int nb = 0;

  for (int i = 0; i < n; ++i) {
    vals[nb] = x(i);
    wts[nb] = std::max(w(i), 0.0);  // R: w <- pmax(w, 0)
    cnts[nb] = 1;
    ++nb;

    // Merge while the last two blocks violate monotonicity
    while (nb > 1 && vals[nb - 2] > vals[nb - 1]) {
      const double w_sum = wts[nb - 2] + wts[nb - 1];
      if (w_sum > 0) {
        vals[nb - 2] = (vals[nb - 2] * wts[nb - 2] + vals[nb - 1] * wts[nb - 1])
          / w_sum;
      } else {
        // Both blocks have zero weight: use simple average
        vals[nb - 2] = (vals[nb - 2] + vals[nb - 1]) / 2.0;
      }
      wts[nb - 2] = w_sum;
      cnts[nb - 2] += cnts[nb - 1];
      --nb;
    }
  }

  arma::vec out(n);
  int pos = 0;
  for (int b = 0; b < nb; ++b) {
    for (int k = 0; k < cnts[b]; ++k) out(pos++) = vals[b];
  }
  return out;
}

static arma::vec pava_decreasing_core(const arma::vec& x, const arma::vec& w) {
  return arma::reverse(pava_increasing_core(arma::reverse(x), arma::reverse(w)));
}

// ---------------------------------------------------------------------------
// Dykstra alternating projections for double monotonicity (mirrors
// dykstra_dm_projection() in R/constraints.R)
// ---------------------------------------------------------------------------

// order0 is the item order converted to 0-based indices; class_weights must
// already be floored at 1e-10 by the caller.
static arma::mat dykstra_core(const arma::mat& item_probs,
                              const arma::uvec& order0,
                              const arma::vec& class_weights,
                              double tol, int max_cycles) {
  const arma::uword n_items = item_probs.n_rows;
  const arma::uword n_classes = item_probs.n_cols;

  arma::mat x = item_probs;
  arma::mat incr_row(n_items, n_classes, arma::fill::zeros);
  arma::mat incr_col(n_items, n_classes, arma::fill::zeros);
  arma::vec unit_w(n_items, arma::fill::ones);

  for (int cycle = 0; cycle < max_cycles; ++cycle) {
    const arma::mat x_old = x;

    // Project onto row-monotone set (weighted PAVA per row)
    arma::mat z = x + incr_row;
    arma::mat y = z;
    for (arma::uword i = 0; i < n_items; ++i) {
      y.row(i) = pava_increasing_core(z.row(i).t(), class_weights).t();
    }
    incr_row = z - y;

    // Project onto column-ordered set (plain PAVA per column in item order)
    z = y + incr_col;
    x = z;
    for (arma::uword c = 0; c < n_classes; ++c) {
      arma::vec sub(n_items);
      for (arma::uword k = 0; k < n_items; ++k) sub(k) = z(order0(k), c);
      sub = pava_decreasing_core(sub, unit_w);
      for (arma::uword k = 0; k < n_items; ++k) x(order0(k), c) = sub(k);
    }
    incr_col = z - x;

    if (arma::abs(x - x_old).max() < tol) break;
  }

  return x;
}

// ---------------------------------------------------------------------------
// Weighted projection onto the constraint space (mirrors
// project_constraints_weighted() in R/constraints.R)
// ---------------------------------------------------------------------------

static arma::mat project_constraints_core(const arma::mat& item_probs,
                                          bool class_monotonicity,
                                          bool item_ordering,
                                          const arma::uvec& order0,
                                          arma::vec class_weights) {
  const arma::uword n_items = item_probs.n_rows;
  const arma::uword n_classes = item_probs.n_cols;

  // Guard against fully collapsed classes (R: pmax(class_weights, 1e-10))
  for (arma::uword c = 0; c < n_classes; ++c) {
    class_weights(c) = std::max(class_weights(c), 1e-10);
  }

  if (class_monotonicity && item_ordering) {
    // Double monotonicity: exact 2D weighted isotonic regression via Dykstra
    return bound_probs_mat(
      dykstra_core(item_probs, order0, class_weights, 1e-10, 500));
  }

  arma::mat proj = item_probs;

  if (class_monotonicity) {
    for (arma::uword i = 0; i < n_items; ++i) {
      proj.row(i) = pava_increasing_core(proj.row(i).t(), class_weights).t();
    }
  }

  if (item_ordering) {
    arma::vec unit_w(n_items, arma::fill::ones);
    for (arma::uword c = 0; c < n_classes; ++c) {
      arma::vec sub(n_items);
      for (arma::uword k = 0; k < n_items; ++k) sub(k) = proj(order0(k), c);
      sub = pava_decreasing_core(sub, unit_w);
      for (arma::uword k = 0; k < n_items; ++k) proj(order0(k), c) = sub(k);
    }
  }

  return bound_probs_mat(proj);
}

// ---------------------------------------------------------------------------
// Convergence check (mirrors check_convergence() in R/utils.R with the
// default min_iter = 3; decreases > tol never count as convergence -- the R
// wrapper re-scans ll_history to emit the corresponding warnings)
// ---------------------------------------------------------------------------

static bool check_convergence_core(const std::vector<double>& ll_history,
                                   double tol) {
  const int n = static_cast<int>(ll_history.size());
  if (n < 3) return false;

  const double change = ll_history[n - 1] - ll_history[n - 2];
  const double rel_change = change / std::fabs(ll_history[n - 2]);

  if (rel_change < -tol) return false;  // GEM decrease: not convergence
  return std::fabs(rel_change) < tol;
}

// Final degenerate-class check: any(colSums(posteriors) < 1)
static bool degenerate_core(const arma::mat& posteriors) {
  for (arma::uword c = 0; c < posteriors.n_cols; ++c) {
    long double s = 0.0L;
    for (arma::uword i = 0; i < posteriors.n_rows; ++i) s += posteriors(i, c);
    if (static_cast<double>(s) < 1.0) return true;
  }
  return false;
}

// Convert a 1-based item order from R into 0-based indices
static arma::uvec order_to_zero_based(const IntegerVector& item_order) {
  arma::uvec order0(item_order.size());
  for (int k = 0; k < item_order.size(); ++k) {
    order0(k) = static_cast<arma::uword>(item_order[k] - 1);
  }
  return order0;
}

// Build the common EM result list
static List em_result_list(const arma::mat& item_probs,
                           const arma::vec& class_probs,
                           const arma::mat& posteriors,
                           const std::vector<double>& ll_history,
                           bool converged, int iterations, bool degenerate) {
  return List::create(
    _["item_probs"] = item_probs,
    _["class_probs"] = NumericVector(class_probs.begin(), class_probs.end()),
    _["posteriors"] = posteriors,
    _["loglik"] = ll_history.back(),
    _["ll_history"] = NumericVector(ll_history.begin(), ll_history.end()),
    _["converged"] = converged,
    _["iterations"] = iterations,
    _["degenerate"] = degenerate);
}

// ---------------------------------------------------------------------------
// Exported functions
// ---------------------------------------------------------------------------

//' C++ E-step for the latent class EM (internal)
//'
//' Log-domain Bernoulli mixture posteriors with log-sum-exp stabilization.
//' Port of \code{e_step()}.
//'
//' @param data Binary data matrix (n x I)
//' @param item_probs Item probability matrix (I x C)
//' @param class_probs Class probability vector (length C)
//' @return List with \code{posteriors} (n x C) and \code{loglik}
//' @noRd
// [[Rcpp::export]]
List cpp_e_step(const arma::mat& data, const arma::mat& item_probs,
                const arma::vec& class_probs) {
  arma::mat one_minus_data = 1.0 - data;
  arma::mat posteriors;
  const double loglik =
    e_step_core(data, one_minus_data, item_probs, class_probs, posteriors);
  return List::create(_["posteriors"] = posteriors, _["loglik"] = loglik);
}

//' C++ unconstrained M-step for the latent class EM (internal)
//'
//' Closed-form class proportions and weighted item means. Port of
//' \code{m_step()}.
//'
//' @param data Binary data matrix (n x I)
//' @param posteriors Posterior class membership matrix (n x C)
//' @return List with \code{item_probs}, \code{class_probs}, \code{degenerate}
//' @noRd
// [[Rcpp::export]]
List cpp_m_step(const arma::mat& data, const arma::mat& posteriors) {
  arma::mat item_probs;
  arma::vec class_probs, class_counts;
  bool degenerate = false;
  m_step_core(data, posteriors, item_probs, class_probs, class_counts,
              degenerate);
  return List::create(
    _["item_probs"] = item_probs,
    _["class_probs"] = NumericVector(class_probs.begin(), class_probs.end()),
    _["degenerate"] = degenerate);
}

//' C++ weighted PAVA (internal)
//'
//' Weighted isotonic regression by block merging. Port of
//' \code{pava_increasing()} / \code{pava_decreasing()}.
//'
//' @param x Numeric vector
//' @param w Optional non-negative weights (same length as x); NULL = unit
//' @param increasing If TRUE (default) non-decreasing fit, else non-increasing
//' @return Isotonic (weighted L2) projection of x
//' @noRd
// [[Rcpp::export]]
NumericVector cpp_weighted_pava(const arma::vec& x,
                                Nullable<NumericVector> w = R_NilValue,
                                bool increasing = true) {
  const arma::uword n = x.n_elem;
  if (n <= 1) return NumericVector(x.begin(), x.end());

  arma::vec wv(n, arma::fill::ones);
  if (w.isNotNull()) {
    NumericVector w_in(w);
    if (static_cast<arma::uword>(w_in.size()) != n) {
      stop("weights must have same length as x");
    }
    for (arma::uword i = 0; i < n; ++i) wv(i) = w_in[i];
  }

  const arma::vec out = increasing ? pava_increasing_core(x, wv)
                                   : pava_decreasing_core(x, wv);
  return NumericVector(out.begin(), out.end());
}

//' C++ Dykstra projection for double monotonicity (internal)
//'
//' Exact weighted L2 projection onto the intersection of the row-monotone
//' and column-ordered constraint sets. Port of
//' \code{dykstra_dm_projection()}.
//'
//' @param item_probs Item probability matrix (I x C)
//' @param item_order Item order (1-based, easiest to hardest)
//' @param class_weights Optional class weights (length C); NULL = unit
//' @param tol Convergence tolerance on the iterate change
//' @param max_cycles Maximum number of projection cycles
//' @return Projected matrix satisfying both constraint sets
//' @noRd
// [[Rcpp::export]]
arma::mat cpp_dykstra_dm(const arma::mat& item_probs,
                         const IntegerVector& item_order,
                         Nullable<NumericVector> class_weights = R_NilValue,
                         double tol = 1e-10, int max_cycles = 500) {
  const arma::uword n_classes = item_probs.n_cols;

  arma::vec cw(n_classes, arma::fill::ones);
  if (class_weights.isNotNull()) {
    NumericVector cw_in(class_weights);
    if (static_cast<arma::uword>(cw_in.size()) != n_classes) {
      stop("class_weights must have length ncol(item_probs)");
    }
    for (arma::uword c = 0; c < n_classes; ++c) cw(c) = cw_in[c];
  }
  for (arma::uword c = 0; c < n_classes; ++c) cw(c) = std::max(cw(c), 1e-10);

  return dykstra_core(item_probs, order_to_zero_based(item_order), cw, tol,
                      max_cycles);
}

//' C++ weighted projection onto the constraint space (internal)
//'
//' Port of \code{project_constraints_weighted()}: weighted PAVA per item for
//' class monotonicity, plain PAVA per class for item ordering, Dykstra for
//' double monotonicity, followed by probability bounding.
//'
//' @param item_probs Item probability matrix (I x C)
//' @param class_monotonicity Enforce class monotonicity?
//' @param item_ordering Enforce invariant item ordering?
//' @param item_order Item order (1-based); required if item_ordering
//' @param class_weights Optional class weights (length C); NULL = unit
//' @return Projected item probability matrix
//' @noRd
// [[Rcpp::export]]
arma::mat cpp_project_constraints(const arma::mat& item_probs,
                                  bool class_monotonicity,
                                  bool item_ordering,
                                  const IntegerVector& item_order,
                                  Nullable<NumericVector> class_weights = R_NilValue) {
  const arma::uword n_classes = item_probs.n_cols;

  arma::vec cw(n_classes, arma::fill::ones);
  if (class_weights.isNotNull()) {
    NumericVector cw_in(class_weights);
    if (static_cast<arma::uword>(cw_in.size()) != n_classes) {
      stop("class_weights must have length ncol(item_probs)");
    }
    for (arma::uword c = 0; c < n_classes; ++c) cw(c) = cw_in[c];
  }

  arma::uvec order0;
  if (item_ordering) {
    if (item_order.size() != static_cast<int>(item_probs.n_rows)) {
      stop("item_order must be specified for item ordering constraints");
    }
    order0 = order_to_zero_based(item_order);
  }

  return project_constraints_core(item_probs, class_monotonicity,
                                  item_ordering, order0, cw);
}

//' C++ EM driver for unconstrained LCA (internal)
//'
//' Full EM loop of \code{em_lca()} given explicit initial values. Convergence
//' semantics match the R implementation exactly (including the final E-step
//' when max_iter is reached without convergence). GEM-decrease and
//' degenerate-class warnings are emitted by the R wrapper.
//'
//' @param data Binary data matrix (n x I)
//' @param init_probs Initial item probabilities (I x C)
//' @param init_class_probs Initial class probabilities (length C)
//' @param max_iter Maximum number of iterations
//' @param tol Convergence tolerance
//' @return List mirroring \code{em_lca()} output
//' @noRd
// [[Rcpp::export]]
List cpp_em_lca(const arma::mat& data, const arma::mat& init_probs,
                const arma::vec& init_class_probs, int max_iter, double tol) {
  arma::mat item_probs = init_probs;
  arma::vec class_probs = init_class_probs;
  const arma::mat one_minus_data = 1.0 - data;

  std::vector<double> ll_history;
  ll_history.reserve(max_iter + 1);
  bool converged = false;
  arma::mat posteriors;
  int iter = 1;

  for (iter = 1; iter <= max_iter; ++iter) {
    // E-step
    ll_history.push_back(
      e_step_core(data, one_minus_data, item_probs, class_probs, posteriors));

    // Check convergence
    if (iter > 1 && check_convergence_core(ll_history, tol)) {
      converged = true;
      break;
    }

    // M-step
    arma::vec class_counts;
    bool degen = false;
    m_step_core(data, posteriors, item_probs, class_probs, class_counts,
                degen);
  }
  if (iter > max_iter) iter = max_iter;  // loop ran to completion

  // If we exited on max_iter, run a final E-step so loglik/posteriors match
  // the returned parameters
  if (!converged) {
    ll_history.push_back(
      e_step_core(data, one_minus_data, item_probs, class_probs, posteriors));
  }

  return em_result_list(item_probs, class_probs, posteriors, ll_history,
                        converged, iter, degenerate_core(posteriors));
}

//' C++ EM driver for constrained LCA with the exact PAVA M-step (internal)
//'
//' Full EM loop of \code{em_constrained(method = "pava")} given explicit
//' initial values: E-step, unconstrained M-step, then exact weighted
//' projection onto the constraint space (weighted PAVA / Dykstra) with the
//' expected class counts as weights, exactly as \code{m_step_exact()}.
//'
//' @param data Binary data matrix (n x I)
//' @param init_probs Initial item probabilities (I x C)
//' @param init_class_probs Initial class probabilities (length C)
//' @param class_monotonicity Enforce class monotonicity?
//' @param item_ordering Enforce invariant item ordering?
//' @param item_order Item order (1-based); required if item_ordering
//' @param max_iter Maximum number of iterations
//' @param tol Convergence tolerance
//' @return List mirroring \code{em_constrained()} output
//' @noRd
// [[Rcpp::export]]
List cpp_em_constrained(const arma::mat& data, const arma::mat& init_probs,
                        const arma::vec& init_class_probs,
                        bool class_monotonicity, bool item_ordering,
                        const IntegerVector& item_order, int max_iter,
                        double tol) {
  arma::uvec order0;
  if (item_ordering) {
    if (item_order.size() != static_cast<int>(data.n_cols)) {
      stop("item_order must be specified for item ordering constraints");
    }
    order0 = order_to_zero_based(item_order);
  }

  arma::mat item_probs = init_probs;
  arma::vec class_probs = init_class_probs;
  const arma::mat one_minus_data = 1.0 - data;

  std::vector<double> ll_history;
  ll_history.reserve(max_iter + 1);
  bool converged = false;
  arma::mat posteriors;
  int iter = 1;

  for (iter = 1; iter <= max_iter; ++iter) {
    // E-step (same as unconstrained)
    ll_history.push_back(
      e_step_core(data, one_minus_data, item_probs, class_probs, posteriors));

    // Check convergence
    if (iter > 1 && check_convergence_core(ll_history, tol)) {
      converged = true;
      break;
    }

    // Exact constrained M-step: unconstrained M-step, then weighted
    // projection with the expected class counts as weights (m_step_exact)
    arma::vec class_counts;
    bool degen = false;
    m_step_core(data, posteriors, item_probs, class_probs, class_counts,
                degen);
    item_probs = project_constraints_core(item_probs, class_monotonicity,
                                          item_ordering, order0, class_counts);
  }
  if (iter > max_iter) iter = max_iter;  // loop ran to completion

  // If we exited on max_iter, run a final E-step so loglik/posteriors match
  // the returned parameters
  if (!converged) {
    ll_history.push_back(
      e_step_core(data, one_minus_data, item_probs, class_probs, posteriors));
  }

  return em_result_list(item_probs, class_probs, posteriors, ll_history,
                        converged, iter, degenerate_core(posteriors));
}

//' C++ objective for the Latent Class Rasch M-step (internal)
//'
//' Negative expected complete-data log-likelihood as a function of the free
//' parameter vector \code{par = c(theta, delta\[-1\])} with the mean(delta) = 0
//' identification constraint (delta\[1\] = -sum(delta\[-1\])). Exact port of the
//' objective closure inside \code{m_step_rasch()}; the surrounding BFGS
//' optimization (\code{stats::optim}, finite-difference gradients) stays in R
//' so the optimizer trajectory is identical to the pure-R path.
//'
//' @param par Parameter vector c(theta, delta\[-1\]) of length C + I - 1
//' @param data Binary data matrix (n x I)
//' @param posteriors Posterior class membership matrix (n x C)
//' @param n_classes Number of latent classes
//' @return Negative expected complete-data log-likelihood (scalar)
//' @noRd
// [[Rcpp::export]]
double cpp_lcr_q(const arma::vec& par, const arma::mat& data,
                 const arma::mat& posteriors, int n_classes) {
  const arma::uword n_obs = data.n_rows;
  const arma::uword n_items = data.n_cols;
  const arma::uword C = static_cast<arma::uword>(n_classes);

  if (par.n_elem != C + n_items - 1) {
    stop("par must have length n_classes + n_items - 1");
  }

  // theta_new <- par[1:C]; delta_new <- c(-sum(delta_free), delta_free)
  const arma::vec theta = par.subvec(0, C - 1);
  arma::vec delta(n_items);
  long double dsum = 0.0L;
  for (arma::uword k = C; k < par.n_elem; ++k) dsum += par(k);
  delta(0) = -static_cast<double>(dsum);
  for (arma::uword j = 1; j < n_items; ++j) delta(j) = par(C + j - 1);

  // item_probs <- bound_probs(compute_rasch_probs(theta_new, delta_new))
  arma::mat item_probs(n_items, C);
  for (arma::uword c = 0; c < C; ++c) {
    for (arma::uword j = 0; j < n_items; ++j) {
      const double p = 1.0 / (1.0 + std::exp(-(theta(c) - delta(j))));
      item_probs(j, c) = bound_prob(p);
    }
  }

  const arma::mat one_minus_data = 1.0 - data;
  arma::vec log_p(n_items), log_1mp(n_items);

  double ll = 0.0;
  for (arma::uword c = 0; c < C; ++c) {
    for (arma::uword j = 0; j < n_items; ++j) {
      log_p(j) = std::log(item_probs(j, c));
      log_1mp(j) = std::log(1.0 - item_probs(j, c));
    }
    // sum(weights * (data %*% log_p + (1 - data) %*% log_1mp))
    const arma::vec contrib = data * log_p + one_minus_data * log_1mp;
    long double s = 0.0L;
    for (arma::uword i = 0; i < n_obs; ++i) s += posteriors(i, c) * contrib(i);
    ll += static_cast<double>(s);
  }

  return -ll;
}
