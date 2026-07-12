// Constrained polytomous M-step solved entirely in C++.
//
// The M-step maximises the weighted multinomial log-likelihood
//   sum_i w_i log p_i
// subject to linear constraints: per-(item,class) simplices (Aeq p = 1) and the
// monotonicity / item-ordering inequalities (B p <= 0). These are exactly the
// constraints of the R reference solver; here we call NLopt's SLSQP (the same
// algorithm nloptr wraps) directly through nloptrAPI.h, with C++ objective,
// gradient and constraint callbacks, so there is no per-iteration R round-trip.

#include <RcppArmadillo.h>
#include <nloptrAPI.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// objective: minimise -sum w log x, gradient -w/x
struct ObjData { const double* w; };
static double poly_obj(unsigned n, const double* x, double* grad, void* data) {
  const double* w = static_cast<ObjData*>(data)->w;
  double f = 0.0;
  for (unsigned i = 0; i < n; ++i) {
    double xi = x[i] < 1e-12 ? 1e-12 : x[i];
    f += w[i] * std::log(xi);
    if (grad) grad[i] = -w[i] / xi;
  }
  return -f;
}

// vector-valued linear constraint c(x) = M x - b (row-major jacobian = M)
struct LinData { arma::mat M; arma::vec b; };
static void lin_con(unsigned m, double* result, unsigned n, const double* x,
                    double* grad, void* data) {
  LinData* d = static_cast<LinData*>(data);
  arma::vec xv(const_cast<double*>(x), n, false, true);
  arma::vec r = d->M * xv - d->b;
  for (unsigned i = 0; i < m; ++i) result[i] = r[i];
  if (grad) {
    for (unsigned i = 0; i < m; ++i)
      for (unsigned j = 0; j < n; ++j) grad[i * n + j] = d->M(i, j);
  }
}

//' Constrained multinomial M-step via NLopt SLSQP (internal)
//'
//' @param w Packed expected-count weights (length = total free cells).
//' @param p0 Packed warm-start probabilities (same length).
//' @param Aeq Equality-constraint matrix (rows sum to 1).
//' @param B Inequality-constraint matrix (B p <= 0); may have zero rows.
//' @param xtol_rel Relative x tolerance for SLSQP.
//' @param maxeval Maximum objective evaluations.
//' @return The solution vector (unnormalised; the caller renormalises rows).
//' @noRd
// [[Rcpp::export]]
NumericVector cpp_poly_mstep_solve(const arma::vec& w, const arma::vec& p0,
                                   const arma::mat& Aeq, const arma::mat& B,
                                   double xtol_rel = 1e-8, int maxeval = 500) {
  unsigned n = p0.n_elem;
  ObjData od{ w.memptr() };
  LinData eqd{ Aeq, arma::ones<arma::vec>(Aeq.n_rows) };
  LinData ind{ B, arma::zeros<arma::vec>(B.n_rows) };

  nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, n);
  std::vector<double> lb(n, 0.0), ub(n, 1.0);
  nlopt_set_lower_bounds(opt, lb.data());
  nlopt_set_upper_bounds(opt, ub.data());
  nlopt_set_min_objective(opt, poly_obj, &od);

  std::vector<double> tol_eq(Aeq.n_rows, 1e-8);
  nlopt_add_equality_mconstraint(opt, Aeq.n_rows, lin_con, &eqd, tol_eq.data());
  std::vector<double> tol_in;
  if (B.n_rows > 0) {
    tol_in.assign(B.n_rows, 1e-8);
    nlopt_add_inequality_mconstraint(opt, B.n_rows, lin_con, &ind, tol_in.data());
  }
  nlopt_set_xtol_rel(opt, xtol_rel);
  nlopt_set_maxeval(opt, maxeval);

  std::vector<double> x(p0.begin(), p0.end());
  double minf;
  nlopt_optimize(opt, x.data(), &minf);
  nlopt_destroy(opt);
  return NumericVector(x.begin(), x.end());
}
