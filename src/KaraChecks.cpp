#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector lsqisotonic(NumericVector x, NumericVector y, NumericVector w) {
  int n = x.size();

  std::vector<int> ord(n);
  for (int i = 0; i < n; i++) ord[i] = i;

  std::sort(ord.begin(), ord.end(), [&](int i, int j) {
    if (x[i] != x[j]) return x[i] < x[j];
    return y[i] < y[j];
  });

  std::vector<double> yhat(n), ww(n);
  for (int i = 0; i < n; i++) {
    yhat[i] = y[ord[i]];
    ww[i] = w[ord[i]];
  }

  std::vector<double> pooled_y, pooled_w;
  std::vector<int> pooled_count;

  for (int i = 0; i < n; i++) {
    pooled_y.push_back(yhat[i]);
    pooled_w.push_back(ww[i]);
    pooled_count.push_back(1);

    while (pooled_y.size() > 1) {
      int m = pooled_y.size();
      if (pooled_y[m-1] >= pooled_y[m-2]) break;

      double sum_wy = pooled_w[m-1] * pooled_y[m-1] + pooled_w[m-2] * pooled_y[m-2];
      double sum_w = pooled_w[m-1] + pooled_w[m-2];
      int sum_count = pooled_count[m-1] + pooled_count[m-2];

      pooled_y.pop_back();
      pooled_w.pop_back();
      pooled_count.pop_back();

      pooled_y.back() = sum_wy / sum_w;
      pooled_w.back() = sum_w;
      pooled_count.back() = sum_count;
    }
  }

  std::vector<double> result_sorted(n);
  int idx = 0;
  for (size_t b = 0; b < pooled_y.size(); b++) {
    for (int c = 0; c < pooled_count[b]; c++) {
      result_sorted[idx++] = pooled_y[b];
    }
  }

  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[ord[i]] = result_sorted[i];
  }

  return result;
}

// Kernel density at a point using Scott's rule with MAD-based sigma (matches MATLAB's default).
// MATLAB uses robust sigma: sigma = MAD / 0.6745, where MAD = median(abs(x - median(x))).
// Bandwidth: h = sigma * (4 / (3*n))^0.2
// [[Rcpp::export]]
double ksdensity(NumericVector data, double point) {
  int n = data.size();

  std::vector<double> sorted_data(data.begin(), data.end());
  std::sort(sorted_data.begin(), sorted_data.end());
  double median_val = (n % 2 == 0)
    ? (sorted_data[n/2 - 1] + sorted_data[n/2]) / 2.0
    : sorted_data[n/2];

  std::vector<double> abs_dev(n);
  for (int i = 0; i < n; i++) {
    abs_dev[i] = std::abs(data[i] - median_val);
  }
  std::sort(abs_dev.begin(), abs_dev.end());
  double mad = (n % 2 == 0)
    ? (abs_dev[n/2 - 1] + abs_dev[n/2]) / 2.0
    : abs_dev[n/2];

  double sigma = mad / 0.6745;
  double h;

  if (sigma == 0 || std::isnan(sigma)) {
    double min_val = *std::min_element(data.begin(), data.end());
    double max_val = *std::max_element(data.begin(), data.end());
    double range = max_val - min_val;
    sigma = (range > 0) ? range : 0.0;
  }

  if (sigma == 0) {
    h = 1.0;
  } else {
    h = sigma * std::pow(4.0 / (3.0 * n), 0.2);
  }

  double sum = 0;
  double sqrt2pi = std::sqrt(2.0 * M_PI);
  for (int i = 0; i < n; i++) {
    double z = (point - data[i]) / h;
    sum += std::exp(-0.5 * z * z) / (h * sqrt2pi);
  }

  return sum / n;
}
