#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

double nw_corner_distance_cpp(
    NumericVector x, NumericVector y,
    bool presorted = false,
    double p = 2.0,
    double eps = 1.490116e-08
) {
  int n = x.size();
  int m = y.size();

  // Sort if needed
  NumericVector xs = presorted ? x : clone(x);
  NumericVector ys = presorted ? y : clone(y);
  if (!presorted) {
    std::sort(xs.begin(), xs.end());
    std::sort(ys.begin(), ys.end());
  }

  // If equal-sized, no need for run-length encoding
  if (n == m) {
    double total = 0.0;
    for (int i = 0; i < n; i++) {
      total += std::pow(std::abs(xs[i] - ys[i]), p);
    }
    return std::pow(total / n, 1.0 / p);
  }

  // Otherwise, use run-length encoding
  std::vector<double> xv, yv, alpha, beta;
  int k = 1;
  for (int i = 1; i < n; i++) {
    if (std::abs(xs[i] - xs[i-1]) > eps) {
      xv.push_back(xs[i-1]);
      alpha.push_back((double)k / n);
      k = 1;
    } else k++;
  }
  xv.push_back(xs[n-1]); alpha.push_back((double)k / n);

  k = 1;
  for (int i = 1; i < m; i++) {
    if (std::abs(ys[i] - ys[i-1]) > eps) {
      yv.push_back(ys[i-1]);
      beta.push_back((double)k / m);
      k = 1;
    } else k++;
  }
  yv.push_back(ys[m-1]); beta.push_back((double)k / m);

  // North-west corner algorithm
  int nx = xv.size(), ny = yv.size();
  double total = 0.0;
  int i = 0, j = 0;
  double a = alpha[0], b = beta[0];
  while (i < nx && j < ny) {
    double t = (a < b) ? a : b;
    total += t * std::pow(std::abs(xv[i] - yv[j]), p);
    a -= t; b -= t;
    if (a <= eps) { i++; if (i < nx) a = alpha[i]; }
    else          { j++; if (j < ny) b = beta[j];  }
  }
  return std::pow(total, 1.0 / p);
}
