// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

// utility function
inline void multModify(double *col, int D, double adj) {
  for (int i = 0; i != D; i++) col[i] *= adj;
}

bool negativeGradient(double* i,
                      double* k,
                      double* holder,
                      const double alpha,
                      const double gamma,
                      const double cap,
                      const int D) {
  const double dist_squared = distAndVector(i, k, holder, D);
  if (dist_squared == 0) return true; // If the two points are in the same place, skip
  const double adk = alpha * dist_squared;
  const double grad = gamma * ((alpha == 0) ?
           ((dist_squared > gamma * gamma) ? 0 : 1 / (1 + exp(dist_squared))) :
           alpha / (adk * (adk + 1)));
  //multModify(holder, D, (grad > gamma / 4) ? gamma / 4 : grad);
  multModify(holder, D, grad);
  for (int d = 0; d != D; d++) {
    if (holder[d] > cap) holder[d] = cap;
    if (holder[d] < - cap) holder[d] = - cap;
  }
  return false;
};

void positiveGradient(double* i, double* j,
                      double* holder,
                      const double alpha,
                      const int D) {
  const double dist_squared = distAndVector(i, j, holder, D);
  const double grad = (alpha == 0) ?
        ((dist_squared > 4) ? -1 : -(exp(dist_squared) / (exp(dist_squared) + 1))) :
        -2 * alpha / (1 + alpha * dist_squared);
  multModify(holder, D, grad);
};

// [[Rcpp::export]]
arma::vec testPositiveGradient(arma::vec i, arma::vec j,
                               NumericVector alpha, NumericVector f) {
  double a = alpha[0];
  vec holder = vec(i.size());
  positiveGradient(i.memptr(), j.memptr(), holder.memptr(), a, i.size());
  return holder;
};
// [[Rcpp::export]]
arma::vec testNegativeGradient(arma::vec i, arma::vec j,
                               NumericVector alpha, NumericVector gamma, NumericVector f) {
  double a = alpha[0];
  double g = gamma[0];
  vec holder = vec(i.size());
  negativeGradient(i.memptr(), j.memptr(), holder.memptr(), a, g, i.size(), g / 4);
  return holder;
};
