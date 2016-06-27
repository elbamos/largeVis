// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

// utility function
inline void multModify(double *col, int d, double adj) {
  for (int i = 0; i < d; i++) col[i] *= adj;
}

bool negativeGradient(double* i,
                      double* k,
                      double* holder,
                      const double alpha,
                      const double gamma,
                      const int D) {
  const double dist_ik = sqrt(distAndVector(i, k, holder, D));
  if (alpha == 0 && dist_ik > 7) { // Address issue of overflow as exp(x^2) approaches inf.
    for (int d = 0; d < D; d++) holder[d] = 0;
    return false;
  }
  const double adk = alpha * dist_ik * dist_ik;
  if (dist_ik == 0) return true; // If the two points are in the same place, skip
  // df/dd has a dist_ij factor in the numerator that cancels with dd/dx
  const double grad = gamma * 2 * ((alpha == 0) ?
           1 / (1 + exp(dist_ik * dist_ik)) :
           alpha / (adk * (adk + 1)));
  multModify(holder, D, (grad > gamma) ? gamma : grad);
  return false;
};

void positiveGradient(double* i, double* j,
                      double* holder,
                      const double alpha,
                      const int D) {
  double dist_ij = sqrt(distAndVector(i, j, holder, D));
  const double powdist = dist_ij * dist_ij;
  // df/dd has a dist_ij factor in the numerator that cancels with dd/dx
  const double grad = (-2) * ((alpha == 0) ?
                                ((dist_ij > 7) ? 1 : (exp(powdist) / (exp(powdist) + 1))) :
                                alpha / (1 + alpha * powdist));
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
  negativeGradient(i.memptr(), j.memptr(), holder.memptr(), a, g, i.size());
  return holder;
};
