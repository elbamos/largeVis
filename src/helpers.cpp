#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
using namespace Rcpp;

double dist(const arma::vec& i, const arma::vec& j) {
  return sqrt(sum(square(i - j)));
}

double cosDist(const arma::vec& i, const arma::vec& j) {
  return 1 - (dot(i,j) / (sqrt(sum(square(i))) * sqrt(sum(square(j)))));
}
