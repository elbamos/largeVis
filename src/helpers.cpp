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

double sparseDist(const arma::sp_mat& i, const arma::sp_mat& j) {
  return arma::as_scalar(sqrt(sum(square(i - j))));
}

double sparseCosDist(const arma::sp_mat& i, const arma::sp_mat& j) {
  return 1 - (arma::as_scalar((dot(i,j)) / arma::as_scalar(arma::norm(i,2) * arma::norm(j,2))));
}
