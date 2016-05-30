#ifndef LARGEVISHELPERS
#define LARGEVISHELPERS
using namespace Rcpp;

double dist(const arma::vec& i, const arma::vec& j);

double cosDist(const arma::vec& i, const arma::vec& j);

double sparseDist(const arma::sp_mat& i, const arma::sp_mat& j);

double sparseCosDist(const arma::sp_mat& i, const arma::sp_mat& j);

#endif
