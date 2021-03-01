#include "largeVis.h"
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;

distancetype dist(const arma::vec& i, const arma::vec& j);
distancetype relDist(const arma::vec& i, const arma::vec& j);
distancetype cosDist(const arma::vec& i, const arma::vec& j);
distancetype sparseDist(const arma::sp_mat& i, const arma::sp_mat& j);
distancetype sparseRelDist(const arma::sp_mat& i, const arma::sp_mat& j);
distancetype sparseCosDist(const arma::sp_mat& i, const arma::sp_mat& j);

// Exported distance functions for high dimensional space
arma::vec fastDistance(const NumericVector is,
                       const NumericVector js,
                       const arma::mat& data,
                       const std::string& distMethod,
                       bool verbose);
arma::vec fastSparseDistance(const arma::vec& is,
                             const arma::vec& js,
                             const arma::sp_mat& data,
                             const std::string& distMethod,
                             bool verbose);
arma::vec fastCDistance(const arma::vec& is,
                        const arma::vec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& p_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose);
arma::vec fastSDistance(const arma::vec& is,
                        const arma::vec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& j_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose);


