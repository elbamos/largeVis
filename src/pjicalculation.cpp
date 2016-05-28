#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp,cpp11)]]
#include "progress.hpp"
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <queue>
#include <vector>
#include <set>
#include "helpers.h"
using namespace Rcpp;
using namespace std;

/*
 * Fast calculation of pairwise distances with the result stored in a pre-allocated vector.
 */
// [[Rcpp::export]]
arma::vec fastDistance(const NumericVector is,
              const NumericVector js,
              const arma::mat& data,
              const std::string& distMethod,
              bool verbose) {

  Progress p(is.size(), verbose);
  arma::vec xs = arma::vec(is.size());
  double (*distanceFunction)(const arma::vec& x_i, const arma::vec& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = dist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = cosDist;

  #pragma omp parallel for shared (xs)
  for (int i=0; i < is.length(); i++) if (p.increment()) xs[i] =
    distanceFunction(data.col(is[i]), data.col(js[i]));
  return xs;
};

arma::vec fastSparseDistance(const arma::vec& is,
                             const arma::vec& js,
                             const arma::sp_mat& data,
                             const std::string& distMethod,
                             bool verbose) {

  Progress p(is.size(), verbose);
  arma::vec xs = arma::vec(is.size());
  double (*distanceFunction)(
      const arma::sp_mat& x_i,
      const arma::sp_mat& x_j);
  if (distMethod.compare(std::string("Euclidean")) == 0) distanceFunction = sparseDist;
  else if (distMethod.compare(std::string("Cosine")) == 0) distanceFunction = sparseCosDist;

#pragma omp parallel for shared (xs)
  for (int i=0; i < is.size(); i++) if (p.increment()) xs[i] =
    distanceFunction(data.col(is[i]), data.col(js[i]));
  return xs;
};

// [[Rcpp::export]]
arma::vec fastCDistance(const arma::vec& is,
                       const arma::vec& js,
                       const arma::uvec& i_locations,
                       const arma::uvec& p_locations,
                       const arma::vec& x,
                       const std::string& distMethod,
                       bool verbose) {
  const int N = p_locations.size() - 1;
  const arma::sp_mat data = arma::sp_mat(i_locations, p_locations, x, N, N);
  return fastSparseDistance(is,js,data,distMethod,verbose);
}

// [[Rcpp::export]]
arma::vec fastSDistance(const arma::vec& is,
                        const arma::vec& js,
                        const arma::uvec& i_locations,
                        const arma::uvec& j_locations,
                        const arma::vec& x,
                        const std::string& distMethod,
                        bool verbose) {
  const arma::umat locations = arma::join_cols(i_locations, j_locations);
  const arma::sp_mat data = arma::sp_mat(locations, x);
  return fastSparseDistance(is,js,data,distMethod,verbose);
}

// Take four vectors (i indices, j indices, edge distances, and sigmas), and calculate
// p(j|i) and then w_{ij}.
// [[Rcpp::export]]
arma::sp_mat distMatrixTowij(
    const NumericVector is,
    const NumericVector js,
    const NumericVector xs,
    const NumericVector sigmas,
    const int N,
    bool verbose
) {

  Progress p(xs.size() * 2, verbose);
  arma::vec rowSums = arma::vec(N);
  arma::vec pjis = arma::vec(is.length());
  for (int idx=0; idx < N; idx++) rowSums[idx] = 0;
  // Compute pji, accumulate rowSums at the same time
  #pragma omp parallel for shared(pjis, rowSums)
  for (int e=0; e < pjis.size(); e++) if (p.increment()){
    const int i = is[e];
    const double pji = exp(- pow(xs[e], 2)) / sigmas[i];
    pjis[e] = pji;
    #pragma omp atomic
    rowSums[i] += pji;
  }
  if (p.check_abort()) return arma::sp_mat(0);
  // Now convert p(j|i) to w_{ij} by symmetrizing.
  // Loop through the edges, and populate a location matrix and value vector for
  // the arma::sp_mat batch insertion constructor.  Put all coordinates in the
  // lower triangle.  The constructor will automatically add duplicates.
  arma::vec values = arma::vec(pjis.size());
  arma::umat locations = arma::umat(2, pjis.size());
  #pragma omp parallel for shared(locations, values)
  for (int e = 0; e < pjis.size(); e++) if (p.increment()) {
    int newi = is[e], newj = js[e];
    if (newi < newj) std::swap(newi, newj);
    values[e] =  ((pjis[e] / rowSums[is[e]]) / (2 * N));
    locations(1,e) = newi;
    locations(0,e) = newj;
  }
  arma::sp_mat wij = arma::sp_mat(
    true, // add_values
    locations,
    values,
    N, N // n_col and n_row
    );
  wij = wij + wij.t();
  return wij;
};


// [[Rcpp::export]]
double sigFunc(const double sigma,
               const NumericVector x_i,
               const double perplexity) {
  const NumericVector xs = exp(- pow(x_i,2) / sigma);
  const NumericVector softxs = xs / sum(xs);
  const double p2 = - sum(log(softxs) / log(2)) / xs.length();
  return pow(perplexity - p2, 2);
};
