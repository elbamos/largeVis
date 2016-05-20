#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp)]]
#include "progress.hpp"
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <queue>
#include <vector>
#include <set>
using namespace Rcpp;
using namespace std;

/*
 * Fast calculation of pairwise euclidean distances with the result stored in a pre-allocated vector.
 */
// [[Rcpp::export]]
arma::vec distance(const NumericVector is,
              const NumericVector js,
              const arma::mat& data,
              bool verbose) {

  Progress p(is.size(), verbose);
  arma::vec xs = arma::vec(is.size());
  #pragma omp parallel for shared (xs)
  for (int i=0; i < is.length(); i++) if (p.increment()) xs[i] = sqrt(sum(pow(data.col(is[i]) - data.col(js[i]), 2)));
  return xs;
};

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
  for (int e=0; e < pjis.size(); e++) if (p.increment()){
    const int i = is[e];
    const double pji = exp(- pow(xs[e], 2)) / sigmas[i];
    pjis[e] = pji;
    rowSums[i] = rowSums[i] + pji;
  }
  if (p.check_abort()) return arma::sp_mat(0);
  // Now convert p(j|i) to w_{ij} by symmetrizing.
  // At this point, if both ij and ji are edges, tehre will be two entries in the data structures.
  // Loop through the structure, and place each in the lower-triangle of the sparse matrix.
  arma::sp_mat wij = arma::sp_mat(N, N);
  for (int e=0; e < pjis.size(); e++) if (p.increment()) {
    int newi = is[e], newj = js[e];
    if (newi < newj) std::swap(newi, newj);
    wij(newi, newj) +=  ((pjis[e] / rowSums[is[e]]) / (2 * N));
  }
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
