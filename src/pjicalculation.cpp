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
