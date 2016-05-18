#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp,cpp11)]]
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
void distance(const NumericVector is,
              const NumericVector js,
              NumericVector xs,
              const NumericMatrix data,
              const Function callback) {
  for (int i=0; i < is.length(); i++) {
    xs[i] = sqrt(sum(pow(data.row(is[i]) - data.row(js[i]), 2)));
    if (i > 0 && i % 1000 == 0) callback(1000);
  }
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
    const Function callback
) {
  arma::vec rowSums = arma::vec(N);
  arma::vec pjis = arma::vec(is.length());
  for (int idx=0; idx < N; idx++) rowSums[idx] = 0;
  // Compute pji, accumulate rowSums at the same time
  #pragma omp parallel for shared(pjis, rowSums)
  for (int e=0; e < pjis.size(); e++) {
    const int i = is[e];
    const double pji = exp(- pow(xs[e], 2)) / sigmas[i];
    pjis[e] = pji;
    #pragma omp atomic
    rowSums[i] = rowSums[i] + pji;
    if (e > 0 && e % 1000 == 0) callback(1000);
  }
  // Now convert p(j|i) to w_{ij} by symmetrizing.
  // At this point, if both ij and ji are edges, tehre will be two entries in the data structures.
  // Loop through the structure, and place each in the lower-triangle of the sparse matrix.
  arma::sp_mat wij = arma::sp_mat(N, N);
  for (int e=0; e < pjis.size(); e++) {
    int newi = is[e], newj = js[e];
    if (newi < newj) std::swap(newi, newj);
    wij(newi, newj) +=  ((pjis[e] / rowSums[is[e]]) / (2 * N));
    if (e % 1000 == 0) callback(1000);
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
