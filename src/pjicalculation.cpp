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
void distance(NumericVector is, NumericVector js, NumericVector xs, NumericMatrix data, Function callback) {
  for (int i=0; i < is.length(); i++) {
    xs[i] = sqrt(sum(pow(data.row(is[i]) - data.row(js[i]), 2)));
    if (i > 0 && i % 1000 == 0) callback(1000);
  }
};

// Take four vectors (i indices, j indices, edge distances, and sigmas), and calculate
// p(j|i) and then w_{ij}.
// [[Rcpp::export]]
arma::sp_mat distMatrixTowij(
    NumericVector is,
    NumericVector js,
    NumericVector xs,
    NumericVector sigmas,
    int N,
    Function callback
) {
  NumericVector rowSums = NumericVector(N);
  NumericVector pjis = NumericVector(is.length());
  for (int idx=0; idx < N; idx++) rowSums[idx] = 0;
  // Compute pji, accumulate rowSums at the same time
  #pragma omp parallel for shared(pjis, rowSums)
  for (int e=0; e < pjis.length(); e++) {
    const int i = is[e];
    const double pji = exp(- pow(xs[e], 2)) / sigmas[i];
    pjis[e] = pji;
    #pragma omp atomic
    rowSums[i] = rowSums[i] + pji;
    if (e > 0 && e % 1000 == 0) callback(1000);
  }
  // Now convert p(j|i) to w_{ij} by symmetrizing.
  arma::sp_mat wij = arma::sp_mat(N, N);
  for (int e=0; e < pjis.length(); e++) {
    int newi = is[e], newj = js[e];
    if (newi < newj) {
      const int t = newi;
      newi = newj;
      newj = t;
    }
    const double oldw = wij(newi, newj);
    wij(newi, newj) = (oldw > 0) ?
    (oldw / 2) + ((pjis[e] / rowSums[is[e]]) / (2 * N )) :
      (pjis[e] / rowSums[is[e]]) / N;
    if (e > 0 && e % 1000 == 0) callback(1000);
  }
  return wij;
};


// [[Rcpp::export]]
double sigFunc(double sigma, NumericVector x_i, double perplexity) {
  NumericVector xs = exp(- pow(x_i,2) / sigma);
  NumericVector softxs = xs / sum(xs);
  double p2 = - sum(log(softxs) / log(2)) / xs.length();
  return pow(perplexity - p2, 2);
};
