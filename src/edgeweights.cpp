// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

// Take four vectors (i indices, j indices, edge distances, and sigmas), and calculate
// p(j|i) and then w_{ij}.
// [[Rcpp::export]]
arma::sp_mat distMatrixTowij( const NumericVector sources,
                              const NumericVector targets,
                              const NumericVector weights,
                              const NumericVector sigmas,
                              const int N,
                              bool verbose) {
  Progress p(weights.size() * 2, verbose);
  vec rowSums = vec(N);
  vec pjis = vec(sources.length());
  for (int idx=0; idx != N; idx++) rowSums[idx] = 0;
  // Compute pji, accumulate rowSums at the same time
#ifdef _OPENMP
#pragma omp parallel for shared(pjis, rowSums)
#endif
  for (int e=0; e < pjis.size(); e++) if (p.increment()){
    const int i = sources[e];
    const double pji = exp(- pow(weights[e], 2)) / sigmas[i];
    pjis[e] = pji;
#ifdef _OPENMP
#pragma omp atomic
#endif
    rowSums[i] += pji;
  }
  if (p.check_abort()) return sp_mat(0);
  // Now convert p(j|i) to w_{ij} by symmetrizing.
  // Loop through the edges, and populate a location matrix and value vector for
  // the sp_mat batch insertion constructor.  Put all coordinates in the
  // lower triangle.  The constructor will automatically add duplicates.
  vec values = vec(pjis.size());
  umat locations = umat(2, pjis.size());
#ifdef _OPENMP
#pragma omp parallel for shared(locations, values)
#endif
  for (int e = 0; e < pjis.size(); e++) if (p.increment()) {
    int newi = sources[e], newj = targets[e];
    if (newi < newj) swap(newi, newj);
    values[e] =  ((pjis[e] / rowSums[sources[e]]) / (2 * N));
    locations(1,e) = newi;
    locations(0,e) = newj;
  }
  sp_mat wij = sp_mat(
    true, // add_values
    locations,
    values,
    N, N // n_col and n_row
  );
  wij = wij + wij.t();
  return wij;
};


// [[Rcpp::export]]
double sigFunc(const double& sigma,
               const NumericVector& x_i,
               const double& perplexity) {
  const NumericVector xs = exp(- pow(x_i,2) / sigma);
  const NumericVector softxs = xs / sum(xs);
  const double p2 = - sum(log(softxs) / log(2)) / xs.length();
  return pow(perplexity - p2, 2);
};
