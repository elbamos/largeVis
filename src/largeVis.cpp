// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

/*
* The stochastic gradient descent function.
*/
// [[Rcpp::export]]
arma::mat sgd(arma::mat coords,
              arma::ivec& targets_i, // vary randomly
              const IntegerVector sources_j, // ordered
              const IntegerVector ps, // N+1 length vector of indices to start of each row j in vector is
              const arma::vec weights, // w{ij}
              const double gamma,
              const double rho,
              const double minRho,
              const bool useWeights,
              const long nBatches,
              const int M,
              const double alpha,
              bool verbose) {

  Progress progress(nBatches, verbose);

  const int D = coords.n_rows;
  if (D > 10) stop("Low dimensional space cannot have more than 10 dimensions.");
  const int N = ps.size() - 1;
  const int E = weights.size();
  double *coordsPtr = coords.memptr();

  double* negProb = new double[N];
  int* negAlias = new int[N];
  makeAliasTable(N, pow(diff(ps), 0.75), negProb, negAlias);
  double* posProb = new double[E];
  int* posAlias = new int[E];
  if (! useWeights) makeAliasTable(E, weights, posProb, posAlias);
  else for (int i = 0; i < E; i++) posAlias[i] = 1;

  const int posSampleLength = ((nBatches > 1000000) ? 1000000 : (int) nBatches);
  mat positiveSamples = randu<mat>(2, posSampleLength);
  double *posRandomPtr = positiveSamples.memptr();
  
  const double cap = gamma / 4;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) shared (coords, positiveSamples)
#endif
  for (long eIdx=0; eIdx < nBatches; eIdx++) if (progress.increment()) {

    const int e_ij = searchAliasTable(posRandomPtr + ((eIdx % posSampleLength) * 2),
                                      posProb,
                                      posAlias,
                                      E);

    const int i = targets_i[e_ij];
    const int j = sources_j[e_ij];

    // mix weight into learning rate
    const double localRho =  ((useWeights) ? weights[e_ij] : 1.0) * (rho - ((rho - minRho) * eIdx / nBatches));

    double *y_i = coordsPtr + (i * D);
    double *y_j = coordsPtr + (j * D);

    double firstholder[10];
    double secondholder[10];

    positiveGradient(y_i, y_j, firstholder, alpha, D);

    for (int d = 0; d < D; d++) y_j[d] -= firstholder[d] * localRho;

    mat negSamples = mat(2, M * 2);
    double *samplesPtr = negSamples.memptr();
    int sampleIdx = 0;
    ivec searchRange = targets_i.subvec(ps[i], ps[i + 1] - 1);
    ivec::iterator searchBegin = searchRange.begin();
    ivec::iterator searchEnd = searchRange.end();
    int m = 0;
    int k;
    while (m < M) {
      if (sampleIdx % (M * 2) == 0) negSamples.randu();
      k = searchAliasTable(samplesPtr + (sampleIdx++ % (M * 2) * 2),
                           negProb,
                           negAlias,
                           N);
      // Check that the draw isn't one of i's edges
      if (k == i ||
          k == j ||
          binary_search( searchBegin,
                         searchEnd,
                         k)) continue;

      double *y_k = coordsPtr + (k * D);

      if (negativeGradient(y_i, y_k, secondholder, 
                           alpha, gamma, cap, D)) continue;

      for (int d = 0; d < D; d++) firstholder[d] += secondholder[d];
      for (int d = 0; d < D; d++) y_k[d] -= secondholder[d] * localRho;

      m++;
      if (sampleIdx > M * 10) stop("Bad sampleidx");
    }
    for (int d = 0; d < D; d++) y_i[d] += firstholder[d] * localRho;

    if (eIdx > 0 &&
        eIdx % posSampleLength == 0) positiveSamples.randu();
  }
  return coords;
};
