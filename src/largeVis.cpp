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

  AliasTable* negAlias = new AliasTable(N, pow(diff(ps), 0.75));
  AliasTable* posAlias = new AliasTable(E, weights);

  const int posSampleLength = ((nBatches > 1000000) ? 1000000 : (int) nBatches);
  mat positiveSamples = randu<mat>(2, posSampleLength);
  double *posRandomPtr = positiveSamples.memptr();

  Gradient* grad;
  if (alpha == 0) grad = new ExpGradient(gamma, D);
  else if (alpha == 1) grad = new AlphaOneGradient(gamma, D);
  else grad = new AlphaGradient(alpha, gamma, D);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
    shared (coords, positiveSamples, posAlias, negAlias)
#endif
  for (long eIdx=0; eIdx < nBatches; eIdx++) if (progress.increment()) {
    const int e_ij = posAlias -> search(posRandomPtr + ((eIdx % posSampleLength) * 2));
    const int i = targets_i[e_ij];
    const int j = sources_j[e_ij];
    double firstholder[10];
    double secondholder[10];
    // mix weight into learning rate
    const double localRho = rho - ((rho - minRho) * eIdx / nBatches);

    double *y_i = coordsPtr + (i * D);
    double *y_j = coordsPtr + (j * D);

    grad -> positiveGradient(y_i, y_j, firstholder);

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
      k = negAlias -> search(samplesPtr + (sampleIdx++ % (M * 2) * 2));
      // Check that the draw isn't one of i's edges
      if (k == i ||
          k == j ||
          binary_search( searchBegin,
                         searchEnd,
                         k)) continue;

      double *y_k = coordsPtr + (k * D);

      grad -> negativeGradient(y_i, y_k, secondholder);

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
