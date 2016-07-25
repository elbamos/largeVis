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
              NumericVector& weights, // w{ij}
              const double gamma,
              const double rho,
              const double minRho,
              const long nBatches,
              const int M,
              const double alpha,
              const bool verbose) {

  Progress progress(nBatches, verbose);

  const int D = coords.n_rows;
  if (D > 10) stop("Low dimensional space cannot have more than 10 dimensions.");
  const int N = ps.size() - 1;
  const long E = weights.size();
  double * const coordsPtr = coords.memptr();
  
  const long blockSize = 1024;
  const long nBlocks = nBatches / blockSize;

  NumericVector pdiffs = pow(diff(ps), 0.75);
  AliasTable<int>* const negAlias = new AliasTable<int>(N, pdiffs);
  AliasTable<long>* const posAlias = new AliasTable<long>(E, weights);
  negAlias -> initRandom();
  posAlias -> initRandom();

  // const int posSampleLength = ((nBatches > 1000000) ? 1000000 : (int) nBatches);
  // mat positiveSamples = randu<mat>(2, posSampleLength);
  // double * const posRandomPtr = positiveSamples.memptr();

  Gradient* grad;
  if (alpha == 0) grad = new ExpGradient(gamma, D);
  else if (alpha == 1) grad = new AlphaOneGradient(gamma, D);
  else grad = new AlphaGradient(alpha, gamma, D);
  
  long long  * targetPointer = targets_i.memptr();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1024) \
    shared (grad, coords, posAlias, negAlias)
#endif
  for (long eIdx=0; eIdx != nBlocks; eIdx++) if (progress.increment(blockSize)) {
    long e_ij;
    int i, j, k, d, m, shortcircuit;
    double firstholder[10], secondholder[10];
    double localRho = rho - ((rho - minRho) * eIdx * blockSize / nBatches);
    const double rhoIncrement = (rho - minRho) / nBatches;
    double * y_i, * y_j, * y_k;
    //ivec searchRange;
    long long * searchBegin, * searchEnd;
    
    for (long b = 0; b != blockSize; b++) {
      e_ij = posAlias -> sample(); // posAlias -> search(posRandomPtr + ((eIdx % posSampleLength) * 2));
      j = targets_i[e_ij];
      i = sources_j[e_ij];

      y_i = coordsPtr + (i * D);
      y_j = coordsPtr + (j * D);

      grad -> positiveGradient(y_i, y_j, firstholder);

      for (d = 0; d != D; d++) y_j[d] -= firstholder[d] * localRho;

      // mat negSamples = mat(2, M * 2);
      // double * const samplesPtr = negSamples.memptr();
      //searchRange = targets_i.subvec(ps[i], ps[i + 1] - 1);
      searchBegin = targetPointer + ps[i]; 
      searchEnd = targetPointer + ps[i + 1];
      shortcircuit = 0;

      while (m != M && shortcircuit != 10) {
        // if (sampleIdx == 0) negSamples.randu();
  //k = negAlias -> search(samplesPtr + (sampleIdx++ % (M * 2) * 2));
  			// k = negAlias -> search(samplesPtr + sampleIdx);
  			k = negAlias -> sample();
  			// sampleIdx = (sampleIdx + 2) % (M * 2);
  			shortcircuit++;
        // Check that the draw isn't one of i's edges
        if (k == i ||
            k == j ||
            binary_search( searchBegin,
                           searchEnd,
                           k)) continue;
  
        y_k = coordsPtr + (k * D);
  
        grad -> negativeGradient(y_i, y_k, secondholder);
  
        for (d = 0; d != D; d++) firstholder[d] += secondholder[d];
        for (d = 0; d != D; d++) y_k[d] -= secondholder[d] * localRho;
  
        m++;
      }
      for (d = 0; d != D; d++) y_i[d] += firstholder[d] * localRho;
      localRho -= rhoIncrement;
  
      // if (eIdx != 0 &&
      //     eIdx % posSampleLength == 0) positiveSamples.randu();
    }
  }
  return coords;
};

