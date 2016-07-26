// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

using namespace Rcpp;
using namespace std;
using namespace arma;


class Visualizer {
protected:
  AliasTable<int>* negAlias;
  AliasTable<long>* posAlias;
  Gradient* grad;
  int D;
  int M;

  double rho;
  double rhoIncrement;

  long long * targetPointer;
  long long * sourcePointer;
  double * coordsPtr;
  IntegerVector ps;

public:
  Visualizer(long long * newSourcePtr,
             long long * newTargetPointer,
             int newD,
             double * newCoordPtr) {
    targetPointer = newTargetPointer;
    D = newD;
    coordsPtr = newCoordPtr;
    sourcePointer = newSourcePtr;
  }

  void initAlias(IntegerVector& newps,
                 const NumericVector& weights) {
    ps = newps;
    NumericVector pdiffs = pow(diff(newps), 0.75);
    negAlias = new AliasTable<int>(pdiffs);
    posAlias = new AliasTable<long>(weights);
    negAlias -> initRandom();
    posAlias -> initRandom();
  }

  void setGradient(Gradient * newGrad,
                   int newM,
                   double newRho,
                   long batches) {
    grad = newGrad;
    M = newM;
    rho = newRho;
    rhoIncrement = newRho / batches;
  }

  void batch(int batchSize) {
    long long e_ij;
    int i, j, k, d, m, shortcircuit, pstart, pstop, example = 0;
    double firstholder[10], secondholder[10];
    double * y_i, * y_j;
    long long * searchBegin, * searchEnd;

    double localRho = rho;
#ifdef _OPENMP
#pragma omp atomic
#endif
    rho -= (rhoIncrement * batchSize);

    while (example++ != batchSize && localRho > 0) {
      e_ij = posAlias -> sample();
      j = targetPointer[e_ij];
      i = sourcePointer[e_ij];

      y_i = coordsPtr + (i * D);
      y_j = coordsPtr + (j * D);
      grad -> positiveGradient(y_i, y_j, firstholder);
      for (d = 0; d != D; d++) y_j[d] -= firstholder[d] * localRho;

      searchBegin = targetPointer + ps[i];
      searchEnd = targetPointer + ps[i + 1];
      shortcircuit = 0; m = 0;

      while (m != M && shortcircuit != 10) {
        k = negAlias -> sample();
        shortcircuit++;
        // Check that the draw isn't one of i's edges
        if (k == i ||
            k == j ||
            binary_search( searchBegin,
                           searchEnd,
                           k)) continue;

        y_j = coordsPtr + (k * D);
        grad -> negativeGradient(y_i, y_j, secondholder);

        for (d = 0; d != D; d++) firstholder[d] += secondholder[d];
        for (d = 0; d != D; d++) y_j[d] -= secondholder[d] * localRho;

        m++;
      }
      for (d = 0; d != D; d++) y_i[d] += firstholder[d] * localRho;
      localRho -= rhoIncrement;
    }
  }
};

// [[Rcpp::export]]
arma::mat sgd(arma::mat coords,
              arma::ivec& targets_i, // vary randomly
              arma::ivec& sources_j, // ordered
              IntegerVector& ps, // N+1 length vector of indices to start of each row j in vector is
              NumericVector& weights, // w{ij}
              const double gamma,
              const double rho,
              const long nBatches,
              const int M,
              const double alpha,
              const bool verbose) {

  Progress progress(nBatches, verbose);
  int D = coords.n_rows;
  if (D > 10) stop("Limit of 10 dimensions for low-dimensional space.");
  Visualizer* v = new Visualizer(sources_j.memptr(),
                                 targets_i.memptr(),
                                 coords.n_rows,
                                 coords.memptr());
  v -> initAlias(ps, weights);

  if (alpha == 0) v -> setGradient(new ExpGradient(gamma, D), M, rho, nBatches);
  else if (alpha == 1) v -> setGradient(new AlphaOneGradient(gamma, D), M, rho, nBatches);
  else v -> setGradient(new AlphaGradient(alpha, gamma, D), M, rho, nBatches);

  const int batchSize = 16384;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (long eIdx = 0; eIdx < nBatches; eIdx += batchSize) if (progress.increment(batchSize)) {
    v -> batch(batchSize);
  }
  return coords;
};

