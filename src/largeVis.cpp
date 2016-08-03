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
  const int D;
  const int M;
  const int M2;

  long long * const targetPointer;
  long long * const sourcePointer;
  double * const coordsPtr;
  const long long n_samples;

  double rho;
  double rhoIncrement;

  AliasTable<int>* negAlias;
  AliasTable<long>* posAlias;
  Gradient* grad;

  IntegerVector ps;

public:
  Visualizer(long long * sourcePtr,
             long long * targetPtr,
             int D,
             double * coordPtr,
             int M,
             double rho,
             long long n_samples) : D{D}, M{M}, M2(M * 2),
                                    targetPointer{targetPtr},
                                    sourcePointer{sourcePtr},
                                    coordsPtr{coordPtr},
                                    n_samples{n_samples},
                                    rho{rho},
                                    rhoIncrement(rho / n_samples) { }

  void initAlias(IntegerVector& newps,
                 const NumericVector& weights) {
    ps = newps;
    NumericVector pdiffs = pow(diff(newps), 0.75);
    negAlias = new AliasTable<int>(pdiffs);
    posAlias = new AliasTable<long>(weights);
    negAlias -> initRandom();
    posAlias -> initRandom();
  }

  void setGradient(Gradient * newGrad) {
    grad = newGrad;
  }

  void batch(long long startSampleIdx, int batchSize) {
    long long e_ij;
    int i, j, k, d, m, shortcircuit, example = 0;
    double firstholder[10], secondholder[10];
    double * y_i, * y_j;
    long long * searchBegin, * searchEnd;

    double localRho = rho;
    while (example++ != batchSize && localRho > 0) {
      // * (1 - (startSampleIdx / n_samples));
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

      while (m != M && shortcircuit != M2) {
        k = negAlias -> sample();
        shortcircuit++;
        // Check that the draw isn't one of i's edges
        if (k == i ||
            k == j ||
            binary_search( searchBegin,
                           searchEnd,
                           k)) continue;
        m++;

        y_j = coordsPtr + (k * D);
        grad -> negativeGradient(y_i, y_j, secondholder);


        for (d = 0; d != D; d++) y_j[d] -= secondholder[d] * localRho;
        for (d = 0; d != D; d++) firstholder[d] += secondholder[d];
      }
      for (d = 0; d != D; d++) y_i[d] += firstholder[d] * localRho;
      localRho -= rhoIncrement;
    }
    rho -= (rhoIncrement * batchSize);
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
              const long long n_samples,
              const int M,
              const double alpha,
              const bool verbose) {

  Progress progress(n_samples, verbose);
  int D = coords.n_rows;
  if (D > 10) stop("Limit of 10 dimensions for low-dimensional space.");
  Visualizer* v = new Visualizer(sources_j.memptr(),
                                 targets_i.memptr(),
                                 coords.n_rows,
                                 coords.memptr(),
                                 M,
                                 rho,
                                 n_samples);
  v -> initAlias(ps, weights);

  if (alpha == 0) v -> setGradient(new ExpGradient(gamma, D));
  else if (alpha == 1) v -> setGradient(new AlphaOneGradient(gamma, D));
  else v -> setGradient(new AlphaGradient(alpha, gamma, D));

  const int batchSize = 8192;
  const long long barrier = (n_samples * .9 < n_samples - coords.n_cols) ? n_samples * .9 : n_samples - coords.n_cols;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (long long eIdx = 0; eIdx < barrier; eIdx += batchSize) if (progress.increment(batchSize)) {
    v -> batch(eIdx, batchSize);
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
  for (long long eIdx = barrier; eIdx < n_samples; eIdx += batchSize) if (progress.increment(batchSize)) v -> batch(eIdx, batchSize);
  return coords;
};

