// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

class Random {
  const gsl_rng_type * gsl_T;
  gsl_rng * gsl_r;
public:
  Random() {
    Random(314159265);
  }
  Random(unsigned long int seed) {
    Rcout << "rand";
    gsl_rng_env_setup();
    Rcout << "a";
    gsl_T = gsl_rng_rand48;
    Rcout << "b";
    gsl_r = gsl_rng_alloc(gsl_T);
    Rcout << "c";
    gsl_rng_set(gsl_r, seed);
    Rcout << "d";
  };
  double getRandom() {
    Rcout << "e";
    double ret = gsl_rng_uniform(gsl_r);
    Rcout << "f";
    return ret;
  }
  int getRandomInt(int range) {
    return floor(gsl_rng_uniform(gsl_r) * (range - 0.1));
  }
};

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
    const int j = targets_i[e_ij];
    const int i = sources_j[e_ij];
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


inline double distAndVector(const double *x_i,
                            const double *x_j,
                            double *output,
                            int D) {
  double cnt = 0;
  for (int d = 0; d < D; d++) {
    double t = x_i[d] - x_j[d];
    output[d] = t;
    cnt += t * t;
  }
  return cnt;
}

inline double clamp(double val) {
  return fmin(fmax(val, -5.0), 5.0);
}

inline void multModify(double *col, int D, double adj) {
  for (int i = 0; i != D; i++) col[i] = clamp(col[i] * adj);
}

// [[Rcpp::export]]
arma::mat referenceSgd(arma::mat coords,
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

 // Random rand = Random();

// #ifdef _OPENMP
// #pragma omp parallel for schedule(static) \
//   shared (coords, positiveSamples, posAlias, negAlias)
// #endif
    for (long eIdx=0; eIdx < nBatches; eIdx++) if (progress.increment()) {
      const double localRho = rho - ((rho - minRho) * eIdx / nBatches);
      vec random = randu(2);
      const int e_ij = posAlias -> search(random[0], random[1]);
      const int i = sources_j[e_ij];
      int j = targets_i[e_ij];
      double firstholder[10];
      double secondholder[10];
      // mix weight into learning rate

      double *y_i = coordsPtr + (i * D);
      for (int d = 0; d < D; ++d) firstholder[i] = 0;
      for (int m = 0; m < M + 1; ++m) {
        if (m > 0) {
          random = randu(2);
          j = negAlias -> search(random[0], random[1]);
          if (j == i) continue;
        }
        double *y_j = coordsPtr + (j * D);
        double dist = distAndVector(y_i, y_j, secondholder, D);
        double grad;
        if (m == 0) grad = -2 / (1 + dist);
        else grad = 2 * gamma / (1 + dist) / (0.1 + dist);
        multModify(secondholder, D, grad);
        for (int d = 0; d < D; d++) {
          y_j[d] -= secondholder[d] * localRho;
          firstholder[d] += secondholder[d];
        }
      }
      for (int d = 0; d < D; d++) y_i[d] += firstholder[d] * localRho;
    }
    return coords;
};

