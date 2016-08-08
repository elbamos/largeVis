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
  const dimidxtype D;
  const int M;
  const int M2;

  vertexidxtype * const targetPointer;
  vertexidxtype * const sourcePointer;
  coordinatetype * const coordsPtr;
  const iterationtype n_samples;

  distancetype rho;
  distancetype initialRho;
  distancetype rhoIncrement;

  AliasTable< vertexidxtype > negAlias;
  AliasTable< edgeidxtype > posAlias;
  Gradient* grad;

  IntegerVector ps;

  long actualIterationCount = 0;

public:
  Visualizer(vertexidxtype * sourcePtr,
             vertexidxtype * targetPtr,
             dimidxtype D,
             coordinatetype * coordPtr,
             int M,
             distancetype rho,
             iterationtype n_samples) : D{D}, M{M}, M2(M * 2),
                                    targetPointer{targetPtr},
                                    sourcePointer{sourcePtr},
                                    coordsPtr{coordPtr},
                                    n_samples{n_samples},
                                    rho{rho},
                                    initialRho{rho},
                                    rhoIncrement(rho / n_samples) { }
	void initAlias(IntegerVector& newps,
                const NumericVector& weights,
                Rcpp::Nullable<Rcpp::NumericVector> seed) {
		ps = newps;
		NumericVector pdiffs = pow(diff(newps), 0.75);
		negAlias.initialize(pdiffs);
		posAlias.initialize(weights);
		if (seed.isNotNull()) {
#ifdef _OPENMP
			omp_set_num_threads(1);
#endif
			long long innerSeed = Rcpp::NumericVector(seed)[0];
			innerSeed = negAlias.initRandom(innerSeed);
			posAlias.initRandom(innerSeed);
		} else {
			negAlias.initRandom();
			posAlias.initRandom();
		}
	}

  void setGradient(Gradient * newGrad) {
    grad = newGrad;
  }

  void operator()(iterationtype startSampleIdx, int batchSize) {
  	edgeidxtype e_ij;
  	int m, shortcircuit, example = 0;
  	vertexidxtype i, j, k;
  	dimidxtype d;
  	vertexidxtype * searchBegin, * searchEnd;
  	coordinatetype firstholder[10], secondholder[10], * y_i, * y_j;

    distancetype localRho = rho;
    while (example++ != batchSize && localRho > 0) {
      // * (1 - (startSampleIdx / n_samples));
      e_ij = posAlias();
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
        k = negAlias();
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
//#ifdef _OPENMP
//#pragma omp critical
//#endif/
//		{/
//	    actualIterationCount += batchSize;
//	    rho = initialRho * (1 - actualIterationCount / n_samples + 1.0);
	    // rho -= (rhoIncrement * batchSize);/
//		}
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
              const Rcpp::Nullable<Rcpp::NumericVector> seed,
              const bool verbose) {
  Progress progress(n_samples, verbose);
	dimidxtype D = coords.n_rows;
  if (D > 10) stop("Limit of 10 dimensions for low-dimensional space.");
  Visualizer v(sources_j.memptr(),
               targets_i.memptr(),
               coords.n_rows,
               (coordinatetype*) coords.memptr(),
               M,
               (distancetype) rho,
               (iterationtype) n_samples);
  v.initAlias(ps, weights, seed);
  if (alpha == 0) v.setGradient(new ExpGradient(gamma, D));
  else if (alpha == 1) v.setGradient(new AlphaOneGradient(gamma, D));
  else v.setGradient(new AlphaGradient(alpha, gamma, D));
  const int batchSize = 8192;
  const iterationtype barrier = (n_samples * .95 < n_samples - coords.n_cols) ? n_samples * .95 : n_samples - coords.n_cols;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (iterationtype eIdx = 0; eIdx < barrier; eIdx += batchSize) if (progress.increment(batchSize)) {
    v(eIdx, batchSize);
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
  for (iterationtype eIdx = barrier; eIdx < n_samples; eIdx += batchSize) if (progress.increment(batchSize)) v(eIdx, batchSize);
  return coords;
};

