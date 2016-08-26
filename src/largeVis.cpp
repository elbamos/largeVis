// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

#ifdef _OPENMP
void checkCRAN(Rcpp::Nullable<Rcpp::NumericVector> threads) {
	if (threads.isNotNull()) {
		int nthreads = NumericVector(threads)[0];
		if (nthreads > 0) omp_set_num_threads(nthreads);
	}
}
#endif

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
  distancetype rhoIncrement;

  AliasTable< vertexidxtype > negAlias;
  AliasTable< edgeidxtype > posAlias;
  std::shared_ptr< Gradient > grad;

  vertexidxtype* ps;

  int storedThreads = 0;

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
                                    rhoIncrement((rho - 0.0001) / n_samples) { }
#ifdef _OPENMP
	~Visualizer() {
		if (storedThreads > 0) omp_set_num_threads(storedThreads);
	}
#endif

	void initAlias(arma::ivec& newps,
                 const arma::vec& weights,
                 const arma::ivec& targets,
                 Rcpp::Nullable<Rcpp::NumericVector> seed) {
		vertexidxtype N = newps.n_elem - 1;
		ps = newps.memptr();
		distancetype* negweights = new distancetype[N];
		for (vertexidxtype n = 0; n < N; n++) negweights[n] = 0;
		for (vertexidxtype p = 0; p < N; p++) {
			for (edgeidxtype e = newps[p];
        	 e != newps[p + 1];
        	 e++) {
				//negweights[targets[e]] += weights[e];
				negweights[p] += weights[e];
			}
		}
		distancetype sum_weight = 0;
		for (vertexidxtype n = 0; n < N; n++) sum_weight += negweights[n] = pow(negweights[n], 0.75);
		negAlias.initialize(negweights, N);
		posAlias.initialize(weights.memptr(), weights.n_elem);

		if (seed.isNotNull()) {
#ifdef _OPENMP
			storedThreads = omp_get_max_threads();
			omp_set_num_threads(1);
#endif
			long long innerSeed = Rcpp::NumericVector(seed)[0];
			innerSeed = negAlias.initRandom(innerSeed);
			posAlias.initRandom(innerSeed);
		} else {
			negAlias.initRandom();
			posAlias.initRandom();
		}
		delete[] negweights;
	}

  void setGradient(double alpha, double gamma, dimidxtype D) {
  	if (alpha == 0) grad = std::shared_ptr< Gradient > (new ExpGradient(gamma, D));
  	else if (alpha == 1) grad = std::shared_ptr< Gradient > (new AlphaOneGradient(gamma, D));
  	else grad = std::shared_ptr< Gradient > (new AlphaGradient(alpha, gamma, D));
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
              arma::ivec& ps, // N+1 length vector of indices to start of each row j in vector is
              arma::vec& weights, // w{ij}
              const double gamma,
              const double rho,
              const long long n_samples,
              const int M,
              const double alpha,
              const Rcpp::Nullable<Rcpp::NumericVector> seed,
              Rcpp::Nullable<Rcpp::NumericVector> threads,
              const bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
  Progress progress(n_samples, verbose);
	dimidxtype D = coords.n_rows;
  if (D > 10) stop("Limit of 10 dimensions for low-dimensional space.");
  Visualizer v( (vertexidxtype*) sources_j.memptr(),
                (vertexidxtype*) targets_i.memptr(),
                coords.n_rows,
                (coordinatetype*) coords.memptr(),
                M,
                (distancetype) rho,
                (iterationtype) n_samples);
  v.initAlias(ps, weights, targets_i, seed);
  v.setGradient(alpha, gamma, D);
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

