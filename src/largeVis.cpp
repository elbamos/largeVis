// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "alias.h"

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

// [[Rcpp::export]]
bool checkBits() {
	if (sizeof( vertexidxtype ) == 8) return true;
	else return false;
}

// [[Rcpp::export]]
bool checkOpenMP() {
#ifdef _OPENMP
	return true;
#else
	return false;
#endif
}

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

  AliasTable< vertexidxtype, coordinatetype, double > negAlias;
  AliasTable< edgeidxtype, coordinatetype, double > posAlias;
  shared_ptr< Gradient > grad;

  vertexidxtype* ps;

  int storedThreads = 0;

  virtual void updateMinus(coordinatetype* from,
                           const vertexidxtype& i,
                           const distancetype& rho) {
  	coordinatetype* to = coordsPtr + (i * D);
  	for (dimidxtype d = 0; d != D; ++d) to[d] -= from[d] * rho;
  }

public:
  Visualizer(vertexidxtype * sourcePtr,
             vertexidxtype * targetPtr,
             dimidxtype D,
             coordinatetype * coordPtr,
             int M,
             distancetype rho,
             vertexidxtype* ps,
             iterationtype n_samples) : D{D}, M{M}, M2(M * 2),
                                    targetPointer{targetPtr},
                                    sourcePointer{sourcePtr},
                                    coordsPtr{coordPtr},
                                    n_samples{n_samples},
                                    rho{rho},
                                    rhoIncrement((rho - 0.0001) / n_samples),
                                    ps{ps} { }
#ifdef _OPENMP
	~Visualizer() {
		if (storedThreads > 0) omp_set_num_threads(storedThreads);
	}
#endif

	void initAlias(const distancetype* posWeights,
                 const distancetype* negWeights,
                 const vertexidxtype& N,
                 const edgeidxtype& E,
                 Rcpp::Nullable<Rcpp::NumericVector> seed) {
		negAlias.initialize(negWeights, N);
		posAlias.initialize(posWeights, E);

		if (seed.isNotNull()) {
#ifdef _OPENMP
			storedThreads = omp_get_max_threads();
			omp_set_num_threads(1);
#endif
			vertexidxtype innerSeed = Rcpp::NumericVector(seed)[0];
			innerSeed = negAlias.initRandom(innerSeed);
			posAlias.initRandom(innerSeed);
		} else {
			negAlias.initRandom();
			posAlias.initRandom();
		}
	}

  void setGradient(double alpha, double gamma, dimidxtype D) {
  	if (alpha == 0) grad = shared_ptr< Gradient > (new ExpGradient(gamma, D));
  	else if (alpha == 1) grad = shared_ptr< Gradient > (new AlphaOneGradient(gamma, D));
  	else grad = shared_ptr< Gradient > (new AlphaGradient(alpha, gamma, D));
  }

  void operator()(iterationtype startSampleIdx, int batchSize) {
  	edgeidxtype e_ij;
  	int example = 0;
  	vertexidxtype i, j, k;

    distancetype localRho = rho;
    while (example++ != batchSize && localRho > 0) {
    	int m, shortcircuit;
    	vertexidxtype * searchBegin, * searchEnd;
    	coordinatetype firstholder[10], secondholder[10], * y_i, * y_j;

      e_ij = posAlias();

      j = targetPointer[e_ij];
      i = sourcePointer[e_ij];

      y_i = coordsPtr + (i * D);
      y_j = coordsPtr + (j * D);
			grad -> positiveGradient(y_i, y_j, firstholder);
			updateMinus(firstholder, j, localRho);

      searchBegin = targetPointer + ps[i];
      searchEnd = targetPointer + ps[i + 1];
      shortcircuit = m = 0;
      while (m != M && shortcircuit != M2) {
        k =  negAlias();
        ++shortcircuit;
        // Check that the draw isn't one of i's edges
        if (k == i ||
            k == j ||
            binary_search( searchBegin,
                           searchEnd,
                           k)) continue;
        ++m;

        y_j = coordsPtr + (k * D);
        grad -> negativeGradient(y_i, y_j, secondholder);

        updateMinus(secondholder, k, localRho);
        for (dimidxtype d = 0; d != D; ++d) firstholder[d] += secondholder[d];
      }
      updateMinus(firstholder, i, - localRho);
      localRho -= rhoIncrement;
    }
    rho -= (rhoIncrement * batchSize);
  }
};

class MomentumVisualizer : public Visualizer {
protected:
	float momentum;
	coordinatetype* momentumarray;

	virtual void updateMinus(coordinatetype* from,
                          const vertexidxtype& i,
                          const distancetype& rho) {
		coordinatetype* moment = momentumarray + (i * D);
		coordinatetype* to = coordsPtr + (i * D);
		for (dimidxtype d = 0; d != D; ++d) to[d] -= moment[d] = (moment[d] * momentum) + (from[d] * rho);
	}

public:
	MomentumVisualizer(vertexidxtype * sourcePtr,
            vertexidxtype * targetPtr,
            dimidxtype D,
            coordinatetype * coordPtr,
            int M,
            distancetype rho,
            vertexidxtype* ps,
            iterationtype n_samples,
            const float& momentum,
            const vertexidxtype& N ) : Visualizer(sourcePtr, targetPtr, D, coordPtr,
            																			M, rho, ps, n_samples) {
		this -> momentum = momentum;
		momentumarray = new coordinatetype[D * N];
		for (vertexidxtype i = 0; i != D*N; ++i) momentumarray[i] = 0;
	}
	~MomentumVisualizer() {
		delete[] momentumarray;
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
              const Rcpp::Nullable<Rcpp::NumericVector> momentum,
              const Rcpp::Nullable<Rcpp::NumericVector> seed,
              const Rcpp::Nullable<Rcpp::NumericVector> threads,
              const bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
	int useDegree = 0;
  Progress progress(n_samples, verbose);
	dimidxtype D = coords.n_rows;
  if (D > 10) stop("Limit of 10 dimensions for low-dimensional space.");
  Visualizer* v;
  if (momentum.isNull()) v = new Visualizer(  sources_j.memptr(),
															                targets_i.memptr(),
															                coords.n_rows,
															                coords.memptr(),
															                M,
															                rho,
															                ps.memptr(),
															                n_samples);
  else {
  	float moment = NumericVector(momentum)[0];
  	if (moment < 0) stop("Momentum cannot be negative.");
  	if (moment > 1) stop("Momentum canot be > 1.");
  	v = new MomentumVisualizer( sources_j.memptr(),
				                        targets_i.memptr(),
				                        coords.n_rows,
				                        coords.memptr(),
				                        M,
				                        rho,
				                        ps.memptr(),
				                        n_samples,
				                        moment,
				                        coords.n_cols);
  }
  // Make negative weights

  vertexidxtype N = coords.n_cols;
  distancetype* negweights = new distancetype[N];
  for (vertexidxtype n = 0; n < N; ++n) negweights[n] = 0;
  if (useDegree) {
  	for (edgeidxtype e = 0; e < targets_i.n_elem; ++e) negweights[targets_i[e]]++;
  } else {
  	for (vertexidxtype p = 0; p < N; ++p) {
  		for (edgeidxtype e = ps[p];
         e != ps[p + 1];
         ++e) {
  			//negweights[targets[e]] += weights[e];
  			negweights[p] += weights[e];
  		}
  	}
  }
  for (vertexidxtype n = 0; n < N; ++n) negweights[n] = pow(negweights[n], 0.75);
  v -> initAlias(weights.memptr(), negweights, N, sources_j.n_elem, seed);
  delete[] negweights;

  v -> setGradient(alpha, gamma, D);

  const int batchSize = 8192;
  const iterationtype barrier = (n_samples * .95 < n_samples - coords.n_cols) ? n_samples * .95 : n_samples - coords.n_cols;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (iterationtype eIdx = 0; eIdx < barrier; eIdx += batchSize) if (progress.increment(batchSize)) {
    (*v)(eIdx, batchSize);
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
  for (iterationtype eIdx = barrier; eIdx < n_samples; eIdx += batchSize) if (progress.increment(batchSize)) (*v)(eIdx, batchSize);
  return coords;
};

