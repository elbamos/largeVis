// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "alias.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

class Visualizer {
protected:
  const dimidxtype D;
  const int M;

  vertexidxtype * const targetPointer;
  vertexidxtype * const sourcePointer;
  coordinatetype * const coordsPtr;
  const iterationtype n_samples;

  distancetype rho;
  const distancetype rhoIncrement;

  AliasTable< vertexidxtype, coordinatetype, double > negAlias;
  AliasTable< edgeidxtype, coordinatetype, double > posAlias;
  Gradient* grad;

  vertexidxtype* ps;

  int storedThreads = 0;

public:
  Visualizer(vertexidxtype * sourcePtr,
             vertexidxtype * targetPtr,
             const dimidxtype& D,
             coordinatetype * coordPtr,
             const int& M,
             const distancetype& rho,
             vertexidxtype* ps,
             const iterationtype& n_samples) : D{D}, M{M},
                                    targetPointer{targetPtr},
                                    sourcePointer{sourcePtr},
                                    coordsPtr{coordPtr},
                                    n_samples{n_samples},
                                    rho{rho},
                                    rhoIncrement((rho - 0.0001) / n_samples),
                                    ps{ps} { }
	virtual ~Visualizer() {
#ifdef _OPENMP
		if (storedThreads > 0) omp_set_num_threads(storedThreads);
#endif
		delete grad;
	}

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

  void setGradient(const double& alpha, const double& gamma) {
  	if (alpha == 0) grad = new ExpGradient(gamma, D);
  	else if (alpha == 1) grad = new AlphaOneGradient(gamma, D);
  	else grad = new AlphaGradient(alpha, gamma, D);
  }

  virtual void operator()(const iterationtype& startSampleIdx, const int& batchSize) {
  	edgeidxtype e_ij;
  	int m, example = 0;
  	vertexidxtype i, j, k;
  	coordinatetype firstholder[10], secondholder[10], * y_i, * y_j;

    const distancetype localRho = rho;
    const distancetype negRho = - localRho;
    while (example++ != batchSize && localRho > 0) {
      e_ij = posAlias();

      j = targetPointer[e_ij];
      i = sourcePointer[e_ij];

      y_i = coordsPtr + (i * D);
      y_j = coordsPtr + (j * D);
			grad -> positiveGradient(y_i, y_j, firstholder);
			for (dimidxtype d = 0; d != D; ++d) y_j[d] -= firstholder[d] * rho;

			m = 0;
      while (m != M) {
        k =  negAlias();

        // Check that the draw isn't one of i's edges
        if (k == i ||
            k == j) continue;
        m++;

        y_j = coordsPtr + (k * D);
        grad -> negativeGradient(y_i, y_j, secondholder);

        for (dimidxtype d = 0; d != D; ++d) y_j[d] -= secondholder[d] * rho;
        for (dimidxtype d = 0; d != D; ++d) firstholder[d] += secondholder[d];
      }
      for (dimidxtype d = 0; d != D; ++d) y_i[d] -= firstholder[d] * rho;
    }
    rho -= (rhoIncrement * batchSize);
  }
};

class MomentumVisualizer : public Visualizer {
protected:
	const float momentum;
	coordinatetype* momentumarray;

	void updateMinus(coordinatetype* from,
                   coordinatetype* to,
                   const vertexidxtype& i,
                   const distancetype& rho) {
		coordinatetype* moment = momentumarray + (i * D);
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
            																			M, rho, ps, n_samples), momentum{momentum} {
		momentumarray = new coordinatetype[D * N];
		for (vertexidxtype i = 0; i != D*N; ++i) momentumarray[i] = 0;
	}

	~MomentumVisualizer() {
		delete[] momentumarray;
	}

	virtual void operator()(const iterationtype& startSampleIdx, const int& batchSize) {
		edgeidxtype e_ij;
		int m, example = 0;
		vertexidxtype i, j, k;
		coordinatetype firstholder[10], secondholder[10], * y_i, * y_j;

		const distancetype localRho = rho;
		const distancetype negRho = - localRho;
		while (example++ != batchSize && localRho > 0) {
			e_ij = posAlias();

			j = targetPointer[e_ij];
			i = sourcePointer[e_ij];

			y_i = coordsPtr + (i * D);
			y_j = coordsPtr + (j * D);
			grad -> positiveGradient(y_i, y_j, firstholder);
			updateMinus(firstholder, y_j, j, localRho);

			m = 0;
			while (m != M) {
				k =  negAlias();

				// Check that the draw isn't one of i's edges
				if (k == i ||
        k == j) continue;
				m++;

				y_j = coordsPtr + (k * D);
				grad -> negativeGradient(y_i, y_j, secondholder);

				updateMinus(secondholder, y_j, k, localRho);
				for (dimidxtype d = 0; d != D; ++d) firstholder[d] += secondholder[d];
			}
			updateMinus(firstholder, y_i, i, negRho);
		}
		rho -= (rhoIncrement * batchSize);
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
  vertexidxtype N = coords.n_cols;
  distancetype* negweights = new distancetype[N];
  for (vertexidxtype n = 0; n < N; ++n) negweights[n] = 0;
  if (useDegree) {
  	for (edgeidxtype e = 0; e < targets_i.n_elem; ++e) negweights[targets_i[e]]++;
  } else {
  	for (vertexidxtype p = 0; p < N; ++p) {
  		for (edgeidxtype e = ps[p]; e != ps[p + 1]; ++e) {
  			negweights[p] += weights[e];
  		}
  	}
  }
  for (vertexidxtype n = 0; n < N; ++n) negweights[n] = pow(negweights[n], 0.75);
  v -> initAlias(weights.memptr(), negweights, N, sources_j.n_elem, seed);
  delete[] negweights;

  v -> setGradient(alpha, gamma);

  const int batchSize = 8192;
  const iterationtype barrier = (n_samples * .99 < n_samples - coords.n_cols) ? n_samples * .99 : n_samples - coords.n_cols;
  iterationtype eIdx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
  for (eIdx = 0; eIdx < barrier; eIdx += batchSize) if (progress.increment(batchSize)) {
#else
  for (eIdx = 0; eIdx < n_samples; eIdx += batchSize) if (progress.increment(batchSize)) {
#endif
    (*v)(eIdx, batchSize);
  }
#ifdef _OPENMP
#pragma omp barrier
  for ( ; eIdx < n_samples; eIdx += batchSize) if (progress.increment(batchSize)) (*v)(eIdx, batchSize);
#endif
  delete v;
  return coords;
};

