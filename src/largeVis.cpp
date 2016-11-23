// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "alias.h"
#include "progress.hpp"
#include "gradients.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

class Visualizer {
private:
	inline void updateMinus(const coordinatetype* from,
                          coordinatetype* to,
                          const distancetype& rho) {
		for (dimidxtype d = 0; d != D; ++d) to[d] -= from[d] * rho;
	}
protected:
	const dimidxtype D;
	const unsigned int M;

	vertexidxtype * const targetPointer;
	vertexidxtype * const sourcePointer;
	coordinatetype * const coordsPtr;

	distancetype rho;
	const distancetype rhoIncrement;

	AliasTable< vertexidxtype, coordinatetype, double > negAlias;
	AliasTable< edgeidxtype, coordinatetype, double > posAlias;
	Gradient* grad;

	unsigned int storedThreads = 0;

public:
	Visualizer(vertexidxtype *sourcePtr,
            vertexidxtype *targetPtr,
            coordinatetype *coordPtr,

            const dimidxtype& D,
            const vertexidxtype& N,
            const edgeidxtype& E,

            distancetype rho,
            const iterationtype& n_samples,

            const unsigned int& M,
            const double& alpha,
            const double& gamma) : D{D}, M{M},
	            targetPointer{targetPtr},
	            sourcePointer{sourcePtr},
	            coordsPtr{coordPtr},
	            rho{rho},
	            rhoIncrement((rho - 0.0001) / n_samples),
	            negAlias(AliasTable< vertexidxtype, coordinatetype, double >(N)),
	            posAlias(AliasTable< edgeidxtype, coordinatetype, double >(E)){
    	if (alpha == 0) grad = new ExpGradient(gamma, D);
    	else if (alpha == 1) grad = new AlphaOneGradient(gamma, D);
    	else grad = new AlphaGradient(alpha, gamma, D);
    }
	virtual ~Visualizer() {
#ifdef _OPENMP
		if (storedThreads > 0) omp_set_num_threads(storedThreads);
#endif
		delete grad;
	}

	void initAlias(const distancetype* posWeights,
                 const distancetype* negWeights,
                Rcpp::Nullable<Rcpp::NumericVector> seed) {
		negAlias.initialize(negWeights);
		posAlias.initialize(posWeights);

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

	void operator()(const iterationtype& startSampleIdx, const unsigned int& batchSize) {
		edgeidxtype e_ij;
		unsigned int example = 0;
		coordinatetype *firstholder = new coordinatetype[D * 2];
		coordinatetype *secondholder = firstholder + D;

		const distancetype localRho = rho;
		if (localRho < 0) return;
		while (example++ != batchSize) {
			e_ij = posAlias();

			const vertexidxtype j = targetPointer[e_ij];
			const vertexidxtype i = sourcePointer[e_ij];

			coordinatetype* y_i = coordsPtr + (i * D);
			coordinatetype* y_j = coordsPtr + (j * D);
			grad -> positiveGradient(y_i, y_j, firstholder);
			updateMinus(firstholder, y_j, localRho);

			unsigned int m = 0;
			while (m != M) {
				const vertexidxtype k =  negAlias();

				// Check that the draw isn't one of i's edges
				if (k == i || k == j) continue;
				m++;

				y_j = coordsPtr + (k * D);
				grad -> negativeGradient(y_i, y_j, secondholder);

				updateMinus(secondholder, y_j, localRho);
				for (dimidxtype d = 0; d != D; ++d) firstholder[d] += secondholder[d];
			}
			updateMinus(firstholder, y_i, - localRho);
		}
		rho -= (rhoIncrement * batchSize);
		delete[] firstholder;
	}
};

class MomentumVisualizer : public Visualizer {
private:
	inline void updateMinus(const coordinatetype* from, const vertexidxtype& i,
                          coordinatetype* to, const distancetype& rho) {
		coordinatetype* moment = momentumarray + (i * D);
		for (dimidxtype d = 0; d != D; ++d) to[d] -= moment[d] = (moment[d] * momentum) + (from[d] * rho);
	}
protected:
	float momentum;
	coordinatetype* momentumarray;

public:
	MomentumVisualizer(vertexidxtype *sourcePtr,
                    vertexidxtype *targetPtr,
                    coordinatetype *coordPtr,

                    const dimidxtype& D,
                    const vertexidxtype& N,
                    const edgeidxtype& E,

                    distancetype rho,
                    const iterationtype& n_samples,
                    const float& momentum,

                    const unsigned int& M,
                    const double& alpha,
                    const double& gamma) : Visualizer(sourcePtr, targetPtr, coordPtr, D,
                    																	N, E, rho, n_samples, M, alpha, gamma) {
		this -> momentum = momentum;
		momentumarray = new coordinatetype[D * N];
		for (vertexidxtype i = 0; i != D*N; ++i) momentumarray[i] = 0;
	}
	~MomentumVisualizer() {
		delete[] momentumarray;
	}
	void operator()(const iterationtype& startSampleIdx, const unsigned int& batchSize) {
		edgeidxtype e_ij;
		unsigned int example = 0;
		coordinatetype *firstholder = new coordinatetype[D * 2];
		coordinatetype *secondholder = firstholder + D;
		coordinatetype * y_i, * y_j;

		const distancetype localRho = rho;
		if (localRho < 0) return;
		while (example++ != batchSize) {
			e_ij = posAlias();

			const vertexidxtype j = targetPointer[e_ij];
			const vertexidxtype i = sourcePointer[e_ij];

			y_i = coordsPtr + (i * D);
			y_j = coordsPtr + (j * D);
			grad -> positiveGradient(y_i, y_j, firstholder);
			updateMinus(firstholder, j, y_j, localRho);

			unsigned int m = 0;
			while (m != M) {
				const vertexidxtype k =  negAlias();

				// Check that the draw isn't one of i's edges
				if (k == i || k == j) continue;
				m++;

				y_j = coordsPtr + (k * D);
				grad -> negativeGradient(y_i, y_j, secondholder);

				updateMinus(secondholder, k, y_j, localRho);
				for (dimidxtype d = 0; d != D; ++d) firstholder[d] += secondholder[d];
			}
			updateMinus(firstholder, i, y_i, - localRho);
		}
		rho -= (rhoIncrement * batchSize);
		delete[] firstholder;
	}
};

// [[Rcpp::export]]
arma::mat sgd(arma::mat& coords,
              arma::ivec& targets_i, // vary randomly
              arma::ivec& sources_j, // ordered
              arma::ivec& ps, // N+1 length vector of indices to start of each row j in vector is
              arma::vec& weights, // w{ij}
              const double& gamma,
              const double& rho,
              const arma::uword& n_samples,
              const int& M,
              const double& alpha,
              const Rcpp::Nullable<Rcpp::NumericVector> momentum,
              const bool& useDegree,
              const Rcpp::Nullable<Rcpp::NumericVector> seed,
              const Rcpp::Nullable<Rcpp::NumericVector> threads,
              const bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
	Progress progress(n_samples, verbose);
	const dimidxtype D = coords.n_rows;
	const vertexidxtype N = coords.n_cols;
	const edgeidxtype E = targets_i.n_elem;

	Visualizer* v;
	if (momentum.isNull()) v = new Visualizer(
			sources_j.memptr(), targets_i.memptr(), coords.memptr(),
     	D, N, E,
     	rho, n_samples,
     	M, alpha, gamma);
	else {
		float moment = NumericVector(momentum)[0];
		if (moment < 0) stop("Momentum cannot be negative.");
		if (moment > 1) stop("Momentum canot be > 1.");
		v = new MomentumVisualizer(
			 sources_j.memptr(), targets_i.memptr(), coords.memptr(),
	     D, N, E,
	     rho, n_samples, moment,
	     M, alpha, gamma);
	}

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
	v -> initAlias(weights.memptr(), negweights, seed);
	delete[] negweights;

	const unsigned int batchSize = 8192;
	const iterationtype barrier = (n_samples * .99 < n_samples - coords.n_cols) ? n_samples * .99 : n_samples - coords.n_cols;
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
	delete v;
	return coords;
};