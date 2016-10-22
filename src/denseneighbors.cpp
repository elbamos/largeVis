// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "neighbors.h"
#include "distance.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

class EuclideanAdder : public DistanceAdder<Mat<double>, Col<double>> {
protected:
	virtual distancetype distanceFunction(const Col<double>& x_i, const Col<double>& x_j) const {
		return relDist(x_i, x_j);
	}
public:
	EuclideanAdder(const Mat<double>& data, const kidxtype& K) : DistanceAdder(data, K) {}
};

class CosineAdder : public DistanceAdder<Mat<double>, Col<double>> {
protected:
	virtual distancetype distanceFunction(const Col<double>& x_i, const Col<double>& x_j) const {
		return cosDist(x_i, x_j);
	}
public:
	CosineAdder(const Mat<double>& data, const kidxtype& K) : DistanceAdder(data, K) {}
};

class DenseAnnoySearch : public AnnoySearch<mat, vec> {
protected:
	virtual vec hyperplane(const ivec& indices) {
		vec direction = vec(indices.size());
		vertexidxtype x1idx, x2idx;
		vec v;
		vec m;
		do {
			x1idx = indices[sample(indices.n_elem)];
			x2idx = indices[sample(indices.n_elem)];
			if (x1idx == x2idx) x2idx = indices[sample(indices.n_elem)];
			const vec x2 = data.col(x2idx);
			const vec x1 = data.col(x1idx);
			// Get hyperplane
			m =  (x1 + x2) / 2; // Base point of hyperplane
			const vec d = x1 - x2;
			v =  d / as_scalar(norm(d, 2)); // unit vector
		} while (x1idx == x2idx);

		for (vertexidxtype i = 0; i != indices.n_elem; i++) {
			const vec X = data.col(indices[i]);
			direction[i] = dot((X - m), v);
		}
		return direction;
	}
public:
	DenseAnnoySearch(const mat& data, Progress& p) : AnnoySearch(data, p) {}
};

// [[Rcpp::export]]
arma::imat searchTrees(const int& threshold,
                       const int& n_trees,
                       const int& K,
                       const int& maxIter,
                       const arma::mat& data,
                       const std::string& distMethod,
                       Rcpp::Nullable< NumericVector > seed,
                       Rcpp::Nullable< NumericVector > threads,
                       bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
  const vertexidxtype N = data.n_cols;
	DistanceAdder<mat, Col<double>>* adder;
	if (distMethod.compare(string("Cosine")) == 0) adder = new CosineAdder(data, K);
	else adder = new EuclideanAdder(data, K);

  Progress p((N * n_trees) + (3 * N) + (N * maxIter), verbose);

	mat dataMat;
	if (distMethod.compare(string("Cosine")) == 0) dataMat = normalise(data, 2, 0);
	else dataMat = data;
	DenseAnnoySearch annoy = DenseAnnoySearch(dataMat, p);
	annoy.setSeed(seed);
	annoy.trees(n_trees, threshold);
	annoy.reduce(K, adder);
	imat ret = annoy.exploreNeighborhood(maxIter, adder);
	delete adder;
	return ret;
}
