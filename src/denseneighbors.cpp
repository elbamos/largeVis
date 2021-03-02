#include "neighbors.h"
#include "distance.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

class DenseAnnoySearch : public AnnoySearch<arma::Mat<double>, arma::Col<double>> {
protected:
	virtual vec hyperplane(const ivec& indices) {
		const vertexidxtype I = indices.n_elem;
		vec direction = vec(I);

		const vertexidxtype idx1 = sample(I);
		vertexidxtype idx2 = sample(I - 1);
		idx2 = (idx2 >= idx1) ? (idx2 + 1) % I : idx2;

		const vec x2 = data.col(indices[idx1]);
		const vec x1 = data.col(indices[idx2]);
			// Get hyperplane
		const vec m =  (x1 + x2) / 2; // Base point of hyperplane
		const vec d = x1 - x2;
		const vec v =  d / as_scalar(norm(d, 2)); // unit vector

		for (vertexidxtype i = 0; i != I; i++) {
			const vec X = data.col(indices[i]);
			direction[i] = dot((X - m), v);
		}
		return direction;
	}
public:
	DenseAnnoySearch(const mat& data, const kidxtype& K, Progress& p) : AnnoySearch(data, K, p) {}
};

class DenseEuclidean : public DenseAnnoySearch {
protected:
	virtual distancetype distanceFunction(const Col<double>& x_i, const Col<double>& x_j) const {
		return relDist(x_i, x_j);
	}
public:
	DenseEuclidean(const Mat<double>& data, const kidxtype& K, Progress& p) : DenseAnnoySearch(data, K, p) {}
};

class DenseCosine : public DenseAnnoySearch {
protected:
	virtual distancetype distanceFunction(const Col<double>& x_i, const Col<double>& x_j) const {
		return cosDist(x_i, x_j);
	}
public:
	DenseCosine(const Mat<double>& data, const kidxtype& K, Progress& p) : DenseAnnoySearch(data, K, p) {}
};


// [[Rcpp::export]]
arma::imat searchTrees(const int& threshold,
                       const int& n_trees,
                       const int& K,
                       const int& maxIter,
                       const arma::mat& data,
                       const std::string& distMethod,
                       Rcpp::Nullable< NumericVector > seed,
                       bool verbose) {
  const vertexidxtype N = data.n_cols;

  Progress p((N * n_trees) + (3 * N) + (N * maxIter), verbose);

  const mat dataMat = (distMethod.compare(string("Cosine")) == 0) ? normalise(data) : mat();

	DenseAnnoySearch* annoy;
	if (distMethod.compare(string("Cosine")) == 0) {
		annoy = new DenseCosine(dataMat, K, p);
	} else {
		annoy = new DenseEuclidean(data, K, p);
	}

	annoy->setSeed(seed);
	annoy->trees(n_trees, threshold);
	annoy->reduce();
	annoy->exploreNeighborhood(maxIter);
	imat ret = annoy->sortAndReturn();
	delete annoy;
	return ret;
}
