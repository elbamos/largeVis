#include "neighbors.h"
#include "distance.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

template<typename Distance>
class DenseAnnoySearch : public AnnoySearch<arma::Mat<double>, arma::Col<double>, Distance> {

public:
	DenseAnnoySearch<Distance>(const mat& data, const kidxtype& K, const bool &verbose, const int &maxIter, const int&n_trees) :
		AnnoySearch<arma::Mat<double>, arma::Col<double>, Distance>(data, K, verbose, maxIter, n_trees) {}
};

class DenseEuclidean : public DenseAnnoySearch<Euclidean> {
public:
	DenseEuclidean(const Mat<double>& data, const kidxtype& K, const bool &verbose, const int &maxIter, const int&n_trees) :
		DenseAnnoySearch(data, K, verbose, maxIter, n_trees) {}
};

class DenseCosine : public DenseAnnoySearch<Angular> {
protected:
	virtual distancetype distanceFunction(const Col<double>& x_i, const Col<double>& x_j) const {
		return cosDist(x_i, x_j);
	}
public:
	DenseCosine(const Mat<double>& data, const kidxtype& K, const bool &verbose, const int &maxIter, const int&n_trees) :
		DenseAnnoySearch(data, K, verbose, maxIter, n_trees) {}
};


// [[Rcpp::export]]
arma::imat searchTrees(
                       const int& n_trees,
                       const int& K,
                       const int& maxIter,
                       const arma::mat& data,
                       const std::string& distMethod,
                       Rcpp::Nullable< NumericVector > seed,
                       bool verbose) {

  const mat dataMat = (distMethod.compare(string("Cosine")) == 0) ? normalise(data) : mat();

	imat ret;
	if (distMethod.compare(string("Cosine")) == 0) {
		DenseAnnoySearch<Angular>* annoy = new DenseCosine(dataMat, K, verbose, maxIter, n_trees);
		annoy->setSeed(seed);
		annoy->trees(n_trees);
		annoy->reduce();
		annoy->exploreNeighborhood(maxIter);
		ret = annoy->sortAndReturn();
		delete annoy;
	} else {
		DenseAnnoySearch<Euclidean>* annoy = new DenseEuclidean(data, K, verbose, maxIter, n_trees);
		annoy->setSeed(seed);
		annoy->trees(n_trees);
		annoy->reduce();
		annoy->exploreNeighborhood(maxIter);
		ret = annoy->sortAndReturn();
		delete annoy;
	}

	return ret;
}
