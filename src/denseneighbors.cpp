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

typedef DenseAnnoySearch<Euclidean> DenseEuclidean ;
typedef DenseAnnoySearch<Angular> DenseCosine ;


// [[Rcpp::export]]
arma::imat searchTrees(
                       const int& n_trees,
                       const int& K,
                       const int& maxIter,
                       const arma::mat& data,
                       const std::string& distMethod,
                       Rcpp::Nullable< Rcpp::String > &saveFile,
                       bool verbose) {

  const mat dataMat = (distMethod.compare(string("Cosine")) == 0) ? normalise(data) : mat();

	imat ret;
	if (distMethod.compare(string("Cosine")) == 0) {
		DenseAnnoySearch<Angular>* annoy = new DenseCosine(dataMat, K, verbose, maxIter, n_trees);
		annoy->trees(n_trees, saveFile);
		annoy->reduce();
		annoy->exploreNeighborhood(maxIter);
		ret = annoy->sortAndReturn();
		delete annoy;
	} else {
		DenseAnnoySearch<Euclidean>* annoy = new DenseEuclidean(data, K, verbose, maxIter, n_trees);
		annoy->trees(n_trees, saveFile);
		annoy->reduce();
		annoy->exploreNeighborhood(maxIter);
		ret = annoy->sortAndReturn();
		delete annoy;
	}

	return ret;
}
