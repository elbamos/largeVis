#include "neighbors.h"
#include "distance.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

template<typename Distance>
class SparseAnnoySearch : public AnnoySearch<arma::SpMat<double>, arma::SpMat<double>, Distance> {
public:
	SparseAnnoySearch<Distance>(const sp_mat& data, const kidxtype& K, const bool &verbose, const int &maxIter, const int&n_trees) :
		AnnoySearch<arma::SpMat<double>, arma::SpMat<double>, Distance>(data, K, verbose, maxIter, n_trees) {}
};


class SparseEuclidean : public SparseAnnoySearch<Euclidean> {
protected:
	virtual distancetype distanceFunction(const SpMat<double> &x_i, const SpMat<double> &x_j) const {
		return sparseRelDist(x_i, x_j);
	}
public:
	SparseEuclidean(const SpMat<double>& data, const kidxtype& K, const bool &verbose, const int &maxIter, const int&n_trees) :
		SparseAnnoySearch(data, K, verbose, maxIter, n_trees) {}
};

class SparseCosine : public SparseAnnoySearch<Angular> {
protected:
	virtual distancetype distanceFunction(const SpMat<double>& x_i, const SpMat<double>& x_j) const {
		return sparseCosDist(x_i, x_j);
	}
public:
	SparseCosine(const SpMat<double>& data, const kidxtype& K, const bool &verbose, const int &maxIter, const int&n_trees) :
		SparseAnnoySearch(data, K, verbose, maxIter, n_trees) {}
};

imat searchTreesSparse(
                        const int& n_trees,
                        const kidxtype& K,
                        const int& maxIter,
                        const sp_mat& data,
                        const string& distMethod,
                        Rcpp::Nullable< NumericVector> seed,
                        bool verbose) {
	const vertexidxtype N = data.n_cols;

	Progress p((N * n_trees) + (3 * N) + (N * maxIter), verbose);

	sp_mat dataMat;


	if (distMethod.compare(string("Cosine")) == 0) {
		dataMat = sp_mat(data);
		for (arma::uword d = 0; d < dataMat.n_cols; d++) dataMat.col(d) /= norm(dataMat.col(d));
		SparseAnnoySearch<Angular>* annoy = new SparseCosine(dataMat, K, verbose, maxIter, n_trees);
		annoy->setSeed(seed);
		annoy->trees(n_trees);
		annoy->reduce();
		annoy->exploreNeighborhood(maxIter);
		imat ret = annoy->sortAndReturn();
		delete annoy;
		return ret;
	} else {
		SparseAnnoySearch<Euclidean>* annoy = new SparseEuclidean(data, K, verbose, maxIter, n_trees);
		annoy->setSeed(seed);
		annoy->trees(n_trees);
		annoy->reduce();
		annoy->exploreNeighborhood(maxIter);
		imat ret = annoy->sortAndReturn();
		delete annoy;
		return ret;
		}


}

// [[Rcpp::export]]
arma::imat searchTreesCSparse(
                             const int& n_trees,
                             const int& K,
                             const int& maxIter,
                             const arma::uvec& i,
                             const arma::uvec& p,
                             const arma::vec& x,
                             const std::string& distMethod,
                             Rcpp::Nullable< Rcpp::NumericVector> seed,
                             bool verbose) {
  const vertexidxtype N = p.size() -1;
  const sp_mat data = sp_mat(i,p,x,N,N);
  return searchTreesSparse(n_trees,K,maxIter,data,distMethod,seed, verbose);
}

// [[Rcpp::export]]
arma::imat searchTreesTSparse(
                             const int& n_trees,
                             const int& K,
                             const int& maxIter,
                             const arma::uvec& i,
                             const arma::uvec& j,
                             const arma::vec& x,
                             const std::string& distMethod,
                             Rcpp::Nullable< NumericVector> seed,
                             bool verbose) {
  const umat locations = join_cols(i,j);
  const sp_mat data = sp_mat(locations,x);
  return searchTreesSparse(n_trees,K,maxIter,data,distMethod,seed, verbose);
}
