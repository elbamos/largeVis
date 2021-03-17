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

typedef SparseAnnoySearch<Euclidean> SparseEuclidean ;
typedef SparseAnnoySearch<Angular> SparseCosine ;

imat searchTreesSparse(
                        const int& n_trees,
                        const kidxtype& K,
                        const int& maxIter,
                        const sp_mat& data,
                        const string& distMethod,
                        const Rcpp::Nullable< Rcpp::String >& saveFile,
                        bool verbose) {
	const vertexidxtype N = data.n_cols;

	Progress p((N * n_trees) + (3 * N) + (N * maxIter), verbose);

	sp_mat dataMat;


	if (distMethod.compare(string("Cosine")) == 0) {
		dataMat = sp_mat(data);
		for (arma::uword d = 0; d < dataMat.n_cols; d++) dataMat.col(d) /= norm(dataMat.col(d));
		SparseCosine* annoy = new SparseCosine(dataMat, K, verbose, maxIter, n_trees);
		annoy->trees(n_trees, saveFile);
		annoy->reduce();
		annoy->exploreNeighborhood(maxIter);
		imat ret = annoy->sortAndReturn();
		delete annoy;
		return ret;
	} else {
		SparseEuclidean* annoy = new SparseEuclidean(data, K, verbose, maxIter, n_trees);
		annoy->trees(n_trees, saveFile);
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
                             const Rcpp::Nullable< Rcpp::String > &saveFile,
                             bool verbose) {
  const vertexidxtype N = p.size() -1;
  const sp_mat data = sp_mat(i,p,x,N,N);
  return searchTreesSparse(n_trees,K,maxIter,data,distMethod,saveFile, verbose);
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
                             const Rcpp::Nullable< Rcpp::String > &saveFile,
                             bool verbose) {
  const umat locations = join_cols(i,j);
  const sp_mat data = sp_mat(locations,x);
  return searchTreesSparse(n_trees,K,maxIter,data,distMethod,saveFile, verbose);
}
