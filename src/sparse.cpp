// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "neighbors.h"
#include "distance.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

class SparseAnnoySearch : public AnnoySearch<arma::SpMat<double>, arma::SpMat<double>> {
protected:
	virtual vec hyperplane(const ivec& indices) {
		const vertexidxtype I = indices.n_elem;
		vec direction = vec(I);
		const vertexidxtype x1idx  = sample(I);
		vertexidxtype x2idx = sample(I - 1);
		x2idx = (x2idx >= x1idx) ? (x2idx + 1) % I : x2idx;

		const sp_mat x2 = data.col(indices[x1idx]);
		const sp_mat x1 = data.col(indices[x2idx]);

		const sp_mat m =  (x1 + x2) / 2;
		const sp_mat d = x1 - x2;
		const distancetype dn = as_scalar(norm(d, 2)) + 1e-5;
		const sp_mat v =  d / dn; // unit vector

		for (vertexidxtype i = 0; i < I; i++) {
			const vertexidxtype I2 = indices[i];
			const sp_mat X = data.col(I2);
			direction[i] = dot((X - m), v);
		}
		return direction;
	}
public:
	SparseAnnoySearch(const sp_mat& data, const kidxtype& K, Progress& p) : AnnoySearch(data, K, p) {}
};


class SparseEuclidean : public SparseAnnoySearch {
protected:
	virtual distancetype distanceFunction(const SpMat<double> &x_i, const SpMat<double> &x_j) const {
		return sparseRelDist(x_i, x_j);
	}
public:
	SparseEuclidean(const SpMat<double>& data, const kidxtype& K, Progress& p) : SparseAnnoySearch(data, K, p) {}
};

class SparseCosine : public SparseAnnoySearch {
protected:
	virtual distancetype distanceFunction(const SpMat<double>& x_i, const SpMat<double>& x_j) const {
		return sparseCosDist(x_i, x_j);
	}
public:
	SparseCosine(const SpMat<double>& data, const kidxtype& K, Progress& p) : SparseAnnoySearch(data, K, p) {}
};

imat searchTreesSparse(const int& threshold,
                            const int& n_trees,
                            const kidxtype& K,
                            const int& maxIter,
                            const sp_mat& data,
                            const string& distMethod,
                            Rcpp::Nullable< NumericVector> seed,
                            Rcpp::Nullable< NumericVector> threads,
                            bool verbose) {
	const vertexidxtype N = data.n_cols;

	Progress p((N * n_trees) + (3 * N) + (N * maxIter), verbose);

	sp_mat dataMat;

	if (distMethod.compare(string("Cosine")) == 0) {
		dataMat = sp_mat(data);
		for (vertexidxtype d = 0; d < dataMat.n_cols; d++) dataMat.col(d) /= norm(dataMat.col(d));
	} else {
		dataMat = data;
	}

	SparseAnnoySearch* annoy;
	if (distMethod.compare(string("Cosine")) == 0) {
		annoy = new SparseCosine(dataMat, K, p);
	} else {
		annoy = new SparseEuclidean(dataMat, K, p);
	}

	annoy->setSeed(seed);
	annoy->trees(n_trees, threshold);
	annoy->reduce();
	annoy->exploreNeighborhood(maxIter);
	imat ret = annoy->sortAndReturn();
	delete annoy;
	return ret;
}

// [[Rcpp::export]]
arma::imat searchTreesCSparse(const int& threshold,
                             const int& n_trees,
                             const int& K,
                             const int& maxIter,
                             const arma::uvec& i,
                             const arma::uvec& p,
                             const arma::vec& x,
                             const std::string& distMethod,
                             Rcpp::Nullable< Rcpp::NumericVector> seed,
                             Rcpp::Nullable< Rcpp::NumericVector> threads,
                             bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
  const vertexidxtype N = p.size() -1;
  const sp_mat data = sp_mat(i,p,x,N,N);
  return searchTreesSparse(threshold,n_trees,K,maxIter,data,distMethod,seed,threads, verbose);
}

// [[Rcpp::export]]
arma::imat searchTreesTSparse(const int& threshold,
                             const int& n_trees,
                             const int& K,
                             const int& maxIter,
                             const arma::uvec& i,
                             const arma::uvec& j,
                             const arma::vec& x,
                             const std::string& distMethod,
                             Rcpp::Nullable< NumericVector> seed,
                             Rcpp::Nullable< NumericVector> threads,
                             bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
  const umat locations = join_cols(i,j);
  const sp_mat data = sp_mat(locations,x);
  return searchTreesSparse(threshold,n_trees,K,maxIter,data,distMethod,seed,threads,verbose);
}
