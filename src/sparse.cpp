// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "neighbors.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

class EuclideanSparseAdder : public DistanceAdder<arma::sp_mat, arma::sp_mat> {
protected:
	virtual distancetype distanceFunction(const arma::sp_mat& x_i, const arma::sp_mat& x_j) {
		return sparseRelDist(x_i, x_j);
	}
public:
	EuclideanSparseAdder(const arma::sp_mat& data, const kidxtype K) : DistanceAdder(data, K) {}
};

class CosineSparseAdder : public DistanceAdder<arma::sp_mat, arma::sp_mat> {
protected:
	virtual distancetype distanceFunction(const arma::sp_mat& x_i, const arma::sp_mat& x_j) {
		return sparseCosDist(x_i, x_j);
	}
public:
	CosineSparseAdder(const arma::sp_mat& data, const kidxtype K) : DistanceAdder(data, K) {}
};

class SparseAnnoySearch : public AnnoySearch<arma::sp_mat, arma::sp_mat> {
protected:
	virtual arma::vec hyperplane(const arma::ivec& indices) {
		const vertexidxtype I = indices.n_elem;
		vec direction = vec(indices.size());
		{
			vertexidxtype x1idx, x2idx;
			sp_mat v;
			sp_mat m;
			do {
				x1idx = indices[sample(I)];
				x2idx = indices[sample(I)];
				if (x1idx == x2idx) x2idx = indices[sample(I)];
				const sp_mat x2 = data.col(x2idx);
				const sp_mat x1 = data.col(x1idx);
				// Get hyperplane
				m =  (x1 + x2) / 2; // Base point of hyperplane
				const sp_mat d = x1 - x2;
				const distancetype dn = as_scalar(norm(d, 2));
				v =  d / dn; // unit vector
			} while (x1idx == x2idx);

			for (vertexidxtype i = 0; i < indices.size(); i++) {
				const vertexidxtype I = indices[i];
				const sp_mat X = data.col(I);
				direction[i] = dot((X - m), v);
			}
		}
		return direction;
	}
public:
	SparseAnnoySearch(const arma::sp_mat& data, Progress& p) : AnnoySearch(data, p) {}
};

arma::imat searchTreesSparse(const int& threshold,
                            const int& n_trees,
                            const kidxtype& K,
                            const int& maxIter,
                            const sp_mat& data,
                            const string& distMethod,
                            Rcpp::Nullable< Rcpp::NumericVector> seed,
                            bool verbose) {
	const vertexidxtype N = data.n_cols;

	shared_ptr< DistanceAdder<arma::sp_mat, arma::sp_mat> >  adder;
	if (distMethod.compare(string("Cosine")) == 0) adder = shared_ptr< DistanceAdder<arma::sp_mat, arma::sp_mat> >(new CosineSparseAdder(data, K));
	else adder = shared_ptr< DistanceAdder<arma::sp_mat, arma::sp_mat> >(new EuclideanSparseAdder(data, K));

	Progress p((N * n_trees) + (3 * N) + (N * maxIter), verbose);

	sp_mat dataMat;
	if (distMethod.compare(string("Cosine")) == 0) {
		dataMat = sp_mat(data);
		for (vertexidxtype d = 0; d < dataMat.n_cols; d++) dataMat.col(d) /= norm(dataMat.col(d));
	} else {
		dataMat = data;
	}
	SparseAnnoySearch annoy = SparseAnnoySearch(dataMat, p);
	annoy.setSeed(seed);
	annoy.trees(n_trees, threshold);
	annoy.reduce(K, adder);
	annoy.convertToMatrix(K);
	return (maxIter == 0) ? annoy.getMatrix(adder) : annoy.exploreNeighborhood(maxIter, adder);
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
  return searchTreesSparse(threshold,n_trees,K,maxIter,data,distMethod,seed,verbose);
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
                             Rcpp::Nullable< Rcpp::NumericVector> seed,
                             Rcpp::Nullable< Rcpp::NumericVector> threads,
                             bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
  const umat locations = join_cols(i,j);
  const sp_mat data = sp_mat(locations,x);
  return searchTreesSparse(threshold,n_trees,K,maxIter,data,distMethod,seed,verbose);
}
