// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "hdbscan.h"
#include "primsalgorithm.h"
//#define DEBUG

// [[Rcpp::export]]
List hdbscanc(const arma::sp_mat& edges,
              const IntegerMatrix& neighbors,
              const int& K,
              const int& minPts,
              const Rcpp::Nullable<Rcpp::NumericVector> threads,
              const bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
	HDBSCAN object = HDBSCAN(edges.n_cols, verbose);
	// 1 N
	IntegerVector tree = object.build(K, edges, minPts, neighbors); // 4N
	NumericMatrix clusters = NumericMatrix(2, edges.n_cols);
	object.condenseAndExtract(minPts, REAL(clusters)); // 3N
	List hierarchy = object.getHierarchy();
	return List::create(Named("clusters") = clusters,
                      Named("tree") = IntegerVector(tree),
                      Named("hierarchy") = hierarchy);
}
