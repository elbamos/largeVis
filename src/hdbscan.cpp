#include "largeVis.h"
#include "hdbscan.h"
#include "primsalgorithm.h"
//#define DEBUG

// [[Rcpp::export]]
List hdbscanc(const arma::sp_mat& edges,
              const IntegerMatrix& neighbors,
              const int& K,
              const int& minPts,
              const bool verbose) {
	HDBSCAN object = HDBSCAN(edges.n_cols, verbose);
	// 1 N
	IntegerVector tree = object.build(K, edges, minPts, neighbors); // 4N

	IntegerVector clusters = IntegerVector(edges.n_cols);
	NumericVector lambdas = NumericVector(edges.n_cols);
	object.condenseAndExtract(minPts, INTEGER(clusters), REAL(lambdas)); // 3N
	List hierarchy = object.getHierarchy();
	return List::create(Named("clusters") = clusters,
                      Named("lambdas") = lambdas,
                      Named("tree") = IntegerVector(tree),
                      Named("hierarchy") = hierarchy);
}
