// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "largeVis.h"
#include "pq.h"

class HDBSCAN : public UF<long long> {
private:
	double* coreDistances;

public:
  HDBSCAN(const int& N,
          const bool& verbose) : UF(N, verbose, true) {
  	coreDistances = new double[N];
  }
	virtual ~HDBSCAN() {
		delete[] coreDistances;
	}

  void makeCoreDistances(const sp_mat& edges,
                         const IntegerMatrix& neighbors,
                         const int& K) {
  	if (neighbors.nrow() < K) stop("Specified K bigger than the number of neighbors in the adjacency matrix.");
 		//if (K < 4) stop("K must be >= 4 when used with neighbors.");
    IntegerVector kthNeighbors = neighbors.row(K - 1);
    for (long long n = 0; n < N; n++) if (p.increment()) {
    	long long q = kthNeighbors[n];
    	if (q == -1 || q == NA_INTEGER) stop("Insufficient neighbors.");
      coreDistances[n] = edges(n, q);
    	if (coreDistances[n] == 0) coreDistances[n] = max(edges(q, n), 1e-5);
    }
  }

  IntegerVector build( const unsigned int& K,
				               const sp_mat& edges,
				               const IntegerMatrix& neighbors) {
  	makeCoreDistances(edges, neighbors, K);
  	PrimsAlgorithm<long long, double> prim = PrimsAlgorithm<long long, double>(N, coreDistances);
  	const long long* minimum_spanning_tree = prim.run(edges, neighbors, p, 0);
  	vector< pair<double, long long> > mergeSequence = prim.getMergeSequence();
  	buildHierarchy(mergeSequence, minimum_spanning_tree); // 2 N
  	vector<long long> treevector(minimum_spanning_tree, minimum_spanning_tree + N);
  	return IntegerVector(treevector.begin(), treevector.end());
  }

	void condenseAndExtract(const unsigned int& minPts, double* clusters) {
    condense(minPts); // 2 N
    determineStability(minPts); // N
    extractClusters(minPts); // N
    getClusters(clusters);
  }

  Rcpp::List reportHierarchy() {
    long long survivingClusterCnt = survivingClusters.size();
    vector<int> parent(survivingClusterCnt);
    vector<int> nodeMembership(N);
    vector<double> stabilities(survivingClusterCnt);
    vector<int> selected(survivingClusterCnt);
    set< long long >::iterator it;
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp single private(it)
{
#endif
    for (it = roots.begin();
         it != roots.end();
         ++it)
#ifdef _OPENMP
#pragma omp task
#endif
{
    	reportAHierarchy(*it, *it, parent, nodeMembership, stabilities, selected, 0, 0);
}
#ifdef _OPENMP
}
}
#endif
    NumericVector lambdas = NumericVector(N);
    // FIXME - need to adjust this to be lambda_p not birth
    for (long long n = 0; n != N; n++) lambdas[n] = lambda_births[n];

    return List::create(Named("nodemembership") = IntegerVector(nodeMembership.begin(), nodeMembership.end()),
                        Named("lambda") = lambdas,
                        Named("parent") = IntegerVector(parent.begin(), parent.end()),
                        Named("stability") = NumericVector(stabilities.begin(), stabilities.end()),
                        Named("selected") = IntegerVector(selected.begin(), selected.end()),
                        Named("coredistances") = wrap(vector<double>(coreDistances, coreDistances + N)));
  }
};

// [[Rcpp::export]]
List hdbscanc(const arma::sp_mat& edges,
              const IntegerMatrix& neighbors,
              const int& K,
              const int& minPts,
              Rcpp::Nullable<Rcpp::NumericVector> threads,
              const bool& verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
  HDBSCAN object = HDBSCAN(edges.n_cols, verbose);
 // 1 N
  IntegerVector tree = object.build(K, edges, neighbors);
  NumericMatrix clusters = NumericMatrix(2, edges.n_cols);
  object.condenseAndExtract(minPts, REAL(clusters));
  List hierarchy = object.reportHierarchy();
  return List::create(Named("clusters") = clusters,
                      Named("tree") = IntegerVector(tree),
                      Named("hierarchy") = hierarchy);
}
