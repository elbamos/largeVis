// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "pq.h"

class HDBSCAN : public UF<long long> {
public:

  HDBSCAN(const int N,
          const long long nedges,
          bool verbose) : UF(N, nedges, verbose) {

  }

  void makeCoreDistances(const arma::sp_mat& edges, const int K) {
    coreDistances = arma::vec(N);
    for (long long n = 0; n < N; n++) //if (p.increment())
    	{
      DistanceSorter srtr = DistanceSorter();
      for (auto it = edges.begin_row(n);
           it != edges.end_row(n);
           it++) srtr.emplace(iddist(it.col(), *it));
      int k = 0;
      for (k = 0; k != K && srtr.size() > 1; k++) srtr.pop();
      if (k != K) stop("Insufficient neighbors.");
      coreDistances[n] = srtr.top().second;
    }
  }

  void makeCoreDistances(const arma::sp_mat& edges,
                         const IntegerMatrix& neighbors,
                         const int K) {
    coreDistances = arma::vec(N);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long long n = 0; n < N; n++) {
      coreDistances[n] = edges(neighbors(K, n), n);
    }
    p.increment(N);
  }

	void primsAlgorithm(const sp_mat& edges) {
		UF<long long>::primsAlgorithm(edges, 0);
	}

  arma::mat process(const int minPts) {
  	buildHierarchy();
    condense(minPts);
    determineStability(minPts);
    extractClusters();
    return getClusters();
  }

  Rcpp::List reportHierarchy() {
    long long survivingClusterCnt = survivingClusters.size();
    IntegerVector parent = IntegerVector(survivingClusterCnt);
    IntegerVector nodeMembership = IntegerVector(N);
    NumericVector stabilities = NumericVector(survivingClusterCnt);
    IntegerVector selected = IntegerVector(survivingClusterCnt);
    std::set< long long >::iterator it;
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp single private(it)
{
#endif
    for (it = roots.begin();
         it != roots.end();
         it++)
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

    return Rcpp::List::create(Rcpp::Named("nodemembership") = nodeMembership,
                              Rcpp::Named("lambda") = lambdas,
                              Rcpp::Named("parent") = parent,
                              Rcpp::Named("stability") = stabilities,
                              Rcpp::Named("selected") = selected);
  }

	long long* getMinimumSpanningTree() {
		return minimum_spanning_tree.get();
	}
};

// [[Rcpp::export]]
List hdbscanc(const arma::sp_mat& edges,
              Rcpp::Nullable< Rcpp::IntegerMatrix > neighbors,
              const int K,
              const int minPts,
              Rcpp::Nullable<Rcpp::NumericVector> threads,
              const bool verbose) {
#ifdef _OPENMP
	checkCRAN(threads);
#endif
  HDBSCAN object = HDBSCAN(edges.n_cols, edges.n_elem, verbose);
  if (neighbors.isNotNull()) {
    IntegerMatrix neigh = IntegerMatrix(neighbors);
    object.makeCoreDistances(edges, neigh, K);
  } else {
    object.makeCoreDistances(edges, K);
  }
  object.primsAlgorithm(edges);
  arma::mat clusters = object.process(minPts);
  arma::ivec tree = arma::ivec(edges.n_cols);
  long long* mst = object.getMinimumSpanningTree();
  for (int n = 0; n != edges.n_cols; n++) {
    tree[n] = mst[n];
  }
  Rcpp::List hierarchy = object.reportHierarchy();
  return Rcpp::List::create(Rcpp::Named("clusters") = clusters,
                            Rcpp::Named("tree") = tree,
                            Rcpp::Named("hierarchy") = hierarchy);
}