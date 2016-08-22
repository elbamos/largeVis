// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "pq.h"

class HDBSCAN : UF<long long> {
public:
  class CompareDist {
  public:
    bool operator()(iddist n1, iddist n2) {
      return n1.second > n2.second;
    }
  };
  typedef std::priority_queue<iddist,
                              std::vector<iddist>,
                              CompareDist> DistanceSorter;

  long long starterIndex;
  arma::sp_mat mrdMatrix;

  std::unique_ptr<long long[]>   minimum_spanning_tree;
  std::unique_ptr<double[]>    minimum_spanning_distances;

  HDBSCAN(const int N,
          const long long nedges,
          bool verbose) : UF(N, nedges, verbose) {
            minimum_spanning_tree = std::unique_ptr<long long[]>(new long long[N]);
            minimum_spanning_distances = std::unique_ptr<double[]>(new double[N]);
          }

  arma::vec coreDistances;

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

  void makeMRD(const arma::sp_mat& edges, const int K) {
    mrdMatrix = arma::sp_mat(edges);
    double bestOverall = INFINITY;
    arma::sp_mat::iterator it;
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp single private(it)
{
#endif
    for (it = mrdMatrix.begin();
         it != mrdMatrix.end();
         it++) if (p.increment())
#pragma omp task
         	{
      long long i = it.row();
      long long j = it.col();
      double d = max(coreDistances[i], coreDistances[j]);
      d = (d > *it) ? d : *it;
      if (d < bestOverall) {
        starterIndex = i;
        bestOverall = d;
      }
      *it = d;
    }
#ifdef _OPENMP
}
}
#endif
  }

  void primsAlgorithm() {
 //   double* Cv = minimum_spanning_distances;
 //   long long* Ev = minimum_spanning_tree;
    MinIndexedPQ<long long, double> Q = MinIndexedPQ<long long, double>(N);
    for (long long n = 0; n != N; n++) {
      minimum_spanning_distances[n] = (n == starterIndex) ? -1 : INFINITY;
      minimum_spanning_tree[n] = -1;
      Q.insert(n, minimum_spanning_distances[n]);
    }

    long long v;
    while (! Q.isEmpty() && p.increment()) {
      v = Q.deleteMin();
      for (auto it = mrdMatrix.begin_row(v);
           it != mrdMatrix.end_row(v);
           it++) {
        long long w = it.col();
        if (Q.contains(w) && *it < minimum_spanning_distances[w]) {
          Q.changeKey(w, *it);
          minimum_spanning_distances[w] = *it;
          minimum_spanning_tree[w] = v;
        }
      }
    }
  }

  void buildHierarchy() {
    setup();
    DistanceSorter srtr = DistanceSorter();
    for (long long n = 0; n != N && p.increment(); n++) {
      double distance = minimum_spanning_distances[n];
      srtr.push(iddist(n, distance));
    }

    while (! srtr.empty() && p.increment()) {
      const iddist ijd = srtr.top();
      srtr.pop();
      agglomerate(ijd.first, minimum_spanning_tree[ijd.first], ijd.second);
    }
  }

  arma::mat process(const int minPts) {
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
    	reportAHierarchy(*it, *it, parent, nodeMembership,
         stabilities, selected, 0, 0);
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
};

// [[Rcpp::export]]
List hdbscanc(const arma::sp_mat& edges,
              Rcpp::Nullable< Rcpp::IntegerMatrix > neighbors,
              const int K,
              const int minPts,
              const bool verbose) {
  HDBSCAN object = HDBSCAN(edges.n_cols, edges.n_elem, verbose);
  if (neighbors.isNotNull()) {
    IntegerMatrix neigh = IntegerMatrix(neighbors);
    object.makeCoreDistances(edges, neigh, K);
  } else {
    object.makeCoreDistances(edges, K);
  }
  object.makeMRD(edges, K);
  object.primsAlgorithm();
  object.buildHierarchy();
  arma::mat clusters = object.process(minPts);
  arma::ivec tree = arma::ivec(edges.n_cols);
  for (int n = 0; n != edges.n_cols; n++) {
    tree[n] = object.minimum_spanning_tree[n];
  }
  Rcpp::List hierarchy = object.reportHierarchy();
  return Rcpp::List::create(Rcpp::Named("clusters") = clusters,
                            Rcpp::Named("tree") = tree,
                            Rcpp::Named("hierarchy") = hierarchy);
}